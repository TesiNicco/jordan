#!/usr/bin/Rscript

# This script should be used to make PRSs
# We do PRS given:
# 1. genotype files (PLINK or VCF)
# 2. a file with snps and relative weights

# Libraries
    library(argparse)
    library(data.table)
    library(stringr)
    library(ggplot2)

# Functions
    # Function to check output file
    checkOutputFile = function(outname){
        # check if directory exists otherwise create it
        if (dir.exists(outname)){
            cat('** Output directory is valid and exists.\n')
            run=TRUE
            res = list(run, outname)            
        } else {
            tryCatch({
                system(paste0('mkdir ', outname))
                cat('** Output directory does not exist. Will create.\n')
                run = TRUE
                res = list(run, outname)
            }, error = function(e) {
                print(paste("** Can not create output directory. Path is invalid.\n", conditionMessage(e)))
                run = FALSE
                res = list(run, outname)
            })            
        }
        return(res)
    }

    # Function to check input genotype file
    checkGenoFile = function(genotype_file, outdir, isdosage){
        # check if file is a vcf/bcf
        if (endsWith(genotype_file, 'vcf') | endsWith(genotype_file, 'vcf.gz') | endsWith(genotype_file, 'bcf') | endsWith(genotype_file, 'bcf.gz')){
            # vcf file
            if (isdosage == TRUE){
                cat('** VCF file found. Converting to PLINK assuming it is imputed data from Minimac-4 (dosage=HDS).\n')
                system(paste0('plink2 --vcf ', genotype_file, ' dosage=HDS --make-pgen --out ', outdir, '/tmp'))
            } else {
                cat('** VCF file found. Converting to PLINK.\n')
                system(paste0('plink2 --vcf ', genotype_file, ' --make-pgen --out ', outdir, '/tmp'))
            }
            run = TRUE
            data_path = paste0(outdir, '/tmp')
            res = list(run, data_path, 'plink2')
        } else if (file.exists(paste0(genotype_file, '.pvar'))){
            cat('** PLINK2 file found.\n')
            run = TRUE
            res = list(run, genotype_file, 'plink2')
        } else if (file.exists(paste0(genotype_file, '.bim'))){
            cat('** PLINK file found.\n')
            run = TRUE
            res = list(run, genotype_file, 'plink')
        } else {
            cat('** Unknown genotype file provided. \n')
            run = FALSE
            res = list(run, genotype_file, NA)
        }
        return(res)
    }

    # Function to check input genotype file
    checkInputSNPs = function(snps_file){
        # check if file exists
        if (file.exists(snps_file)){
            # open file
            snps_data = data.table::fread(snps_file, h=T, stringsAsFactors=F, sep="\t")
            # check column names
            if ('CHROM' %in% colnames(snps_data)){
                if ('POS' %in% colnames(snps_data)){
                    if ('EFFECT_ALLELE' %in% colnames(snps_data)){
                        if ('OTHER_ALLELE' %in% colnames(snps_data)){
                            if ('BETA' %in% colnames(snps_data)){
                                run = TRUE
                                res = list(run, snps_data)
                                cat('** Required columns found\n')
                            } else if ('OR' %in% colnames(snps_data)){
                                # convert to BETA
                                snps_data$BETA = log(as.numeric(snps_data$OR))
                                run = TRUE
                                res = list(run, snps_data)
                                cat('** Required columns found\n')            
                            } else {
                                cat('** No BETA column found in SNP data.\n')
                                run = FALSE
                                res = list(run, NA)
                            }
                        } else {
                            cat('** No OTHER_ALLELE column found in SNP data.\n')
                            run = FALSE
                            res = list(run, NA)
                        }
                    } else {
                        cat('** No EFFECT_ALLELE column found in SNP data.\n')
                        run = FALSE
                        res = list(run, NA)
                    }
                } else {
                    cat('** No POS column found in SNP data.\n')
                    run = FALSE
                    res = list(run, NA)
                }
            } else {
                cat('** No CHROM column found in SNP data.\n')
                run = FALSE
                res = list(run, NA)
            }
        } else {
            cat('** SNP data file does not exist.\n')
            run = FALSE
            res = list(run, NA)
        }
        return(res)
    }

    # Function to guide PRS
    makePRS = function(outdir, genotype_path, snps_data, genotype_type){
        # match ids in the plink file and extract dosages/genotypes
        res = matchIDs(genotype_path, snps_data, genotype_type, outdir)
        dosages = res[[1]]
        mappingSnp = res[[2]]
        # snps info flip to risk allele
        snps_data$risk_allele = snps_data$EFFECT_ALLELE
        snps_data$risk_beta = abs(snps_data$BETA)
        snps_data$risk_allele[which(snps_data$BETA <0)] = snps_data$OTHER_ALLELE[which(snps_data$BETA <0)]
        # then do the prs
        res = prs(snps_data, dosages, mappingSnp)
        return(res)
    }

    # Function to match snp IDs and extract genotypes
    matchIDs = function(genotype_path, snps_data, genotype_type, outdir){
        # read bim/pvar file
        if (genotype_type == 'plink'){
            snpsinfo = data.table::fread(paste0(genotype_path, '.bim'), h=F, stringsAsFactors=F)
            colnames(snpsinfo) = c('chr', 'id', 'na', 'pos', 'ref', 'alt')
        } else {
            snpsinfo = data.table::fread(paste0(genotype_path, '.pvar'), h=T, stringsAsFactors=F)
            colnames(snpsinfo)[1:5] = c('chr', 'pos', 'id', 'ref', 'alt')
        }
        # create identifier with chrom:pos
        snpsinfo$unique_id = paste(snpsinfo$chr, snpsinfo$pos, sep=":")
        # check which snps of interest are in there
        snps_data$id = paste(snps_data$CHROM, snps_data$POS, sep=":")
        # first match snps based on their chromosome:position
        matchingsnps_lv1 = snpsinfo[which(snpsinfo$unique_id %in% snps_data$id),]
        # then match snps based on their alleles
        matchingsnps_lv1$allele_match = NA
        for (i in 1:nrow(matchingsnps_lv1)){
            # get variant info from data
            tmp_info = as.character(matchingsnps_lv1[i, c('ref', 'alt')])
            # get original snp info
            ori_info = as.character(snps_data[which(snps_data$id == matchingsnps_lv1$unique_id[i]), c('EFFECT_ALLELE', 'OTHER_ALLELE')])
            # check if they are the same
            if (length(table(tmp_info %in% ori_info)) == 1){
                # then ok
                matchingsnps_lv1$allele_match[i] = 'yes'
            }
        }
        # filter out snps where the allele doesn't match
        matchingsnps = matchingsnps_lv1[which(matchingsnps_lv1$allele_match == 'yes'),]
        # write this output
        if (genotype_type == 'plink'){
            write.table(matchingsnps[, c('chr', 'id', 'na', 'pos', 'ref', 'alt')], paste0(outdir, '/snpsInterest.txt'), quote=F, row.names=F, col.names=F, sep="\t")
        } else {
            write.table(matchingsnps[, c('chr', 'pos', 'id', 'ref', 'alt')], paste0(outdir, '/snpsInterest.txt'), quote=F, row.names=F, col.names=c('#CHROM', 'POS', 'ID', 'REF', 'ALT'), sep="\t")
        }
        # extract these snps
        if (genotype_type == 'plink'){
            system(paste0('plink --bfile ', genotype_path, ' --extract ', outdir, '/snpsInterest.txt --export A include-alt --out ', outdir, '/dosages'))
        } else {
            system(paste0('plink2 --pfile ', genotype_path, ' --extract ', outdir, '/snpsInterest.txt --export A include-alt --out ', outdir, '/dosages'))
        }
        # read dosages
        dos = data.table::fread(paste0(outdir, '/dosages.raw'), h=T, stringsAsFactors=F)
        res = list(dos, matchingsnps)
        return(res)
    }

    # Function to match risk alleles
    prs = function(snps_data, dosages, mappingSnp){
        prs_df = data.frame(iid = dosages$IID, PRS = 0)
        included_snps = data.frame()
        dosages = data.frame(dosages, check.names=F)
        # iterate over dosages columns
        for (i in seq(7, ncol(dosages))){
            # extract data
            temp = dosages[, c(1, 2, 3, 4, 5, 6, i)]
            # check alleles
            temp_name = stringr::str_split_fixed(colnames(dosages)[i], '_', 2)[, 1]
            temp_alleles = stringr::str_split_fixed(colnames(dosages)[i], '_', 2)[, 2]
            temp_effect = stringr::str_split_fixed(temp_alleles, '\\(', 2)[, 1]
            temp_other = stringr::str_replace_all(stringr::str_replace_all(stringr::str_split_fixed(temp_alleles, '\\(', 2)[, 2], '/', ''), '\\)', '')
            # get the info from the snpdata
            temp_snpinfo = mappingSnp[which(mappingSnp$id == temp_name),]
            temp_snpdata = snps_data[which(snps_data$CHROM == temp_snpinfo$chr & snps_data$POS == temp_snpinfo$pos),]
            # check if the alleles are ok
            if (temp_effect %in% c(temp_snpdata$EFFECT_ALLELE, temp_snpdata$OTHER_ALLELE) & temp_other %in% c(temp_snpdata$EFFECT_ALLELE, temp_snpdata$OTHER_ALLELE)){
                # flip to risk
                if (temp_snpdata$risk_allele == temp_effect){
                    temp$score = temp[, 7] * temp_snpdata$risk_beta
                } else {
                    temp$score = (2 - temp[, 7]) * temp_snpdata$risk_beta
                }
                # add to main df
                prs_df$PRS = prs_df$PRS + temp$score
                # also save the included snps
                included_snps = rbind(included_snps, data.frame(SNP = temp_name, BETA = temp_snpdata$risk_beta, ALLELE = temp_snpdata$risk_allele, OTHER_ALLELE = temp_snpdata$OTHER_ALLELE, TYPE = 'Included', CHROM = temp_snpinfo$chr, POS = temp_snpinfo$pos))
            } else {
                cat('** Alleles do not match for snp ', temp_name, '. Skipping.\n')
                included_snps = rbind(included_snps, data.frame(SNP = temp_name, BETA = temp_snpdata$risk_beta, ALLELE = temp_snpdata$risk_allele, OTHER_ALLELE = temp_snpdata$OTHER_ALLELE, TYPE = 'Excluded', CHROM = temp_snpinfo$chr, POS = temp_snpinfo$pos))
            }
        }
        res = list(prs_df, included_snps)
        return(res)
    }

    # Function to draw plot
    makePlot = function(prs_df, included_snps, outdir){
        # plot density
        pdf(paste0(outdir, '/PRS_density.pdf'), height = 7, width = 7)
        plt_density = ggplot(data = prs_df, aes(x = PRS, fill = 'red')) + geom_density() + theme(legend.position = 'none') + xlab('Polygenic Risk Score') + ylab('Density') + ggtitle(paste0('PRS based on ', nrow(included_snps), ' SNPs'))
        print(plt_density)
        dev.off()
        # plot snps
        pdf(paste0(outdir, '/PRS_SNPs.pdf'), height = 7, width = 7)
        plt_snps = ggplot(data = included_snps[which(included_snps$TYPE == 'Included'),], aes(y = SNP, x = BETA)) + geom_point(stat = 'identity') + xlab('Beta') + ylab('SNP') + ggtitle('Beta of SNPs included in PRS')
        print(plt_snps)
        dev.off()     
        return('Plots are done!')
    }

# Parse arguments
    # Add required arguments
        parser <- ArgumentParser(description = "makePRS: a script to make PRS in R")
        parser$add_argument("--genotype", help = "Path to the genotype file. If PLINK, do not provide file extension. If VCF, please provide the file extension.")
        parser$add_argument("--snplist", help = "Path to the SNPs file. This file will define the SNPs to include in the PRS and the relative weights. By default, required column names are CHROM, POS, EFFECT_ALLELE, OTHER_ALLELE and BETA, (or OR)")
        parser$add_argument("--outname", help = "Path to output PRS file. A new directory can be specified. In that case, the directory will be created and the file with be written to the new directory.")
        parser$add_argument("--isdosage", help="Whether Input data is imputed or genotyped. The information is used to read genotypes in PLINK.", default=FALSE)
        parser$add_argument("--plot", help="Whether to plot (default = FALSE) or not the PRS densities.", default=FALSE)
    # Read arguments
        args <- parser$parse_args()
        genotype_file <- args$genotype
        snps_file <- args$snplist
        outfile <- args$outname
        isdosage <- args$isdosage
        plt <- args$plot
    # Print arguments on screen
        cat("\nGenotype file: ", genotype_file)
        cat("\nSNPs file: ", snps_file)
        cat("\nOutput file: ", outfile)
        cat("\nDosage: ", isdosage)
        cat("\nPlot: ", plt, '\n\n')

# Check inputs
    # Check output directory
        res = checkOutputFile(outfile)
        run1 = res[[1]]
        outdir = res[[2]]
    # Check input genotype file
        res = checkGenoFile(genotype_file, outdir, isdosage)
        run2 = res[[1]]
        genotype_path = res[[2]]
        genotype_type = res[[3]]
    # Check input snplist
        res = checkInputSNPs(snps_file)
        run3 = res[[1]]
        snps_data = res[[2]]
    # Check before running
        if (run1 == FALSE | run2 == FALSE | run3 == FALSE){
            stop("** Inputs are not valid. Check above messages and try again.\n\n")
        } else {
            cat('** Inputs are valid. Starting the script.\n\n')
            res = makePRS(outdir, genotype_path, snps_data, genotype_type)
            prs_df = res[[1]]
            included_snps = res[[2]]
            # write prs output
            write.table(prs_df, paste0(outdir, '/PRS_table.txt'), quote=F, row.names=F, sep="\t")
            # fill included snps file with those that are excluded as well
            #included_snps$POS = stringr::str_split_fixed(included_snps$SNP, ':', 4)[, 2]
            #included_snps$CHROM = stringr::str_replace_all(stringr::str_split_fixed(included_snps$SNP, ':', 4)[, 1], 'chr', '')
            included_snps$ID = paste(included_snps$CHROM, included_snps$POS, sep=":")
            snps_data$ID = paste(stringr::str_replace_all(snps_data$CHROM, 'chr', ''), snps_data$POS, sep=":")
            missing = snps_data[which(!(snps_data$ID %in% included_snps$ID)),]
            if (nrow(missing) >0){ included_snps = rbind(included_snps, data.frame(SNP = NA, BETA = missing$BETA, ALLELE = missing$EFFECT_ALLELE, OTHER_ALLELE = missing$OTHER_ALLELE, TYPE = 'Excluded', POS = missing$POS, CHROM = missing$CHROM, ID = missing$ID)) }
            # write these outputs
            write.table(included_snps, paste0(outdir, '/SNPs_included_PRS.txt'), quote=F, row.names=F, sep="\t")
            # clean temporary data
            system(paste0('rm ', outdir, '/tmp*'))
            system(paste0('rm ', outdir, '/dosage*'))
            # check if plot needs to be done
            if (plt == TRUE){
                makePlot(prs_df, included_snps, outdir)
            }
            # end message
            cat('\n\n** Analysis over. Ciao! \n\n')
        }
