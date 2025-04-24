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
    checkGenoFile = function(genotype_file, outdir, isdosage, multiple){
        # check if file is a vcf/bcf
        if (endsWith(genotype_file, 'vcf') | endsWith(genotype_file, 'vcf.gz') | endsWith(genotype_file, 'bcf') | endsWith(genotype_file, 'bcf.gz')){
            # vcf file
            if (isdosage == TRUE){
                cat('** VCF file found. Converting to PLINK assuming it is imputed data from Minimac-4 (dosage=HDS).\n')
                system(paste0('plink2 --vcf ', genotype_file, ' dosage=HDS --make-pgen --out ', outdir, '/tmp > /dev/null 2>&1'))
            } else {
                cat('** VCF file found. Converting to PLINK.\n')
                system(paste0('plink2 --vcf ', genotype_file, ' --make-pgen --out ', outdir, '/tmp > /dev/null 2>&1'))
            }
            run = TRUE
            data_path = paste0(outdir, '/tmp')
            res = list(run, data_path, 'plink2')
        } else if (file.exists(paste0(genotype_file, '.pvar'))){
            genotype_file = paste0(genotype_file, '.pvar')
            cat('** PLINK2 file found.\n')
            run = TRUE
            # check other chromosomes if present
            if (multiple == TRUE){
                all_files = system(paste0('ls ', dirname(genotype_file), '/*pvar'), intern=TRUE)
                genotype_file = all_files
            }
            res = list(run, genotype_file, 'plink2')
        } else if (file.exists(paste0(genotype_file, '.bim'))){
            genotype_file = paste0(genotype_file, '.bim')
            cat('** PLINK file found.\n')
            run = TRUE
            if (multiple == TRUE){
                all_files = system(paste0('ls ', dirname(genotype_file), '/*bim'), intern=TRUE)
                genotype_file = all_files
            }
            res = list(run, genotype_file, 'plink')
        } else {
            cat('** Unknown genotype file provided. \n')
            run = FALSE
            res = list(run, genotype_file, NA)
        }
        return(res)
    }

    # Function to check input genotype file
    checkInputSNPs = function(snps_file, addWeight){
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
                                if (addWeight != FALSE){
                                    if (addWeight %in% colnames(snps_data)){
                                        run = TRUE
                                        res = list(run, snps_data)
                                        cat('** Required columns found\n')
                                    } else {
                                        cat('** Additional weight was selected, but column was not found in SNP data.\n')
                                        run = FALSE
                                        res = list(run, NA)
                                    }
                                } else {
                                    run = TRUE
                                    res = list(run, snps_data)
                                    cat('** Required columns found\n')
                                }
                            } else if ('OR' %in% colnames(snps_data)){
                                # convert to BETA
                                snps_data$BETA = log(as.numeric(snps_data$OR))
                                if (addWeight != FALSE){
                                    if (addWeight %in% colnames(snps_data)){
                                        run = TRUE
                                        res = list(run, snps_data)
                                        cat('** Required columns found\n')
                                    } else {
                                        cat('** Additional weight was selected, but column was not found in SNP data.\n')
                                        run = FALSE
                                        res = list(run, NA)
                                    }
                                } else {
                                    run = TRUE
                                    res = list(run, snps_data)
                                    cat('** Required columns found\n')
                                }
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
    makePRS = function(outdir, genotype_path, snps_data, genotype_type, multiple, excludeAPOE, maf, fliprisk, keepDos, addWeight, freq){
        # match ids in the plink file and extract dosages/genotypes
        cat('**** Matching SNPs and extracting dosages.\n')
        # Match variants of interest and get the dosages
        res = matchIDs_multiple(genotype_path, snps_data, genotype_type, outdir, maf, freq)
        dosages = res[[1]]
        mappingSnp = res[[2]]
	    # if dosages need to be outputted, write them now
	    if (keepDos == TRUE){
		    write.table(dosages, paste0(outdir, '/chrAll_dosages.txt'), quote=F, row.names=F, sep="\t")
	    }
        # check for nas and in case exclude them as well as update the mappingSNP
        dosages = excludeNAs(dosages)
        # snps info flip to risk allele if this was requested
        snps_data$risk_allele = snps_data$EFFECT_ALLELE
	    if (fliprisk == FALSE){
		    snps_data$risk_beta = abs(snps_data$BETA)
	        snps_data$risk_allele[which(snps_data$BETA <0)] = snps_data$OTHER_ALLELE[which(snps_data$BETA <0)]
	    } else {
		    snps_data$risk_beta = snps_data$BETA
	    }
        # then do the prs
        cat('**** Calculating PRS.\n')
        res = prs(snps_data, dosages, mappingSnp, addWeight)
        # if requested, do without apoe as well
        if (excludeAPOE != FALSE){
            cat('**** Removing APOE SNPs and re-calculating PRS.\n')
            snps_data_noAPOE = snps_data[which(!(snps_data$POS %in% c(44908684, 44908822))),]
            res_noapoe = prs(snps_data_noAPOE, dosages, mappingSnp, addWeight)
            return(list(res, res_noapoe))
        } else {
            return(res)
        }
    }

    # Function to check for NAs
    excludeNAs = function(dosages){
        # check for NAs
        if (any(is.na(dosages))){
            cat('** There are NAs in the dosages. Data is likely from WGS. Excluding them.\n')
            # make sure dosages is a dataframe
            dosages = data.frame(dosages, check.names=F)
            # remove the columns with NAs
            df_clean <- dosages[, colSums(is.na(dosages)) == 0]
            return(df_clean)
        } else {
            return(dosages)
        }
    }

    matchIDs_multiple = function(genotype_path, snps_data, genotype_type, outdir, maf, freq){
        # container for all dosages and matching snps and frequencies
        matchingsnps_all = data.frame()
        all_dos = data.frame()
        all_freq = data.frame()
        # write positions of interest
        write.table(snps_data$POS, paste0(outdir, '/tmp_positions.txt'), quote=F, row.names=F, col.names=F)
        # iterate over all files
        for (f in genotype_path){
            tryCatch({
                # read bim/pvar file
                if (genotype_type == 'plink'){
                    #snpsinfo = data.table::fread(f, h=F, stringsAsFactors=F)
                    snpsinfo = data.frame(str_split_fixed(system(paste0('grep -f ', outdir, '/tmp_positions.txt ', f), intern=T, ignore.stderr=T), '\t', 7))
                    snpsinfo <- snpsinfo[, apply(snpsinfo, 2, function(x) any(x != "" & !is.na(x)))]
                    colnames(snpsinfo) = c('chr', 'id', 'na', 'pos', 'ref', 'alt')
                } else {
                    #snpsinfo = data.table::fread(f, h=T, stringsAsFactors=F)
                    snpsinfo = suppressWarnings(data.frame(str_split_fixed(system(paste0('grep -f ', outdir, '/tmp_positions.txt ', f), intern=T), '\t', 7)))
                    snpsinfo <- snpsinfo[, apply(snpsinfo, 2, function(x) any(x != "" & !is.na(x)))]
                    colnames(snpsinfo)[1:5] = c('chr', 'pos', 'id', 'ref', 'alt')
                }
                # first rough filter on the position
                snpsinfo = snpsinfo[which(snpsinfo$pos %in% snps_data$POS),]
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
                    tmp_info = tmp_info[order(tmp_info)]
                    # get original snp info
                    ori_info = as.character(snps_data[which(snps_data$id == matchingsnps_lv1$unique_id[i]), c('EFFECT_ALLELE', 'OTHER_ALLELE')])
                    ori_info = ori_info[order(ori_info)]
                    # check if they are the same
                    if (identical(tmp_info, ori_info)){
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
                    if (maf == TRUE){
                        # Extract dosages
                        system(paste0('plink --bfile ', str_replace_all(f, '.bim', ''), ' --extract ', outdir, '/snpsInterest.txt --maf ', maf, ' --export A include-alt --out ', outdir, '/dosages > /dev/null 2>&1'))
                        # Calculate frequencies if requested
                        if (freq == TRUE){
                            system(paste0('plink --bfile ', str_replace_all(f, '.bim', ''), ' --extract ', outdir, '/snpsInterest.txt --maf ', maf, ' --freq --out ', outdir, '/frequencies > /dev/null 2>&1'))
                            # Read frequencies
                            freq_file = data.table::fread(paste0(outdir, '/frequencies.frq'), h=T, stringsAsFactors=F)
                            # Add to dataframe
                            all_freq = rbind(all_freq, freq_file)
                        }
                    } else {
                        # Extract dosages
                        system(paste0('plink --bfile ', str_replace_all(f, '.bim', ''), ' --extract ', outdir, '/snpsInterest.txt --export A include-alt --out ', outdir, '/dosages > /dev/null 2>&1'))
                        # Calculate frequencies if requested
                        if (freq == TRUE){
                            system(paste0('plink --bfile ', str_replace_all(f, '.bim', ''), ' --extract ', outdir, '/snpsInterest.txt --freq --out ', outdir, '/frequencies > /dev/null 2>&1'))
                            # Read frequencies
                            freq_file = data.table::fread(paste0(outdir, '/frequencies.frq'), h=T, stringsAsFactors=F)
                            # Add to dataframe
                            all_freq = rbind(all_freq, freq_file)
                        }
                    }
                } else {
                    if (maf == TRUE){
                        # Extract dosages
                        system(paste0('plink2 --pfile ', str_replace_all(f, '.pvar', ''), ' --extract ', outdir, '/snpsInterest.txt --maf ', maf, ' --export A include-alt --out ', outdir, '/dosages > /dev/null 2>&1'))
                        # Calculate frequencies if requested
                        if (freq == TRUE){
                            system(paste0('plink2 --pfile ', str_replace_all(f, '.pvar', ''), ' --extract ', outdir, '/snpsInterest.txt --maf ', maf, ' --freq --out ', outdir, '/frequencies > /dev/null 2>&1'))
                            # Read frequencies
                            freq_file = data.table::fread(paste0(outdir, '/frequencies.frq'), h=T, stringsAsFactors=F)
                            # Add to dataframe
                            all_freq = rbind(all_freq, freq_file)
                        }
                    } else {
                        # Extract dosages
                        system(paste0('plink2 --pfile ', str_replace_all(f, '.pvar', ''), ' --extract ', outdir, '/snpsInterest.txt --export A include-alt --out ', outdir, '/dosages > /dev/null 2>&1'))
                        # Calculate frequencies if requested
                        if (freq == TRUE){
                            system(paste0('plink2 --pfile ', str_replace_all(f, '.pvar', ''), ' --freq --out ', outdir, '/frequencies > /dev/null 2>&1'))
                            # Read frequencies
                            freq_file = data.table::fread(paste0(outdir, '/frequencies.frq'), h=T, stringsAsFactors=F)
                            # Add to dataframe
                            all_freq = rbind(all_freq, freq_file)
                        }
                    }
                }
                # read dosages
                dos = data.table::fread(paste0(outdir, '/dosages.raw'), h=T, stringsAsFactors=F)
                dos$FID = NULL
                dos$PAT = NULL
                dos$MAT = NULL
                dos$SEX = NULL
                dos$PHENOTYPE = NULL
                # add to all dosages
                if (nrow(all_dos) == 0){
                    all_dos = dos
                } else {
                    all_dos = merge(all_dos, dos, by = 'IID') 
                }
                matchingsnps_all = rbind(matchingsnps_all, matchingsnps)
                # remove temporary files
                system(paste0('rm ', outdir, '/snpsInterest.txt'), ignore.stderr=T, ignore.stdout=T)
                system(paste0('rm ', outdir, '/dosages.raw'), ignore.stderr=T, ignore.stdout=T)
                system(paste0('rm ', outdir, '/frequencies.*'), ignore.stderr=T, ignore.stdout=T)
            },
            error = function(e){
                cat(paste0('****** An error occurred with file ', f, ': Skipping this chromosome.\n'))
            })
        }
        # Write frequencies
        if (nrow(all_freq) > 0){
            write.table(all_freq, paste0(outdir, '/frequencies.txt'), quote=F, row.names=F, sep="\t")
        }
        res = list(all_dos, matchingsnps_all)
        return(res)
    }

    # Function to match risk alleles
    prs = function(snps_data, dosages, mappingSnp, addWeight){
        prs_df = data.frame(iid = dosages$IID, PRS = 0)
        included_snps = data.frame()
        dosages = data.frame(dosages, check.names=F)
        # iterate over dosages columns
        for (i in seq(2, ncol(dosages))){
            # extract data
            temp = dosages[, c(1, i)]
            # check alleles
            temp_name = stringr::str_split_fixed(colnames(dosages)[i], '_', 2)[, 1]
            temp_alleles = stringr::str_split_fixed(colnames(dosages)[i], '_', 2)[, 2]
            temp_effect = stringr::str_split_fixed(temp_alleles, '\\(', 2)[, 1]
            temp_other = stringr::str_replace_all(stringr::str_replace_all(stringr::str_split_fixed(temp_alleles, '\\(', 2)[, 2], '/', ''), '\\)', '')
            # get the info from the snpdata
            temp_snpinfo = mappingSnp[which(mappingSnp$id == temp_name),]
            temp_snpdata = snps_data[which(snps_data$CHROM == temp_snpinfo$chr & snps_data$POS == temp_snpinfo$pos),]
            # get additional weight if requested
            if (addWeight != FALSE){
                addWeight_num = as.numeric(temp_snpdata[, ..addWeight])
            } else {
                addWeight_num = 1
            }
            # check if there are lines
            if (nrow(temp_snpdata) >0){
                # check if the alleles are ok
                if (temp_effect %in% c(temp_snpdata$EFFECT_ALLELE, temp_snpdata$OTHER_ALLELE) & temp_other %in% c(temp_snpdata$EFFECT_ALLELE, temp_snpdata$OTHER_ALLELE)){
                    # flip to risk
                    if (temp_snpdata$risk_allele == temp_effect){
                        temp$score = temp[, 2] * temp_snpdata$risk_beta * addWeight_num
                    } else {
                        temp$score = (2 - temp[, 2]) * temp_snpdata$risk_beta * addWeight_num
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

    # Function to write log file
    writeLog = function(outdir, genotype_file, snps_file, outfile, isdosage, plt, maf, multiple, excludeAPOE, fliprisk, keepDos, addWeight, freq){
        # Define output name
        outname = paste0(outdir, '/run_info.log')
        # Create log info
        info = paste0("Genotype file: ", genotype_file, "\nMultiple files: ", multiple, "\nSNPs file: ", snps_file, "\nOutput file: ", outfile, "\nDosage: ", isdosage, "\nMAF: ", maf, "\nWith and Without APOE: ", excludeAPOE, "\nDirect effects (Risk and Protective): ", fliprisk, "\nKeep dosages: ", keepDos, "\nAdditional weight: ", addWeight, "\nCalculate frequency: ", freq, "\nPlot: ", plt, '\n\n')
        # Write log file
        writeLines(info, outname)
        return(info)
    }