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
            res = outname
        } else {
            tryCatch({
                system(paste0('mkdir ', outname))
                cat('** Output directory does not exist. Will create.\n')
                res = outname
            }, error = function(e) {
                stop("** Error! Can not create output directory. Path is invalid.\n\n", call. = FALSE)
            })            
        }
        return(res)
    }

    # Function to check input genotype file
    checkGenoFile = function(genotype_file, outdir, isdosage, multiple, run){
        # check if file is a vcf/bcf
        if (file.exists(genotype_file) && (endsWith(genotype_file, 'vcf') | endsWith(genotype_file, 'vcf.gz') | endsWith(genotype_file, 'bcf') | endsWith(genotype_file, 'bcf.gz'))){
            # check if the file is a vcf/bcf
            vcftype = ifelse(endsWith(genotype_file, 'vcf') | endsWith(genotype_file, 'vcf.gz'), 'vcf', 'bcf')
            # vcf file
            if (isdosage == TRUE){
                cat('** VCF file found. Converting to PLINK assuming it is imputed data from Minimac-4 (dosage=HDS).\n')
                system(paste0('plink2 --', vcftype, ' ', genotype_file, ' dosage=HDS --make-pgen --out ', outdir, '/tmp > /dev/null 2>&1'))
            } else {
                cat('** VCF file found. Converting to PLINK.\n')
                system(paste0('plink2 --', vcftype, ' ', genotype_file, ' --make-pgen --out ', outdir, '/tmp > /dev/null 2>&1'))
            }
            data_path = paste0(outdir, '/tmp.pvar')
            res = list(data_path, 'plink2')
        } else if (file.exists(paste0(genotype_file, '.pvar'))){
            genotype_file = paste0(genotype_file, '.pvar')
            cat('** PLINK2 file found.\n')
            # check other chromosomes if present
            if (multiple == TRUE){
                all_files = system(paste0('ls ', dirname(genotype_file), '/*pvar'), intern=TRUE)
                genotype_file = all_files
            }
            res = list(genotype_file, 'plink2')
        } else if (file.exists(paste0(genotype_file, '.bim'))){
            genotype_file = paste0(genotype_file, '.bim')
            cat('** PLINK file found.\n')
            if (multiple == TRUE){
                all_files = system(paste0('ls ', dirname(genotype_file), '/*bim'), intern=TRUE)
                genotype_file = all_files
            }
            res = list(genotype_file, 'plink')
        } else {
            stop('** Genotype file not found. Please provide a valid file.\n\n', call. = FALSE)
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
            if ('CHROM' %in% toupper(colnames(snps_data))){
                # adjust the column by removing 'chr' if present
                if (any(grepl('chr', snps_data$CHROM))){ snps_data$CHROM = gsub('chr', '', snps_data$CHROM) }
                if ('POS' %in% toupper(colnames(snps_data))){
                    if ('EFFECT_ALLELE' %in% toupper(colnames(snps_data))){
                        if ('OTHER_ALLELE' %in% toupper(colnames(snps_data))){
                            if ('BETA' %in% toupper(colnames(snps_data))){
                                if (addWeight != FALSE){
                                    if (toupper(addWeight) %in% toupper(colnames(snps_data))){
                                        res = snps_data
                                        cat('** Required columns found\n')
                                    } else {
                                        stop('** Additional weight was selected, but column was not found in SNP data.\n\n', call. = FALSE)
                                    }
                                } else {
                                    res = snps_data
                                    cat('** Required columns found\n')
                                }
                            } else if ('OR' %in% toupper(colnames(snps_data))){
                                # convert to BETA
                                snps_data$BETA = log(as.numeric(snps_data$OR))
                                if (addWeight != FALSE){
                                    if (addWeight %in% colnames(snps_data)){
                                        res = snps_data
                                        cat('** Required columns found\n')
                                    } else {
                                        stop('** Additional weight was selected, but column was not found in SNP data.\n\n', call. = FALSE)
                                    }
                                } else {
                                    res = snps_data
                                    cat('** Required columns found\n')
                                }
                            } else {
                                stop('** No BETA or OR column found in SNP data.\n\n', call. = FALSE)
                            }
                        } else {
                            stop('** No OTHER_ALLELE column found in SNP data.\n\n', call. = FALSE)
                        }
                    } else {
                        stop('** No EFFECT_ALLELE column found in SNP data.\n\n', call. = FALSE)
                    }
                } else {
                    cat('** No POS column found in SNP data.\n\n', call. = FALSE)
                }
            } else {
                stop('** No CHROM column found in SNP data.\n\n', call. = FALSE)
            }
        } else {
            stop('** SNP data file does not exist.\n\n', call. = FALSE)
        }
        return(res)
    }

    # Function to define frequency mode
    defineFreqMode = function(assoc_file, assoc_info, freq){
        # check if frequency need to be done
        if (freq == TRUE){
            # then check if association need to be done
            if (assoc_file != FALSE){
                # extract association info
                assoc_vars = assoc_info[[3]]
                # case-control frequency only for binary traits --> model should be binomial
                if (any(grepl('binomial', assoc_vars))){
                    # then frequencies need to be calculated
                    cc_freq_var = paste(assoc_vars$variable[which(assoc_vars$mode == 'binomial')], collapse = ',')
                    return(cc_freq_var)
                } else {
                    # otherwise not
                    return(FALSE)
                }
            } else {
                return(FALSE)
            }
        } else {
            return(FALSE)
        }
    }

    # Function to guide PRS
    makePRS = function(outdir, genotype_path, snps_data, genotype_type, multiple, excludeAPOE, maf, fliprisk, keepDos, addWeight, freq, assoc_file, assoc_info){
        # decide whether frequencies in cases and controls need to be calculated
        freq_mode = defineFreqMode(assoc_file, assoc_info, freq)
        # match ids in the plink file and extract dosages/genotypes
        cat('**** Matching SNPs and extracting dosages.\n')
        # Match variants of interest and get the dosages
        res = matchIDs_multiple(genotype_path, snps_data, genotype_type, outdir, maf, freq, freq_mode)
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
        cat('\n**** Calculating PRS.\n')
        res = prs(snps_data, dosages, mappingSnp, addWeight)
        # if requested, do without apoe as well
        if (excludeAPOE != FALSE){
            cat('**** Removing APOE SNPs and re-calculating PRS.\n')
            snps_data_noAPOE = snps_data[which(!(snps_data$POS %in% c(44908684, 44908822))),]
            res_noapoe = prs(snps_data_noAPOE, dosages, mappingSnp, addWeight)
            return(list(res, res_noapoe, dosages))
        } else {
            return(list(res, dosages))
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

    matchIDs_multiple = function(genotype_path, snps_data, genotype_type, outdir, maf, freq, freq_mode){
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
                    snpsinfo = suppressWarnings(data.frame(str_split_fixed(system(paste0('grep -f ', outdir, '/tmp_positions.txt ', f), intern=T, ignore.stderr=T), '\t', 7)))
                    snpsinfo <- snpsinfo[, apply(snpsinfo, 2, function(x) any(x != "" & !is.na(x)))]
                    colnames(snpsinfo) = c('chr', 'id', 'na', 'pos', 'ref', 'alt')
                } else {
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
                        system(paste0('plink2 --bfile ', str_replace_all(f, '.bim', ''), ' --extract ', outdir, '/snpsInterest.txt --maf ', maf, ' --export A include-alt --out ', outdir, '/dosages > /dev/null 2>&1'))
                        # Calculate frequencies if requested
                        if (freq == TRUE){
                            system(paste0('plink2 --bfile ', str_replace_all(f, '.bim', ''), ' --extract ', outdir, '/snpsInterest.txt --maf ', maf, ' --freq --out ', outdir, '/frequencies > /dev/null 2>&1'))
                            # Read frequencies
                            freq_file = data.table::fread(paste0(outdir, '/frequencies.afreq'), h=T, stringsAsFactors=F)
                            # Check if case-control frequency should be done as well
                            if (freq_mode != FALSE){
                                freq_file = CaseControlFreq(freq_file, freq_mode, f, outdir, assoc_info)
                            }
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
                            freq_file = data.table::fread(paste0(outdir, '/frequencies.afreq'), h=T, stringsAsFactors=F)
                            # Check if case-control frequency should be done as well
                            if (freq_mode != FALSE){
                                freq_file = CaseControlFreq(freq_file, freq_mode, f, outdir, assoc_info)
                            }
                            # Add to dataframe
                            all_freq = rbind(all_freq, freq_file)
                        }
                    }
                } else if (genotype_type == 'plink2'){
                    if (maf == TRUE){
                        # Extract dosages
                        system(paste0('plink2 --pfile ', str_replace_all(f, '.pvar', ''), ' --extract ', outdir, '/snpsInterest.txt --maf ', maf, ' --export A include-alt --out ', outdir, '/dosages > /dev/null 2>&1'))
                        # Calculate frequencies if requested
                        if (freq == TRUE){
                            system(paste0('plink2 --pfile ', str_replace_all(f, '.pvar', ''), ' --extract ', outdir, '/snpsInterest.txt --maf ', maf, ' --freq --out ', outdir, '/frequencies > /dev/null 2>&1'))
                            # Read frequencies
                            freq_file = data.table::fread(paste0(outdir, '/frequencies.afreq'), h=T, stringsAsFactors=F)
                            # Check if case-control frequency should be done as well
                            if (freq_mode != FALSE){
                                freq_file = CaseControlFreq(freq_file, freq_mode, f, outdir, assoc_info)
                            }
                            # Add to dataframe
                            all_freq = rbind(all_freq, freq_file)
                        }
                    } else {
                        # Extract dosages
                        system(paste0('plink2 --pfile ', str_replace_all(f, '.pvar', ''), ' --extract ', outdir, '/snpsInterest.txt --export A include-alt --out ', outdir, '/dosages > /dev/null 2>&1'))
                        # Calculate frequencies if requested
                        if (freq == TRUE){
                            system(paste0('plink2 --pfile ', str_replace_all(f, '.pvar', ''), ' --extract ', outdir, '/snpsInterest.txt --freq --out ', outdir, '/frequencies > /dev/null 2>&1'))
                            # Read frequencies
                            freq_file = data.table::fread(paste0(outdir, '/frequencies.afreq'), h=T, stringsAsFactors=F)
                            # Check if case-control frequency should be done as well
                            if (freq_mode != FALSE){
                                freq_file = CaseControlFreq(freq_file, freq_mode, f, outdir, assoc_info)
                            }
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
                system(paste0('rm ', outdir, '/frequen*.afreq'), ignore.stderr=T, ignore.stdout=T)
                system(paste0('rm ', outdir, '/frequen*.log'), ignore.stderr=T, ignore.stdout=T)
            },
            error = function(e){
                cat(paste0('****** Skipping file ', basename(f), ' as no variants were found.\n'))
            })
        }
        # Write frequencies
        if (nrow(all_freq) > 0){
            write.table(all_freq, paste0(outdir, '/frequencies.txt'), quote=F, row.names=F, sep="\t")
        }
        res = list(all_dos, matchingsnps_all)
        return(res)
    }

    # Function to calculate case-control frequencies
    CaseControlFreq = function(freq_file, freq_mode, f, outdir, assoc_info){
        # split the variables of interest
        var_interest = unlist(str_split(freq_mode, ','))
        # iterate over variables
        for (v in var_interest){
            # extract phenotypes of interest
            ph_sub = assoc_info[[2]][, c(assoc_info[[1]], v)]
            # extract cases and controls based on the mapping stats
            pairs <- strsplit(assoc_info[[3]]$mapping, ";\\s*")[[1]]
            kv <- strsplit(pairs, " -> ")
            # Extract group labels based on RHS values
            rhs <- sapply(kv, function(x) as.integer(x[2]))
            lhs <- sapply(kv, function(x) as.integer(x[1]))
            controls <- lhs[rhs == 0]
            cases <- lhs[rhs == 1]
            # get the ids
            controls_ids = ph_sub[which(ph_sub[, v] == controls), assoc_info[[1]]]
            cases_ids = ph_sub[which(ph_sub[, v] == cases), assoc_info[[1]]]
            # write the ids
            write.table(controls_ids, paste0(outdir, '/tmp_controls.txt'), quote=F, row.names=F, col.names=F)
            write.table(cases_ids, paste0(outdir, '/tmp_cases.txt'), quote=F, row.names=F, col.names=F)
            # then calculate frequency
            system(paste0('plink2 --pfile ', str_replace_all(f, '.pvar', ''), ' --extract ', outdir, '/snpsInterest.txt --keep ', outdir, '/tmp_controls.txt --freq --out ', outdir, '/frequencies_controls > /dev/null 2>&1'))
            system(paste0('plink2 --pfile ', str_replace_all(f, '.pvar', ''), ' --extract ', outdir, '/snpsInterest.txt --keep ', outdir, '/tmp_cases.txt --freq --out ', outdir, '/frequencies_cases > /dev/null 2>&1'))
            # read frequencies
            freq_file_controls = data.table::fread(paste0(outdir, '/frequencies_controls.afreq'), h=T, stringsAsFactors=F)
            freq_file_cases = data.table::fread(paste0(outdir, '/frequencies_cases.afreq'), h=T, stringsAsFactors=F)
            # subset of columns
            freq_file_controls = freq_file_controls[, c('ID', 'ALT_FREQS', 'OBS_CT')]
            freq_file_cases = freq_file_cases[, c('ID', 'ALT_FREQS', 'OBS_CT')]
            # rename columns
            colnames(freq_file_controls) = c('ID', paste0('ALT_FREQS_', v, '_controls'), paste0('ALT_FREQS_', v, '_controls_obs'))
            colnames(freq_file_cases) = c('ID', paste0('ALT_FREQS_', v, '_cases'), paste0('ALT_FREQS_', v, '_cases_obs'))
            # combine files
            freq_cc_merged = merge(freq_file_controls, freq_file_cases, by = 'ID', all = T)
            # combine with all frequencies
            freq_file = merge(freq_file, freq_cc_merged, by = 'ID', all = T)
        }
        return(freq_file)
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
    makePlot = function(prs_df, included_snps, suffix, outdir, snps_data){
        # define title of the plots
        density_title = ifelse(suffix == "", paste0('PRS based on ', nrow(included_snps), ' SNPs'), paste0('PRS based on ', nrow(included_snps), ' SNPs (APOE excluded)'))
        snps_title = ifelse(suffix == "", 'Beta of SNPs', 'Beta of SNPs (APOE excluded)')
        # plot density
        plt_density = ggplot(data = prs_df, aes(x = PRS, fill = 'red')) + geom_density() + theme(legend.position = 'none') + xlab('Polygenic Risk Score') + ylab('Density') + ggtitle(density_title) + theme_bw() + theme(axis.title = element_text(size = 16), axis.text = element_text(size = 14), plot.title = element_text(size = 18), legend.position = 'none')
        pdf(paste0(outdir, '/PRS_density', suffix, '.pdf'), height = 7, width = 7)
        print(plt_density)
        invisible(dev.off())
        # plot snps
        # add the missing snps to the included snps
        # construct ID fields
        included_snps$ID <- paste(included_snps$CHROM, included_snps$POS, sep = ":")
        snps_data$ID <- paste(stringr::str_replace_all(snps_data$CHROM, 'chr', ''), snps_data$POS, sep = ":")
        snps_data$SNP <- paste(snps_data$CHROM, snps_data$POS, sep = ":")
        # identify missing SNPs
        missing <- snps_data[!(snps_data$ID %in% included_snps$ID), ]
        if (nrow(missing) > 0) {
            excluded <- data.frame(SNP = missing$SNP, BETA = missing$BETA, ALLELE = missing$EFFECT_ALLELE, OTHER_ALLELE = missing$OTHER_ALLELE, TYPE = 'Excluded', POS = missing$POS, CHROM = missing$CHROM, ID = missing$ID)
            included_snps <- rbind(included_snps, excluded)
        }
        included_snps$TYPE = factor(included_snps$TYPE, levels = c('Included', 'Excluded'))
        # order
        included_snps = included_snps[order(included_snps$TYPE, included_snps$BETA),]
        # set order
        included_snps$SNP = factor(included_snps$SNP, levels = unique(included_snps$SNP))
        # define dynamically the height of the plot
        base_height = 5
        height_per_row = 0.125
        plot_height = max(5, min(60, base_height + nrow(included_snps) * height_per_row))
        # then the plot
        plt_snps = ggplot(data = included_snps, aes(y = SNP, x = BETA, color = TYPE)) + geom_point(stat = 'identity', size = 2) + xlab('Beta') + ylab('SNP') + ggtitle(snps_title) + theme_bw() + theme(axis.title = element_text(size = 16), axis.text = element_text(size = 14), plot.title = element_text(size = 18), legend.position = 'top', legend.text = element_text(size = 14), legend.title = element_text(size = 14)) + scale_color_manual(values = c('Included' = 'navy', 'Excluded' = 'red')) + labs(color = "Type") + scale_x_continuous(expand = expansion(mult = c(0.05, 0.20)))
        pdf(paste0(outdir, '/PRS_SNPs', suffix, '.pdf'), height = plot_height, width = 7)
        print(plt_snps)
        invisible(dev.off())
    }

    # Function to write log file
    writeLog = function(outdir, genotype_file, snps_file, outfile, isdosage, plt, maf, multiple, excludeAPOE, fliprisk, keepDos, addWeight, freq, assoc, assoc_var, assoc_cov){
        # Define output name
        outname = paste0(outdir, '/run_info.log')
        # Create log info
        info = paste0("Genotype file: ", genotype_file, "\nMultiple files: ", multiple, "\nSNPs file: ", snps_file, "\nOutput file: ", outfile, "\nDosage: ", isdosage, "\nMAF: ", maf, "\nWith and Without APOE: ", excludeAPOE, "\nDirect effects (Risk and Protective): ", fliprisk, "\nKeep dosages: ", keepDos, "\nAdditional weight: ", addWeight, "\nCalculate frequency: ", freq, "\nPlot: ", plt, '\n\n', 'Association file: ', assoc, '\nAssociation variables: ', assoc_var, '\nAssociation covariates: ', assoc_cov, '\n\n')
        # Write log file
        writeLines(info, outname)
        return(info)
    }

    # Function to write output
    write_prs_outputs = function(prs_df, included_snps, snps_data, suffix, outdir) {
        # Write PRS table
        write.table(prs_df, paste0(outdir, '/PRS_table', suffix, '.txt'), quote = FALSE, row.names = FALSE, sep = "\t")
        # Construct ID fields
        included_snps$ID <- paste(included_snps$CHROM, included_snps$POS, sep = ":")
        snps_data$ID <- paste(stringr::str_replace_all(snps_data$CHROM, 'chr', ''), snps_data$POS, sep = ":")
        # Identify missing SNPs
        missing <- snps_data[!(snps_data$ID %in% included_snps$ID), ]
        if (nrow(missing) > 0) {
            excluded <- data.frame(SNP = NA, BETA = missing$BETA, ALLELE = missing$EFFECT_ALLELE, OTHER_ALLELE = missing$OTHER_ALLELE, TYPE = 'Excluded', POS = missing$POS, CHROM = missing$CHROM, ID = missing$ID)
            included_snps <- rbind(included_snps, excluded)
        }
        # Write SNPs included (and excluded) file
        write.table(included_snps, paste0(outdir, '/SNPs_included_PRS', suffix, '.txt'), quote = FALSE, row.names = FALSE, sep = "\t")
    }

    # Function to check association variables
    checkAssocVar = function(assoc_var, assoc_data){
        # check if variables were provided
        if (assoc_var[1] == FALSE){
            stop('** No association variables provided.\n\n', call. = FALSE)
        }
        # unlist the variables
        assoc_var = unlist(strsplit(assoc_var, ','))
        # check if variables are in the data
        assoc_var = assoc_var[which(toupper(assoc_var) %in% toupper(colnames(assoc_data)))]
        if (length(assoc_var) == 0){
            stop('** No association variables found in the data. Did you mispelled the variable names?.\n\n', call. = FALSE)
        }
        # Define a dataframe with the association variables
        df_variables = data.frame()
        # Check if the variables are numeric or categorical
        for (i in 1:length(assoc_var)){
            if (length(table(assoc_data[, assoc_var[i]])) == 2){
                # check if the variable is numeric
                if (is.numeric(assoc_data[, assoc_var[i]])){
                    cat('**** Variable ', assoc_var[i], ': numeric with 2 values. Logistic regression will be used.\n')
                    # make sure the smaller is 0 and the larger is 1
                    min_value = min(assoc_data[, assoc_var[i]], na.rm=T)
                    # store the mapping before updating the data
                    var_mapping = paste0(min_value, ' -> 0; ', max(assoc_data[, assoc_var[i]], na.rm=T), ' -> 1')
                    # update values
                    assoc_data[, assoc_var[i]] = ifelse(assoc_data[, assoc_var[i]] == min_value, 0, 1)
                    # update dataframe
                    df_variables = rbind(df_variables, data.frame(variable = assoc_var[i], type = 'numeric', model = 'binomial', mapping = var_mapping))
                } else {
                    cat('**** Variable ', assoc_var[i], ': categorical with 2 values. Logistic regression will be used.\n')
                    # make sure the smaller is 0 and the larger is 1
                    min_value = min(assoc_data[, assoc_var[i]], na.rm=T)
                    # store the mapping before updating the data
                    var_mapping = paste0(min_value, ' -> 0; ', max(assoc_data[, assoc_var[i]], na.rm=T), ' -> 1')
                    # update values
                    assoc_data[, assoc_var[i]] = ifelse(assoc_data[, assoc_var[i]] == min_value, 0, 1)
                    # update dataframe
                    df_variables = rbind(df_variables, data.frame(variable = assoc_var[i], type = 'categorical', model = 'binomial', mapping = var_mapping))
                }
            } else {
                # check if the variable is numeric
                if (is.numeric(assoc_data[, assoc_var[i]])){
                    cat('**** Variable ', assoc_var[i], ': numeric with more than 2 values. Linear regression will be used.\n')
                    # store the mapping before updating the data
                    var_mapping = paste0('No mapping needed')
                    # update dataframe
                    df_variables = rbind(df_variables, data.frame(variable = assoc_var[i], type = 'numeric', model = 'gaussian', mapping = var_mapping))
                } else {
                    cat('**** Variable ', assoc_var[i], ': categorical with more than 2 values. Linear regression will be used.\n')
                    # store the mapping before updating the data
                    var_mapping = paste0('No mapping needed')
                    # update dataframe
                    df_variables = rbind(df_variables, data.frame(variable = assoc_var[i], type = 'categorical', model = 'gaussian', mapping = var_mapping))
                }
            }
        }
        return(list(df_variables, assoc_data))
    }

    # Function to check covariates
    checkAssocCov = function(assoc_cov, assoc_data){
        # check if covariates were provided
        if (assoc_cov[1] == FALSE){
            cat('**** No association covariates provided.\n\n')
            return(NA)
        } else {
            # unlist the variables
            assoc_cov = unlist(strsplit(assoc_cov, ','))
            # check if variables are in the data
            assoc_cov = assoc_cov[which(toupper(assoc_cov) %in% toupper(colnames(assoc_data)))]
            if (length(assoc_cov) == 0){
                stop('** No association covariates found in the data. Did you mispelled the variable names?.\n\n', call. = FALSE)
            }
            cat('**** Covariates found: ', paste(assoc_cov, collapse = ', '), '\n')
            return(assoc_cov)
        }
    }

    # Function to check association mode
    checkAssocMode = function(assoc_mode){
        # check if the mode is valid
        if (toupper(assoc_mode) %in% c('PRS', 'SINGLE', 'BOTH')){
            assoc_mode = ifelse(toupper(assoc_mode) == 'PRS', 'prs', ifelse(toupper(assoc_mode) == 'SINGLE', 'single', 'both'))
            cat('** Association mode: ', toupper(assoc_mode), '\n')
            return(assoc_mode)
        } else {
            stop('** Association mode not valid. Please use prs, single or both.\n\n', call. = FALSE)
        }
    }

    # Function to check association file
    checkAssocFile = function(assoc, assoc_var, assoc_cov){
        # check if file exists
        if (file.exists(assoc)){
            cat('** Association data file found.\n')
            # open file
            assoc_data = fread(assoc, h=T, stringsAsFactors=F, sep="\t")
            # check if the file is empty
            if (nrow(assoc_data) == 0){
                stop('** Association data file is empty.\n\n', call. = FALSE)
            }
            # convert to dataframe
            assoc_data = data.frame(assoc_data, check.names=F)
            # convert nas to NA
            assoc_data[assoc_data == ''] = NA
            assoc_data[toupper(assoc_data) == 'NA'] = NA
            assoc_data[toupper(assoc_data) == 'N/A'] = NA
            assoc_data[toupper(assoc_data) == 'NAN'] = NA
            # check if IID is in the columns
            if ('IID' %in% toupper(colnames(assoc_data))){
                # get the name of the IID column
                idname = colnames(assoc_data)[which(toupper(colnames(assoc_data)) == 'IID')]
                # check variable to associate
                outcome_info = checkAssocVar(assoc_var, assoc_data)
                df_variables = outcome_info[[1]]
                assoc_data = outcome_info[[2]]
                # check covariates
                covar_names = checkAssocCov(assoc_cov, assoc_data)
            } else {
                stop('** No IID column found in association data.\n\n', call. = FALSE)
            }
        } else {
            stop('** Association data file does not exist.\n\n', call. = FALSE)
        }
        return(list(idname, assoc_data, df_variables, covar_names))
    }

    # Function to run association
    assoc_test = function(prs_df, assoc_info, outdir, suffix, assoc_mode, dosages){
        # extract association info
        assoc_idname = assoc_info[[1]]
        assoc_data = assoc_info[[2]]
        assoc_variables = assoc_info[[3]]
        assoc_covariates = assoc_info[[4]]
        # merge the prs with association info
        prs_df_pheno = merge(prs_df, assoc_data, by.x = 'iid', by.y = assoc_idname)
        # check if single-variant association should be done
        if (assoc_mode %in% c('single', 'both')){
            # get snp names first
            snp_names = colnames(dosages)[-1]
            # if so, merge single-variant data with association info
            dosages_pheno = merge(dosages, assoc_data, by.x = 'IID', by.y = assoc_idname)
            # check if the merge was successful for single-variant -- stop only for single mode
            if (nrow(dosages_pheno) == 0 & assoc_mode == 'single'){
                stop('** No matching individuals between single-variant and association data.\n\n', call. = FALSE)
            }
        }
        # check if the merge was successful for PRS
        if (nrow(prs_df_pheno) == 0 & assoc_mode %in% c('prs', 'both')){
            stop('** No matching individuals between PRS and association data.\n\n', call. = FALSE)
        }
        # check which associations need to be done and do them
        if (assoc_mode == 'prs'){
            cat('**** Running association for PRS.\n')
            assoc_results_prs = prs_assoc(prs_df_pheno, assoc_variables, assoc_covariates, suffix)
            # write the results
            write.table(assoc_results_prs, paste0(outdir, '/association_results_PRS', suffix, '.txt'), quote = F, row.names = F, sep = "\t")
        } else if (assoc_mode == 'single'){
            cat('**** Running association for single-variant.\n')
            assoc_results_single = singleVar_assoc(dosages_pheno, assoc_variables, assoc_covariates, suffix, snp_names)
            # write the results
            write.table(assoc_results_single, paste0(outdir, '/association_results_single', suffix, '.txt'), quote = F, row.names = F, sep = "\t")
        } else {
            cat('**** Running association for both PRS and single-variant.\n')
            assoc_results_prs = prs_assoc(prs_df_pheno, assoc_variables, assoc_covariates, suffix)
            assoc_results_single = singleVar_assoc(dosages_pheno, assoc_variables, assoc_covariates, suffix, snp_names)
            # write the results
            write.table(assoc_results_prs, paste0(outdir, '/association_results_PRS', suffix, '.txt'), quote = F, row.names = F, sep = "\t")
            write.table(assoc_results_single, paste0(outdir, '/association_results_single', suffix, '.txt'), quote = F, row.names = F, sep = "\t")
        }
    }

    # Function to run association for the PRS
    prs_assoc = function(prs_df_pheno, assoc_variables, assoc_covariates, suffix){
        # Scale PRS before association
        prs_df_pheno$PRS = scale(prs_df_pheno$PRS)
        # Define container for the results
        assoc_results = data.frame()
        # iterate over the variables
        for (i in 1:nrow(assoc_variables)){
            cat(paste0('****** Running PRS association for variable: ', assoc_variables$variable[i], ' ', str_replace_all(suffix, '_', ''), '\n'))
            # get the variable name
            var_name = assoc_variables$variable[i]
            # get the model type
            model_type = assoc_variables$model[i]
            # get the covariates
            covar_names = assoc_covariates
            # get the mapping
            var_mapping = assoc_variables$mapping[i]
            # check if there are covariates
            if (length(na.omit(covar_names)) > 0){
                # create formula
                formula = as.formula(paste0(var_name, ' ~ PRS + ', paste(covar_names, collapse = ' + ')))
            } else {
                formula = as.formula(paste0(var_name, ' ~ PRS'))
            }
            # calculate number of cases and controls for logistic regression
            if (model_type == 'binomial'){
                n_cases = nrow(prs_df_pheno[which(prs_df_pheno[, var_name] == 1),])
                n_controls = nrow(prs_df_pheno[which(prs_df_pheno[, var_name] == 0),])
            } else {
                n_cases = NA
                n_controls = NA
            }
            # perform association test
            model = suppressWarnings(glm(formula, data = prs_df_pheno, family = model_type))
            # add to the results
            assoc_results = rbind(assoc_results, data.frame(Predictor = 'PRS', Outcome = var_name, Covariates = paste(covar_names, collapse = ','), Beta_PRS = summary(model)$coefficients[2, 1], SE_PRS = summary(model)$coefficients[2, 2], P_PRS = summary(model)$coefficients[2, 4], Model = model_type, N_tot = nrow(prs_df_pheno), N_missing = sum(is.na(prs_df_pheno[, var_name])), Mapping = var_mapping, Model_converged = model$converged, N_effective = nrow(prs_df_pheno) - sum(is.na(prs_df_pheno[, var_name])), N_cases = n_cases, N_controls = n_controls, stringsAsFactors = F))
        }
        return(assoc_results)
    }

    # Function to run association for the single-variants
    singleVar_assoc = function(dosages_pheno, assoc_variables, assoc_covariates, suffix, snp_names){
        # make sure data is a dataframe
        dosages_pheno = data.frame(dosages_pheno, check.names=F)
        # Define container for the results
        assoc_results = data.frame()
        # iterate over the variables
        for (i in 1:nrow(assoc_variables)){
            cat(paste0('****** Running single-variant association for variable: ', assoc_variables$variable[i], ' ', str_replace_all(suffix, '_', ''), '\n'))
            # get the variable name
            var_name = assoc_variables$variable[i]
            # get the model type
            model_type = assoc_variables$model[i]
            # get the covariates
            covar_names = assoc_covariates
            # get the mapping
            var_mapping = assoc_variables$mapping[i]
            # then iterate over the snps
            for (snp in snp_names){
                tryCatch({
                    # update the snp name as otherwise there are bad characters in the name (e.g. /)
                    colnames(dosages_pheno)[which(colnames(dosages_pheno) == snp)] = 'SNP'
                    # check if there are covariates
                    if (length(na.omit(covar_names)) > 0){
                        # create formula
                        formula = as.formula(paste0(var_name, ' ~ SNP + ', paste(covar_names, collapse = ' + ')))
                    } else {
                        formula = as.formula(paste0(var_name, ' ~ SNP'))
                    }
                    # calculate number of cases and controls for logistic regression
                    if (model_type == 'binomial'){
                        n_cases = nrow(dosages_pheno[which(dosages_pheno[, var_name] == 1),])
                        n_controls = nrow(dosages_pheno[which(dosages_pheno[, var_name] == 0),])
                    } else {
                        n_cases = NA
                        n_controls = NA
                    }
                    # perform association test
                    model = suppressWarnings(glm(formula, data = dosages_pheno, family = model_type))
                    # add to the results
                    assoc_results = rbind(assoc_results, data.frame(Predictor = snp, Outcome = var_name, Covariates = paste(covar_names, collapse = ','), Beta_PRS = summary(model)$coefficients[2, 1], SE_PRS = summary(model)$coefficients[2, 2], P_PRS = summary(model)$coefficients[2, 4], Model = model_type, N_tot = nrow(dosages_pheno), N_missing = sum(is.na(dosages_pheno[, var_name])), Mapping = var_mapping, Model_converged = model$converged, N_effective = nrow(dosages_pheno) - sum(is.na(dosages_pheno[, var_name])), N_cases = n_cases, N_controls = n_controls, stringsAsFactors = F))
                    # restore the name
                    colnames(dosages_pheno)[which(colnames(dosages_pheno) == 'SNP')] = snp
                },
                error = function(e){
                    # add to the results
                    assoc_results = rbind(assoc_results, data.frame(Predictor = snp, Outcome = var_name, Covariates = paste(covar_names, collapse = ','), Beta_PRS = NA, SE_PRS = NA, P_PRS = NA, Model = model_type, N_tot = nrow(dosages_pheno), N_missing = sum(is.na(dosages_pheno[, var_name])), Mapping = var_mapping, Model_converged = NA, N_effective = nrow(dosages_pheno) - sum(is.na(dosages_pheno[, var_name])), N_cases = NA, N_controls = NA, stringsAsFactors = F))
                })
            }
        }
        return(assoc_results)
    }