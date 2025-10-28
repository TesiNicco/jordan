# Libraries
    suppressWarnings({
        suppressMessages({
            library(argparse)
            library(data.table)
            library(stringr)
            library(survival)
            library(ggplot2)
            args <- commandArgs(trailingOnly = FALSE)
            library(survminer)
            library(ggpubr)
        })
    })

# Functions
    # Function to detect the system and adapt plink executables
    detect_system <- function(script_path) {
        os <- Sys.info()[["sysname"]]
        arch <- .Machine$sizeof.pointer
        arch_str <- if (arch == 8) "64-bit" else "32-bit"
        
        if (os == "Linux") {
            estimate = 'Linux (Unknown CPU)'
            info <- readLines("/proc/cpuinfo")
            flags <- grep("flags", info, value = TRUE)[1]
            vendor <- grep("vendor_id", info, value = TRUE)[1]
            is_avx2 <- grepl("avx2", flags)
            is_amd <- grepl("AuthenticAMD", vendor)
            is_intel <- grepl("GenuineIntel", vendor)
            
            if (is_avx2 && is_intel) estimate = "Linux AVX2 Intel"
            if (is_avx2 && is_amd) estimate = "Linux AVX2 AMD"
            if (arch == 8 && is_intel) estimate = "Linux 64-bit Intel"
            if (arch == 4) estimate = "Linux 32-bit"
        } else if (os == "Darwin") {
            estimate = 'macOS AVX2'
            chip <- system("sysctl -n machdep.cpu.brand_string", intern = TRUE)
            is_m1 <- grepl("Apple", chip)
            
            if (is_m1) estimate = "macOS M1"
        } else {
            # For other systems, we can only return a generic message and stop as it's not supported
            stop("Unsupported operating system. Please use Linux or macOS.")
        }
        # The define the plink executables based on the system
        if (estimate == 'Linux AVX2 Intel'){
            plink_path = file.path(script_path, 'plink_executables', 'plink_linux_64bit')
            plink2_path = file.path(script_path, 'plink_executables', 'plink2_linux_avx2_intel')
        } else if (estimate == 'Linux AVX2 AMD'){
            plink_path = file.path(script_path, 'plink_executables', 'plink_linux_64bit')
            plink2_path = file.path(script_path, 'plink_executables', 'plink2_linux_avx2_amd')
        } else if (estimate == 'Linux 64-bit Intel' | estimate == 'Linux (Unknown CPU)'){
            plink_path = file.path(script_path, 'plink_executables', 'plink_linux_64bit')
            plink2_path = file.path(script_path, 'plink_executables', 'plink2_linux_intel')
        } else if (estimate == 'Linux 32-bit'){
            plink_path = file.path(script_path, 'plink_executables', 'plink_linux_32bit')
            plink2_path = file.path(script_path, 'plink_executables', 'plink2_linux_32bit')
        } else if (estimate == 'macOS M1'){
            plink_path = file.path(script_path, 'plink_executables', 'plink_macos')
            plink2_path = file.path(script_path, 'plink_executables', 'plink2_macos_m1')
        } else if (estimate == 'macOS AVX2'){
            plink_path = file.path(script_path, 'plink_executables', 'plink_macos')
            plink2_path = file.path(script_path, 'plink_executables', 'plink2_macos_avx2')
        }
        return(list(plink_path, plink2_path, estimate))
    }

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
    checkGenoFile = function(genotype_file, outdir, isdosage, multiple, script_path, plink_path, plink2_path){
        # check if file is a vcf/bcf
        if (file.exists(genotype_file) && (endsWith(genotype_file, 'vcf') | endsWith(genotype_file, 'vcf.gz') | endsWith(genotype_file, 'bcf') | endsWith(genotype_file, 'bcf.gz'))){
            # check if the file is a vcf/bcf
            vcftype = ifelse(endsWith(genotype_file, 'vcf') | endsWith(genotype_file, 'vcf.gz'), 'vcf', 'bcf')
            # vcf file
            if (isdosage == TRUE){
                if (multiple == TRUE){
                    cat('** Multiple VCF files found. Converting to PLINK assuming it is imputed data from Minimac-4 (dosage=HDS).\n')
                    all_files = c(system(paste0("find ", dirname(genotype_file), " -maxdepth 1 -type f -regex '.*\\.", vcftype, "$'"), intern = TRUE), system(paste0("find ", dirname(genotype_file), " -maxdepth 1 -type f -regex '.*\\.", vcftype, ".gz$'"), intern = TRUE))
                    # create temporary ids
                    tmp_ids = paste0(outdir, '/tmp_', seq_along(all_files))
                    # convert all files to plink
                    for (i in seq_along(all_files)){
                        system(paste0(plink2_path, ' --', vcftype, ' ', all_files[i], ' dosage=HDS --make-pgen --out ', tmp_ids[i], ' > /dev/null 2>&1'))
                    }
                    data_path = system(paste0("find ", outdir, " -maxdepth 1 -type f -name 'tmp_*.pvar'"), intern = TRUE)
                } else {
                    cat('** VCF file found. Converting to PLINK assuming it is imputed data from Minimac-4 (dosage=HDS).\n')
                    system(paste0(plink2_path, ' --', vcftype, ' ', genotype_file, ' dosage=HDS --make-pgen --out ', outdir, '/tmp > /dev/null 2>&1'))
                    data_path = paste0(outdir, '/tmp.pvar')
                }
            } else {
                if (multiple == TRUE){
                    cat('** Multiple VCF files found. Converting to PLINK.\n')
                    all_files = c(system(paste0("find ", dirname(genotype_file), " -maxdepth 1 -type f -regex '.*\\.", vcftype, "$'"), intern = TRUE), system(paste0("find ", dirname(genotype_file), " -maxdepth 1 -type f -regex '.*\\.", vcftype, ".gz$'"), intern = TRUE))
                    # create temporary ids
                    tmp_ids = paste0(outdir, '/tmp_', seq_along(all_files))
                    # convert all files to plink
                    for (i in seq_along(all_files)){
                        system(paste0(plink2_path, ' --', vcftype, ' ', all_files[i], ' --make-pgen --out ', tmp_ids[i], ' > /dev/null 2>&1'))
                    }
                    data_path = system(paste0("find ", outdir, " -maxdepth 1 -type f -name 'tmp_*.pvar'"), intern = TRUE)
                } else {
                    cat('** VCF file found. Converting to PLINK.\n')
                    system(paste0(plink2_path, ' --', vcftype, ' ', genotype_file, ' --make-pgen --out ', outdir, '/tmp > /dev/null 2>&1'))
                    data_path = paste0(outdir, '/tmp.pvar')
                }
            }
            res = list(data_path, 'plink2')
        } else if (file.exists(paste0(genotype_file, '.pvar'))){
            genotype_file = paste0(genotype_file, '.pvar')
            cat('** PLINK2 file found.\n')
            # check other chromosomes if present
            if (multiple == TRUE){
                all_files <- system(paste0("find ", dirname(genotype_file), " -maxdepth 1 -type f -regex '.*\\.pvar$'"), intern = TRUE)
                genotype_file = all_files
            }
            res = list(genotype_file, 'plink2')
        } else if (file.exists(paste0(genotype_file, '.bim'))){
            genotype_file = paste0(genotype_file, '.bim')
            cat('** PLINK file found.\n')
            if (multiple == TRUE){
                all_files <- system(paste0("find ", dirname(genotype_file), " -maxdepth 1 -type f -regex '.*\\.bim$'"), intern = TRUE)
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
                    cc_freq_var = paste(unique(assoc_vars$variable[which(assoc_vars$mode == 'binomial')]), collapse = ',')
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
    makePRS = function(outdir, genotype_path, snps_data, genotype_type, multiple, excludeAPOE, maf, fliprisk, keepDos, addWeight, freq, assoc_file, assoc_info, script_path, sex_strata, plink_path, plink2_path){
        # decide whether frequencies in cases and controls need to be calculated
        freq_mode = defineFreqMode(assoc_file, assoc_info, freq)
        # match ids in the plink file and extract dosages/genotypes
        cat('**** Matching SNPs and extracting dosages.\n')
        # Match variants of interest and get the dosages
        res = matchIDs_multiple(genotype_path, snps_data, genotype_type, outdir, maf, freq, freq_mode, script_path, plink_path, plink2_path)
        dosages = res[[1]]
        mappingSnp = res[[2]]
        all_freq = res[[3]]
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
            return(list(res, res_noapoe, dosages, all_freq))
        } else {
            return(list(res, dosages, all_freq))
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

    matchIDs_multiple = function(genotype_path, snps_data, genotype_type, outdir, maf, freq, freq_mode, script_path, plink_path, plink2_path){
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
                    } else{
                        matchingsnps_lv1$allele_match[i] = 'no'
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
                        system(paste0(plink2_path, ' --bfile ', str_replace_all(f, '.bim', ''), ' --extract ', outdir, '/snpsInterest.txt --maf ', maf, ' --export A include-alt --out ', outdir, '/dosages > /dev/null 2>&1'))
                        # Calculate frequencies if requested
                        if (freq == TRUE){
                            system(paste0(plink2_path, ' --bfile ', str_replace_all(f, '.bim', ''), ' --extract ', outdir, '/snpsInterest.txt --maf ', maf, ' --freq cols=+pos --out ', outdir, '/frequencies > /dev/null 2>&1'))
                            # Read frequencies
                            freq_file = data.table::fread(paste0(outdir, '/frequencies.afreq'), h=T, stringsAsFactors=F)
                            # Add unique identifier
                            freq_file$UNIQUE_ID = paste(freq_file$ID, freq_file$REF, freq_file$ALT, sep = "_")
                            # Check if case-control frequency should be done as well
                            if (freq_mode != FALSE){
                                freq_file = CaseControlFreq(freq_file, freq_mode, f, outdir, assoc_info, genotype_type, script_path, plink2_path)
                            }
                            # Add to dataframe
                            all_freq = rbind(all_freq, freq_file)
                        }
                    } else {
                        # Extract dosages
                        system(paste0(plink_path, ' --bfile ', str_replace_all(f, '.bim', ''), ' --extract ', outdir, '/snpsInterest.txt --export A include-alt --out ', outdir, '/dosages > /dev/null 2>&1'))
                        # Calculate frequencies if requested
                        if (freq == TRUE){
                            system(paste0(plink2_path, ' --bfile ', str_replace_all(f, '.bim', ''), ' --extract ', outdir, '/snpsInterest.txt --freq cols=+pos --out ', outdir, '/frequencies > /dev/null 2>&1'))
                            # Read frequencies
                            freq_file = data.table::fread(paste0(outdir, '/frequencies.afreq'), h=T, stringsAsFactors=F)
                            # Add unique identifier
                            freq_file$UNIQUE_ID = paste(freq_file$ID, freq_file$REF, freq_file$ALT, sep = "_")
                            # Check if case-control frequency should be done as well
                            if (freq_mode != FALSE){
                                freq_file = CaseControlFreq(freq_file, freq_mode, f, outdir, assoc_info, genotype_type, script_path, plink2_path)
                            }
                            # Add to dataframe
                            all_freq = rbind(all_freq, freq_file)
                        }
                    }
                } else if (genotype_type == 'plink2'){
                    if (maf == TRUE){
                        # Extract dosages
                        system(paste0(plink2_path, ' --pfile ', str_replace_all(f, '.pvar', ''), ' --extract ', outdir, '/snpsInterest.txt --maf ', maf, ' --export A include-alt --out ', outdir, '/dosages > /dev/null 2>&1'))
                        # Calculate frequencies if requested
                        if (freq == TRUE){
                            system(paste0(plink2_path, ' --pfile ', str_replace_all(f, '.pvar', ''), ' --extract ', outdir, '/snpsInterest.txt --maf ', maf, ' --freq cols=+pos --out ', outdir, '/frequencies > /dev/null 2>&1'))
                            # Read frequencies
                            freq_file = data.table::fread(paste0(outdir, '/frequencies.afreq'), h=T, stringsAsFactors=F)
                            # Add unique identifier
                            freq_file$UNIQUE_ID = paste(freq_file$ID, freq_file$REF, freq_file$ALT, sep = "_")
                            # Check if case-control frequency should be done as well
                            if (freq_mode != FALSE){
                                freq_file = CaseControlFreq(freq_file, freq_mode, f, outdir, assoc_info, genotype_type, script_path, plink2_path)
                            }
                            # Add to dataframe
                            all_freq = rbind(all_freq, freq_file)
                        }
                    } else {
                        # Extract dosages
                        system(paste0(plink2_path, ' --pfile ', str_replace_all(f, '.pvar', ''), ' --extract ', outdir, '/snpsInterest.txt --export A include-alt --out ', outdir, '/dosages > /dev/null 2>&1'))
                        # Calculate frequencies if requested
                        if (freq == TRUE){
                            system(paste0(plink2_path, ' --pfile ', str_replace_all(f, '.pvar', ''), ' --extract ', outdir, '/snpsInterest.txt --freq cols=+pos --out ', outdir, '/frequencies > /dev/null 2>&1'))
                            # Read frequencies
                            freq_file = data.table::fread(paste0(outdir, '/frequencies.afreq'), h=T, stringsAsFactors=F)
                            # Add unique identifier
                            freq_file$UNIQUE_ID = paste(freq_file$ID, freq_file$REF, freq_file$ALT, sep = "_")
                            # Check if case-control frequency should be done as well
                            if (freq_mode != FALSE){
                                freq_file = CaseControlFreq(freq_file, freq_mode, f, outdir, assoc_info, genotype_type, script_path, plink2_path)
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
        res = list(all_dos, matchingsnps_all, all_freq)
        return(res)
    }

    # Function to calculate case-control frequencies
    CaseControlFreq = function(freq_file, freq_mode, f, outdir, assoc_info, genotype_type, script_path, plink2_path){
        # split the variables of interest
        var_interest = unlist(str_split(freq_mode, ','))
        # iterate over variables
        for (v in var_interest){
            # extract phenotypes of interest
            ph_sub = assoc_info[[2]][, c(assoc_info[[1]], v)]
            # extract cases and controls based on the mapping stats
            #pairs <- strsplit(assoc_info[[3]]$mapping, ";\\s*")[[1]]
            #kv <- strsplit(pairs, " -> ")
            # Extract group labels based on RHS values
            #rhs <- sapply(kv, function(x) as.integer(x[2]))
            #lhs <- sapply(kv, function(x) as.integer(x[1]))
            #controls <- lhs[rhs == 0]
            #cases <- lhs[rhs == 1]
            # get the ids
            controls_ids = ph_sub[which(ph_sub[, v] == 0), assoc_info[[1]]]
            cases_ids = ph_sub[which(ph_sub[, v] == 1), assoc_info[[1]]]
            # write the ids
            write.table(controls_ids, paste0(outdir, '/tmp_controls.txt'), quote=F, row.names=F, col.names=F)
            write.table(cases_ids, paste0(outdir, '/tmp_cases.txt'), quote=F, row.names=F, col.names=F)
            # define input file based on the type of genotype file
            input_file = ifelse(genotype_type == 'plink', paste0('--bfile ', str_replace_all(f, '.bim', '')), paste0('--pfile ', str_replace_all(f, '.pvar', '')))
            # then calculate frequency
            system(paste0(plink2_path, ' ', input_file, ' --extract ', outdir, '/snpsInterest.txt --keep ', outdir, '/tmp_controls.txt --freq cols=+pos --out ', outdir, '/frequencies_controls > /dev/null 2>&1'))
            system(paste0(plink2_path, ' ', input_file, ' --extract ', outdir, '/snpsInterest.txt --keep ', outdir, '/tmp_cases.txt --freq cols=+pos --out ', outdir, '/frequencies_cases > /dev/null 2>&1'))
            # read frequencies
            freq_file_controls = data.table::fread(paste0(outdir, '/frequencies_controls.afreq'), h=T, stringsAsFactors=F)
            freq_file_cases = data.table::fread(paste0(outdir, '/frequencies_cases.afreq'), h=T, stringsAsFactors=F)
            # create unique identifier
            freq_file_controls$UNIQUE_ID = paste(freq_file_controls$ID, freq_file_controls$REF, freq_file_controls$ALT, sep = "_")
            freq_file_cases$UNIQUE_ID = paste(freq_file_cases$ID, freq_file_cases$REF, freq_file_cases$ALT, sep = "_")
            # subset of columns
            freq_file_controls = freq_file_controls[, c('ALT_FREQS', 'OBS_CT', 'UNIQUE_ID')]
            freq_file_cases = freq_file_cases[, c('ALT_FREQS', 'OBS_CT', 'UNIQUE_ID')]
            # rename columns
            colnames(freq_file_controls) = c(paste0('ALT_FREQS_', v, '_controls'), paste0('ALT_FREQS_', v, '_controls_obs'), 'UNIQUE_ID')
            colnames(freq_file_cases) = c(paste0('ALT_FREQS_', v, '_cases'), paste0('ALT_FREQS_', v, '_cases_obs'), 'UNIQUE_ID')
            # combine files
            freq_cc_merged = merge(freq_file_controls, freq_file_cases, by = 'UNIQUE_ID', all = T)
            # combine with all frequencies
            freq_file = merge(freq_file, freq_cc_merged, by = 'UNIQUE_ID', all = T)
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
            # check if the beta of the snp is NA otherwise put 0
            temp_snpdata$risk_beta[is.na(temp_snpdata$risk_beta)] = 0
            # get additional weight if requested
            if (addWeight != FALSE){
                addWeight_num = as.numeric(temp_snpdata[, ..addWeight])
                addWeight_num[is.na(addWeight_num)] = 0
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
                    included_snps = rbind(included_snps, data.frame(SNP = temp_name, BETA = temp_snpdata$risk_beta, ALLELE = temp_snpdata$risk_allele, OTHER_ALLELE = temp_snpdata$OTHER_ALLELE, TYPE = 'Included', CHROM = temp_snpinfo$chr, POS = temp_snpinfo$pos, SNPID_DATA = colnames(dosages)[i], REASON = NA))
                } else {
                    cat(paste0('** Alleles do not match for snp ', temp_name, ': Expected alleles: ', temp_snpdata$EFFECT_ALLELE, '/', temp_snpdata$OTHER_ALLELE, ' - Found alleles: ', temp_effect, '/', temp_other, '. Skipping.\n'))
                    included_snps = rbind(included_snps, data.frame(SNP = temp_name, BETA = temp_snpdata$risk_beta, ALLELE = temp_snpdata$risk_allele, OTHER_ALLELE = temp_snpdata$OTHER_ALLELE, TYPE = 'Excluded', CHROM = temp_snpinfo$chr, POS = temp_snpinfo$pos, SNPID_DATA = colnames(dosages)[i], REASON = 'Alleles_do_not_match'))
                }
            }
        }
        res = list(prs_df, included_snps)
        return(res)
    }

    # Function to draw plot
    makePlot = function(prs_df, included_snps, suffix, outdir, snps_data, assoc_info, all_freq, freq, assoc_file, sex_strata, tiles_prs_df, split_info_df, na_be){
        suppressWarnings({
            # Plot PRS distribution independently from phenotype -- eventually sex-strata
            # define title of the plots
            density_title = ifelse(suffix == "", paste0('PRS ~ ', nrow(included_snps), ' SNPs ~ All'), paste0('PRS ~ ', nrow(included_snps), ' SNPs (APOE excluded) ~ All'))
            snps_title = ifelse(suffix == "", 'Beta of SNPs', 'Beta of SNPs (APOE excluded)')
            # plot density
            tryCatch({
                cat('**** Plotting PRS distribution.\n')
                suppressMessages({
                    if (length(assoc_info) >1 & sex_strata == TRUE){
                        # merge assoc with prs
                        prs_df_pheno = merge(prs_df[, c('iid', 'PRS')], assoc_info[[2]][, c(assoc_info[[1]], 'SEX')], by.x = 'iid', by.y = assoc_info[[1]], all = T)
                        prs_df_pheno = rbind(prs_df_pheno, data.frame(iid = prs_df_pheno$iid, PRS = prs_df_pheno$PRS, SEX = 'all'))
                        # check if nas need to be excluded
                        if (na_beh == 'exclude'){
                            tmp = prs_df_pheno[!is.na(prs_df_pheno$SEX),]
                        } else {
                            tmp = prs_df_pheno
                        }
                        plt_density = ggplot(data = tmp, aes(x = PRS, fill = SEX)) + geom_density(alpha = 0.5) + theme(legend.position = 'none') + xlab('Polygenic Risk Score') + ylab('Density') + ggtitle(density_title) + theme_bw() + theme(axis.title = element_text(size = 16), axis.text = element_text(size = 14), plot.title = element_text(size = 18), legend.position = 'top')
                    } else {
                        plt_density = ggplot(data = prs_df, aes(x = PRS, fill = 'red')) + geom_density() + theme(legend.position = 'none') + xlab('Polygenic Risk Score') + ylab('Density') + ggtitle(density_title) + theme_bw() + theme(axis.title = element_text(size = 16), axis.text = element_text(size = 14), plot.title = element_text(size = 18), legend.position = 'none')
                    }
                    pdf(paste0(outdir, '/PRS_density', suffix, '.pdf'), height = 7, width = 7)
                    print(plt_density)
                    invisible(dev.off())
                })
            }, error = function(e){
                cat('** Error in plotting PRS distribution. Skipping this plot.\n')
            })

            # if assoc_info is not FALSE and there is a binary trait, then plot the PRS vs the binary trait
            if (assoc_file != FALSE){
                # check if the trait is binary
                tryCatch({
                    if (any(grepl('binomial', assoc_info[[3]]$model))){
                        cat('**** Plotting PRS distribution vs. phenotype.\n')
                        # get the variable of interest
                        var_interest = assoc_info[[3]]$variable[which(assoc_info[[3]]$mode == 'binomial')]
                        # iterate over the variables of interest
                        for (v in var_interest){
                            # get the phenotype of interest
                            ph_sub = assoc_info[[2]][, c(assoc_info[[1]], v)]
                            # change the value in the v column from 0/1 based on the mapping
                            pairs <- strsplit(assoc_info[[3]][which(assoc_info[[3]]$variable == v), "mapping"], ";\\s*")[[1]]
                            kv <- strsplit(pairs, " -> ")
                            # Change the values
                            tmp1 = ph_sub[which(ph_sub[, v] == kv[[1]][2]), ]
                            tmp2 = ph_sub[which(ph_sub[, v] == kv[[2]][2]), ]
                            tmp1[which(tmp1[, v] == kv[[1]][2]), v] = kv[[1]][1]
                            tmp2[which(tmp2[, v] == kv[[2]][2]), v] = kv[[2]][1]
                            ph_sub = rbind(tmp1, tmp2)
                            # merge with the prs_df
                            ph_sub_prs = merge(ph_sub, prs_df, by.x = assoc_info[[1]], by.y = 'iid', all = T)
                            # check if nas need to be excluded
                            if (na_beh == 'exclude'){
                                ph_sub_prs <- ph_sub_prs[!is.na(ph_sub_prs[[v]]), ]
                            }
                            prs_title = ifelse(suffix == "", 'PRS', 'PRS (APOE excluded)')
                            # plot
                            suppressMessages({
                                plt_pheno = ggplot(data = ph_sub_prs, aes_string(x = "PRS", fill = v)) + geom_density(alpha = 0.6) + xlab(prs_title) + ylab('Density') + ggtitle(v) + theme_bw() + theme(axis.title = element_text(size = 16), axis.text = element_text(size = 14), plot.title = element_text(size = 18), legend.position = 'top', legend.text = element_text(size = 14), legend.title = element_text(size = 14))
                                pdf(paste0(outdir, '/PRS_vs_', v, suffix, '.pdf'), height = 7, width = 7)
                                print(plt_pheno)
                                invisible(dev.off())
                            })
                        }
                    }
                    if (any(grepl('gaussian', assoc_info[[3]]$model))){
                        # get the variable of interest
                        var_interest = assoc_info[[3]]$variable[which(assoc_info[[3]]$mode == 'gaussian')]
                        # iterate over the variables of interest
                        for (v in var_interest){
                            # get the phenotype of interest
                            ph_sub = assoc_info[[2]][, c(assoc_info[[1]], v)]
                            # merge with the prs_df
                            ph_sub_prs = merge(ph_sub, prs_df, by.x = assoc_info[[1]], by.y = 'iid', all = T)
                            # check if nas need to be excluded
                            if (na_beh == 'exclude'){
                                ph_sub_prs <- ph_sub_prs[!is.na(ph_sub_prs[[v]]), ]
                            }
                            # compute spearman correlation
                            cor_test <- cor.test(ph_sub_prs$PRS, ph_sub_prs[[v]], method = "spearman", na.rm = TRUE)
                            cor_text <- paste0("Spearman=", round(cor_test$estimate, 2), ", p=", signif(cor_test$p.value, 3))
                            prs_title = ifelse(suffix == "", 'PRS', 'PRS (APOE excluded)')
                            # plot
                            suppressMessages({
                                plt_pheno = ggplot(data = ph_sub_prs, aes_string(x = "PRS", y = v)) + geom_point(alpha = 0.6, size = 2, color = '#008080') + geom_smooth(method = 'lm', col = 'red') + xlab(prs_title) + ylab(paste0('Phenotype (', v, ')')) + ggtitle(paste0('PRS vs ', v)) + theme_bw() + theme(axis.title = element_text(size = 16), axis.text = element_text(size = 14), plot.title = element_text(size = 18), legend.position = 'top', legend.text = element_text(size = 14), legend.title = element_text(size = 14)) + annotate("text", x = min(ph_sub_prs$PRS, na.rm = TRUE), y = min(ph_sub_prs[[v]], na.rm = TRUE), label = cor_text, hjust = 0, vjust = 0, size = 5)
                                # add marginal density plots using ggMarginal
                                pdf(paste0(outdir, '/PRS_vs_', v, suffix, '.pdf'), height = 7, width = 7)
                                print(ggExtra::ggMarginal(plt_pheno, type = "density", fill = '#008080', size = 5))
                                invisible(dev.off())
                            })
                        }
                    }
                }, error = function(e){
                    cat('** Error in plotting PRS vs phenotype. Skipping this plot.\n')
                })
            }

            # plot snps
            # add the missing snps to the included snps
            # construct ID fields
            tryCatch({
                cat('**** SNP information.\n')
                included_snps$ID <- paste(included_snps$CHROM, included_snps$POS, sep = ":")
                snps_data$ID <- paste(stringr::str_replace_all(snps_data$CHROM, 'chr', ''), snps_data$POS, sep = ":")
                snps_data$SNP <- paste(snps_data$CHROM, snps_data$POS, sep = ":")
                # identify missing SNPs
                missing <- snps_data[!(snps_data$ID %in% included_snps$ID), ]
                if (nrow(missing) > 0){
                    excluded <- data.frame(SNP = missing$SNP, BETA = missing$BETA, ALLELE = missing$EFFECT_ALLELE, OTHER_ALLELE = missing$OTHER_ALLELE, TYPE = 'Excluded', POS = missing$POS, CHROM = missing$CHROM, ID = missing$ID)
                    included_snps <- rbind(included_snps[, c('SNP', 'BETA', 'ALLELE', 'OTHER_ALLELE', 'TYPE', 'POS', 'CHROM', 'ID')], excluded)
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
                suppressMessages({
                    plt_snps = ggplot(data = included_snps, aes(y = SNP, x = BETA, color = TYPE)) + geom_point(stat = 'identity', size = 2) + xlab('Beta') + ylab('SNP') + ggtitle(snps_title) + theme_bw() + theme(axis.title = element_text(size = 16), axis.text = element_text(size = 14), plot.title = element_text(size = 18), legend.position = 'top', legend.text = element_text(size = 14), legend.title = element_text(size = 14)) + scale_color_manual(values = c('Included' = 'navy', 'Excluded' = 'red')) + labs(color = "Type") + scale_x_continuous(expand = expansion(mult = c(0.05, 0.20)))
                    pdf(paste0(outdir, '/PRS_SNPs', suffix, '.pdf'), height = plot_height, width = 7)
                    print(plt_snps)
                    invisible(dev.off())
                })
            }, error = function(e){
                cat('** Error in plotting SNPs. Skipping this plot.\n')
            })

            # plot frequencies if these were calculated
            if (freq & nrow(all_freq) > 0){
                # check if association labels are present
                if (assoc_file != FALSE){
                    tryCatch({
                        cat('**** Plotting SNP frequencies vs. phenotype.\n')
                        # get the variable of interest
                        var_interest = assoc_info[[3]]$variable[which(assoc_info[[3]]$mode == 'binomial')]
                        # iterate over the variables of interest
                        for (v in var_interest){
                            # get data by grepping v _controls or v_cases and end line
                            v_controls = paste0(v, '_controls$')
                            v_cases = paste0(v, '_cases$')
                            # get the columns of interest
                            v_controls_cols = colnames(all_freq)[grep(v_controls, colnames(all_freq))]
                            v_cases_cols = colnames(all_freq)[grep(v_cases, colnames(all_freq))]
                            # get the variable of interest
                            all_freq = data.frame(all_freq, check.names=F)
                            ph_sub = all_freq[, c('ID', '#CHROM', 'POS', 'REF', 'ALT', v_controls_cols, v_cases_cols)]
                            # rename columns
                            colnames(ph_sub) = c('ID', 'CHROM', 'POS', 'REF', 'ALT', paste0(v, '_controls'), paste0(v, '_cases'))
                            # check if frequency columns are < 0.5 for both v_controls and v_cases, for each line
                            ph_sub$MINOR_ALLELE = ph_sub$ALT
                            ph_sub[which(ph_sub[, paste0(v, '_controls')] > 0.5 & ph_sub[, paste0(v, '_cases')] > 0.5), 'MINOR_ALLELE'] = ph_sub$REF[which(ph_sub[, paste0(v, '_controls')] > 0.5 & ph_sub[, paste0(v, '_cases')] > 0.5)]
                            ph_sub[which(ph_sub$ALT != ph_sub$MINOR_ALLELE), paste0(v, '_controls')] = 1 - ph_sub[which(ph_sub$ALT != ph_sub$MINOR_ALLELE), paste0(v, '_controls')]
                            ph_sub[which(ph_sub$ALT != ph_sub$MINOR_ALLELE), paste0(v, '_cases')] = 1 - ph_sub[which(ph_sub$ALT != ph_sub$MINOR_ALLELE), paste0(v, '_cases')]
                            # define final_ID as ID (MINOR_ALLELE)
                            ph_sub$SNP = paste0(ph_sub$ID, ' (', ph_sub$MINOR_ALLELE, ")")
                            # Add risk/protective allele
                            ph_sub$SNPID = paste(ph_sub$CHROM, ph_sub$POS, sep = ":")
                            ph_sub_info = merge(ph_sub, snps_data[, c('SNP', 'EFFECT_ALLELE', 'BETA')], by.x = 'SNPID', by.y = 'SNP')
                            ph_sub_info$Effect = ifelse(ph_sub_info$MINOR_ALLELE == ph_sub_info$EFFECT_ALLELE & ph_sub_info$BETA >0, 'Risk', ifelse(ph_sub_info$MINOR_ALLELE == ph_sub_info$EFFECT_ALLELE & ph_sub_info$BETA <0, 'Protective', ifelse(ph_sub_info$MINOR_ALLELE != ph_sub_info$EFFECT_ALLELE & ph_sub_info$BETA >0, 'Protective', 'Risk')))
                            # reshape to long format -- the AD_controls and AD_cases columns should be combined
                            # define the ID columns
                            id_cols <- c('SNPID', 'ID', 'SNP', 'CHROM', 'POS', 'REF', 'ALT', 'MINOR_ALLELE', 'EFFECT_ALLELE', 'BETA', 'Effect')
                            # create long format
                            frequencies_long <- reshape(ph_sub_info, varying = setdiff(names(ph_sub_info), id_cols), v.names = "MAF", timevar = "Group", times = setdiff(names(ph_sub_info), id_cols), direction = "long")
                            # reset rownames if needed
                            rownames(frequencies_long) <- NULL
                            # sort by frequency
                            frequencies_long = frequencies_long[order(frequencies_long$MAF),]
                            # set the order of the SNPs
                            frequencies_long$SNP = factor(frequencies_long$SNP, levels = unique(frequencies_long$SNP))
                            # plot
                            suppressMessages({
                                p4 = ggplot(frequencies_long, aes(x=SNP, y=MAF, color=Group, shape=Group)) + geom_point(size = 3, stat = 'identity') + labs(title="Case/Control Frequencies", x="SNPs", y="Minor Allele Frequency", caption="**<span style='color:navy;'>Blue</span>** labels indicate risk-increasing alleles, **<span style='color:orange;'>Orange</span>** labels indicate protective alleles") + theme_bw() + theme(title = element_text(size = 18), axis.text.x = ggtext::element_markdown(size = 13, angle = 60, hjust = 1, vjust = 1, color = ifelse(frequencies_long$Effect == "Risk", "navy", "orange")), axis.text.y = element_text(size = 14), axis.title = element_text(size = 16), legend.position = 'top', legend.text = element_text(size = 14), legend.title = element_text(size = 16), plot.caption = ggtext::element_markdown(size = 12, hjust = 0.5)) + guides(color = guide_legend(override.aes = list(size = 4)), shape = guide_legend(override.aes = list(size = 4)))
                                # set the width of the plot dynamically
                                base_width = 5
                                width_per_snp = 0.125
                                plot_width = max(5, min(60, base_width + nrow(frequencies_long) * width_per_snp))
                                pdf(paste0(outdir, '/PRS_Frequencies_', v, '.pdf'), height = 7, width = plot_width)
                                print(p4)
                                invisible(dev.off())
                            })
                        }
                    }, error = function(e){
                        cat('** Error in plotting frequencies for case-control. Skipping this plot.\n')
                    })
                } else {
                    tryCatch({
                        cat('**** Plotting SNP frequencies.\n')
                        # subset of data
                        all_freq = data.frame(all_freq, check.names=F)
                        ph_sub = all_freq[, c('ID', '#CHROM', 'POS', 'REF', 'ALT', 'ALT_FREQS')]
                        # rename columns
                        colnames(ph_sub) = c('ID', 'CHROM', 'POS', 'REF', 'ALT', 'ALT_FREQS')
                        # add minor allele
                        ph_sub$MINOR_ALLELE = ph_sub$ALT
                        ph_sub[which(ph_sub$ALT_FREQS > 0.5), 'MINOR_ALLELE'] = ph_sub$REF[which(ph_sub$ALT_FREQS > 0.5)]
                        ph_sub$MAF = ifelse(ph_sub$ALT_FREQS > 0.5, 1 - ph_sub$ALT_FREQS, ph_sub$ALT_FREQS)
                        # define final_ID as ID (MINOR_ALLELE)
                        ph_sub$SNP = paste0(ph_sub$ID, ' (', ph_sub$MINOR_ALLELE, ")")
                        # Add risk/protective allele
                        ph_sub$SNPID = paste(ph_sub$CHROM, ph_sub$POS, sep = ":")
                        ph_sub_info = merge(ph_sub, snps_data[, c('SNP', 'EFFECT_ALLELE', 'BETA')], by.x = 'SNPID', by.y = 'SNP')
                        ph_sub_info$Effect = ifelse(ph_sub_info$MINOR_ALLELE == ph_sub_info$EFFECT_ALLELE & ph_sub_info$BETA >0, 'Risk', ifelse(ph_sub_info$MINOR_ALLELE == ph_sub_info$EFFECT_ALLELE & ph_sub_info$BETA <0, 'Protective', ifelse(ph_sub_info$MINOR_ALLELE != ph_sub_info$EFFECT_ALLELE & ph_sub_info$BETA >0, 'Protective', 'Risk')))
                        # sort
                        ph_sub_info = ph_sub_info[order(ph_sub_info$MAF),]
                        # set the order of the SNPs
                        ph_sub_info$SNP = factor(ph_sub_info$SNP, levels = unique(ph_sub_info$SNP))
                        # plot
                        suppressMessages({
                            p4 = ggplot(data = ph_sub_info, aes(x=SNP, y=MAF)) + geom_point(size = 3, col = 'coral') + labs(title="All samples frequencies", x="SNPs", y="Minor Allele Frequency", caption="**<span style='color:navy;'>Blue</span>** labels indicate risk-increasing alleles, **<span style='color:orange;'>Orange</span>** labels indicate protective alleles") + theme_bw() + theme(title = element_text(size = 18), axis.text.x = ggtext::element_markdown(size = 13, angle = 60, hjust = 1, vjust = 1, color = ifelse(ph_sub_info$Effect == "Risk", "navy", "orange")), axis.text.y = element_text(size = 14), axis.title = element_text(size = 16), legend.position = 'none', plot.caption = ggtext::element_markdown(size = 12, hjust = 0.5)) + guides(color = guide_legend(override.aes = list(size = 4)))
                            # set the width of the plot dynamically
                            base_width = 5
                            width_per_snp = 0.15
                            plot_width = max(5, min(60, base_width + nrow(ph_sub_info) * width_per_snp))
                            pdf(paste0(outdir, '/PRS_Frequencies.pdf'), height = 7, width = plot_width)
                            print(p4)
                            invisible(dev.off())
                        })
                    }, error = function(e){
                        cat('** Error in plotting frequencies. Skipping this plot.\n')
                    })
                }
            }

            # plot tiles
            if (nrow(tiles_prs_df) >0){
                tryCatch({
                    cat('**** Plotting PRS tiles.\n')
                    # add tile names
                    tiles_prs_df$tile_name = paste0(tiles_prs_df$n_tiles, '_Tiles_', tiles_prs_df$variable, '_', tiles_prs_df$group, '_', tiles_prs_df$sex)
                    # density plot lists
                    plotlist = list()
                    proplist = list()
                    # iterate over the tiles
                    for (i in 1:nrow(tiles_prs_df)){
                        # get data of interest
                        tmp_data = prs_df[, c('iid', 'PRS', tiles_prs_df$tile_name[i])]
                        # add phenotype of interest
                        tmp_pheno = assoc_info[[2]][, c(assoc_info[[1]], tiles_prs_df$variable[i])]
                        # merge data
                        tmp_data = merge(tmp_data, tmp_pheno, by.x = 'iid', by.y = assoc_info[[1]], all = T)
                        # rename columns
                        colnames(tmp_data) = c('iid', 'PRS', 'tile', 'variable')
                        # density title
                        tmp_title = paste0('Tiles: ', tiles_prs_df$variable[i], ' (Sex=', tiles_prs_df$sex[i], ')')
                        prs_title = ifelse(suffix == '', paste0('PRS ~ ', tiles_prs_df$n_tiles[i], ' tiles'), paste0('PRS (APOE excluded) ~ ', tiles_prs_df$n_tiles[i], ' tiles'))
                        tmp_data$variable = factor(tmp_data$variable)
                        # find deciles
                        prs_deciles = c()
                        for (x in 1:max(as.numeric(as.character(tmp_data$tile)), na.rm = TRUE)){
                            tmp = tmp_data[which(tmp_data$tile == x), 'PRS']
                            prs_deciles = c(prs_deciles, max(tmp))
                        }
                        # find proportions
                        prs_summary = data.frame(table(tmp_data$tile, tmp_data$variable))
                        colnames(prs_summary) = c('Tile', 'Variable', 'Count')
                        totals <- tapply(prs_summary$Count, prs_summary$Variable, sum)
                        prs_summary$Proportion <- prs_summary$Count / totals[as.character(prs_summary$Variable)]
                        # check if nas need to be excluded
                            if (na_beh == 'exclude'){
                                tmp_data <- tmp_data[!is.na(tmp_data$variable), ]
                            }
                        # plot
                        suppressMessages({
                            # density plot with dashed lines indicating tiles
                            plt_density = ggplot(tmp_data[!is.na(tmp_data$tile),], aes(x=PRS, fill=variable)) + geom_density(alpha=0.5, bw = 0.15) + scale_fill_manual(values=c("#FF9999", "#66B3FF")) + labs(title=tmp_title, x=prs_title, y="Density") + theme_bw() + theme(title = element_text(size = 18), axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), axis.title = element_text(size = 16), legend.position = 'top', legend.text = element_text(size = 14), legend.title = element_text(size = 16)) + geom_vline(xintercept = as.numeric(prs_deciles), linetype = "dashed", color = "red")
                            # proportion plot
                            prp_plot = ggplot(prs_summary, aes(x=Tile, y=Proportion, fill=Variable)) + geom_bar(stat='identity', position='dodge', width=0.7, alpha = 0.7) + labs(title="Proportions per decile and group", x="Decile", y="Proportion") + scale_fill_manual(values=c("#FF9999", "#66B3FF")) + theme_bw() + theme(title = element_text(size = 18), axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), axis.title = element_text(size = 16), legend.position = 'top', legend.text = element_text(size = 14), legend.title = element_text(size = 16))
                            # add the plots to the list
                            plotlist[[i]] = plt_density
                            proplist[[i]] = prp_plot
                        })
                    }
                    # save the plots
                    pdf(paste0(outdir, '/PRS_Tiles', suffix, '.pdf'), height = 14, width = 7*length(plotlist))
                    allplots = c(plotlist, proplist)
                    suppressMessages({
                        combined = ggarrange(plotlist = allplots, ncol = length(plotlist), nrow = 2, common.legend = FALSE)
                        print(combined)
                        invisible(dev.off())
                    })                    
                }, error = function(e){
                    cat('** Error in plotting tiles. Skipping this plot.\n')
                })
            }

            # plot splits
            if (nrow(split_info_df) >0){
                tryCatch({
                    cat('**** Plotting splits.\n')
                    # add split names
                    split_info_df$split_name = paste0('Split_', split_info_df$variable, '_', split_info_df$thresholds, '_', split_info_df$sex)
                    # define container for plots
                    plotlist = list()
                    # iterate over the splits
                    for (i in 1:nrow(split_info_df)){
                        # get data of interest
                        tmp_data = prs_df[, c('iid', 'PRS', split_info_df$split_name[i])]
                        # add phenotype of interest
                        tmp_pheno = assoc_info[[2]][, c(assoc_info[[1]], split_info_df$variable[i])]
                        # merge data
                        tmp_data = merge(tmp_data, tmp_pheno, by.x = 'iid', by.y = assoc_info[[1]], all = T)
                        # rename columns
                        colnames(tmp_data) = c('iid', 'PRS', 'split', 'variable')
                        # density title
                        tmp_title = paste0('PRS across ', split_info_df$n_split[i], ' splits of ', split_info_df$variable[i], ' (Sex=', split_info_df$sex[i], ')', sep = '')
                        # prs title
                        prs_title = ifelse(suffix == '', 'PRS', 'PRS (APOE excluded)')
                        # plot -- violin and boxplots
                        suppressMessages({
                            plt = ggplot(tmp_data[!is.na(tmp_data$split),], aes(x = split, y = PRS, fill = split)) + geom_violin(alpha = 0.5) + geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5) + scale_fill_brewer(palette = "Set3") + labs(title = tmp_title, x = "Splits", y = prs_title) + theme_bw() + theme(title = element_text(size = 18), axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), axis.title = element_text(size = 16), legend.position = 'none')
                        })
                        # add the plot to the list
                        plotlist[[i]] = plt
                    }
                    # save the plots
                    pdf(paste0(outdir, '/PRS_Splits', suffix, '.pdf'), height = 7, width = 7*length(plotlist))
                    suppressMessages({
                        combined = ggarrange(plotlist = plotlist, ncol = length(plotlist), nrow = 1, common.legend = FALSE)
                        print(combined)
                        invisible(dev.off())
                    })
                }, error = function(e){
                    cat('** Error in plotting splits. Skipping this plot.\n')
                })
            }
        })
    }

    # Function to write log file
    writeLog = function(outdir, genotype_file, snps_file, outfile, isdosage, plt, maf, multiple, excludeAPOE, fliprisk, keepDos, addWeight, freq, assoc, assoc_var, assoc_cov, assoc_survival, sex_strata){
        # Define output name
        outname = paste0(outdir, '/run_info.log')
        # Create log info
        info = paste0("Genotype file: ", genotype_file, "\nMultiple files: ", multiple, "\nSNPs file: ", snps_file, "\nOutput file: ", outfile, "\nDosage: ", isdosage, "\nMAF: ", maf, "\nWith and Without APOE: ", excludeAPOE, "\nDirect effects (Risk and Protective): ", fliprisk, "\nKeep dosages: ", keepDos, "\nAdditional weight: ", addWeight, "\nCalculate frequency: ", freq, "\nPlot: ", plt, '\n\n', 'Association file: ', assoc, '\nAssociation variables: ', assoc_var, '\nAssociation covariates: ', assoc_cov, '\nAssociation with survival: ', assoc_survival, '\nSex-stratified: ', sex_strata, '\n\n')
        # Write log file
        writeLines(info, outname)
        return(outname)
    }

    # Function to write output
    write_prs_outputs = function(prs_df, included_snps, snps_data, suffix, outdir, assoc_info, addPheno) {
        # Check if assoc_info is provided
        if (length(assoc_info) >1){
            if (addPheno == TRUE){
                assoc_idname = assoc_info[[1]]
                assoc_data = assoc_info[[2]]
                # merge the prs with association info
                prs_df = merge(prs_df, assoc_data, by.x = 'iid', by.y = assoc_idname)
            }
        }
        # Write PRS table
        write.table(prs_df, paste0(outdir, '/PRS_table', suffix, '.txt'), quote = FALSE, row.names = FALSE, sep = "\t")
        # Construct ID fields
        included_snps$ID <- paste(included_snps$CHROM, included_snps$POS, sep = ":")
        snps_data$ID <- paste(stringr::str_replace_all(snps_data$CHROM, 'chr', ''), snps_data$POS, sep = ":")
        # Identify missing SNPs
        missing <- snps_data[!(snps_data$ID %in% included_snps$ID), ]
        if (nrow(missing) > 0) {
            excluded <- data.frame(SNP = NA, BETA = missing$BETA, ALLELE = missing$EFFECT_ALLELE, OTHER_ALLELE = missing$OTHER_ALLELE, TYPE = 'Excluded', CHROM = missing$CHROM, POS = missing$POS, SNPID_DATA = NA, REASON = 'Variant_not_found', ID = missing$ID)
            included_snps <- rbind(included_snps, excluded)
        }
        # Write SNPs included (and excluded) file
        write.table(included_snps, paste0(outdir, '/SNPs_included_PRS', suffix, '.txt'), quote = FALSE, row.names = FALSE, sep = "\t")
    }

    # Function to check association variables
    checkAssocVar = function(assoc_var, assoc_data, assoc_survival, sex_strata){
        # check if variables were provided
        if (assoc_var[1] == FALSE){
            stop('** No association variables provided.\n\n', call. = FALSE)
        }
        
        # unlist the variables
        assoc_var = unlist(strsplit(assoc_var, ','))
        # check if variables are in the data
        assoc_var = toupper(assoc_var)
        assoc_var = assoc_var[which(assoc_var %in% colnames(assoc_data))]
        if (length(assoc_var) == 0){
            stop('** No association variables found in the data. Did you mispelled the variable names?.\n\n', call. = FALSE)
        }
        
        # check if survival analysis is requested
        if (assoc_survival != FALSE){
            # check if survival variable is among the variables
            assoc_survival = toupper(unlist(strsplit(assoc_survival, ',')))
            # check if all assoc_survival variables are in the assoc_var
            if (!all(assoc_survival %in% toupper(assoc_var))){
                stop('** Survival variables not found in the association variables. Please check the variable names.\n\n', call. = FALSE)
            } else {
                cat('**** Survival analysis requested for variables: ', paste(assoc_survival, collapse = ', '), '\n')
            }
        }
        
        # Define a dataframe with the association variables
        df_variables = data.frame()
        # Check if the variables are numeric or categorical
        for (i in 1:length(assoc_var)){
            # check if the variable is in the data
            if (!(toupper(assoc_var[i]) %in% toupper(colnames(assoc_data)))){
                stop(paste0('** Variable ', assoc_var[i], ' not found in the association data. Please check the variable names.\n\n'), call. = FALSE)
            }
            # check if the variable is survival
            if (toupper(assoc_var[i]) %in% toupper(assoc_survival)){
                cat('**** Variable ', assoc_var[i], ': survival variable. Cox regression will be used.\n')
                # check if there is an event variable
                if (toupper(paste0(assoc_var[i], '_event')) %in% toupper(colnames(assoc_data))){
                    # store the mapping before updating the data
                    cat('**** Event variable found: ', paste0(assoc_var[i], '_EVENT'), '\n')
                    var_mapping = paste0('No mapping needed;', assoc_var[i], '_EVENT')
                } else {
                    cat('**** No event variable found for survival analysis. Using default event variable.\n')
                    var_mapping = paste0('No mapping needed;')
                }
                # update dataframe
                df_variables = rbind(df_variables, data.frame(variable = assoc_var[i], type = 'survival', model = 'cox', mapping = var_mapping, sex = 'all'))
            # check if the variable has 2 values (binary, numerical or categorical)
            } else if (length(table(assoc_data[, assoc_var[i]])) == 2){
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
                    df_variables = rbind(df_variables, data.frame(variable = assoc_var[i], type = 'numeric', model = 'binomial', mapping = var_mapping, sex = 'all'))
                } else {
                    cat('**** Variable ', assoc_var[i], ': categorical with 2 values. Logistic regression will be used.\n')
                    # make sure the smaller is 0 and the larger is 1
                    min_value = min(assoc_data[, assoc_var[i]], na.rm=T)
                    # store the mapping before updating the data
                    var_mapping = paste0(min_value, ' -> 0; ', max(assoc_data[, assoc_var[i]], na.rm=T), ' -> 1')
                    # update values
                    assoc_data[, assoc_var[i]] = ifelse(assoc_data[, assoc_var[i]] == min_value, 0, 1)
                    # update dataframe
                    df_variables = rbind(df_variables, data.frame(variable = assoc_var[i], type = 'categorical', model = 'binomial', mapping = var_mapping, sex = 'all'))
                }
            } else {
                # check if the variable is numeric
                if (is.numeric(assoc_data[, assoc_var[i]])){
                    cat('**** Variable ', assoc_var[i], ': numeric with more than 2 values. Linear regression will be used.\n')
                    # store the mapping before updating the data
                    var_mapping = paste0('No mapping needed')
                    # update dataframe
                    df_variables = rbind(df_variables, data.frame(variable = assoc_var[i], type = 'numeric', model = 'gaussian', mapping = var_mapping, sex = 'all'))
                } else {
                    cat('**** Variable ', assoc_var[i], ': categorical with', length(table(assoc_data[, assoc_var[i]])), 'values. Pairwise logistc regressions will be used.\n')
                    # collect the values of the variable
                    unique_values = unique(assoc_data[, assoc_var[i]])
                    # exclude NA values
                    unique_values = unique_values[!is.na(unique_values)]
                    # sort alphabetically
                    unique_values = sort(unique_values, na.last = TRUE)
                    # get combinations of the values
                    pairs <- combn(unique_values, 2, simplify = FALSE)
                    # iterate over the pairs
                    for (p in pairs){
                        # store the mapping before updating the data
                        var_mapping = paste0(p[1], ' -> 0; ', p[2], ' -> 1')
                        # update values
                        assoc_data[, paste(p, collapse = '_vs_')] = ifelse(assoc_data[, assoc_var[i]] == p[1], 0, ifelse(assoc_data[, assoc_var[i]] == p[2], 1, NA))
                        # update dataframe
                        df_variables = rbind(df_variables, data.frame(variable = paste(p, collapse = '_vs_'), type = 'categorical', model = 'binomial', mapping = var_mapping, sex = 'all'))
                    }
                }
            }
        }
        # check if sex-stratafied analysis is requested
        if (sex_strata){
            # collect the values of SEX
            values_sex = unique(assoc_data[, 'SEX'], na.rm=TRUE)
            # iterate over the values
            for (s in values_sex){
                tmp = df_variables[which(df_variables$sex == 'all'), ]
                # update sex
                tmp$sex = s
                # add to dataframe
                df_variables = rbind(df_variables, tmp)
                df_variables$sex = as.character(df_variables$sex)
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
    checkAssocFile = function(assoc, assoc_var, assoc_cov, assoc_survival, sex_strata){
        # check if file exists
        if (file.exists(assoc)){
            cat('** Association data file found.\n')
            # open file
            assoc_data = fread(assoc, h=T, stringsAsFactors=F, sep="\t")
            colnames(assoc_data) = toupper(colnames(assoc_data))  # convert column names to uppercase for consistency
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
            if (sex_strata){
                # check if SEX variable is present
                if ('SEX' %in% toupper(colnames(assoc_data))){
                    cat('**** Sex stratified analysis requested. SEX column found.\n')
                    # update the sex column to be SEX
                    colnames(assoc_data)[which(toupper(colnames(assoc_data)) == 'SEX')] = 'SEX'
                } else {
                    cat('**** !!! Sex stratified analysis requested but SEX column not found. Association will be run in the whole sample.\n')
                    sex_strata = FALSE
                }
            }
            # check if IID is in the columns
            if ('IID' %in% toupper(colnames(assoc_data))){
                # get the name of the IID column
                idname = colnames(assoc_data)[which(toupper(colnames(assoc_data)) == 'IID')]
                # check variable to associate
                outcome_info = checkAssocVar(assoc_var, assoc_data, assoc_survival, sex_strata)
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
    assoc_test = function(prs_df, assoc_info, outdir, suffix, assoc_mode, dosages, sex_strata, plt){
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
            assoc_results_prs = prs_assoc(prs_df_pheno, assoc_variables, assoc_covariates, suffix, plt, outdir)
            # write the results
            write.table(assoc_results_prs, paste0(outdir, '/association_results_PRS', suffix, '.txt'), quote = F, row.names = F, sep = "\t")
        
        } else if (assoc_mode == 'single'){
            cat('**** Running association for single-variant.\n')
            assoc_results_single = singleVar_assoc(dosages_pheno, assoc_variables, assoc_covariates, suffix, snp_names)
            # write the results
            write.table(assoc_results_single, paste0(outdir, '/association_results_single', suffix, '.txt'), quote = F, row.names = F, sep = "\t")
        
        } else {
            cat('**** Running association for both PRS and single-variant.\n')
            assoc_results_prs = prs_assoc(prs_df_pheno, assoc_variables, assoc_covariates, suffix, plt, outdir)
            assoc_results_single = singleVar_assoc(dosages_pheno, assoc_variables, assoc_covariates, suffix, snp_names)
            # write the results
            write.table(assoc_results_prs, paste0(outdir, '/association_results_PRS', suffix, '.txt'), quote = F, row.names = F, sep = "\t")
            write.table(assoc_results_single, paste0(outdir, '/association_results_single', suffix, '.txt'), quote = F, row.names = F, sep = "\t")
        }
    }

    # Function to run association for the PRS
    prs_assoc = function(prs_df_pheno, assoc_variables, assoc_covariates, suffix, plt, outdir){
        # Scale PRS before association
        prs_df_pheno$PRS = scale(prs_df_pheno$PRS)
        # Define container for the results
        assoc_results = data.frame()
        # iterate over the variables
        for (i in 1:nrow(assoc_variables)){
            # check sex first
            if (assoc_variables$sex[i] != 'all'){
                prs_df_pheno_sex = prs_df_pheno[which(prs_df_pheno$SEX == assoc_variables$sex[i]), ]
                tmp_sex = assoc_variables$sex[i]
            } else {
                prs_df_pheno_sex = prs_df_pheno
                tmp_sex = 'all'
            }
            cat(paste0('****** Running PRS association in SEX=', assoc_variables$sex[i], ' for variable: ', assoc_variables$variable[i], ' ', str_replace_all(suffix, '_', ''), '\n'))
            if (nrow(prs_df_pheno_sex) == 0){
                cat('******** No individuals left. Skipping!\n')
                assoc_results = rbind(assoc_results, data.frame(Predictor = 'PRS', Outcome = paste(var_name, event_variable, sep=";"), Covariates = paste(covar_names, collapse = ','), Beta_PRS = NA, SE_PRS = NA, P_PRS = NA, Model = model_type, N_tot = nrow(prs_df_pheno_sex), N_missing = sum(is.na(prs_df_pheno_sex[, var_name])), Mapping = var_mapping, Model_converged = 'NO', N_effective = 0, N_cases = NA, N_controls = NA, sex = assoc_variables$sex[i], stringsAsFactors = F))
                next
            }
            # get the variable name
            var_name = assoc_variables$variable[i]
            # get the model type
            model_type = assoc_variables$model[i]
            # get the covariates
            covar_names = assoc_covariates
            # get the mapping
            var_mapping = assoc_variables$mapping[i]
            # initialize the survival object
            surv_obj = NULL
            # check if it is survival analysis
            if (model_type == 'cox'){
                # take the event if present
                event_variable = str_split_fixed(assoc_variables$mapping[i], ';', 2)[2]
                if (event_variable == ''){
                    event_variable = NA
                    prs_df_pheno_sex$event = 1
                    # create survival object
                    surv_obj = Surv(time = prs_df_pheno_sex[, var_name], event = prs_df_pheno_sex$event)
                } else {
                    surv_obj = Surv(time = prs_df_pheno_sex[, var_name], event = prs_df_pheno_sex[, event_variable])
                }
            }
            # check if there are covariates
            if (length(na.omit(covar_names)) > 0){
                # remove sex from covariates if it is not 'all'
                if (tmp_sex != 'all'){ covar_names = covar_names[which(covar_names != 'SEX')]}
                # create formula -- check if it is survival analysis
                if (model_type == 'cox'){
                    formula = as.formula(paste0('surv_obj ~ PRS + ', paste(covar_names, collapse = ' + ')))
                } else {
                    formula = as.formula(paste0(var_name, ' ~ PRS + ', paste(covar_names, collapse = ' + ')))
                }
            } else {
                # create formula -- check if it is survival analysis
                if (model_type == 'cox'){
                    formula = as.formula(paste0('surv_obj ~ PRS'))
                } else {
                    formula = as.formula(paste0(var_name, ' ~ PRS'))
                }
            }
            # calculate number of cases and controls for logistic regression
            if (model_type == 'binomial'){
                n_cases = nrow(prs_df_pheno_sex[which(prs_df_pheno_sex[, var_name] == 1),])
                n_controls = nrow(prs_df_pheno_sex[which(prs_df_pheno_sex[, var_name] == 0),])
            } else {
                n_cases = NA
                n_controls = NA
            }
            # perform association test -- check if it is survival analysis
            if (model_type == 'cox'){
                model = suppressWarnings(coxph(formula, data = prs_df_pheno_sex))
            } else {
                model = suppressWarnings(glm(formula, data = prs_df_pheno_sex, family = model_type))
            }
            # add to the results -- check if it is survival analysis
            if (model_type == 'cox'){
                assoc_results = rbind(assoc_results, data.frame(Predictor = 'PRS', Outcome = paste(var_name, event_variable, sep=";"), Covariates = paste(covar_names, collapse = ','), Beta_PRS = summary(model)$coefficients[1, 1], SE_PRS = summary(model)$coefficients[1, 3], P_PRS = summary(model)$coefficients[1, 5], Model = model_type, N_tot = nrow(prs_df_pheno_sex), N_missing = sum(is.na(prs_df_pheno_sex[, var_name])), Mapping = var_mapping, Model_converged = NA, N_effective = nrow(prs_df_pheno_sex) - sum(is.na(prs_df_pheno_sex[, var_name])), N_cases = n_cases, N_controls = n_controls, sex = assoc_variables$sex[i], stringsAsFactors = F))
                # if plot was requested, make the survival plot now
                if (!is.null(plt)){
                    makeSurvPlot(prs_df_pheno_sex, suffix, outdir, assoc_variables$variable[i], assoc_variables$sex[i], surv_obj)
                }
            } else {
                assoc_results = rbind(assoc_results, data.frame(Predictor = 'PRS', Outcome = var_name, Covariates = paste(covar_names, collapse = ','), Beta_PRS = summary(model)$coefficients[2, 1], SE_PRS = summary(model)$coefficients[2, 2], P_PRS = summary(model)$coefficients[2, 4], Model = model_type, N_tot = nrow(prs_df_pheno_sex), N_missing = sum(is.na(prs_df_pheno_sex[, var_name])), Mapping = var_mapping, Model_converged = model$converged, N_effective = nrow(prs_df_pheno_sex) - sum(is.na(prs_df_pheno_sex[, var_name])), N_cases = n_cases, N_controls = n_controls, sex = assoc_variables$sex[i], stringsAsFactors = F))
            }
        }
        return(assoc_results)
    }

    # Function to make survival plot
    makeSurvPlot = function(prs_df_pheno_sex, suffix, outdir, variable, sexinfo, surv_obj){
        suppressMessages({
            suppressWarnings({
                # split by median value of PRS
                prs_df_pheno_sex$PRS_group <- ifelse(prs_df_pheno_sex$PRS >= median(prs_df_pheno_sex$PRS, na.rm = TRUE), "High", "Low")
                prs_df_pheno_sex$PRS_group <- factor(prs_df_pheno_sex$PRS_group, levels = c("Low", "High"))
                # bind survival object as column
                prs_df_pheno_sex$.__surv__ <- surv_obj
                # Fit model using the injected survival column
                fit <- survfit(.__surv__ ~ PRS_group, data = prs_df_pheno_sex)
                # Custom theme
                custom_theme <- theme_minimal() + theme(legend.position = "bottom", axis.text = element_text(size = 14), axis.title = element_text(size = 16), legend.text = element_text(size = 14), legend.title = element_text(size = 16))
                # Plot
                plt_title = ifelse(suffix == '', paste0('Survival: PRS~', variable, ' (Sex=', sexinfo, ')'), paste0('Survival: PRS~', variable, ' ~ APOE excluded (Sex=', sexinfo, ')'))
                plt = ggsurvplot(fit, data = prs_df_pheno_sex, risk.table = TRUE, pval = TRUE, conf.int = TRUE, risk.table.height = 0.25, ggtheme = custom_theme, title = plt_title)
                # Save the plot
                pdf(paste0(outdir, '/survival_plot_', variable, '_', sexinfo, suffix, '.pdf'), width = 8, height = 8)
                print(plt, newpage = FALSE)
                dev.off()
            })
        })
    }

    # Function to run association for the single-variants
    singleVar_assoc = function(dosages_pheno, assoc_variables, assoc_covariates, suffix, snp_names){
        # make sure data is a dataframe
        dosages_pheno = data.frame(dosages_pheno, check.names=F)
        # Define container for the results
        assoc_results = data.frame()
        # iterate over the variables
        for (i in 1:nrow(assoc_variables)){
            # check sex first
            if (assoc_variables$sex[i] != 'all'){
                dosages_pheno_sex = dosages_pheno[which(dosages_pheno$SEX == assoc_variables$sex[i]), ]
                tmp_sex = assoc_variables$sex[i]
            } else {
                dosages_pheno_sex = dosages_pheno
                tmp_sex = 'all'
            }
            cat(paste0('****** Running single-variant analysis in SEX=', assoc_variables$sex[i], ' for variable: ', assoc_variables$variable[i], ' ', str_replace_all(suffix, '_', ''), '\n'))
            if (nrow(dosages_pheno_sex) == 0){
                cat('******** No individuals left. Skipping!\n')
                assoc_results = rbind(assoc_results, data.frame(Predictor = 'PRS', Outcome = paste(var_name, event_variable, sep=";"), Covariates = paste(covar_names, collapse = ','), Beta_PRS = NA, SE_PRS = NA, P_PRS = NA, Model = model_type, N_tot = nrow(dosages_pheno_sex), N_missing = sum(is.na(dosages_pheno_sex[, var_name])), Mapping = var_mapping, Model_converged = 'NO', N_effective = 0, N_cases = NA, N_controls = NA, sex = assoc_variables$sex[i], stringsAsFactors = F))
                next
            }
            # get the variable name
            var_name = assoc_variables$variable[i]
            # get the model type
            model_type = assoc_variables$model[i]
            # get the covariates
            covar_names = assoc_covariates
            # get the mapping
            var_mapping = assoc_variables$mapping[i]
            # check if it is survival analysis
            if (model_type == 'cox'){
                # take the event if present
                event_variable = str_split_fixed(assoc_variables$mapping[i], ';', 2)[2]
                if (event_variable == ''){
                    event_variable = NA
                    dosages_pheno_sex$event = 1
                    # create survival object
                    surv_obj = Surv(time = dosages_pheno_sex[, var_name], event = dosages_pheno_sex$event)
                } else {
                    surv_obj = Surv(time = dosages_pheno_sex[, var_name], event = dosages_pheno_sex[, event_variable])
                }
            }
            # then iterate over the snps
            for (snp in snp_names){
                tryCatch({
                    # update the snp name as otherwise there are bad characters in the name (e.g. /)
                    colnames(dosages_pheno_sex)[which(colnames(dosages_pheno_sex) == snp)] = 'SNP'
                    # check if there are covariates
                    if (length(na.omit(covar_names)) > 0){
                        # remove sex from covariates if it is not 'all'
                        if (tmp_sex != 'all'){ covar_names = covar_names[which(covar_names != 'SEX')]}
                        # create formula -- check if it is survival analysis
                        if (model_type == 'cox'){
                            formula = as.formula(paste0('surv_obj ~ SNP + ', paste(covar_names, collapse = ' + ')))
                        } else {
                            formula = as.formula(paste0(var_name, ' ~ SNP + ', paste(covar_names, collapse = ' + ')))
                        }
                    } else {
                        # create formula -- check if it is survival analysis
                        if (model_type == 'cox'){
                            formula = as.formula(paste0('surv_obj ~ SNP'))
                        } else {
                            formula = as.formula(paste0(var_name, ' ~ SNP'))
                        }
                    }
                    # calculate number of cases and controls for logistic regression
                    if (model_type == 'binomial'){
                        n_cases = nrow(dosages_pheno_sex[which(dosages_pheno_sex[, var_name] == 1),])
                        n_controls = nrow(dosages_pheno_sex[which(dosages_pheno_sex[, var_name] == 0),])
                    } else {
                        n_cases = NA
                        n_controls = NA
                    }
                    # perform association test -- check if it is survival analysis
                    if (model_type == 'cox'){
                        model = suppressWarnings(coxph(formula, data = dosages_pheno_sex))
                    } else {
                        model = suppressWarnings(glm(formula, data = dosages_pheno_sex, family = model_type))
                    }
                    # add to the results -- check if it is survival analysis
                    if (model_type == 'cox'){
                        assoc_results = rbind(assoc_results, data.frame(Predictor = snp, Outcome = paste(var_name, event_variable, sep=";"), Covariates = paste(covar_names, collapse = ','), Beta_PRS = summary(model)$coefficients[1, 1], SE_PRS = summary(model)$coefficients[1, 3], P_PRS = summary(model)$coefficients[1, 5], Model = model_type, N_tot = nrow(dosages_pheno_sex), N_missing = sum(is.na(dosages_pheno_sex[, var_name])), Mapping = var_mapping, Model_converged = NA, N_effective = nrow(dosages_pheno_sex) - sum(is.na(dosages_pheno_sex[, var_name])), N_cases = n_cases, N_controls = n_controls, sex = tmp_sex, stringsAsFactors = F))
                    } else {
                        assoc_results = rbind(assoc_results, data.frame(Predictor = snp, Outcome = var_name, Covariates = paste(covar_names, collapse = ','), Beta_PRS = summary(model)$coefficients[2, 1], SE_PRS = summary(model)$coefficients[2, 2], P_PRS = summary(model)$coefficients[2, 4], Model = model_type, N_tot = nrow(dosages_pheno_sex), N_missing = sum(is.na(dosages_pheno_sex[, var_name])), Mapping = var_mapping, Model_converged = model$converged, N_effective = nrow(dosages_pheno_sex) - sum(is.na(dosages_pheno_sex[, var_name])), N_cases = n_cases, N_controls = n_controls, sex = tmp_sex, stringsAsFactors = F))
                    }
                    # restore the name
                    colnames(dosages_pheno_sex)[which(colnames(dosages_pheno_sex) == 'SNP')] = snp
                },
                error = function(e){
                    # add to the results
                    assoc_results = rbind(assoc_results, data.frame(Predictor = snp, Outcome = var_name, Covariates = paste(covar_names, collapse = ','), Beta_PRS = NA, SE_PRS = NA, P_PRS = NA, Model = model_type, N_tot = nrow(dosages_pheno_sex), N_missing = sum(is.na(dosages_pheno_sex[, var_name])), Mapping = var_mapping, Model_converged = NA, N_effective = nrow(dosages_pheno_sex) - sum(is.na(dosages_pheno_sex[, var_name])), N_cases = NA, N_controls = NA, sex = tmp_sex, stringsAsFactors = F))
                })
            }
        }
        return(assoc_results)
    }

    # Function to check tiles
    checkTiles = function(tiles_prs, assoc_info, sex_strata){
        # split the tiles by comma
        tiles_prs = unlist(strsplit(tiles_prs, ','))
        # remove spaces
        tiles_prs = gsub(" ", "", tiles_prs)
        # check if we have tiles
        if (length(tiles_prs) == 0 || all(tiles_prs == '')){
            stop('**** No tiles provided for PRS calculation. Please check examples and try again.\n\n', call. = FALSE)
        }
        # dataframe to store the tiles
        tiles_info_df = data.frame()
        # iterate over the tiles
        for (tile in tiles_prs){
            # split the tile by colon
            tile_info = unlist(strsplit(tile, ';'))
            # check if the tile has 3 parts
            if (length(tile_info) != 3){
                stop(paste0('**** Tile ', tile, ' is not valid. Please provide tiles in the format n_tiles;variable;group.\n\n'), call. = FALSE)
            }
            # check if the variable is in the association data
            if (toupper(tile_info[2]) %in% colnames(assoc_info[[2]])){
                # check for sex stratification
                if (sex_strata){
                    # find sexes code
                    sex_codes = unique(assoc_info[[3]]$sex)
                    tiles_info_df = rbind(tiles_info_df, data.frame(n_tiles = rep(as.numeric(tile_info[1]), length(sex_codes)), variable = rep(toupper(tile_info[2]), length(sex_codes)), group = rep(tile_info[3], length(sex_codes)), sex = sex_codes, stringsAsFactors = FALSE))
                } else {
                    tiles_info_df = rbind(tiles_info_df, data.frame(n_tiles = as.numeric(tile_info[1]), variable = toupper(tile_info[2]), group = tile_info[3], sex = 'all', stringsAsFactors = FALSE))
                }
            } else {
                stop(paste0('**** Variable ', tile_info[2], ' not found in the association data. Please check the variable names.\n\n'), call. = FALSE)
            }
        }
        return(tiles_info_df)
    }

    # Function to make tiles -- to fix from here
    makeTiles = function(res_prs, tiles_prs, assoc_info){
        # read mapping info from the association info
        mapping_info = assoc_info[[3]]
        # define output dataframe
        res_prs_with_tiles = data.frame()
        # iterate over tiles
        for (tile in 1:nrow(tiles_prs)){
            # get the tile info
            n_tiles = tiles_prs$n_tiles[tile]
            variable = tiles_prs$variable[tile]
            group = tiles_prs$group[tile]
            orig_group = tiles_prs$group[tile]
            sex_info = tiles_prs$sex[tile]
            # take the variable from the association info and add it to the PRS results
            if (sex_info == 'all'){
                tmp_assoc = assoc_info[[2]][, c(assoc_info[[1]], variable)]
                tmp_assoc$SEX = 'all'
            } else {
                tmp_assoc = assoc_info[[2]][, c(assoc_info[[1]], variable, 'SEX')]
            }
            tmp_assoc_prs = merge(res_prs, tmp_assoc, by.x = 'iid', by.y = assoc_info[[1]], all.x = TRUE)
            # Check if the reference group is present or not
            if (group == 'NA'){
                # Compute decile cutoffs from the whole sample
                tmp_prs_tiles = quantile(tmp_assoc_prs$PRS[which(tmp_assoc_prs$SEX == sex_info)], probs = seq(0, 1, by = 1/n_tiles), na.rm = TRUE)
            } else {
                # Find the control group
                tmp_assoc_info = mapping_info[which(mapping_info$variable == variable & mapping_info$sex == sex_info),]
                if (tmp_assoc_info$model == 'binomial'){
                    tmp_mapping_orig_controls = str_split_fixed(str_split_fixed(tmp_assoc_info$mapping, ';', 2)[1], ' ', 3)[1]
                    tmp_mapping_new_controls = str_split_fixed(str_split_fixed(tmp_assoc_info$mapping, ';', 2)[1], ' ', 3)[3]
                    tmp_mapping_orig_cases = str_split_fixed(str_split_fixed(tmp_assoc_info$mapping, ';', 2)[2], ' ', 4)[2]
                    tmp_mapping_new_cases = str_split_fixed(str_split_fixed(tmp_assoc_info$mapping, ';', 2)[2], ' ', 4)[4]
                    # check if the group is
                    group = ifelse(group == tmp_mapping_orig_controls, tmp_mapping_new_controls, tmp_mapping_new_cases)
                }
                # Compute decile cutoffs from the control group
                tmp_prs_tiles = quantile(tmp_assoc_prs$PRS[which((tmp_assoc_prs[, variable] == group) & (tmp_assoc_prs[, 'SEX'] == sex_info))], probs = seq(0, 1, by = 1/n_tiles), na.rm = TRUE)
            }
            # Create a new column in the data frame to assign deciles
            tmp_assoc_prs[, paste0(n_tiles, '_Tiles_', variable, '_', orig_group, '_', sex_info)] = cut(tmp_assoc_prs$PRS, breaks = tmp_prs_tiles, labels = 1:n_tiles, include.lowest = TRUE, right = TRUE)
            # Check for sex stratification -- if so we need to remove the other sex
            if (sex_info != 'all'){
                tmp_assoc_prs[which(tmp_assoc_prs$SEX != sex_info), paste0(n_tiles, '_Tiles_', variable, '_', orig_group, '_', sex_info)] = NA
                tmp_assoc_prs[is.na(tmp_assoc_prs$SEX), paste0(n_tiles, '_Tiles_', variable, '_', orig_group, '_', sex_info)] = NA
            }
            # add to the results
            if (nrow(res_prs_with_tiles) == 0){
                res_prs_with_tiles = tmp_assoc_prs[, c('iid', 'PRS', paste0(n_tiles, '_Tiles_', variable, '_', orig_group, '_', sex_info))]
            } else {
                # merge with the previous results
                res_prs_with_tiles = merge(res_prs_with_tiles, tmp_assoc_prs[, c('iid', paste0(n_tiles, '_Tiles_', variable, '_', orig_group, '_', sex_info))], by = 'iid', all.x = TRUE)
            }
        }
        # remove sex column
        res_prs_with_tiles$SEX <- NULL
        return(res_prs_with_tiles)
    }

    # Function to check split info
    checkSplit = function(split_info, assoc_info, sex_strata){
        # split by comma
        split_prs = unlist(strsplit(split_info, ','))
        # remove spaces
        split_prs = gsub(" ", "", split_prs)
        # check if we have splits
        if (length(split_prs) == 0 || all(split_prs == '')){
            stop('**** No split information provided for PRS calculation. Please check examples and try again.\n\n', call. = FALSE)
        }
        # dataframe to store the tiles
        split_info_df = data.frame()
        # iterate over the tiles
        for (spl in split_prs){
            # split the tile by colon
            split_info = unlist(strsplit(spl, ';'))
            # check if the tile has 2 parts
            if (length(split_info) != 2){
                stop(paste0('**** Split info ', spl, ' is not valid. Please provide splits in the format variable:threshold.\n\n'), call. = FALSE)
            }
            # check if the variable is in the association data
            if (toupper(split_info[1]) %in% toupper(assoc_info[[3]]$variable)){
                if (sex_strata){
                    sex_codes = unique(assoc_info[[3]]$sex)
                    split_info_df = rbind(split_info_df, data.frame(n_split = rep(length(strsplit(split_info[2], '-')[[1]])+1, length(sex_codes)), variable = rep(toupper(split_info[1]), length(sex_codes)), thresholds = rep(split_info[2], length(sex_codes)), sex = sex_codes, stringsAsFactors = FALSE))
                } else {
                    split_info_df = rbind(split_info_df, data.frame(n_split = length(strsplit(split_info[2], '-')[[1]])+1, variable = toupper(split_info[1]), thresholds = split_info[2], sex = 'all', stringsAsFactors = FALSE))
                }
            } else {
                stop(paste0('**** Variable ', tile_info[2], ' not found in the association data. Please check the variable names.\n\n'), call. = FALSE)
            }
        }
        return(split_info_df)
    }

    # Function to make splits
    makeSplit = function(res_prs, split_info, assoc_info){
        # define output dataframe
        res_prs_with_split = data.frame()
        # iterate over splits
        for (spl in 1:nrow(split_info)){
            # get the split info
            variable = split_info$variable[spl]
            thresholds = c(-Inf, as.numeric(unlist(strsplit(split_info$thresholds[spl], '-'))), Inf)
            sex_info = split_info$sex[spl]
            # take the variable from the association info and add it to the PRS results
            if (sex_info == 'all'){
                tmp_assoc = assoc_info[[2]][, c(assoc_info[[1]], variable)]
                tmp_assoc$SEX = 'all'
            } else {
                tmp_assoc = assoc_info[[2]][, c(assoc_info[[1]], variable, 'SEX')]
            }
            tmp_assoc_prs = merge(res_prs, tmp_assoc, by.x = 'iid', by.y = assoc_info[[1]], all.x = TRUE)
            # Create a new column in the data frame to assign splits
            tmp_assoc_prs[, paste0('Split_', variable, '_', split_info$thresholds[spl], '_', sex_info)] = cut(tmp_assoc_prs[, variable], breaks = thresholds, labels = 1:split_info$n_split[spl], include.lowest = TRUE, right = TRUE)
            # check for sex stratification -- if so we need to remove the other sex
            if (sex_info != 'all'){
                tmp_assoc_prs[which(tmp_assoc$SEX != sex_info), paste0('Split_', variable, '_', split_info$thresholds[spl], '_', sex_info)] = NA
                tmp_assoc_prs[is.na(tmp_assoc$SEX), paste0('Split_', variable, '_', split_info$thresholds[spl], '_', sex_info)] = NA
            }
            # add to the results
            if (nrow(res_prs_with_split) == 0){
                tmp_assoc_prs <- tmp_assoc_prs[ , !(names(tmp_assoc_prs) %in% variable)]
                res_prs_with_split = tmp_assoc_prs
            } else {
                # merge with the previous results
                res_prs_with_split = merge(res_prs_with_split, tmp_assoc_prs[, c('iid', paste0('Split_', variable, '_', split_info$thresholds[spl], '_', sex_info))], by = 'iid', all.x = TRUE)
            }
        }
        # remove sex column
        res_prs_with_split$SEX <- NULL
        return(res_prs_with_split)
    }

    # Function to do association of the tiles/splits
    assoc_split_tiles_test = function(prs_df, assoc_info, outdir, suffix, split_info_df, tiles_prs_df){
        # extract association info
        assoc_idname = assoc_info[[1]]
        assoc_data = assoc_info[[2]]
        assoc_variables = assoc_info[[3]]
        assoc_covariates = assoc_info[[4]]

        # merge the prs with association info
        prs_df_pheno = merge(prs_df, assoc_data, by.x = 'iid', by.y = assoc_idname)

        # test tiles
        if (nrow(tiles_prs_df) > 0){
            tiles_results = test_tiles(tiles_prs_df, prs_df_pheno, assoc_covariates, assoc_variables, assoc_idname)
        } else {
            tiles_results = data.frame()
        }
        # test splits
        if (nrow(split_info_df) > 0){
            splits_results = test_splits(split_info_df, prs_df_pheno, assoc_covariates, assoc_variables, assoc_idname)
        } else {
            splits_results = data.frame()
        } 

        # write the results
        write.table(tiles_results, paste0(outdir, '/PRS_results_with_tiles', suffix, '.txt'), quote = F, row.names = F, sep = "\t")
        write.table(splits_results, paste0(outdir, '/PRS_results_with_splits', suffix, '.txt'), quote = F, row.names = F, sep = "\t")
        return(list(tiles_results, splits_results))
    }

    # Function to test tiles
    test_tiles = function(tiles_prs_df, prs_df_pheno, assoc_covariates, assoc_variables, assoc_idname){
        # iterate over the tiles
        res_tiles = data.frame()
        for (tile in 1:nrow(tiles_prs_df)){
            # get variable name
            tmp_var_name = paste0(tiles_prs_df$n_tiles[tile], '_Tiles_', tiles_prs_df$variable[tile], '_', tiles_prs_df$group[tile], '_', tiles_prs_df$sex[tile])
            # get variable info from association
            tmp_var_info = assoc_variables[which(assoc_variables$variable == tiles_prs_df$variable[tile] & assoc_variables$sex == tiles_prs_df$sex[tile]), ]
            # subset the data for the tile of interest
            tmp_data = prs_df_pheno[, c('iid', tmp_var_name, tiles_prs_df$variable[tile], assoc_covariates)]
            colnames(tmp_data)[2] = 'Decile'
            # make sure decile is a factor
            tmp_data$Decile = factor(tmp_data$Decile, levels = rev(1:tiles_prs_df$n_tiles[tile]))
            # define the formula for the model depending on the model type -- check for sex type
            if (tiles_prs_df$sex[tile] != 'all'){
                formula = as.formula(paste0(tmp_var_info$variable, ' ~ ', 'Decile + ', paste(assoc_covariates[which(assoc_covariates != 'SEX')], collapse = ' + ')))
            } else {
                formula = as.formula(paste0(tmp_var_info$variable, ' ~ ', 'Decile + ', paste(assoc_covariates, collapse = ' + ')))
            }
            # define the model 
            model = suppressWarnings(glm(formula, data = tmp_data, family = tmp_var_info$model))
            # save summary of the model
            model_summary = data.frame(summary(model)$coefficients)
            # clean info 
            model_summary$Term = rownames(model_summary)
            model_summary = model_summary[which(model_summary$Term != '(Intercept)'), ]
            colnames(model_summary) = c('Beta', 'SE', 'Statistic', 'P', 'Term')
            model_summary$Variable = tmp_var_name
            rownames(model_summary) = NULL
            model_summary$Sex = tiles_prs_df$sex[tile]
            res_tiles = rbind(res_tiles, model_summary)
        }
        return(res_tiles)
    }

    # Function to test splits
    test_splits = function(split_info_df, prs_df_pheno, assoc_covariates, assoc_variables, assoc_idname){
        # define container for the results
        res_splits = data.frame()
        # iterate over the splits
        for (spl in 1:nrow(split_info_df)){
            # get variable name
            tmp_var_name = paste0('Split_', split_info_df$variable[spl], '_', split_info_df$thresholds[spl], '_', split_info_df$sex[spl])
            # subset the data for the split of interest
            tmp_data = prs_df_pheno[, c('iid', tmp_var_name, 'PRS', assoc_covariates)]
            colnames(tmp_data)[2] = 'Split'
            # make sure split is a factor
            tmp_data$Split = factor(tmp_data$Split, levels = 1:split_info_df$n_split[spl])
            # make formula -- take sex into account
            if (split_info_df$sex[spl] != 'all'){
                formula = as.formula(paste0('PRS ~ Split + ', paste(assoc_covariates[which(assoc_covariates != 'SEX')], collapse = ' + ')))
            } else {
                formula = as.formula(paste0('PRS ~ Split + ', paste(assoc_covariates, collapse = ' + ')))
            }
            # define the model
            model = suppressWarnings(lm(formula, data = tmp_data))
            # save summary of the model
            model_summary = data.frame(summary(model)$coefficients)
            # clean info 
            model_summary$Term = rownames(model_summary)
            model_summary = model_summary[which(model_summary$Term != '(Intercept)'), ]
            colnames(model_summary) = c('Beta', 'SE', 'Statistic', 'P', 'Term')
            model_summary$Variable = tmp_var_name
            rownames(model_summary) = NULL
            model_summary$Sex = split_info_df$sex[spl]
            res_splits = rbind(res_splits, model_summary)
        }
        return(res_splits)
    }


    # Function to compile phenotype report
    reportPhenoStats = function(res_prs, assoc_info, outdir, assoc_mode, assoc_var, assoc_cov, assoc_survival, suffix){
        # Define output file
        report_file = paste0(outdir, '/phenotype_report.txt')
        # Check if the file already exists
        if (!file.exists(report_file)){
            # Write the header
            cat(paste0("Phenotype Report\n", "========================================\n\n"), file = report_file, append = TRUE)
        }

        # Combine phenotypes with prs
        combined_data = merge(assoc_info[[2]], res_prs, by.x = "IID", by.y = 'iid', all.x = TRUE)
        # Isolate the association info
        assoc_info_df = assoc_info[[3]]
        # exclude sex info as it is already handled separately
        assoc_info_df = assoc_info_df[which(assoc_info_df$sex == 'all'), ]

        # All individuals
            # total number of individuals
            total_individuals = nrow(res_prs)
            cat(paste0("Total number of individuals: ", total_individuals, "\n"), file = report_file, append = TRUE)
            
            # total number of individuals with PRS
            total_individuals_prs = nrow(res_prs[!is.na(res_prs$PRS), ])
            cat(paste0("Total number of individuals with PRS: ", total_individuals_prs, "\n\n"), file = report_file, append = TRUE)
            
            # age if present
            if ('AGE' %in% toupper(colnames(combined_data))){
                # Derive mean, median, sd, min, max, and IQR
                mean_age = round(mean(combined_data$AGE, na.rm = TRUE), 3)
                median_age = round(median(combined_data$AGE, na.rm = TRUE), 3)
                sd_age = round(sd(combined_data$AGE, na.rm = TRUE), 3)
                min_age = round(min(combined_data$AGE, na.rm = TRUE), 3)
                max_age = round(max(combined_data$AGE, na.rm = TRUE), 3)
                iqr_age = paste(round(quantile(combined_data$AGE, probs = c(0.25, 0.75)), 3), collapse = '-')
                cat(paste0("Age Statistics ", suffix, "\n"), file = report_file, append = TRUE)
                cat("Mean\tMedian\tSD\tMin\tMax\tIQR\n", file = report_file, append = TRUE)
                cat(paste0(mean_age, "\t", median_age, "\t", sd_age, "\t", min_age, "\t", max_age, "\t", iqr_age), "\n\n", file = report_file, append = TRUE)
            }

            # PRS stats -- all individuals
            mean_prs = round(mean(res_prs$PRS, na.rm = TRUE), 3)
            median_prs = round(median(res_prs$PRS, na.rm = TRUE), 3)
            sd_prs = round(sd(res_prs$PRS, na.rm = TRUE), 3)
            min_prs = round(min(res_prs$PRS, na.rm = TRUE), 3)
            max_prs = round(max(res_prs$PRS, na.rm = TRUE), 3)
            iqr_prs = paste(round(quantile(res_prs$PRS, probs = c(0.25, 0.75)), 3), collapse = '-')
            cat(paste0("PRS Statistics - All Individuals ", suffix, "\n"), file = report_file, append = TRUE)
            cat("Mean\tMedian\tSD\tMin\tMax\tIQR\n", file = report_file, append = TRUE)
            cat(paste0(mean_prs, "\t", median_prs, "\t", sd_prs, "\t", min_prs, "\t", max_prs, "\t", iqr_prs), "\n\n", file = report_file, append = TRUE)
    
            # sex-stratified
            if ('SEX' %in% toupper(colnames(combined_data))){
                # sex table
                tmp_table = table(combined_data$SEX)
                cat('Sex information on all samples:\n', file = report_file, append = TRUE)
                cat(paste0(names(tmp_table), ": ", as.numeric(tmp_table), collapse = ', '), "\n", file = report_file, append = TRUE)
                cat("\n", file = report_file, append = TRUE)
                # sex table of individuals with PRS
                tmp_table_prs = table(combined_data$SEX[!is.na(combined_data$PRS)])
                cat('Sex information on individuals with PRS:\n', file = report_file, append = TRUE)
                cat(paste0(names(tmp_table_prs), ": ", as.numeric(tmp_table_prs), collapse = ', '), "\n\n", file = report_file, append = TRUE)

                # PRS stats by sex
                tmp_sexes = unique(combined_data$SEX)
                for (sex in tmp_sexes){
                    # age if present
                    if ('AGE' %in% toupper(colnames(combined_data))){
                        # Derive mean, median, sd, min, max, and IQR
                        mean_age = round(mean(combined_data$AGE[combined_data$SEX == sex], na.rm = TRUE), 3)
                        median_age = round(median(combined_data$AGE[combined_data$SEX == sex], na.rm = TRUE), 3)
                        sd_age = round(sd(combined_data$AGE[combined_data$SEX == sex], na.rm = TRUE), 3)
                        min_age = round(min(combined_data$AGE[combined_data$SEX == sex], na.rm = TRUE), 3)
                        max_age = round(max(combined_data$AGE[combined_data$SEX == sex], na.rm = TRUE), 3)
                        iqr_age = paste(round(quantile(combined_data$AGE[combined_data$SEX == sex], probs = c(0.25, 0.75)), 3), collapse = '-')
                        # Write the stats to the file
                        cat(paste0("Age Statistics by Sex: ", sex, ' ', suffix, "\n"), file = report_file, append = TRUE)
                        cat("Mean\tMedian\tSD\tMin\tMax\tIQR\n", file = report_file, append = TRUE)
                        cat(paste0(mean_age, "\t", median_age, "\t", sd_age, "\t", min_age, "\t", max_age, "\t", iqr_age), "\n\n", file = report_file, append = TRUE)
                    }

                    cat(paste0("PRS Statistics - Sex:", sex, " ", suffix, "\n"), file = report_file, append = TRUE)
                    cat("Mean\tMedian\tSD\tMin\tMax\tIQR\n", file = report_file, append = TRUE)
                    # Derive mean, median, sd, min, max, and IQR
                    mn = round(mean(res_prs$PRS[combined_data$SEX == sex], na.rm = TRUE), 3)
                    md = round(median(res_prs$PRS[combined_data$SEX == sex], na.rm = TRUE), 3)
                    sd_val = round(sd(res_prs$PRS[combined_data$SEX == sex], na.rm = TRUE), 3)
                    min_val = round(min(res_prs$PRS[combined_data$SEX == sex], na.rm = TRUE), 3)
                    max_val = round(max(res_prs$PRS[combined_data$SEX == sex], na.rm = TRUE), 3)
                    iqr_val = paste(round(quantile(res_prs$PRS[combined_data$SEX == sex], probs = c(0.25, 0.75)), 3), collapse = '-')
                    # Write the stats to the file
                    cat(paste0(mn, "\t", md, "\t", sd_val, "\t", min_val, "\t", max_val, "\t", iqr_val), "\n\n", file = report_file, append = TRUE)
                }
            } else {
                tmp_sexes = ''
            }
            cat("========================================\n", file = report_file, append = TRUE)
            cat("\n\n", file = report_file, append = TRUE)

        # Per-variable
        # Iterate over the association variables
        for (i in 1:nrow(assoc_info_df)){
            # Get the variable names and info
            tmp_var = assoc_info_df$variable[i]
            tmp_var_type = assoc_info_df$type[i]
            tmp_var_model = assoc_info_df$model[i]
            tmp_var_mapping = assoc_info_df$mapping[i]

            # Number of individuals with phenotype
                # Write the header for the variable
                cat(paste0("Phenotype Statistics: ", tmp_var, " ", suffix, "\n"), file = report_file, append = TRUE)
                # Number of individuals with the phenotype
                n_individuals = nrow(combined_data[!is.na(combined_data[, tmp_var]), ])
                cat(paste0("Number of individuals with phenotype: ", tmp_var, ": ", n_individuals, "\n"), file = report_file, append = TRUE)
                # Number of individuals with the phenotype and PRS
                n_individuals_prs = nrow(combined_data[!is.na(combined_data[, tmp_var]) & !is.na(combined_data$PRS), ])
                cat(paste0("Number of individuals with phenotype: ", tmp_var, " and PRS: ", n_individuals_prs, "\n"), file = report_file, append = TRUE)
                # sex strata
                if ('SEX' %in% toupper(colnames(combined_data))){
                    # table of sex and phenotype
                    tmp_table = table(combined_data$SEX[!is.na(combined_data[, tmp_var])])
                    # Write the table to the file
                    cat(paste0("Number of individuals with phenotype: ", tmp_var, " by SEX:\n"), file = report_file, append = TRUE)
                    cat(paste0(names(tmp_table), ": ", as.numeric(tmp_table), collapse = ', '), "\n", file = report_file, append = TRUE)
                    # table of sex and phenotype with PRS
                    tmp_table_prs = table(combined_data$SEX[!is.na(combined_data[, tmp_var]) & !is.na(combined_data$PRS)])
                    # Write the table to the file
                    cat(paste0("Number of individuals with phenotype: ", tmp_var, " and PRS by SEX:\n"), file = report_file, append = TRUE)
                    cat(paste0(names(tmp_table_prs), ": ", as.numeric(tmp_table_prs), collapse = ', '), "\n", file = report_file, append = TRUE)
                }

            # Calculate statistics based on the variable type
            if (tmp_var_model == 'gaussian'){
                cat("\n", file = report_file, append = TRUE)
                # Calculate mean, median, sd, min, max, and IQR of the variable
                mean_val = round(mean(combined_data[!is.na(combined_data[, tmp_var]), tmp_var], na.rm = TRUE), 3)
                median_val = round(median(combined_data[!is.na(combined_data[, tmp_var]), tmp_var], na.rm = TRUE), 3)
                sd_val = round(sd(combined_data[!is.na(combined_data[, tmp_var]), tmp_var], na.rm = TRUE), 3)
                min_val = round(min(combined_data[!is.na(combined_data[, tmp_var]), tmp_var], na.rm = TRUE), 3)
                max_val = round(max(combined_data[!is.na(combined_data[, tmp_var]), tmp_var], na.rm = TRUE), 3)
                iqr_val = paste(round(quantile(combined_data[!is.na(combined_data[, tmp_var]), tmp_var], probs = c(0.25, 0.75), na.rm = TRUE), 3), collapse = '-')
                # Write the stats to the file
                cat(paste0(tmp_var, " Statistics - Phenotype: ", tmp_var, " ", suffix, "\n"), file = report_file, append = TRUE)
                cat("Mean\tMedian\tSD\tMin\tMax\tIQR\n", file = report_file, append = TRUE)
                cat(paste0(mean_val, "\t", median_val, "\t", sd_val, "\t", min_val, "\t", max_val, "\t", iqr_val), "\n", file = report_file, append = TRUE)
                cat("\n", file = report_file, append = TRUE)

                # The related to the PRS
                # Calculate mean, median, sd, min, max, and IQR of the variable
                mean_val = round(mean(combined_data[!is.na(combined_data[, tmp_var]), "PRS"], na.rm = TRUE), 3)
                median_val = round(median(combined_data[!is.na(combined_data[, tmp_var]), "PRS"], na.rm = TRUE), 3)
                sd_val = round(sd(combined_data[!is.na(combined_data[, tmp_var]), "PRS"], na.rm = TRUE), 3)
                min_val = round(min(combined_data[!is.na(combined_data[, tmp_var]), "PRS"], na.rm = TRUE), 3)
                max_val = round(max(combined_data[!is.na(combined_data[, tmp_var]), "PRS"], na.rm = TRUE), 3)
                iqr_val = paste(round(quantile(combined_data[!is.na(combined_data[, tmp_var]), "PRS"], probs = c(0.25, 0.75), na.rm = TRUE), 3), collapse = '-')
                # Write the stats to the file
                cat(paste0("PRS Statistics - phenotype: ", tmp_var, " ", suffix, "\n"), file = report_file, append = TRUE)
                cat("Mean\tMedian\tSD\tMin\tMax\tIQR\n", file = report_file, append = TRUE)
                cat(paste0(mean_val, "\t", median_val, "\t", sd_val, "\t", min_val, "\t", max_val, "\t", iqr_val), "\n", file = report_file, append = TRUE)
                cat("\n", file = report_file, append = TRUE)

                # Age if present
                if ('AGE' %in% toupper(colnames(combined_data))){
                    # Calculate mean, median, sd, min, max, and IQR of the age
                    mean_val = round(mean(combined_data[!is.na(combined_data[, tmp_var]), "AGE"], na.rm = TRUE), 3)
                    median_val = round(median(combined_data[!is.na(combined_data[, tmp_var]), "AGE"], na.rm = TRUE), 3)
                    sd_val = round(sd(combined_data[!is.na(combined_data[, tmp_var]), "AGE"], na.rm = TRUE), 3)
                    min_val = round(min(combined_data[!is.na(combined_data[, tmp_var]), "AGE"], na.rm = TRUE), 3)
                    max_val = round(max(combined_data[!is.na(combined_data[, tmp_var]), "AGE"], na.rm = TRUE), 3)
                    iqr_val = paste(round(quantile(combined_data[!is.na(combined_data[, tmp_var]), "AGE"], probs = c(0.25, 0.75), na.rm = TRUE), 3), collapse = '-')
                    # Write the stats to the file
                    cat(paste0("Age Statistics - phenotype: ", tmp_var, " ", suffix, "\n"), file = report_file, append = TRUE)
                    cat("Mean\tMedian\tSD\tMin\tMax\tIQR\n", file = report_file, append = TRUE)
                    cat(paste0(mean_val, "\t", median_val, "\t", sd_val, "\t", min_val, "\t", max_val, "\t", iqr_val), "\n", file = report_file, append = TRUE)
                    cat("\n", file = report_file, append = TRUE)
                }

                # If present, sex stratified as well
                if (length(tmp_sexes) > 0){
                    for (sex in tmp_sexes){
                        # Calculate mean, median, sd, min, max, and IQR for the different sexes
                        mean_val = round(mean(combined_data[!is.na(combined_data[, tmp_var]) & combined_data$SEX == sex, tmp_var], na.rm = TRUE), 3)
                        median_val = round(median(combined_data[!is.na(combined_data[, tmp_var]) & combined_data$SEX == sex, tmp_var], na.rm = TRUE), 3)
                        sd_val = round(sd(combined_data[!is.na(combined_data[, tmp_var]) & combined_data$SEX == sex, tmp_var], na.rm = TRUE), 3)
                        min_val = round(min(combined_data[!is.na(combined_data[, tmp_var]) & combined_data$SEX == sex, tmp_var], na.rm = TRUE), 3)
                        max_val = round(max(combined_data[!is.na(combined_data[, tmp_var]) & combined_data$SEX == sex, tmp_var], na.rm = TRUE), 3)
                        iqr_val = paste(round(quantile(combined_data[!is.na(combined_data[, tmp_var]) & combined_data$SEX == sex, tmp_var], probs = c(0.25, 0.75), na.rm = TRUE), 3), collapse = '-')
                        # Write the stats to the file
                        cat(paste0(tmp_var, " Statistics - phenotype: ", tmp_var, " ~ SEX: ", sex, " ", suffix, "\n"), file = report_file, append = TRUE)
                        cat("Mean\tMedian\tSD\tMin\tMax\tIQR\n", file = report_file, append = TRUE)
                        cat(paste0(mean_val, "\t", median_val, "\t", sd_val, "\t", min_val, "\t", max_val, "\t", iqr_val), "\n", file = report_file, append = TRUE)
                        cat("\n", file = report_file, append = TRUE)

                        # The related to the PRS
                        # Calculate mean, median, sd, min, max, and IQR of the variable
                        mean_val = round(mean(combined_data[!is.na(combined_data[, tmp_var]) & combined_data$SEX == sex, "PRS"], na.rm = TRUE), 3)
                        median_val = round(median(combined_data[!is.na(combined_data[, tmp_var]) & combined_data$SEX == sex, "PRS"], na.rm = TRUE), 3)
                        sd_val = round(sd(combined_data[!is.na(combined_data[, tmp_var]) & combined_data$SEX == sex, "PRS"], na.rm = TRUE), 3)
                        min_val = round(min(combined_data[!is.na(combined_data[, tmp_var]) & combined_data$SEX == sex, "PRS"], na.rm = TRUE), 3)
                        max_val = round(max(combined_data[!is.na(combined_data[, tmp_var]) & combined_data$SEX == sex, "PRS"], na.rm = TRUE), 3)
                        iqr_val = paste(round(quantile(combined_data[!is.na(combined_data[, tmp_var]) & combined_data$SEX == sex, "PRS"], probs = c(0.25, 0.75), na.rm = TRUE), 3), collapse = '-')
                        # Write the stats to the file
                        cat(paste0("PRS Statistics - phenotype: ", tmp_var, " ~ SEX: ", sex, " ", suffix, "\n"), file = report_file, append = TRUE)
                        cat("Mean\tMedian\tSD\tMin\tMax\tIQR\n", file = report_file, append = TRUE)
                        cat(paste0(mean_val, "\t", median_val, "\t", sd_val, "\t", min_val, "\t", max_val, "\t", iqr_val), "\n", file = report_file, append = TRUE)
                        cat("\n", file = report_file, append = TRUE)

                        # Age if present
                        if ('AGE' %in% toupper(colnames(combined_data))){
                            # Calculate mean, median, sd, min, max, and IQR of the age for the different sex
                            mean_val = round(mean(combined_data[!is.na(combined_data[, tmp_var]) & combined_data$SEX == sex, 'AGE'], na.rm = TRUE), 3)
                            median_val = round(median(combined_data[!is.na(combined_data[, tmp_var]) & combined_data$SEX == sex, 'AGE'], na.rm = TRUE), 3)
                            sd_val = round(sd(combined_data[!is.na(combined_data[, tmp_var]) & combined_data$SEX == sex, 'AGE'], na.rm = TRUE), 3)
                            min_val = round(min(combined_data[!is.na(combined_data[, tmp_var]) & combined_data$SEX == sex, 'AGE'], na.rm = TRUE), 3)
                            max_val = round(max(combined_data[!is.na(combined_data[, tmp_var]) & combined_data$SEX == sex, 'AGE'], na.rm = TRUE), 3)
                            iqr_val = paste(round(quantile(combined_data[!is.na(combined_data[, tmp_var]) & combined_data$SEX == sex, 'AGE'], probs = c(0.25, 0.75), na.rm = TRUE), 3), collapse = '-')
                            # Write the stats to the file
                            cat(paste0("Age Statistics - phenotype: ", tmp_var, " ~ SEX: ", sex, " ", suffix, "\n"), file = report_file, append = TRUE)
                            cat("Mean\tMedian\tSD\tMin\tMax\tIQR\n", file = report_file, append = TRUE)
                            cat(paste0(mean_val, "\t", median_val, "\t", sd_val, "\t", min_val, "\t", max_val, "\t", iqr_val), "\n", file = report_file, append = TRUE)
                            cat("\n", file = report_file, append = TRUE)
                        }
                    }
                }
            } else if (tmp_var_model == 'binomial'){
                # Number of individuals with phenotype
                    # Derive original variables and newer
                    values_df = data.frame(orig_values = c(str_split_fixed(str_split_fixed(tmp_var_mapping, ';', 2)[1], " ", 3)[1], str_split_fixed(str_split_fixed(tmp_var_mapping, ';', 2)[2], " ", 4)[2]), new_values = c(str_split_fixed(str_split_fixed(tmp_var_mapping, ';', 2)[1], " ", 3)[3], str_split_fixed(str_split_fixed(tmp_var_mapping, ';', 2)[2], " ", 4)[4]))
                    values_df = values_df[order(values_df$new_values),]
                    # Number of cases and controls
                    n_controls = nrow(combined_data[!is.na(combined_data[, tmp_var]) & combined_data[, tmp_var] == values_df$new_values[1], ])
                    n_cases = nrow(combined_data[!is.na(combined_data[, tmp_var]) & combined_data[, tmp_var] == values_df$new_values[2], ])
                    # Write the number of cases and controls to the file
                    cat(paste0("Number of controls (phenotype=", values_df$orig_values[1], "): ", n_controls, "\n"), file = report_file, append = TRUE)
                    cat(paste0("Number of cases (phenotype=", values_df$orig_values[2], "): ", n_cases, "\n"), file = report_file, append = TRUE)
                    cat("\n", file = report_file, append = TRUE)
        
                    # Number of cases and controls by sex
                    tmp_table_controls = table(combined_data$SEX[!is.na(combined_data[, tmp_var]) & combined_data[, tmp_var] == values_df$new_values[1]])
                    tmp_table_cases = table(combined_data$SEX[!is.na(combined_data[, tmp_var]) & combined_data[, tmp_var] == values_df$new_values[2]])
                    # Write the number of cases and controls by sex
                    cat(paste0("Number of controls (phenotype=", values_df$orig_values[1], ") by sex:\n"), file = report_file, append = TRUE)
                    cat(paste0(names(tmp_table_controls), ": ", as.numeric(tmp_table_controls), collapse = ', '), "\n\n", file = report_file, append = TRUE)
                    cat(paste0("Number of cases (phenotype=", values_df$orig_values[2], ") by sex:\n"), file = report_file, append = TRUE)
                    cat(paste0(names(tmp_table_cases), ": ", as.numeric(tmp_table_cases), collapse = ', '), "\n\n", file = report_file, append = TRUE)
                    cat("\n", file = report_file, append = TRUE)

                # The related to the PRS -- iterate over the phenotypes
                for (ph in values_df$new_values){
                    # get the original value
                    ph_orig = values_df$orig_values[which(values_df$new_values == ph)]
                    # Calculate mean, median, sd, min, max, and IQR of the variable
                    mean_val = round(mean(combined_data[!is.na(combined_data[, tmp_var]) & combined_data[, tmp_var] == ph, "PRS"], na.rm = TRUE), 3)
                    median_val = round(median(combined_data[!is.na(combined_data[, tmp_var]) & combined_data[, tmp_var] == ph, "PRS"], na.rm = TRUE), 3)
                    sd_val = round(sd(combined_data[!is.na(combined_data[, tmp_var]) & combined_data[, tmp_var] == ph, "PRS"], na.rm = TRUE), 3)
                    min_val = round(min(combined_data[!is.na(combined_data[, tmp_var]) & combined_data[, tmp_var] == ph, "PRS"], na.rm = TRUE), 3)
                    max_val = round(max(combined_data[!is.na(combined_data[, tmp_var]) & combined_data[, tmp_var] == ph, "PRS"], na.rm = TRUE), 3)
                    iqr_val = paste(round(quantile(combined_data[!is.na(combined_data[, tmp_var]) & combined_data[, tmp_var] == ph, "PRS"], probs = c(0.25, 0.75), na.rm = TRUE), 3), collapse = '-')
                    # Write the stats to the file
                    cat(paste0("PRS Statistics - phenotype: ", tmp_var, " ~ ", ph_orig, " ", suffix, "\n"), file = report_file, append = TRUE)
                    cat("Mean\tMedian\tSD\tMin\tMax\tIQR\n", file = report_file, append = TRUE)
                    cat(paste0(mean_val, "\t", median_val, "\t", sd_val, "\t", min_val, "\t", max_val, "\t", iqr_val), "\n", file = report_file, append = TRUE)
                    cat("\n", file = report_file, append = TRUE)

                    # Age if present
                    if ('AGE' %in% toupper(colnames(combined_data))){
                        # Calculate mean, median, sd, min, max, and IQR of the age
                        mean_val = round(mean(combined_data[!is.na(combined_data[, tmp_var]) & combined_data[, tmp_var] == ph, "AGE"], na.rm = TRUE), 3)
                        median_val = round(median(combined_data[!is.na(combined_data[, tmp_var]) & combined_data[, tmp_var] == ph, "AGE"], na.rm = TRUE), 3)
                        sd_val = round(sd(combined_data[!is.na(combined_data[, tmp_var]) & combined_data[, tmp_var] == ph, "AGE"], na.rm = TRUE), 3)
                        min_val = round(min(combined_data[!is.na(combined_data[, tmp_var]) & combined_data[, tmp_var] == ph, "AGE"], na.rm = TRUE), 3)
                        max_val = round(max(combined_data[!is.na(combined_data[, tmp_var]) & combined_data[, tmp_var] == ph, "AGE"], na.rm = TRUE), 3)
                        iqr_val = paste(round(quantile(combined_data[!is.na(combined_data[, tmp_var]) & combined_data[, tmp_var] == ph, "AGE"], probs = c(0.25, 0.75), na.rm = TRUE), 3), collapse = '-')
                        # Write the stats to the file
                        cat(paste0("Age Statistics - phenotype: ", tmp_var, " ~ ", ph_orig, " ", suffix, "\n"), file = report_file, append = TRUE)
                        cat("Mean\tMedian\tSD\tMin\tMax\tIQR\n", file = report_file, append = TRUE)
                        cat(paste0(mean_val, "\t", median_val, "\t", sd_val, "\t", min_val, "\t", max_val, "\t", iqr_val), "\n", file = report_file, append = TRUE)
                        cat("\n", file = report_file, append = TRUE)
                    }

                    # If present, sex stratified as well
                    if (length(tmp_sexes) > 0){
                        for (sex in tmp_sexes){
                            # The related to the PRS
                            # Calculate mean, median, sd, min, max, and IQR of the variable
                            mean_val = round(mean(combined_data[!is.na(combined_data[, tmp_var]) & combined_data$SEX == sex, "PRS"], na.rm = TRUE), 3)
                            median_val = round(median(combined_data[!is.na(combined_data[, tmp_var]) & combined_data$SEX == sex, "PRS"], na.rm = TRUE), 3)
                            sd_val = round(sd(combined_data[!is.na(combined_data[, tmp_var]) & combined_data$SEX == sex, "PRS"], na.rm = TRUE), 3)
                            min_val = round(min(combined_data[!is.na(combined_data[, tmp_var]) & combined_data$SEX == sex, "PRS"], na.rm = TRUE), 3)
                            max_val = round(max(combined_data[!is.na(combined_data[, tmp_var]) & combined_data$SEX == sex, "PRS"], na.rm = TRUE), 3)
                            iqr_val = paste(round(quantile(combined_data[!is.na(combined_data[, tmp_var]) & combined_data$SEX == sex, "PRS"], probs = c(0.25, 0.75), na.rm = TRUE), 3), collapse = '-')
                            # Write the stats to the file
                            cat(paste0("PRS Statistics - phenotype: ", tmp_var, " ~ SEX: ", sex, " ", suffix, "\n"), file = report_file, append = TRUE)
                            cat("Mean\tMedian\tSD\tMin\tMax\tIQR\n", file = report_file, append = TRUE)
                            cat(paste0(mean_val, "\t", median_val, "\t", sd_val, "\t", min_val, "\t", max_val, "\t", iqr_val), "\n", file = report_file, append = TRUE)
                            cat("\n", file = report_file, append = TRUE)

                            # Age if present
                            if ('AGE' %in% toupper(colnames(combined_data))){
                                # Calculate mean, median, sd, min, max, and IQR of the age for the different sex
                                mean_val = round(mean(combined_data[!is.na(combined_data[, tmp_var]) & combined_data$SEX == sex, 'AGE'], na.rm = TRUE), 3)
                                median_val = round(median(combined_data[!is.na(combined_data[, tmp_var]) & combined_data$SEX == sex, 'AGE'], na.rm = TRUE), 3)
                                sd_val = round(sd(combined_data[!is.na(combined_data[, tmp_var]) & combined_data$SEX == sex, 'AGE'], na.rm = TRUE), 3)
                                min_val = round(min(combined_data[!is.na(combined_data[, tmp_var]) & combined_data$SEX == sex, 'AGE'], na.rm = TRUE), 3)
                                max_val = round(max(combined_data[!is.na(combined_data[, tmp_var]) & combined_data$SEX == sex, 'AGE'], na.rm = TRUE), 3)
                                iqr_val = paste(round(quantile(combined_data[!is.na(combined_data[, tmp_var]) & combined_data$SEX == sex, 'AGE'], probs = c(0.25, 0.75), na.rm = TRUE), 3), collapse = '-')
                                # Write the stats to the file
                                cat(paste0("Age Statistics - phenotype: ", tmp_var, " ~ SEX: ", sex, " ", suffix, "\n"), file = report_file, append = TRUE)
                                cat("Mean\tMedian\tSD\tMin\tMax\tIQR\n", file = report_file, append = TRUE)
                                cat(paste0(mean_val, "\t", median_val, "\t", sd_val, "\t", min_val, "\t", max_val, "\t", iqr_val), "\n", file = report_file, append = TRUE)
                                cat("\n", file = report_file, append = TRUE)
                            }
                        }
                    }
                }
            }
            cat("========================================\n", file = report_file, append = TRUE)
            cat("\n\n", file = report_file, append = TRUE)
        }
    }
