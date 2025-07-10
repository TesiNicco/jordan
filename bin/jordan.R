#!/usr/bin/env Rscript
cat("\n** Jordan: a pipeline to make PRS and PRS analyses in R **\n")

# Libraries
    # Check if the required libraries are installed
    suppressWarnings({
        suppressMessages({
            tryCatch({
                library(argparse, quietly = TRUE)
                library(data.table, quietly = TRUE)
                library(stringr, quietly = TRUE)
                library(survival, quietly = TRUE)
                library(ggplot2, quietly = TRUE)
                library(survminer, quietly = TRUE)
                library(ggpubr, quietly = TRUE)
            }, error = function(e) {
                cat('**** Required packages are not installed. Installing now...\n')
                required_packages <- c("argparse", "data.table", "stringr", "survival", "ggplot2", "survminer", "ggpubr")
                installed_packages <- rownames(installed.packages())
                missing_packages <- setdiff(required_packages, installed_packages)
                if (length(missing_packages) > 0) {
                    install.packages(missing_packages, repos = "https://cloud.r-project.org/")
                }
                # Load the libraries again after installation
                library(argparse, quietly = TRUE)
                library(data.table, quietly = TRUE)
                library(stringr, quietly = TRUE)
                library(survival, quietly = TRUE)
                library(ggplot2, quietly = TRUE)
                library(survminer, quietly = TRUE)
                library(ggpubr, quietly = TRUE)
            })
        })
    })

# Parse arguments
    # Required arguments
        # Create parser
        parser <- ArgumentParser(description = "jordan: a pipeline to make PRS and PRS analyses in R")
        # Genotype file
        parser$add_argument("--genotype", help = "Path to the genotype file. If PLINK, do NOT provide file extension. If VCF, please provide the file extension. If multiple PLINK files should be used, input the path of the first file, and all PLINK files in the directory where the input files are will be used.")
        # SNPs file
        parser$add_argument("--snplist", help = "Path to the SNPs file. This file will define the SNPs to include in the PRS and the relative weights. By default, required column names are CHROM, POS, EFFECT_ALLELE, OTHER_ALLELE and BETA, (or OR)")
        # Output directory name
        parser$add_argument("--outname", help = "Path to outputs. A new directory can be specified. In that case, the directory will be created and the file with be written to the new directory.")
    # Flag-like arguments
        # Dosage
        parser$add_argument("--dosage", help="When present, it is assumed data to be imputed. The information is used to read genotypes in PLINK.", action = "store_true", default = FALSE)
        # Multiple input files
        parser$add_argument("--multiple", help="When present, multiple PLINK files will be used. Jordan will then look at all files with the same extension in the same input directory. CURRENTLY, ONLY PLINK FILES ARE SUPPORTED.", default = FALSE, action = "store_true")
        # Exclude APOE
        parser$add_argument("--exclude", help="When present, PRS will be also made excluding APOE SNPs. These by default are APOE e2 and APOE e4 variants.", default = FALSE, action = "store_true")
        # Direct effects
	    parser$add_argument("--directEffects", help="When present, PRS will be calculated following the --snplist file. By default, PRS is calculated based on the risk allele, i.e., always flip BETA or OR to risk (BETA>0, OR>1) and adjust the alleles accordingly.", default = FALSE, action = "store_true")
	    # Keep dosages
        parser$add_argument("--keepDosage", help="When present, dosages will be kept.", default = FALSE, action = "store_true")
        # Calculate frequency
        parser$add_argument("--freq", help="When present, allele frequencies will be calculated. This is useful to check the MAF of the SNPs in the input file.", default = FALSE, action = "store_true")
        # Sex-stratified analysis
        parser$add_argument("--sex-strata", help="When present, the association analyses will be conducted in males, females, and the combined group. Make sure that a variable SEX or sex is present in the association file.", default = FALSE, action = "store_true")
    # Optional arguments
        # Plot
        parser$add_argument("--plot", nargs = "?", const = "default", default = NULL, help = "If specified, saves a plot. Optionally takes 'exclude_NA' to exclude NAs from the plots (Default: include all).")        
        # Additional weight for the prs
        parser$add_argument("--addWeight", help="Additional weight to be applied on top of the BETA, for each SNP. This is useful when, for example, foing pathway-specific PRS. The SNP will be weighted by the BETA and by the optional weight (PRS_snp = SNP_dosage * SNP_beta * SNP_additionalWeight). To enable, insert the name of the column to be used. The column must be present in the snplist file.", default = FALSE)
        # Minor allele frequency filter
        parser$add_argument("--maf", help="When present, a filtering based on Minor Allele Frequency (MAF) will be applied. Usage: 0.01 for MAF>1%%.", default = FALSE)
        # Association file
        parser$add_argument("--assoc", nargs=2, metavar=c("mode", 'file'), help="When present, association testing will be conducted. By default, PRS association will be conducted (--assoc prs). Alternatively, single variant association can be conducted as well (--assoc single), or both (--assoc both). The models are pairwise logistic regression models (if categorical or binary phenotypes are detected), and linear regression models (if numerical phenotypes are detected). Please provide a file with labels to associate: the file should contain IID column with the sample names, and the variables to associate. Use argument --assoc-var to define variables to associate and --assoc-cov to define covariates.", default = c("prs", FALSE))
        # Association variables
        parser$add_argument("--assoc-var", help="Comma-separated list of the variables to associate.", default = FALSE)
        # Association covariates
        parser$add_argument("--assoc-cov", help="Comma-separated list of the covariates to include in the association testing.", default = FALSE)
        # Survival analysis
        parser$add_argument("--survival-var", help="When present, the comma-separated list of variables will be subject to survival analysis. This will be a Cox regression analysis, with the same covariates as in the association testing. The variables should be in the format: 'time, event', where time is the time to event and event is the event status (1 for event, 0 for censored). For example, event should be defined as {variable}_EVENT. If event is not present, will assume the event is the same for all individuals.", default = FALSE)
        # Decile based analysis, accepts a comma-separated list of variables
        parser$add_argument("--tiles", help="When present, the individuals will be divided into n-tiles based on the PRS and the user preference. The tiles can be calculated on a specific subset of individuals, and applied to all other individuals. The variables should be in the format: 'n_tiles;reference_group_variable;reference_group_name, n_tiles;reference_group_variable;reference_group_name, ...'. See example in the README and example scripts.", default = FALSE)
        # Split individuals based on defined thresholds
        parser$add_argument("--split", help="When present, the individuals will be split into groups based on the defined thresholds. The thresholds should be provided as a comma-separated list of values. The values should be numeric and in ascending order. For example, --split 'reference_group_variable;threshold1-threshold2-threshold3, ...'.", default = FALSE)
        # Association of the split/tiles
        parser$add_argument("--assoc-split-tiles", help="When present, the association testing will be conducted on the split/tiles. The association will be conducted on the split/tiles defined by the --split or --tiles arguments.", default = FALSE, action = "store_true")
    # Read arguments
        args <- parser$parse_args()
        genotype_file <- args$genotype
        snps_file <- args$snplist
        outfile <- args$outname
        dosage <- args$dosage
        plt <- args$plot
        maf <- args$maf
        multiple = args$multiple
        excludeAPOE = args$exclude
	    fliprisk = args$directEffects
	    keepDos = args$keepDosage
        addWeight = args$addWeight
        freq = args$freq
        assoc_mode = args$assoc[1]
        assoc_file = args$assoc[2]
        assoc_var = args$assoc_var
        assoc_cov = args$assoc_cov
        assoc_survival = args$survival_var
        sex_strata = args$sex_strata
        tiles_prs = args$tiles
        split_info = args$split
        assoc_split = args$assoc_split_tiles

    # Print arguments on screen
        cat("\nSystem: ", system_config)
        cat("\nGenotype file: ", genotype_file)
        cat("\nMultiple files: ", multiple)
        cat("\nSNPs file: ", snps_file)
        cat("\nOutput file: ", outfile)
        cat("\nDosage: ", dosage)
        cat("\nMAF: ", maf)
        cat("\nWith and Without APOE: ", excludeAPOE)
	    cat("\nUse direct effects (Risk and Protective): ", fliprisk)
	    cat("\nKeep dosages: ", keepDos)
        cat("\nAdditional weight: ", addWeight)
        cat("\nCalculate frequency: ", freq)
        cat("\nAssociation testing: ", assoc_file)
        if (assoc_file != FALSE){
            cat("\nAssociation mode: ", assoc_mode)
        }
        cat("\nAssociation variables: ", assoc_var)
        cat("\nAssociation covariates: ", assoc_cov)
        cat("\nAssociation survival: ", assoc_survival)
        cat("\nSex-stratified analysis: ", sex_strata)
        cat("\nTile-based analysis: ", tiles_prs)
        cat("\nSplit individuals: ", split_info)
        cat("\nAssociation of the split/tiles: ", assoc_split)
        cat("\nPlot: ", plt, '\n\n')

# Import functions from jordan_functions.R
    # Derive directory of the script
    script_path <- dirname(sub("^--file=", "", args[grep("^--file=", args)]))
    # Load functions    
    source(file.path(script_path, "jordan_functions.R"))
    # Detect system and adapt plink executables
    system_info = detect_system(script_path)
    plink_path = system_info[[1]]
    plink2_path = system_info[[2]]
    system_config = system_info[[3]]
    

# Check inputs
    # Check output directory
        outdir = checkOutputFile(outfile)

    # Check input genotype file
        res = checkGenoFile(genotype_file, outdir, dosage, multiple, script_path, plink_path, plink2_path)
        genotype_path = res[[1]]
        genotype_type = res[[2]]

    # Check input snplist
        snps_data = checkInputSNPs(snps_file, addWeight)

    # Check maf
        if (maf != FALSE){
            maf = suppressWarnings(as.numeric(maf))
            if (is.na(maf)){ stop("** The MAF value is not numeric. It should be, for example, 0.01 for MAF=1%.\n\n", call. = FALSE) }
        }
    
    # Check association file
        if (assoc_file != FALSE){
            cat('\n** Association testing requested. Checking parameters.\n')
            # Check association mode
            assoc_mode = checkAssocMode(assoc_mode)
            # Check association file
            assoc_info = checkAssocFile(assoc_file, assoc_var, assoc_cov, assoc_survival, sex_strata)
        } else {
            assoc_info = FALSE
        }

    # Final check before running
        cat('\n** Inputs are valid. Starting the script.\n\n')
            
        # Add log file with run info to the output folder
        log = writeLog(outdir, genotype_file, snps_file, outfile, dosage, plt, maf, multiple, excludeAPOE, fliprisk, keepDos, addWeight, freq, assoc_file, assoc_var, assoc_cov, assoc_survival, sex_strata)
            
        # Calculate PRS
        res = makePRS(outdir, genotype_path, snps_data, genotype_type, multiple, excludeAPOE, maf, fliprisk, keepDos, addWeight, freq, assoc_file, assoc_info, script_path, sex_strata, plink_path, plink2_path)

        # Check if tile-based analysis is requested
        if (tiles_prs != FALSE){
            tryCatch({
                cat('\n**** Tile-based analysis requested.\n')
                # Check tiles_prs format
                cat('****** Checking tiles preferences and making tiles.\n')
                tiles_prs_df = checkTiles(tiles_prs, assoc_info, sex_strata)
                # check if APOE is considered
                if (excludeAPOE){
                    res_prs = res[[1]][[1]]
                    res_noapoe = res[[2]][[1]]
                    # Make tiles for PRS with APOE
                    res_prs_with_tiles = makeTiles(res_prs, tiles_prs_df, assoc_info)
                    res[[1]][[1]] = res_prs_with_tiles
                    # Make tiles for PRS without APOE
                    res_noapoe_with_tiles = makeTiles(res_noapoe, tiles_prs_df, assoc_info)
                    res[[2]][[1]] = res_noapoe_with_tiles
                } else {
                    res_prs = res[[1]][[1]]
                    res_prs_with_tiles = makeTiles(res_prs, tiles_prs_df, assoc_info)
                    res[[1]][[1]] = res_prs_with_tiles
                }
            }, error = function(e) {
                cat('\n**** Tile-based analysis requested, but failed. Make sure you defined all arguments correctly.\n')
                tiles_prs_df = data.frame(n_tiles = c(), variable = c(), group = c())
            })
        } else {
            tiles_prs_df = data.frame(n_tiles = c(), variable = c(), group = c())
        }

        # Check if split is requested
        if (split_info != FALSE){
            tryCatch({
                cat('\n**** Split individuals requested.\n')
                # Check split_info format
                cat('****** Checking split preferences and making splits.\n')
                split_info_df = checkSplit(split_info, assoc_info, sex_strata)
                # check if APOE is considered
                if (excludeAPOE){
                    # Make split for PRS with APOE
                    res_prs = res[[1]][[1]]
                    res_noapoe = res[[2]][[1]]
                    res_prs_with_split = makeSplit(res_prs, split_info_df, assoc_info)
                    res_noapoe_with_split = makeSplit(res_noapoe, split_info_df, assoc_info)
                    res[[1]][[1]] = res_prs_with_split
                    res[[2]][[1]] = res_noapoe_with_split
                } else {
                    res_prs = res[[1]][[1]]
                    # Make split for PRS
                    res_prs_with_split = makeSplit(res_prs, split_info_df, assoc_info)
                    res[[1]][[1]] = res_prs_with_split
                }
            }, error = function(e) {
                cat('\n**** Split individuals requested, but failed. Make sure you defined all arguments correctly.\n')
                split_info_df = data.frame(n_split = c(), variable = c(), thresholds = c())
            })
        } else {
            split_info_df = data.frame(n_split = c(), variable = c(), thresholds = c())
        }

        # Write outputs
        cat('\n**** Writing outputs.\n')
        if (excludeAPOE){
            res_apoe = res[[1]]
            res_noapoe = res[[2]]
            dosages = res[[3]]
            all_freq = res[[4]]
            # Write PRS with and without APOE
            write_prs_outputs(res_apoe[[1]], res_apoe[[2]], snps_data, "", outdir, assoc_info)
            write_prs_outputs(res_noapoe[[1]], res_noapoe[[2]], snps_data, "_noAPOE", outdir, assoc_info)
        } else {
            res_prs = res[[1]]
            dosages = res[[2]]
            all_freq = res[[3]]
            write_prs_outputs(res_prs[[1]], res_prs[[2]], snps_data, "", outdir, assoc_info)
        }

        # Association testing
        if (assoc_file != FALSE){
            cat('\n**** Association testing.\n')
            if (excludeAPOE){
                res_apoe = res[[1]]
                res_noapoe = res[[2]]
                assoc_test(res_apoe[[1]], assoc_info, outdir, "", assoc_mode, dosages, sex_strata, plt)
                assoc_test(res_noapoe[[1]], assoc_info, outdir, "_noAPOE", 'prs', dosages, sex_strata, plt)
            } else {
                assoc_test(res_prs[[1]], assoc_info, outdir, "", assoc_mode, dosages, sex_strata, plt)
            }
        }

        # Association of the split/tiles
        if (assoc_split){
            tryCatch({
                cat('\n**** Association of the PRS with the split/tiles.\n')
                if (excludeAPOE){
                    res_apoe = res[[1]]
                    res_noapoe = res[[2]]
                    test_info = assoc_split_tiles_test(res_apoe[[1]], assoc_info, outdir, "", split_info_df, tiles_prs_df)
                    test_info_noapoe = assoc_split_tiles_test(res_noapoe[[1]], assoc_info, outdir, "_noAPOE", split_info_df, tiles_prs_df)
                } else {
                    res_prs = res[[1]]
                    test_info = assoc_split_tiles_test(res_prs[[1]], assoc_info, outdir, "", split_info_df, tiles_prs_df)
                }
            }, error = function(e){
                cat('\n**** Association of the PRS with the split/tiles requested, but failed. Make sure you defined all arguments correctly.\n')
            })
        }

        # Clean temporary data
        cat('\n**** Cleaning.\n')
        system(paste0('rm ', outdir, '/tmp*'), ignore.stdout = TRUE, ignore.stderr = TRUE)
        system(paste0('rm ', outdir, '/dosage*'), ignore.stdout = TRUE, ignore.stderr = TRUE)
            
        # Check if plot needs to be done
        if (!is.null(plt)) {
            if (plt == "exclude_NA") {
                cat("**** Making plots and excluding NAs.\n")
                na_beh = "exclude"
            } else {
                cat("**** Making plots including all samples.\n")
                na_beh = "include"
            }
            if (excludeAPOE){
                res_apoe = res[[1]]
                res_noapoe = res[[2]]
                all_freq = res[[4]]
                # Make plots
                makePlot(res_apoe[[1]], res_apoe[[2]], "", outdir, snps_data, assoc_info, all_freq, freq, assoc_file, sex_strata, tiles_prs_df, split_info_df, na_beh)
                makePlot(res_noapoe[[1]], res_noapoe[[2]], "_noAPOE", outdir, snps_data, assoc_info, all_freq, freq, assoc_file, sex_strata, tiles_prs_df, split_info_df, na_beh)
            } else {
                res_prs = res[[1]]
                all_freq = res[[3]]
                # Make plots
                makePlot(res_prs[[1]], res_prs[[2]], "", outdir, snps_data, assoc_info, all_freq, freq, assoc_file, sex_strata, tiles_prs_df, split_info_df, na_beh)
            }
        }

        # end message
        cat('\n\n** Analysis over. Ciao! \n\n')
    