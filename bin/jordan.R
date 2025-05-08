#!/usr/bin/Rscript

# Libraries
    library(argparse)
    library(data.table)
    library(stringr)
    library(ggplot2)
    args <- commandArgs(trailingOnly = FALSE)

# Functions: import functions from jordan_functions.R
    # Derive directory of the script
    script_path <- dirname(sub("^--file=", "", args[grep("^--file=", args)]))
    # Load functions    
    source(file.path(script_path, "jordan_functions.R"))
    
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
        # Plot
        parser$add_argument("--plot", help="When present, density plot is saved.", default = FALSE, action = "store_true")
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
    # Optional arguments
        parser$add_argument("--addWeight", help="Additional weight to be applied on top of the BETA, for each SNP. This is useful when, for example, foing pathway-specific PRS. The SNP will be weighted by the BETA and by the optional weight (PRS_snp = SNP_dosage * SNP_beta * SNP_additionalWeight). To enable, insert the name of the column to be used. The column must be present in the snplist file.", default = FALSE)
        # Minor allele frequency filter
        parser$add_argument("--maf", help="When present, a filtering based on Minor Allele Frequency (MAF) will be applied. Usage: 0.01 for MAF>1%%.", default = FALSE)
        # Association file
        parser$add_argument("--assoc", nargs=2, metavar=c("mode", 'file'), help="When present, association testing will be conducted. By default, PRS association will be conducted (--assoc prs). Alternatively, single variant association can be conducted as well (--assoc single), or both (--assoc both). The models are pairwise logistic regression models (if categorical or binary phenotypes are detected), and linear regression models (if numerical phenotypes are detected). Please provide a file with labels to associate: the file should contain IID column with the sample names, and the variables to associate. Use argument --assoc-var to define variables to associate and --assoc-cov to define covariates.", default = c("prs", FALSE))
        # Association variables
        parser$add_argument("--assoc-var", help="Comma-separated list of the variables to associate.", default = FALSE)
        # Association covariates
        parser$add_argument("--assoc-cov", help="Comma-separated list of the covariates to include in the association testing.", default = FALSE)

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

    # Print arguments on screen
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
        cat("\nPlot: ", plt, '\n\n')
    
# Check inputs
    # Check output directory
        outdir = checkOutputFile(outfile)

    # Check input genotype file
        res = checkGenoFile(genotype_file, outdir, dosage, multiple)
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
            assoc_info = checkAssocFile(assoc_file, assoc_var, assoc_cov)
        } else {
            assoc_info = FALSE
        }

    # Final check before running
        cat('\n** Inputs are valid. Starting the script.\n\n')
            
        # Add log file with run info to the output folder
        log = writeLog(outdir, genotype_file, snps_file, outfile, dosage, plt, maf, multiple, excludeAPOE, fliprisk, keepDos, addWeight, freq, assoc_file, assoc_var, assoc_cov)
            
        # Calculate PRS
        res = makePRS(outdir, genotype_path, snps_data, genotype_type, multiple, excludeAPOE, maf, fliprisk, keepDos, addWeight, freq, assoc_file, assoc_info)

        # Write outputs
        cat('\n**** Writing outputs.\n')
        if (excludeAPOE){
            res_apoe = res[[1]]
            res_noapoe = res[[2]]
            dosages = res[[3]]
            # Write PRS with and without APOE
            write_prs_outputs(res_apoe[[1]], res_apoe[[2]], snps_data, "", outdir)
            write_prs_outputs(res_noapoe[[1]], res_noapoe[[2]], snps_data, "_noAPOE", outdir)
        } else {
            res_prs = res[[1]]
            dosages = res[[2]]
            write_prs_outputs(res_prs[[1]], res_prs[[2]], snps_data, "", outdir)
        }

        # Association testing
        if (assoc_file != FALSE){
            cat('\n**** Association testing.\n')
            if (excludeAPOE){
                res_apoe = res[[1]]
                res_noapoe = res[[2]]
                assoc_test(res_apoe[[1]], assoc_info, outdir, "", assoc_mode, dosages)
                assoc_test(res_noapoe[[1]], assoc_info, outdir, "_noAPOE", 'prs', dosages)
            } else {
                assoc_test(res_prs[[1]], assoc_info, outdir, "", assoc_mode, dosages)
            }
        }

        # Clean temporary data
        cat('\n**** Cleaning.\n')
        system(paste0('rm ', outdir, '/tmp*'), ignore.stdout = TRUE, ignore.stderr = TRUE)
        system(paste0('rm ', outdir, '/dosage*'), ignore.stdout = TRUE, ignore.stderr = TRUE)
            
        # Check if plot needs to be done
        if (plt == TRUE){
            cat('**** Making plots.\n')
            if (excludeAPOE){
                res_apoe = res[[1]]
                res_noapoe = res[[2]]
                makePlot(res_apoe[[1]], res_apoe[[2]], "", outdir, snps_data)
                makePlot(res_noapoe[[1]], res_noapoe[[2]], "_noAPOE", outdir, snps_data)
            } else {
                makePlot(res[[1]], res[[2]], "", outdir, snps_data)
            }
        }
        
        # end message
        cat('\n\n** Analysis over. Ciao! \n\n')
    