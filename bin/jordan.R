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
        cat("\nPlot: ", plt, '\n\n')
    
# Check inputs
    # Check output directory
        res = checkOutputFile(outfile)
        run1 = res[[1]]
        outdir = res[[2]]

    # Check input genotype file
        res = checkGenoFile(genotype_file, outdir, dosage, multiple)
        run2 = res[[1]]
        genotype_path = res[[2]]
        genotype_type = res[[3]]

    # Check input snplist
        res = checkInputSNPs(snps_file, addWeight)
        run3 = res[[1]]
        snps_data = res[[2]]

    # Check maf
        if (maf != FALSE){
            tryCatch({
                maf = as.numeric(maf)
                run4 = TRUE
            }, error = function(e) {
                print(paste("** The MAF value is not numeric. It should be, for example, 0.01 for MAF=1%.\n", conditionMessage(e)))
                run4 = FALSE
            })  
        } else {
            maf = FALSE
            run4 = TRUE
        }

    # Final check before running
        if (run1 == FALSE | run2 == FALSE | run3 == FALSE | run4 == FALSE){
            stop("** Inputs are not valid. Check above messages and try again.\n\n")
        } else {
            cat('** Inputs are valid. Starting the script.\n\n')
            
            # Add log file with run info to the output folder
            log = writeLog(outdir, genotype_file, snps_file, outfile, dosage, plt, maf, multiple, excludeAPOE, fliprisk, keepDos, addWeight, freq)
            
            # Calculate PRS
            res = makePRS(outdir, genotype_path, snps_data, genotype_type, multiple, excludeAPOE, maf, fliprisk, keepDos, addWeight, freq)

            # Write outputs
            cat('**** Writing outputs.\n')
            if (excludeAPOE == TRUE){
                # Split results
                res_apoe = res[[1]]
                res_noapoe = res[[2]]
                # Get PRS dataframes and snp information
                prs_df_apoe = res_apoe[[1]]
                included_snps_apoe = res_apoe[[2]]
                prs_df_noapoe = res_noapoe[[1]]
                included_snps_noapoe = res_noapoe[[2]]
                # Write prs output
                write.table(prs_df_apoe, paste0(outdir, '/PRS_table.txt'), quote=F, row.names=F, sep="\t")
                write.table(prs_df_noapoe, paste0(outdir, '/PRS_table_noAPOE.txt'), quote=F, row.names=F, sep="\t")
                # Fill included snps file with those that are excluded as well
                included_snps_apoe$ID = paste(included_snps_apoe$CHROM, included_snps_apoe$POS, sep=":")
                snps_data$ID = paste(stringr::str_replace_all(snps_data$CHROM, 'chr', ''), snps_data$POS, sep=":")
                missing = snps_data[which(!(snps_data$ID %in% included_snps_apoe$ID)),]
                if (nrow(missing) >0){ included_snps_apoe = rbind(included_snps_apoe, data.frame(SNP = NA, BETA = missing$BETA, ALLELE = missing$EFFECT_ALLELE, OTHER_ALLELE = missing$OTHER_ALLELE, TYPE = 'Excluded', POS = missing$POS, CHROM = missing$CHROM, ID = missing$ID)) }
                # Write these outputs
                write.table(included_snps_apoe, paste0(outdir, '/SNPs_included_PRS.txt'), quote=F, row.names=F, sep="\t")
                # Same for set without APOE
                included_snps_noapoe$ID = paste(included_snps_noapoe$CHROM, included_snps_noapoe$POS, sep=":")
                snps_data$ID = paste(stringr::str_replace_all(snps_data$CHROM, 'chr', ''), snps_data$POS, sep=":")
                missing = snps_data[which(!(snps_data$ID %in% included_snps_noapoe$ID)),]
                if (nrow(missing) >0){ included_snps_noapoe = rbind(included_snps_noapoe, data.frame(SNP = NA, BETA = missing$BETA, ALLELE = missing$EFFECT_ALLELE, OTHER_ALLELE = missing$OTHER_ALLELE, TYPE = 'Excluded', POS = missing$POS, CHROM = missing$CHROM, ID = missing$ID)) }
                # Write these outputs
                write.table(included_snps_noapoe, paste0(outdir, '/SNPs_included_PRS_noAPOE.txt'), quote=F, row.names=F, sep="\t")
            } else {
                # Get PRS dataframe and snp information
                prs_df = res[[1]]
                included_snps = res[[2]]
                # Write prs output
                write.table(prs_df, paste0(outdir, '/PRS_table.txt'), quote=F, row.names=F, sep="\t")
                # Fill included snps file with those that are excluded as well
                included_snps$ID = paste(included_snps$CHROM, included_snps$POS, sep=":")
                snps_data$ID = paste(stringr::str_replace_all(snps_data$CHROM, 'chr', ''), snps_data$POS, sep=":")
                missing = snps_data[which(!(snps_data$ID %in% included_snps$ID)),]
                if (nrow(missing) >0){ included_snps = rbind(included_snps, data.frame(SNP = NA, BETA = missing$BETA, ALLELE = missing$EFFECT_ALLELE, OTHER_ALLELE = missing$OTHER_ALLELE, TYPE = 'Excluded', POS = missing$POS, CHROM = missing$CHROM, ID = missing$ID)) }
                # Write these outputs
                write.table(included_snps, paste0(outdir, '/SNPs_included_PRS.txt'), quote=F, row.names=F, sep="\t")
            }

            # Clean temporary data
            cat('**** Cleaning.\n')
            system(paste0('rm ', outdir, '/tmp*'), ignore.stdout = TRUE, ignore.stderr = TRUE)
            system(paste0('rm ', outdir, '/dosage*'), ignore.stdout = TRUE, ignore.stderr = TRUE)
            
            # Check if plot needs to be done
            if (plt == TRUE){
                makePlot(prs_df, included_snps, outdir)
            }
            # end message
            cat('\n\n** Analysis over. Ciao! \n\n')
        }
