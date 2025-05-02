# The following are example scripts for jordan. To replicate, move to jordan/example_data directory

# Calculate PRS using a single PLINK data file (most basic use case)
jordan.R --genotype example_data_plink --snplist AD_snps.txt --outname prs

# Calculate PRS using a multiple PLINK data file (eg. chr1_plink, chr2_plink, etc.) [not warking as no individual chromosomes are provided]
jordan.R --genotype example_data_plink_chr1 --multiple --snplist AD_snps.txt --outname prs

# Calculate PRS using a single VCF file (also working with gzipped VCF files)
jordan.R --genotype example_data.[vcf/vcf.gz/bcf/bcf.gz] --snplist AD_snps.txt --outname prs_vcf

# Calculate PRS using a single PLINK data, defining dosages, and MAF>1%
jordan.R --genotype example_data_plink --dosage --snplist AD_snps.txt --maf 0.01 --outname prs_dos_maf

# Calculate PRS using a single PLINK data, defining dosages, MAF>1%, request frequencies and single-variant dosages
jordan.R --genotype example_data_plink --dosage --snplist AD_snps.txt --maf 0.01 --freq --keepDosage --outname prs_dos_maf_freq

# Calculate PRS using a single PLINK data, defining dosages, MAF>1%, request frequencies and single-variant dosages. Perform association analysis of AD status with PRS, and SEX as covariate
jordan.R --genotype example_data_plink --dosage --snplist AD_snps.txt --maf 0.01 --freq --keepDosage --outname prs_dos_maf_freq_assoc --assoc example_assoc_file.txt --assoc-var AD --assoc-cov SEX

# Calculate PRS using a single PLINK data, defining dosages, MAF>1%, request frequencies and single-variant dosages. Perform association analysis of AD and MMSE status with PRS, and SEX, PC1 and PC2 as covariates
jordan.R --genotype example_data_plink --dosage --snplist AD_snps.txt --maf 0.01 --freq --keepDosage --outname prs_dos_maf_freq_assocMulti --assoc example_assoc_file.txt --assoc-var AD,MMSE --assoc-cov SEX,PC1,PC2



