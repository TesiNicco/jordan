# jordan
R script to calculate Polygenic Risk Score from genotype data, and eventually plot them.

## Install
Clone the repository locally:
`git clone https://github.com/TesiNicco/jordan.git`  
You may need to make the main file executable by typing:
```console
chmod +x ./jordan/bin/jordan.R
```

Make sure you have the required R packages in place:  
`argparse`  
`data.table`  
`stringr`  
`ggplot2`  
In addition to R packages, you need to have [https://www.cog-genomics.org/plink/2.0/](PLINK2) and [https://www.cog-genomics.org/plink/1.9/](PLINK) installed in your system.  
If all this is OK, then you should be able to run the script.  

## How to use
By running:  
`./bin/jordan.R -h`  
will display `jordan` options. To run, `jordan` requires:  
- genotype data (PLINK or VCF format)  
- snplist (list of SNPs to use for the PRS, including CHROM, POS, EFFECT_ALLELE, OTHER_ALLELE, and BETA (or OR)  
- outname (the name of the output folder)  
- isdosage: binary flag whether data is genotyped or imputed (TRUE/FALSE)  
- plot: binary flag to plot densities of not (TRUE/FALSE)  
- multiple: whether multiple input (PLINK) files should be used  
- exclude: whether PRS including and excluding APOE variants (e2 and e4 alleles) should be made  
- maf: minor allele frequency (MAF) threshold to include variants in the PRS. It should be, for example, 0.01 for MAF>1%  

## Example usage
Assuming a single PLINK file storing genotype data as dosages, no plots, PRS with and without APOE variants, with MAF>5%:  
`Rscript ./bin/jordan.R --genotype genotype/example_data_plink --snplist Snps_interest.txt --outname test_output --isdosage TRUE --exclude TRUE --maf 0.05`  
Assuming multiple PLINK files with genotype data as dosages, no plots, PRS with and without APOE variants, with MAF>5%:  
`Rscript ./bin/jordan.R --genotype genotype_data_path/chr1.example_data_plink --snplist Snps_interest.txt --outname test_output --isdosage TRUE --exclude TRUE --multiple TRUE --maf 0.05`  
Assuming a single VCF file with genotype data as dosages, no plots, PRS with and without APOE variants, with MAF>5%:  
`Rscript ./bin/jordan.R --genotype genotype_data_path/example_data.vcf --snplist Snps_interest.txt --outname test_output --isdosage TRUE --exclude TRUE --maf 0.05`  

## Test
To test everything is allright, you can use the following script alongside the example datasets provided in the repository. The example data will make a PRS from 85 SNPs associated with Alzheimer's Disease in individuals from the 1000Genome Project:  
`./bin/jordan.R --genotype example_data/example_data_plink --snplist example_data/AD_snps.txt --outname test_output --isdosage TRUE --plot TRUE`  
This should run successfully, and will display results in the output directory specified:  
- PRS_density.pdf: density of the PRS across individuals  
- PRS_SNPs.pdf: SNPs included/excluded from the PRS, and relative effect sizes  
- PRS_table.txt: tab-delimited table with samples in rows and PRS in column  
- snpsInterest.txt: tab-delimited table with information about included/excluded SNPs.  

## Contact
For comments, feedbacks, or question, feel free to reach me at n.tesi@amsterdamumc.nl or open an issue.



