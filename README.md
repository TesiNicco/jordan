# jordan
R script to calculate Polygenic Risk Score from genotype data, do PRS as well as single-variant associations, and plot.

## Install
Clone the repository locally:
```console
git clone https://github.com/TesiNicco/jordan.git
```

You may need to make the main file executable by typing:
```console
chmod +x ./jordan/bin/jordan.R
```

Make sure you have these required `R` packages in place: `argparse`, `data.table`, `stringr`, `ggplot2`.
In case these are not installed, you can install with the following code:
```console
install.packages(c("argparse", "data.table", "stringr", "ggplot2"))
```
In addition to R packages, you need to have [PLINK2](https://www.cog-genomics.org/plink/2.0/) and [PLINK](https://www.cog-genomics.org/plink/1.9/) installed in your system.  
If all this is OK, then you should be able to run the script. These are also provided in this package. 

## Application
We provide also an application wrapper around `jordan.R`, which allows the user to run the script

## How to use
By running:  
```console
./bin/jordan.R -h
```
will display the help message. `jordan` parameters are:  
- `--genotype` (**Mandatory**): path to genotype data in PLINK/VCF/BCF format. For PLINK data, <u>do not</u> specify file extension. For VCF/BCF data, <u>do</u> specify file extension.  
- `--snplist` (**Mandatory**): path to file containing the list of SNPs to use for the PRS. The file must include at least CHROM, POS, EFFECT_ALLELE, OTHER_ALLELE, and BETA (or OR).  
- `--outname` (**Mandatory**): path to the output folder. If a directory exists, output files will be placed there, otherwise, the directory will be created.  
- `--dosage` (*Optional, flag*): flag to indicate whether data is genotyped or imputed.  
- `--plot` (*Optional, flag*): flag to indicate whether plots should be drawn.  
- `--multiple` (*Optional, flag*): flag to indicate whether multiple input (PLINK) files should be used. In this case, all files with the same extension (PLINK files) in the input data directory will be used.  
- `--exclude` (*Optional, flag*): flag to indicate whether PRS including and excluding APOE variants (e2 and e4 alleles) should be made.  
- `--directEffects` (*Optional, flag*): by default, PRS is calculated with respect to the risk allele (BETA>0, OR>1). When this flag is indicated, PRS will be calculated with respect to the EFFECT_ALLELE as specified in the --snplist file.  
- `--keepDosage` (*Optional, flag*): flag to indicate whether the raw dosages of the variants specified in `--snplist` should be outputted.  
- `--freq` (*Optional, flag*): flag to indicate whether allele frequencies of the variants specified in `--snplist` should be outputted.  
- `--addWeight` (*Optional, default: False*): whether an additional weight (on top of BETA/OR) should be included in the PRS. If desired, add here the column name of the input `--snplist` file that reports the additional weight to use. See *PRS calculation* section below for more information on how this additional weight is used.  
- `--maf` (*Optional, default: False*): minor allele frequency (MAF) threshold to include variants in the PRS. It should be, for example, 0.01 for MAF>1%.  
- `--assoc` (*Optional, default: True*): when present, association analysis will be performed. Association can be done at the PRS level (`--assoc prs`), at the single variant level (`--assoc single`), or both (`--assoc both`). When association is requested, a file including the phenotypes of the individuals to compare should be provided. Also, you will need to specify the `--assoc-var` and `--assoc-cov` parameters.  
- `--assoc-var` (*Optional, default: True*): when association is requested, please specify here a comma-separated list of the outcome variable(s) you want to associate with PRS or single variants.  
- `--assoc-cov` (*Optional, default: True*): when association is requested, please specify here a comma-separated list of the covariate(s) you want to include in the association model with PRS or single variants.  

## PRS and weights
There can be different types of PRS: (i) *unweighted*, (ii) *weighted*, and (iii) *multiple-weighted*. *Unweighted* PRS are simply the sum of trait-associated alleles, with all variants having the same weight. To do this in jordan, simply set the BETA to 1 in the `--snplist` file to all variants. *Weighted* PRS are the most commonly used PRS, defined as the weighted sum of trait-associated alleles, weighted by an effect size, typically originating from a GWAS study. In `jordan`, PRS calculation is implemented as:  
$PRS_{sample} = \sum_{snp}^{SNPs} \alpha_{snp} \cdot \beta_{snp}$  
Finally, *multiple-weighted* PRS can be handy for calculating pathway-specific PRS or cell-specific PRS. In these cases, the additional weight (on top of the BETA, quantifying the effect on the trait of interest) will quantify the effect of each variant on a specific pathway, or cell-type. When an additional weight is selected with `--addWeight` option, for example for the calculation of pathway-specific PRS, then the formula will adapt accordingly to:  
$PRS_{sample} = \sum_{snp}^{SNPs} \alpha_{snp} \cdot \beta_{snp} \cdot w_{snp}$  
We previously used [PRS](https://alz-journals.onlinelibrary.wiley.com/doi/epdf/10.1002/alz.13810) and [pathway-PRS](https://pmc.ncbi.nlm.nih.gov/articles/PMC7524800/#:~:text=Immune%20response%20and%20endocytosis%20pathways,resilience%20against%20Alzheimer's%20disease%20%2D%20PMC) for the associations with Alzheimer's Disease previously.

## Example usage
The `example_scripts.sh` file in `example_data` folder reports multiple example commands. Here below it is reported its content:
Calculate PRS using a single PLINK data file (most basic use case)
```console
jordan.R --genotype example_data_plink --snplist AD_snps.txt --outname prs
```  
Calculate PRS using a multiple PLINK data file (eg. chr1_plink, chr2_plink, etc.) [not warking as no individual chromosomes are provided]
```console
jordan.R --genotype example_data_plink_chr1 --multiple --snplist AD_snps.txt --outname prs
```  
Calculate PRS using a single VCF file (also working with gzipped VCF files)
```console
jordan.R --genotype example_data.[vcf/vcf.gz/bcf/bcf.gz] --snplist AD_snps.txt --outname prs_vcf
```  
Calculate PRS using a single PLINK data, defining dosages, and MAF>1%
```console
jordan.R --genotype example_data_plink --dosage --snplist AD_snps.txt --maf 0.01 --outname prs_dos_maf
```  
Calculate PRS using a single PLINK data, defining dosages, MAF>1%, request frequencies and single-variant dosages
```console
jordan.R --genotype example_data_plink --dosage --snplist AD_snps.txt --maf 0.01 --freq --keepDosage --outname prs_dos_maf_freq
```  
Calculate PRS using a single PLINK data, defining dosages, MAF>1%, request frequencies and single-variant dosages. Perform association analysis of AD status with PRS, and SEX as covariate
```console
jordan.R --genotype example_data_plink --dosage --snplist AD_snps.txt --maf 0.01 --freq --keepDosage --outname prs_dos_maf_freq_assoc --assoc example_assoc_file.txt --assoc-var AD --assoc-cov SEX
```  
Calculate PRS using a single PLINK data, defining dosages, MAF>1%, request frequencies and single-variant dosages. Perform association analysis of AD and MMSE status with PRS, and SEX, PC1 and PC2 as covariates
```console
jordan.R --genotype example_data_plink --dosage --snplist AD_snps.txt --maf 0.01 --freq --keepDosage --outname prs_dos_maf_freq_assocMulti --assoc prs example_assoc_file.txt --assoc-var AD,MMSE --assoc-cov SEX,PC1,PC2
```  
Calculate PRS using a single PLINK data, defining dosages, MAF>1%, request frequencies and single-variant dosages. Perform association analysis of AD and MMSE status with PRS as well as single-variant analysis, and SEX, PC1 and PC2 as covariates
```console
jordan.R --genotype example_data_plink --dosage --snplist AD_snps.txt --maf 0.01 --freq --keepDosage --outname prs_dos_maf_freq_assocMulti --assoc both example_assoc_file.txt --assoc-var AD,MMSE --assoc-cov SEX,PC1,PC2
```  

## Test
To test everything is allright, you can use the following script alongside the example datasets provided in the repository. The example data will make a PRS from 85 SNPs associated with Alzheimer's Disease in individuals from the 1000Genome Project.  

## Contact
For comments, feedback, or questions, feel free to reach me at [n.tesi@amsterdamumc.nl](mailto:n.tesi@amsterdamumc.nl) or open an issue.



