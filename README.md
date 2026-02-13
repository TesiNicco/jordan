<img src="assets/jordan.png" width="180" alt="Jordan icon">

This repository provides a command-line R script and a Shiny-based graphical interface to (1) compute Polygenic Risk Scores (PRS) from genotype data (VCF/BCF/PLINK formats), (2) perform single-variant and PRS association testing, and (3) generate comprehensive diagnostic and association plots.  
The R script (`jordan.R`) can be executed standalone via the command line, while the Shiny app (`shiny_app/`) offers a user-friendly interface that internally wraps and executes the same script using dynamic input.  

`jordan` is currently under active development.

## Install
Clone the repository locally:
```console
git clone https://github.com/TesiNicco/jordan.git
```

You may need to make the main file executable by typing:
```console
chmod +x ./jordan/bin/jordan.R
```

`jordan` requires `R` and few libraries to be installed. You can either use the provided `conda` environment, or install the packages manually.
To install a conda environment with all the required packages, you can use the following code:
```console
conda env create -f jordan.yml
```
If you prefer not to use `conda`, make sure you have these required `R` packages in place: `argparse`, `data.table`, `stringr`, `ggplot2`, `survival`, `survminer`, and `ggpubr`.
In case these are not installed, you can install with the following code:
```console
R  
install.packages(c("argparse", "data.table", "stringr", "ggplot2", "survival", "survminer", "ggpubr"))
```
In addition to R packages, you need to have [PLINK2](https://www.cog-genomics.org/plink/2.0/) and [PLINK](https://www.cog-genomics.org/plink/1.9/) installed in your system.  
PLINK and PLINK2 executables are also provided in this package. 

## Shiny application
The application wrapper around `jordan.R` allows the user to run the script with a graphic interface using dynamic input. In order to use the application, few additional `R` packages need to be installed: `shiny`, `tools`, `shinyjs`, `processx`. These can be installed with the following command:
```console
R
install.packages(c("argparse", "data.table", "stringr", "ggplot2", "shiny", "tools", "shinyjs", "processx"))
```

## How to use -- jordan command line tool
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
- `--survival-var` (*Optional, default: False*): when survival analysis, should be performed, please specify here the variable name. Note that the variable should also be included in `--assoc-var`. The variable is intended to be the *time to event*. If *event* should be considered, a variable with `{var_name}_EVENT` should be included in the phenotype file.  
- `--tiles` (*Optional, default: False*): when present, tile-based analysis will be performed. This analysis divides the PRS distribution into tiles. The number of tiles can be defined. Often, tiles are made relative to a specific subset of samples (e.g healthy controls), and then applied to all individuals. The variable in the phenotype file on which to base the tile-analysis can be defined, as well as the reference group to be used. The format should be `"n_tiles;variable1;reference, n_tiles;variable2;reference"`. For example, indicating `"10;AD;0"` will make 10-tiles PRS based on the PRS distribution of the individuals labeled as 0 in AD column.  
- `--split` (*Optional, default: False*): when present, a split-based analysis will be performed. This analysis groups individuals based on defined grouping parameters. The format should be `"variable1;threshold1-threshold2-threshold3, variable2;threshold1-threshold2"`. For example, indicating `"MMSE;20-26"` will group individuals based on the MMSE column in three groups: value lower than 20, value between 20 and 26, and value larger than 26.  
- `--assoc-split-tiles` (*Optional, default: False*): when present, tile-based and/or split-based analyses will be followed by association analysis. This consist in linear regression analysis of the variable(s) indicated in the `--tiles` analysis, using as predictors the tiles grouping. Associations will be corrected for the covariates defined in `--assoc-cov`. For splits, the PRS will be compared across splits in a linear regression model, accounting for the covariates defined in `--assoc-cov`.  
- `--sex-strata` (*Optional, default: False*): when present, all analyses will be performed in males, females, and the combined set of individuals. This analysis requires that a variable SEX is included in the phenotype data.  
- `--impute-missing` (*Optional, default: False*): when present, missing genotypes will be imputed using 2*MAF dosage across all individuals. For this to be effective, a `EFFECT_ALLELE_FREQUENCY` column should be present in the SNP file.	

## PRS and weights
There can be different types of PRS: (i) *unweighted*, (ii) *weighted*, and (iii) *multiple-weighted*. *Unweighted* PRS are simply the sum of trait-associated alleles, with all variants having the same weight. To do this in jordan, simply set the BETA to 1 in the `--snplist` file to all variants. *Weighted* PRS are the most commonly used PRS, defined as the weighted sum of trait-associated alleles, weighted by an effect size, typically originating from a GWAS study. In `jordan`, PRS calculation is implemented as:  
$PRS_{sample} = \sum_{snp}^{SNPs} \alpha_{snp} \cdot \beta_{snp}$  
Finally, *multiple-weighted* PRS can be handy for calculating pathway-specific PRS or cell-specific PRS. In these cases, the additional weight (on top of the BETA, quantifying the effect on the trait of interest) is needed. This additional weight will quantify the effect of each variant on a specific pathway, or cell-type. When an additional weight is selected with `--addWeight` option, for example for the calculation of pathway-specific PRS, then the formula will adapt accordingly to:  
$PRS_{sample} = \sum_{snp}^{SNPs} \alpha_{snp} \cdot \beta_{snp} \cdot w_{snp}$  
We previously used [PRS](https://alz-journals.onlinelibrary.wiley.com/doi/epdf/10.1002/alz.13810) and [pathway-PRS](https://pmc.ncbi.nlm.nih.gov/articles/PMC7524800/#:~:text=Immune%20response%20and%20endocytosis%20pathways,resilience%20against%20Alzheimer's%20disease%20%2D%20PMC) for the associations with Alzheimer's Disease.

## Example usage
The `example_scripts.sh` file in `example_data` folder reports multiple example commands. Here below it is reported its content:
Calculate PRS using a single PLINK data file (most basic use case)
```console
jordan.R --genotype example_data_plink --snplist ADsnps_cojo_5e-08_list_HG38.txt --impute-missing --outname prs
```  
Calculate PRS using a multiple PLINK data file (eg. chr1_plink, chr2_plink, etc.) [not warking as no individual chromosomes are provided]
```console
jordan.R --genotype example_data_plink_chr1 --multiple --snplist ADsnps_cojo_5e-08_list_HG38.txt --impute-missing --outname prs
```  
Calculate PRS using a single VCF file (also working with gzipped VCF files)
```console
jordan.R --genotype example_data.[vcf/vcf.gz/bcf/bcf.gz] --snplist ADsnps_cojo_5e-08_list_HG38.txt --impute-missing --outname prs_vcf
```  
Calculate PRS using a single PLINK data, defining dosages, and MAF>1%
```console
jordan.R --genotype example_data_plink --dosage --snplist ADsnps_cojo_5e-08_list_HG38.txt --impute-missing --maf 0.01 --outname prs_dos_maf
```  
Calculate PRS using a single PLINK data, defining dosages, MAF>1%, request frequencies and single-variant dosages
```console
jordan.R --genotype example_data_plink --dosage --snplist ADsnps_cojo_5e-08_list_HG38.txt --impute-missing --maf 0.01 --freq --keepDosage --outname prs_dos_maf_freq
```  
Calculate PRS using a single PLINK data, defining dosages, MAF>1%, request frequencies and single-variant dosages. Perform association analysis of AD status with PRS, and SEX as covariate
```console
jordan.R --genotype example_data_plink --dosage --snplist ADsnps_cojo_5e-08_list_HG38.txt --impute-missing --maf 0.01 --freq --keepDosage --outname prs_dos_maf_freq_assoc --assoc example_assoc_file.txt --assoc-var AD --assoc-cov SEX
```  
Calculate PRS using a single PLINK data, defining dosages, MAF>1%, request frequencies and single-variant dosages. Perform association analysis of AD and MMSE status with PRS, and SEX, PC1 and PC2 as covariates
```console
jordan.R --genotype example_data_plink --dosage --snplist ADsnps_cojo_5e-08_list_HG38.txt --impute-missing --maf 0.01 --freq --keepDosage --outname prs_dos_maf_freq_assocMulti --assoc prs example_assoc_file.txt --assoc-var AD,MMSE --assoc-cov SEX,PC1,PC2
```  
Calculate PRS using a single PLINK data, defining dosages, MAF>1%, request frequencies and single-variant dosages. Perform association analysis of AD and MMSE status with PRS as well as single-variant analysis, and SEX, PC1 and PC2 as covariates
```console
jordan.R --genotype example_data_plink --dosage --snplist ADsnps_cojo_5e-08_list_HG38.txt --impute-missing --maf 0.01 --freq --keepDosage --outname prs_dos_maf_freq_assocMulti --assoc both example_assoc_file.txt --assoc-var AD,MMSE --assoc-cov SEX,PC1,PC2
```

## Test
To test everything is allright, you can use the following script alongside the example datasets provided in the repository. The example data will make a PRS from 120 SNPs associated with Alzheimer's Disease in individuals from the 1000Genome Project. Some SNPs are missing in the genotype file and will be imputed.

## Application wrapper
The application can be run from anywhere in your system with the following code:
```console
Rscript -e "shiny::runApp('/path/to/jordan/bin/shiny_app/app.R', port = 4525, host = '127.0.0.1')"
```  
The server normally exposes the graphic interface on a port (in the second comman it is specified as 4525, but can be changed). The user can then user any browser to open the connection, for example by using the URL http://127.0.0.1:4525.

## Centenarian Collaboration
We are conducting single-variant and PRS analyses related to Alzheimer's Disease risk in Centenarians as opposed to Healthy younger individuals, as well as in the age-continuum. To conduct analyses and ensure replicability across cohorts, we recommend using the following code for the analysis including Sex as covariate:
```console
jordan.R --genotype path/to/genotype/chr1.dosage --dosage --multiple --exclude --snplist ADsnps_cojo_5e-08_list_HG38.txt --impute-missing --outname /path/to/output/folder --freq --sex-strata --assoc both /path/to/phenotype/data --assoc-var MMSE,SURVIVAL,AD_Control,AD_Centenarian,Centenarian_Control --assoc-cov SEX,PC1,PC2,PC3,PC4,PC5 --survival-var SURVIVAL --tiles "10;Centenarian_Control;1" --split "MMSE;20-26" --assoc-split-tiles
```
And the following code for conducting analyses excluding Sex as covariate:
```console
jordan.R --genotype path/to/genotype/chr1.dosage --dosage --multiple --exclude --snplist ADsnps_cojo_5e-08_list_HG38.txt --impute-missing --outname /path/to/output/folder --freq --sex-strata --assoc both /path/to/phenotype/data --assoc-var MMSE,SURVIVAL,AD_Control,AD_Centenarian,Centenarian_Control --assoc-cov PC1,PC2,PC3,PC4,PC5 --survival-var SURVIVAL --tiles "10;Centenarian_Control;1" --split "MMSE;20-26" --assoc-split-tiles
```
*Note: These lines of code assume one file per chromosome (PLINK or VCF format). In addition to genetics data, a phenotype file should be present, and formatted as described above. This code will (1) generate PRS with the variants of interest (including and excluding APOE variants), (2) will perform single-variant analyses of the variants of interest, (3) will report on allelic frequency, (4) will perform PRS analysis across the defined groups (correcting for defined covariates), (5) will perform a survival analysis, (6) a percentile-based analysis, and (7) will compare the PRS between groups split based on MMSE. Not all analyses are mandatory, but if data is available, we encourage collaborators to perform all of them.*

## Contact
For comments, feedback, or questions, feel free to reach me at [n.tesi@amsterdamumc.nl](mailto:n.tesi@amsterdamumc.nl) or open an issue.



