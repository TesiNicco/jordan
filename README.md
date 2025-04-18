# jordan
R script to calculate Polygenic Risk Score from genotype data, and eventually plot them.

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
If all this is OK, then you should be able to run the script.  

## How to use
By running:  
```console
./bin/jordan.R -h
```
will display the help message. `jordan` parameters are:  
- `--genotype` (**Mandatory**): genotype data in PLINK or VCF format  
- `--snplist` (**Mandatory**): file containing the list of SNPs to use for the PRS, including at least CHROM, POS, EFFECT_ALLELE, OTHER_ALLELE, and BETA (or OR).  
- `--outname` (**Mandatory**): the name of the output folder. If a directory exists, output files will be placed there, otherwise, the directory will be created.  
- `--isdosage` (*Optional, default: False*): binary flag whether data is genotyped or imputed (TRUE/FALSE)  
- `--multiple` (*Optional, default: False*): whether multiple input (PLINK) files should be used. In this case, all files with the same extension (PLINK files) will be used    
- `--exclude` (*Optional, default: False*): whether PRS including and excluding APOE variants (e2 and e4 alleles) should be made  
- `--maf` (*Optional, default: False*): minor allele frequency (MAF) threshold to include variants in the PRS. It should be, for example, 0.01 for MAF>1%  
- `--risk` (*Optional, default: True*): force the PRS calculation wrt to the risk allele. In this case, all SNPs will be flipped to the risk allele based on BETA (or OR) and PRS will be calculated accordingly  
- `--keepDosage` (*Optional, default: False*): the raw dosages of all the SNPs of interest will be written as a table. The format is `--export A` from [PLINK/PLINK2](https://www.cog-genomics.org/plink/2.0/data#export)  
- `--addWeight` (*Optional, default: False*): whether an additional weight should be included in the PRS. If desired, add the column name of the input SNP-list file that report the additional weight to use. See *PRS calculation* section below for more information  
- `--plot` (*Optional, default: False*): binary flag to plot densities of not (TRUE/FALSE)  

## PRS calculation
PRS are weighted sum of trait-associated alleles, weighted by an effect size, typically originating from a GWAS study. In `jordan`, PRS calculation is implemented as:  
$PRS_{sample} = \sum_{snp}^{SNPs} \alpha_{snp} \cdot \beta_{snp}$  
When an additional weight is selected with `--addWeight` option, for example for the calculation of pathway-specific PRS, then the formula will adapt accordingly to:  
$PRS_{sample} = \sum_{snp}^{SNPs} \alpha_{snp} \cdot \beta_{snp} \cdot w_{snp}$  

## Example usage
Single PLINK file storing genotype data as dosages, no plots, PRS with and without APOE variants, with MAF>5%:  
```console
./bin/jordan.R --genotype genotype/example_data_plink --snplist Snps_interest.txt --outname test_output --isdosage TRUE --exclude TRUE --maf 0.05
```
Multiple PLINK files with genotype data as dosages, no plots, PRS with and without APOE variants, with MAF>5%:  
```console
./bin/jordan.R --genotype genotype_data_path/chr1.example_data_plink --snplist Snps_interest.txt --outname test_output --isdosage TRUE --exclude TRUE --multiple TRUE --maf 0.05
```
Multiple PLINK files with genotype data as dosages, no plots, PRS with and without APOE variants, with MAF>5%, keep dosages and use an additional weight:  
```console
./bin/jordan.R --genotype genotype_data_path/chr1.example_data_plink --snplist Snps_interest.txt --outname test_output --isdosage TRUE --exclude TRUE --multiple TRUE --keepDosage TRUE --addWeight WEIGHT --maf 0.05
```
Single VCF file with genotype data as dosages, no plots, PRS with and without APOE variants, with MAF>5%:  
```console
./bin/jordan.R --genotype genotype_data_path/example_data.vcf --snplist Snps_interest.txt --outname test_output --isdosage TRUE --exclude TRUE --maf 0.05
```

## Test
To test everything is allright, you can use the following script alongside the example datasets provided in the repository. The example data will make a PRS from 85 SNPs associated with Alzheimer's Disease in individuals from the 1000Genome Project:  
Multiple PLINK files with genotype data as dosages, no plots, PRS with and without APOE variants, with MAF>5%:  
```console
./bin/jordan.R --genotype example_data/example_data_plink --snplist example_data/AD_snps.txt --outname test_output --isdosage TRUE`  
```
This should run successfully, and will display results in the output directory specified:  
- PRS_table.txt: tab-delimited table with samples in rows and PRS in column  
- snpsInterest.txt: tab-delimited table with information about included/excluded SNPs.  

## Contact
For comments, feedback, or questions, feel free to reach me at [n.tesi@amsterdamumc.nl](mailto:n.tesi@amsterdamumc.nl) or open an issue.



