# Test run the survival with eamples. 

mkdir temp
# extract vcf.tar.gz in temp
tar -xvzf example/genotype/example.vcf.tar.gz -C temp

VFILE=temp/chr19.vcf
FILE=temp/chr19
N=2


# Get the PASS
bcftools view -f '.,PASS' ${VFILE} \
                  -Oz -o ${FILE}.vcf.gz --threads ${N} 
# Split
bcftools norm -m-both ${FILE}.vcf.gz \
                  -Oz -o ${FILE}_split.vcf.gz --threads ${N} # split
# Create Plink binary (only keep autosome + par)
plink2 --threads ${N} \
           --vcf ${FILE}_split.vcf.gz --make-pgen --allow-extra-chr --autosome-par --out ${FILE}_split
# Rename
plink2 --threads ${N} \
           --pfile ${FILE}_split --make-pgen \
           --set-all-var-ids "chr@:#:\$r:\$a" \
           --new-id-max-allele-len 999 truncate --out ${FILE}_split_temp_renamed
# Remove duplicates
plink2 --threads ${N} \
           --pfile ${FILE}_split_temp_renamed --make-pgen --rm-dup exclude-all --out ${FILE}_split_temp_renamed_uniq
# mac2 snps only
plink2 --threads ${N} \
            --pfile ${FILE}_split_temp_renamed_uniq --make-pgen --snps-only just-acgt --mac 2 --out ${FILE}_split_temp_renamed_uniq_snps
# geno 0.05
plink2 --threads ${N} \
            --pfile ${FILE}_split_temp_renamed_uniq_snps --make-pgen --geno 0.05 dosage --out ${FILE}_split_temp_renamed_uniq_snps_geno
# export A
plink2 --threads ${N} \
           --pfile ${FILE}_split_temp_renamed_uniq_snps_geno --export A --out ${FILE}_split_temp_renamed_uniq_snps_geno_export
# run survival analysis
Rscript --vanilla bin/survival.R\
 --covar example/covariates.tsv\
 --pheno example/phenotype.surv.tsv\
 --covar-name "age_at_baseline SEX"\
 --pheno-name surv_y\
 --rawfile ${FILE}_split_temp_renamed_uniq_snps_geno_export.raw\
 --out ${FILE}_split_temp_renamed_uniq_snps_geno_export_surv.tsv
