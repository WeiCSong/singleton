chr=$1

path_to_500kwes_helper_files="/mnt/project/Bulk/Exome sequences/Population level exome OQFE variants, PLINK format - final release/helper_files"
geno_prefix="/mnt/project/Bulk/Exome sequences/Population level exome OQFE variants, BGEN format - final release/"

regenie \
    --step 2 \
    --ignore-pred \
    --bgen "${geno_prefix}ukb23159_c${chr}_b0_v1.bgen" \
    --ref-first \
    --sample "${geno_prefix}ukb23159_c${chr}_b0_v1.sample" \
    --covarFile british.regeniebase.covar \
    --phenoFile common.regenie.pheno \
    --bsize 200 \
    --force-qt \
    --set-list WES.autosome.sets.txt.gz \
    --anno-file "${path_to_500kwes_helper_files}/ukb23158_500k_OQFE.annotations.txt.gz" \
    --mask-def WES.mask \
    --aaf-bins 0.001 \
    --vc-tests skato-acat \
    --skip-sbat \
    --joint acat \
    --build-mask sum \
    --threads 8 \
    --out commonburden.$chr
    
regenie \
    --step 2 \
    --ignore-pred \
    --bgen "${geno_prefix}ukb23159_c${chr}_b0_v1.bgen" \
    --ref-first \
    --sample "${geno_prefix}ukb23159_c${chr}_b0_v1.sample" \
    --covarFile british.st.covar \
    --phenoFile common.regenie.pheno \
    --bsize 200 \
    --force-qt \
    --set-list WES.autosome.sets.txt.gz \
    --anno-file "${path_to_500kwes_helper_files}/ukb23158_500k_OQFE.annotations.txt.gz" \
    --mask-def WES.mask \
    --aaf-bins 0.001 \
    --vc-tests skato-acat \
    --skip-sbat \
    --joint acat \
    --build-mask sum \
    --threads 8 \
    --out commonburden_st.$chr

#for chr in 1:22 ;do cmd="bash commonburden.sh "$chr; tmp=$(dx run swiss-army-knife -iin="commonburden.sh" -iin="common.regenie.pheno" -iin="WES.autosome.sets.txt.gz" -iin=WES.mask -iin="british.regeniebase.covar" -iin="british.st.covar" -icmd="$cmd" --instance-type mem3_ssd1_v2_x8 --priority normal --destination /songweichen ) ; done
