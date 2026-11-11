file=$1
chr=$(echo "$file" | grep -Eo '[0-9]+')

bgen="/mnt/project/Bulk/bgen/ukb24309_c${chr}_b0_v1.bgen"
sample="/mnt/project/Bulk/bgen/ukb24309_c${chr}_b0_v1.sample"
regenie \
  --step 2 \
  --bgen $bgen \
  --sample $sample \
  --ref-first \
  --extract $file \
  --covarFile wgs.regenie.covar \
  --phenoFile British.normburden.pheno \
  --bsize 200 \
  --pred British.normburden_pred.list \
  --out British.normburden.$chr \
  --minMAC 10 \
  --threads 4
  


#for file in chr*;do cmd="bash stregenie.sh "$file; tmp=$(dx run swiss-army-knife -iin="stregenie.sh" -iin="$file" -iin="British.normburden.pheno" -iin="British.normburden_1.loco" -iin="British.normburden_pred.list" -iin="wgs.regenie.covar" -icmd="$cmd" --instance-type mem3_ssd1_v2_x4 --priority normal --destination /songweichen ) ; done
