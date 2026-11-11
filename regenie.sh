chr=$1

filename="/mnt/project/songweichen/plink/ukb24308_c${chr}_b0_v1"
annoname="/mnt/project/songweichen/burdenmask/c${chr}.regenie.annot"
setname="/mnt/project/songweichen/c${chr}.coding.regenie.setlist"
maskname="/mnt/project/songweichen/masks"
regenie \
  --step 2 \
  --pgen $filename \
  --covarFile british.regeniebase.covar \
  --phenoFile common.regenie.pheno \
  --bsize 200 \
  --pred common_pred.list \
  --out common.burden.$chr \
  --force-qt \
  --anno-file $annoname \
  --set-list $setname \
  --mask-def $maskname \
  --aaf-bins 0.01,0.001 \
  --vc-tests skato-acat \
  --skip-sbat \
  --joint acat \
  --build-mask sum \
  --threads 8
  
regenie \
  --step 2 \
  --pgen $filename \
  --covarFile british.st.covar \
  --phenoFile common.regenie.pheno \
  --bsize 200 \
  --pred common_pred.list \
  --out common.burden.withst.$chr \
  --force-qt \
  --anno-file $annoname \
  --set-list $setname \
  --mask-def $maskname \
  --aaf-bins 0.01,0.001 \
  --vc-tests skato-acat \
  --skip-sbat \
  --joint acat \
  --build-mask sum \
  --threads 8

#for chr in 1:22 ;do cmd="bash commonregenie.sh "$chr; tmp=$(dx run swiss-army-knife -iin="commonregenie.sh" -iin="common.regenie.pheno" -iin="BMI_1.loco" -iin="chronotype_1.loco"  -iin="height_1.loco" -iin="HDL_1.loco" -iin="tanning_1.loco" -iin="RBC.loco" -iin="common_pred.list" -iin="british.regeniebase.covar" -iin="british.st.covar" -icmd="$cmd" --instance-type mem3_ssd1_v2_x8 --priority normal --destination /songweichen ) ; done
