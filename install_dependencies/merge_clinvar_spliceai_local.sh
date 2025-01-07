## From HerediVar see /data/script/merge_clinvar_spliceai.sh

basedir=/home/katzkean

root=$basedir
data=$root/databases/Clinvar
ngsbits=$root/ngs-bits/bin


clinvar_dat=$data/clinvar_snv.vcf.gz
spliceai_snv_dat=$data/spliceai_scores.masked.snv.hg38.vcf.gz
spliceai_indel_dat=$data/spliceai_scores.masked.indel.hg38.vcf.gz


$ngsbits/VcfAnnotateFromVcf -in $clinvar_dat -source $spliceai_snv_dat -threads 5 -info_keys SpliceAI -out $data/clinvar_annotated_snv.vcf
$ngsbits/VcfAnnotateFromVcf -in $data/clinvar_annotated_snv.vcf -source $spliceai_indel_dat -threads 5 -info_keys SpliceAI -out $data/clinvar_annotated_snv_indel.vcf



egrep "^#|;SpliceAI=" $data/clinvar_annotated_snv_indel.vcf > $data/clinvar_only_annotated.vcf
$ngsbits/VcfSubstract -in $data/clinvar_annotated_snv.vcf -in2 $data/clinvar_only_annotated.vcf -out $data/variants_missing_annotation.vcf


bgzip -f $data/variants_missing_annotation.vcf
tabix -f -p vcf $data/variants_missing_annotation.vcf.gz
spliceai -I $data/variants_missing_annotation.vcf.gz -O $data/variants_added_annotation.vcf -R $grch38 -A grch38 -M 1


$ngsbits/VcfAdd -in $data/variants_added_annotation.vcf -in2 $data/clinvar_only_annotated.vcf -out $data/clinvar_spliceai_all.vcf

$ngsbits/VcfSort -in $data/clinvar_spliceai_all.vcf -out $data/clinvar_spliceai_all_sorted.vcf

bgzip clinvar_spliceai_all_sorted.vcf
