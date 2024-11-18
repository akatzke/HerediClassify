## From HerediVar see /data/script/merge_clinvar_spliceai.sh

basedir=/home/user/

root=$basedir
data=$root/data/Clinvar
ngsbits=$data/ngs-bits/bin
grch38=$data/genomes/GRCh38.fa


clinvar_dat=$data/clinvar_snv.vcf.gz
spliceai_snv_dat=$data/spliceai_scores.masked.snv.hg38.vcf.gz
spliceai_indel_dat=$data/spliceai_scores.masked.indel.hg38.vcf.gz

cd $data

# 1. Annotate ClinVar vcf using preannotated SpliceAI scores
$ngsbits/VcfAnnotateFromVcf -in $clinvar_dat -source $spliceai_snv_dat -threads 5 -info_keys SpliceAI -out clinvar_annotated_snv.vcf
$ngsbits/VcfAnnotateFromVcf -in clinvar_annotated_snv.vcf -source $spliceai_indel_dat -threads 5 -info_keys SpliceAI -out clinvar_annotated_snv_indel.vcf


# 2. Select ClinVar entries that don't have a SpliceAI annotation yet
# This step is only needed in case the missing SpliceAI scores will be calculated using SpliceAI
# If not calculate missing SpliceAI scores, go to 5.
egrep "^#|;SpliceAI=" clinvar_annotated_snv_indel.vcf > clinvar_only_annotated.vcf
$ngsbits/VcfSubstract -in clinvar_annotated_snv_indel.vcf -in2 clinvar_only_annotated.vcf -out variants_missing_annotation.vcf

# 3. Calculate missing SpliceAI scores
bgzip -f variants_missing_annotation.vcf
tabix -f -p vcf variants_missing_annotation.vcf.gz
spliceai -I variants_missing_annotation.vcf.gz -O variants_added_annotation.vcf -R $grch38 -A grch38 -M 1

# 4. Add calcualted SpliceAI scores to annotate SpliceAI file
$ngsbits/VcfAdd -in variants_added_annotation.vcf -in2 clinvar_only_annotated.vcf -out clinvar_spliceai_all.vcf

# 5. Sort and zip file
$ngsbits/VcfSort -in clinvar_spliceai_all.vcf -out clinvar_spliceai_all_sorted.vcf

bgzip clinvar_spliceai_all_sorted.vcf
