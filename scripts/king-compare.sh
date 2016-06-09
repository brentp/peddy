module load plink
module load king/1.4
set -ex

VCF=$1
prefix=$(basename $VCF .vcf.gz)
/usr/bin/time plink --const-fid --allow-extra-chr --vcf $VCF --make-bed --out $prefix --biallelic-only --geno 0.95 --vcf-half-call m
/usr/bin/time king --ibs -b ${prefix}.bed --kinship --prefix $prefix
