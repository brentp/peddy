if [[ -x $(which module) ]]; then
	module load plink
	module load king/1.4
fi
set -ex

VCF=$1
prefix=$(basename $VCF .vcf.gz)
echo $prefix
/usr/bin/time plink2 --const-fid --allow-extra-chr --vcf $VCF --make-bed --out $prefix --biallelic-only --geno 0.95 --vcf-half-call m
/usr/bin/time king --ibs -b ${prefix}.bed --kinship --prefix $prefix
