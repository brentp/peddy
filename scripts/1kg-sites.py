import sys
from cyvcf2 import VCF

kg_vcf = sys.argv[1]

pops = ["EUR", "AFR", "AMR", "SAS"]

good = 0
for v in VCF(kg_vcf):
    info = v.INFO
    if 'OLD_MULTIALLELIC' in info: continue
    if info['VT'] != 'SNP': continue
    if info['NS'] < 2500: continue
    if info['AF'] < 0.04: continue
    if info['AF'] > 0.95: continue
    try:
        info['EX_TARGET']
    except KeyError:
        continue

    if not all(info[p + "_AF"] > 0.04 for p in pops):
        continue

    good += 1
    print "%s:%d" % (v.CHROM, v.POS)
