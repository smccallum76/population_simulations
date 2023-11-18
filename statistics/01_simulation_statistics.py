'''
Using output from the simulations calculate statistics using scikit-allel.
- XP-EHH - cross population extended haplotype homozygosity [sa.xpehh]
- iHS - integrated haplotype score (essentially the integral of EHH) [sa.ihs]
- DDAF - difference in derived allele frequency [abs(sa.sfs(pop1_alleles) - sa.sfs(pop2_alleles))]
- Fst - fixation index - ts kit has this option [treeSequence.Fst()] [sa.hudson_fst()]

'''

import allel as sa
import matplotlib.pyplot as plt
import numpy as np

afr_s = 0  # Start of AFR pop
afr_e = 100  # End of AFR pop
eur_s = afr_e  # Start of EUR pop
eur_e = 200  # End of EUR pop

neutral = sa.read_vcf("../simulations/output/ts_neutral.vcf")
sweep = sa.read_vcf("../simulations/output/ts_sweep.vcf")
# convenient function for capturing the genotype array from the vcf import
gt_sweep_AFR = sa.GenotypeArray(sweep['calldata/GT'][:, afr_s:afr_e, :])
gt_sweep_EUR = sa.GenotypeArray(sweep['calldata/GT'][:, eur_s:eur_e, :])
# allele counts at each position
ac_afr = gt_sweep_AFR.count_alleles()
ac_eur = gt_sweep_EUR.count_alleles()
'''
XP-EHH
- Should be comparing AFR to EUR (these are the two populations)
- The VCF file should be sliced according to the positions of the AFR and EUR populations
- Why are there negative values and values > 1? I thought this was supposed to be a probability
        - Disregard - The standardized value can be neg or > 1
- We can use the 'pos' argument instead of the 'map_pos' argument. The 'variants/POS' is the correct indices to use
for this argument (as is done below). 
'''
# include edges ensures the calculation is performed at the beginning and end of the array
xpehh = sa.xpehh(sweep['calldata/GT'][:, afr_s:afr_e, 0], sweep['calldata/GT'][:, eur_s:eur_e, 0],
                 pos=sweep['variants/POS'], include_edges=False, use_threads=True)

plt.plot(np.arange(0, len(xpehh), 1), xpehh, label="AFR/EUR")
plt.axhline(y=0, color='black')
plt.axvline(x=len(xpehh)/2, color='black', linestyle='--', label='Sweep Location')
plt.xlabel("Genomic Position")
plt.ylabel("XP-EHH")
plt.title("XP-EHH using AFR and EUR populations")
plt.legend(loc="upper left")
plt.show()
'''
iHS
- it is normal that iHS returns a lot of nans
- this should only be applied to one population at a time (AFR and then EUR)
'''
ihs_afr = sa.ihs(sweep['calldata/GT'][:, afr_s:afr_e, 0], pos=sweep['variants/POS'], include_edges=False, use_threads=True)
ihs_eur = sa.ihs(sweep['calldata/GT'][:, eur_s:eur_e, 0], pos=sweep['variants/POS'], include_edges=False, use_threads=True)

plt.plot(np.arange(0, len(ihs_afr), 1), ihs_afr, label='AFR')
plt.plot(np.arange(0, len(ihs_eur), 1), ihs_eur, label='EUR')
plt.axhline(y=0, color='black')
plt.axvline(x=len(ihs_afr)/2, color='black', linestyle='--', label="Sweep Location")
plt.xlabel("Genomic Position")
plt.ylabel("iHS")
plt.title("iHS using AFR and EUR populations")
plt.legend(loc="upper left")
plt.show()
'''
Fst 
see example for hudson_fst here:
https://scikit-allel.readthedocs.io/en/stable/stats/fst.html
'''
num, den = sa.hudson_fst(ac_afr, ac_eur)
fst = num / den  # fst for each variant individually
fst_overall = np.sum(num) / np.sum(den)  # fst averaging over all variants

plt.plot(np.arange(0, len(fst), 1), fst, label="Fst")
plt.axvline(x=len(fst)/2, color='black', linestyle='--', label='Sweep Location')
plt.xlabel("Genomic Position")
plt.ylabel("Fst")
plt.title("Fst using AFR and EUR populations")
plt.legend(loc="upper left")
plt.show()
print('done')