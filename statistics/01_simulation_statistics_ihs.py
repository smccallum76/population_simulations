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

gt_neutral_AFR = sa.GenotypeArray(neutral['calldata/GT'][:, afr_s:afr_e, :])
gt_neutral_EUR = sa.GenotypeArray(neutral['calldata/GT'][:, eur_s:eur_e, :])
# allele counts at each position
ac_afr = gt_sweep_AFR.count_alleles()
ac_eur = gt_sweep_EUR.count_alleles()

ac_afrN = gt_neutral_AFR.count_alleles()
ac_eurN = gt_neutral_EUR.count_alleles()

'''
iHS - Unstandardized
- it is normal that iHS returns a lot of nans
- this should only be applied to one population at a time (AFR and then EUR)
'''
ihs_afr = sa.ihs(sweep['calldata/GT'][:, afr_s:afr_e, 0], pos=sweep['variants/POS'], include_edges=False, use_threads=True)
ihs_eur = sa.ihs(sweep['calldata/GT'][:, eur_s:eur_e, 0], pos=sweep['variants/POS'], include_edges=False, use_threads=True)
ihs_afrN = sa.ihs(neutral['calldata/GT'][:, afr_s:afr_e, 0], pos=neutral['variants/POS'], include_edges=False, use_threads=True)
ihs_eurN = sa.ihs(neutral['calldata/GT'][:, eur_s:eur_e, 0], pos=neutral['variants/POS'], include_edges=False, use_threads=True)

'''
iHS - Standardized
'''
ihs_afr_std, ac_afr_bins = sa.standardize_by_allele_count(ihs_afr, ac_afr[:, 1], diagnostics=False)
ihs_eur_std, ac_eur_bins = sa.standardize_by_allele_count(ihs_eur, ac_eur[:, 1], diagnostics=False)
ihs_afrN_std, ac_afrN_bins = sa.standardize_by_allele_count(ihs_afrN, ac_afrN[:, 1], diagnostics=False)
ihs_eurN_std, ac_eurN_bins = sa.standardize_by_allele_count(ihs_eurN, ac_eurN[:, 1], diagnostics=False)

'''SNP Position Plots - iHS'''
fig, ax = plt.subplots(nrows=2, ncols=1, sharey=True, figsize=(15, 8))
ax[0].scatter(sweep['variants/POS'], abs(ihs_afr_std), label='AFR Sweep', s=3, color='blue')
ax[0].scatter(sweep['variants/POS'], abs(ihs_eur_std)*-1, label='EUR Sweep', s=3, color='red')
ax[0].axvline(x=np.max(sweep['variants/POS'])/2, color='black', linestyle='--', label="Sweep Location")
ax[0].set_ylabel("|iHS| or |iHS|*-1")
ax[0].set_title("Sweep Sims -- iHS using AFR and EUR populations")
ax[0].legend(loc="upper left")

ax[1].scatter(neutral['variants/POS'], abs(ihs_afrN_std), label='AFR Neutral', s=3, color='blue')
ax[1].scatter(neutral['variants/POS'], abs(ihs_eurN_std)*-1, label='EUR Neutral', s=3, color='red')
ax[1].axvline(x=np.max(sweep['variants/POS'])/2, color='black', linestyle='--', label="Sweep Location")
ax[1].set_xlabel("SNP Number")
ax[1].set_ylabel("|iHS| or |iHS|*-1")
ax[1].set_title("Neutral Sims -- iHS using AFR and EUR populations")
ax[1].legend(loc="upper left")
plt.show()

'''Histogram Plots - iHS '''
n_bins = 100
fig, ax = plt.subplots(nrows=2, ncols=2, sharey=True, figsize=(12, 8))
ax[0, 0].hist(ihs_afr_std, bins=n_bins, color='blue')
ax[0, 1].hist(ihs_eur_std, bins=n_bins, color='blue')
ax[1, 0].hist(ihs_afrN_std, bins=n_bins, color='red')
ax[1, 1].hist(ihs_eurN_std, bins=n_bins, color='red')
ax[1, 0].set_xlabel('iHS Value')
ax[0, 0].set_title('Histogram of iHS AFR Sweep')
ax[0, 1].set_title('Histogram of iHS EUR with Sweep in AFR')
ax[1, 0].set_title('Histogram of iHS AFR Neutral')
ax[1, 1].set_title('Histogram of iHS EUR Neutral')
plt.show()
