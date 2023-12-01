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

neutral = sa.read_vcf("../simulations/output/ts_neutral_noScaling.vcf")
sweep = sa.read_vcf("../simulations/output/ts_sweep_noScaling.vcf")
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
Plot of allele count -- qc to ensure that we are seeing a bunch of 1's near the sweep
'''
allel_pos = 1
fig, ax = plt.subplots(nrows=2, ncols=1, sharey=False, sharex=True, figsize=(18, 8))
ax[0].plot(sweep['variants/POS'], ac_afr[:, allel_pos], label="AFR ac Sweep", alpha=0.7)
ax[0].plot(sweep['variants/POS'], ac_eur[:, allel_pos]*-1, label="EUR ac Sweep", alpha=0.7, color='red')
ax[1].plot(neutral['variants/POS'], ac_afrN[:, allel_pos], label="AFR ac Neutral", alpha=0.7)
ax[1].plot(neutral['variants/POS'], ac_eurN[:, allel_pos]*-1, label="EUR ac Neutral", alpha=0.7, color='red')
ax[0].axvline(x=np.max(sweep['variants/POS'])/2, color='black', linestyle='--', label='Sweep Location')
ax[0].set_xlabel("SNP Number")
ax[0].set_title("Sweep Sims -- Allele Frequency by SNP")
ax[0].legend(loc="upper left")
ax[1].set_xlabel("SNP Number")
ax[1].axvline(x=np.max(sweep['variants/POS'])/2, color='black', linestyle='--', label='Sweep Location')
ax[1].set_title("Neutral Sims -- Allele Frequency by SNP")
ax[1].legend(loc="upper left")
plt.show()

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
xpehhS = sa.xpehh(sweep['calldata/GT'][:, afr_s:afr_e, 0], sweep['calldata/GT'][:, eur_s:eur_e, 0],
                 pos=sweep['variants/POS'], include_edges=False, use_threads=True)
xpehhN = sa.xpehh(neutral['calldata/GT'][:, afr_s:afr_e, 1], neutral['calldata/GT'][:, eur_s:eur_e, 1],
                 pos=neutral['variants/POS'], include_edges=False, use_threads=True)

'''Genotype Position Plots - XP-EHH'''
plt.figure(figsize=(18, 5))
plt.plot(np.arange(0, len(xpehhS), 1), xpehhS, label="AFR/EUR Sweep", alpha=0.7)
plt.plot(np.arange(0, len(xpehhN), 1), xpehhN, label="AFR/EUR Neutral", alpha=0.7, color='red')
plt.axhline(y=0, color='black')
plt.axvline(x=len(xpehhS)/2, color='black', linestyle='--', label='Sweep Location')
plt.xlabel("Genomic Position")
plt.ylabel("XP-EHH")
plt.title("XP-EHH using AFR and EUR populations")
plt.legend(loc="upper left")
plt.show()

'''Histogram Plots - XP-EHH '''
n_bins = 100
fig, ax = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(12, 5))
ax[0].hist(xpehhS, bins=n_bins, color='blue')
ax[1].hist(xpehhN, bins=n_bins, color='red')
ax[0].set_xlabel('XP-EHH Value')
ax[1].set_xlabel('XP-EHH Value')
ax[0].set_ylabel('Frequency')
ax[0].set_title('Histogram of XP-EHH Sweep')
ax[1].set_title('Histogram of XP-EHH Neutral')
plt.show()

'''
iHS - Unstandardized
- it is normal that iHS returns a lot of nans
- this should only be applied to one population at a time (AFR and then EUR)
'''
ihs_afr = sa.ihs(sweep['calldata/GT'][:, afr_s:afr_e, 0], pos=sweep['variants/POS'], include_edges=False, use_threads=True)
ihs_eur = sa.ihs(sweep['calldata/GT'][:, eur_s:eur_e, 0], pos=sweep['variants/POS'], include_edges=False, use_threads=True)
ihs_afrN = sa.ihs(neutral['calldata/GT'][:, afr_s:afr_e, 0], pos=neutral['variants/POS'], include_edges=False, use_threads=True)
ihs_eurN = sa.ihs(neutral['calldata/GT'][:, eur_s:eur_e, 0], pos=neutral['variants/POS'], include_edges=False, use_threads=True)

# Standardize the iHS scores
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

'''
Fst 
see example for hudson_fst here:
https://scikit-allel.readthedocs.io/en/stable/stats/fst.html
'''
num, den = sa.hudson_fst(ac_afr, ac_eur)
fst = num / den  # fst for each variant individually
fst_overall = np.sum(num) / np.sum(den)  # fst averaging over all variants

num2, den2 = sa.hudson_fst(ac_afrN, ac_eurN)
fstN = num2 / den2  # fst for each variant individually
fst_overall2 = np.sum(num) / np.sum(den)  # fst averaging over all variants

plt.figure(figsize=(18, 5))
plt.plot(np.arange(0, len(fst), 1), fst, label="Fst")
plt.axvline(x=len(fst)/2, color='black', linestyle='--', label='Sweep Location')
plt.xlabel("Genomic Position")
plt.ylabel("Fst")
plt.title("Fst using AFR and EUR populations")
plt.legend(loc="upper left")
plt.show()
print('done')

'''Histogram Plots - Fst '''
n_bins = 100
fig, ax = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(12, 5))
ax[0].hist(fst, bins=n_bins, color='blue', log=True)
ax[1].hist(fstN, bins=n_bins, color='red', log=True)
ax[0].set_xlabel('Fst Value')
ax[1].set_xlabel('Fst Value')
ax[0].set_ylabel('Frequency[log scale]')
ax[0].set_title('Histogram of Fst Sweep')
ax[1].set_title('Histogram of Fst Neutral')
plt.show()