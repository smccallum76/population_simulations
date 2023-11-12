'''
Using output from the simulations calculate statistics using scikit-allel.
- XP-EHH - cross population extended haplotype homozygosity [sa.xpehh]
- iHS - integrated haplotype score (essentially the integral of EHH) [sa.ihs]
- DDAF - difference in derived allele frequency [abs(sa.sfs(pop1_alleles) - sa.sfs(pop2_alleles))]
- Fst - fixation index - ts kit has this option [treeSequence.Fst()] [sa.hudson_fst()]

'''

import allel as sa

neutral = sa.read_vcf("../simulations/output/ts_neutral.vcf")
sweep = sa.read_vcf("../simulations/output/ts_sweep.vcf")
# convenient function for capturing the genotype array from the vcf import
gt_neutral = sa.GenotypeArray(neutral['calldata/GT'])
gt_sweep = sa.GenotypeArray(sweep['calldata/GT'])
# allele counts at each position
ac_neutral = gt_neutral.count_alleles()
ac_sweep = gt_sweep.count_alleles()
'''
XP-EHH
- Are the two populations, in this example, AFR and EUR?
- If so, then is the comparison below correct? Specifically, should I be slicing the such that I'm comparing
SNPs 0:100 to SNPs 100:200?
- Why are there negative values and values > 1? I thought this was supposed to be a probability
        - Disregard - The standardized value can be neg or > 1
- Is the 'pos' argument being correctly considered? 
- Main issue --> the output from the simulation is a single population that represents both AFR and EUR. All of this 
is included in the single VCF file. The first half of the individuals are AFR and the second half are EUR.

'''
# include edges ensures the calculation is performed at the beginning and end of the array
xpehh = sa.xpehh(sweep['calldata/GT'][:, 0:9, 0], sweep['calldata/GT'][:, 10:19, 0],
                 pos=sweep['variants/POS'], include_edges=True, use_threads=True)
'''
iHS
- Seems to be a lot of nans, not sure if I'm doing this correctly. 
'''
ihs = sa.ihs(sweep['calldata/GT'][:, :, 0], pos=sweep['variants/POS'], include_edges=True, use_threads=True)


'''
Fst
- similar questions as above. Need to ensure that I correctly understand the data and what qualifies as a 
population. Is AFR considered one population and EUR the other? 

'''
# fst = sa.hudson_fst(allel_counts_AFR, allel_counts_EUR)

print('done')