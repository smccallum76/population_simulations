'''
Using output from the simulations calculate statistics using scikit-allel.
- XP-EHH - cross population extended haplotype homozygosity [sa.xpehh]
- iHS - integrated haplotype score (essentially the integral of EHH) [sa.ihs]
- DDAF - difference in derived allele frequency
- Fst - fixation index - ts kit has this option [treeSequence.Fst()]

'''

import allel as sa

neutral = sa.read_vcf("../simulations/output/ts_neutral.vcf")
sweep = sa.read_vcf("../simulations/output/ts_sweep.vcf")

# is this statistic really valuable given my data?
# - my data is the product of two populations, but is essentially a single population
xpehh = sa.xpehh(sweep['calldata/GT'][:,0:100,0], sweep['calldata/GT'][:,100:200,0],
                 pos=sweep['variants/POS'], include_edges=True )

ihs = sa.ihs(sweep['calldata/GT'][:,:,0], pos=sweep['variants/POS'], include_edges=True)

print('done')