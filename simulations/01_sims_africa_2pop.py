'''
Simulation of a selective sweep by a single beneficial mutation using the Out of Africa 2 population demographic
https://popsim-consortium.github.io/stdpopsim-docs/stable/catalog.html#sec_catalog_homsap_models_outofafrica_2t12
'''
# sloppy, but the 'warnings' module ignores a couple future warnings in stdpopsim
import warnings
import pandas as pd
warnings.filterwarnings("ignore")
import stdpopsim
import matplotlib.pyplot as plt
import numpy as np
import time

# timing 1e6 no scaling = 5 mins
# timing 1e7 no scaling = 8.6 mins

start = time.time()

species = stdpopsim.get_species("HomSap")  # Homo sapiens
model = species.get_demographic_model("OutOfAfrica_2T12")
samples = {"AFR": 100, "EUR": 100}
# the contig is fully neutral, the sweeping mutation will be inserted later
# note, 'right' implies the length of chr22 to use. The full chromosome is 51 million SNPS in length, but
# we will just use a subset of this chromosome.
contig = species.get_contig("chr22", right=5e6)

'''
The next session will provide settings for the sweep. This includes:
 - randomly selecting a chromosome in a population of our choice
 - The position will be specific within the contig
'''

locus_id = "hard_sweep"
coordinate = round(contig.length / 2)  # this is where the sweep will be located
contig.add_single_site(
    id=locus_id,
    coordinate=coordinate
)

'''
Next the "extended events" are set up, which will modify the demography.
The stdpopsim.ext.selective_sweep() represents a general model for a mutation
that is beneficial within a single population. 
'''

extended_events = stdpopsim.ext.selective_sweep(
    single_site_id=locus_id,
     # this is where the mutation starts
    population="AFR",
    selection_coeff=0.05,  # selection coefficient for the mutation [0.05 for humans]
    mutation_generation_ago=300,  # mutation originates 400 gens ago in AFR pop [5k-30k years, or 10k for one]
    min_freq_at_end=0.3  # mutation frequency at present day [look into this, set to 0.4 originally]
)

'''
Now we run the simulation using SLiM. For comparison a neutral simulation will 
also be run

slim_scaling_factor:
    Reduces population size by a scaling factor (Q) 
    Rescales time by Q 
    If model selection is used it would rescale the selection coefficient by Q 
    In general, this results in indistinguishable output. However, this may NOT be true in simulations with large 
    amount of selection. In my case, I'm assuming we do not have "large" amounts of selection since there is on a single 
    selective sweep. 
    Link --> https://popsim-consortium.github.io/stdpopsim-docs/stable/tutorial.html#slim-scaling-factor 
    
slim_burn_in:
    Number of generations before SLiM kicks in. 
    MSPRIME is used for the prior generations 
    The value of 10 is the default and considered to be "safe" 
    Link --> https://popsim-consortium.github.io/stdpopsim-docs/stable/tutorial.html#the-slim-burn-in 
    
min_freq_at_end 
    Frequency of the mutation in the present data (e.g., at the end of the sweep) should be greater than this value.  
    Link --> https://popsim-consortium.github.io/stdpopsim-docs/stable/tutorial.html#selective-sweeps 
    (scroll down to min_freq_at_end) 
'''

engine = stdpopsim.get_engine("slim")
ts_sweep = engine.simulate(
    model,
    contig,
    samples,
    seed=123,
    extended_events=extended_events,
    slim_scaling_factor=9,
    slim_burn_in=10  # number of gens before SLiM begins (these gens are msprime). 10 is safe and default
)

ts_neutral = engine.simulate(
    model,
    contig,
    samples,
    seed=123,
    # no extended events
    slim_scaling_factor=9,
    slim_burn_in=10  # number of gens before SLiM begins (these gens are msprime). 10 is safe and default
)

'''
Plot the data to compare nucleotide diversity in 10Kb windows for neutral and sweep
'''
# windows = [w for w in range(0, int(ts_neutral.sequence_length), 10000)]
# windows.append(int(ts_neutral.sequence_length))
# neutral_pi = ts_neutral.diversity(windows=windows)
# sweep_pi = ts_sweep.diversity(windows=windows)
# plt.plot(neutral_pi, "b", label="neutral")
# plt.plot(sweep_pi, "r", label="sweep")
# plt.axvline(len(neutral_pi) / 2, color="black", linestyle="dashed")
# plt.legend()
# plt.xlabel("Genomic Window")
# plt.ylabel("Diversity")
# plt.show()

'''
export this data
'''
print("Saving VCF")
# samp_ids = ts_sweep.samples()
# individuals = [ts_sweep.node(i).individual for i in samp_ids]  # diploids - so each individual takes two slots
# geno_array = np.empty((0, ts_sweep.sample_size))  # length is all SNPs and width is all individuals (diploids)
# for v in ts_sweep.variants():
#     # print(f"Site {v.site.id}: {v.genotypes}")
#     vars = v.genotypes.reshape((1, ts_sweep.sample_size))
#     geno_array = np.append(geno_array, vars, axis=0)  # transpose this

# write the simulations in vcf format
with open("output/ts_sweep_Scaling9_300x40_05.vcf", "w") as vcf_file:
    ts_sweep.write_vcf(vcf_file)

with open("output/ts_neutral_Scaling9_300x40_05.vcf", "w") as vcf_file:
    ts_neutral.write_vcf(vcf_file)

end = time.time()
delta = round((end - start) / 60, 2)
print(delta, " minutes")