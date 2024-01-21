'''
Simulation of a selective sweep by a single beneficial mutation using the Out of Africa 2 population demographic
https://popsim-consortium.github.io/stdpopsim-docs/stable/catalog.html#sec_catalog_homsap_models_outofafrica_2t12
'''
# sloppy, but the 'warnings' module ignores a couple future warnings in stdpopsim
import warnings
import pandas as pd
warnings.filterwarnings("ignore")
import stdpopsim
import time

iterations = 100

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
    selection_coeff=0.05,  # selection coefficient for the mutation [0.1 for humans]
    mutation_generation_ago=300,  # mutation originates 300 gens ago in AFR pop [5k-30k years, or 10k for one]
    min_freq_at_end=0.3  # mutation frequency at present day
)

'''
Now we run the simulation using SLiM. For comparison a neutral simulation will 
also be run
'''
engine = stdpopsim.get_engine("slim")
for i in range(iterations):
    # sweep data simulations
    ts_sweep = engine.simulate(
        model,
        contig,
        samples,
        extended_events=extended_events,
        slim_scaling_factor=9,
        slim_burn_in=10
    )
    # neutral data simulations
    ts_neutral = engine.simulate(
        model,
        contig,
        samples,
        # no extended events
        slim_scaling_factor=9,
        slim_burn_in=10
    )

    '''
    export this data
    '''
    print("Saving VCF")
    # write the simulations in vcf format
    with open("output/sweep/ts_sweep_" + str(i) + ".vcf", "w") as vcf_file:
        ts_sweep.write_vcf(vcf_file)

    with open("output/neutral/ts_neutral_" + str(i) + ".vcf", "w") as vcf_file:
        ts_neutral.write_vcf(vcf_file)

    end = time.time()
    delta = round((end - start) / 60, 2)
    print("Simulation #: ", i, " out of ", iterations, " | Cumulative Time: ", delta, " minutes")