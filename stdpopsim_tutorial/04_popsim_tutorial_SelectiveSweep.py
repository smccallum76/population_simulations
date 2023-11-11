'''
Simulation of a selective sweep by a single beneficial mutation
'''
# sloppy, but the 'warnings' module ignores a couple future warnings in stdpopsim
import warnings
import pandas as pd
warnings.filterwarnings("ignore")
import stdpopsim
import matplotlib.pyplot as plt
import numpy as np

# species = stdpopsim.get_species("DroMel")  # Drosophila melanogaster
species = stdpopsim.get_species("HomSap")  # Homo sapiens
# model = stdpopsim.PiecewiseConstantSize(100000)
model = species.get_demographic_model("OutOfAfrica_2T12")
# samples = {"pop_0": 50}
samples = {"AFR":100, "EUR":100}
# the contig is fully neutral, the sweeping mutation will be inserted later
# contig = species.get_contig("2L", right=1e6)
contig = species.get_contig("chr22", right=1e6) # look up chr22 to learn what it implies

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
    # population="pop_0",  # this is where the mutation starts
    population="AFR",
    selection_coeff=0.1,  # selection coefficient for the mutation [0.1 for humans]
    mutation_generation_ago=1000,  # mutation originates 1000 gens ago in pop_0 [5k-30k years, or 10k for one]
    min_freq_at_end=0.8  # mutation frequency at present day [look into this]
)

'''
Now we run the simulation using SLiM. For comparison a neutral simulation will 
also be run
'''

engine = stdpopsim.get_engine("slim")
ts_sweep = engine.simulate(
    model,
    contig,
    samples,
    seed=123,
    extended_events=extended_events,
    slim_scaling_factor=10,
    slim_burn_in=0.1
)

ts_neutral = engine.simulate(
    model,
    contig,
    samples,
    seed=123,
    # no extended events
    slim_scaling_factor=10,
    slim_burn_in=0.1
)

'''
Plot the data to compare nucleotide diversity in 10Kb windows for neutral and sweep
'''
windows = [w for w in range(0, int(ts_neutral.sequence_length), 10000)]
windows.append(int(ts_neutral.sequence_length))
neutral_pi = ts_neutral.diversity(windows=windows)
sweep_pi = ts_sweep.diversity(windows=windows)
plt.plot(neutral_pi, "b", label="neutral")
plt.plot(sweep_pi, "r", label="sweep")
plt.axvline(len(neutral_pi) / 2, color="black", linestyle="dashed")
plt.legend()
plt.xlabel("Genomic Window")
plt.ylabel("Diversity")
plt.show()

# geno_df = pd.DataFrame(columns=['site', ])
# samp_ids = ts_sweep.samples()
# print(" ID of diploid individual: ", " ".join([f"{ts_sweep.node(s).individual:3}" for s in samp_ids]))
# print("      ID of (sample) node: ", " ".join([f"{s:3}" for s in samp_ids]))
# for v in ts_sweep.variants():
#     site = v.site
#     alleles = np.array(v.alleles)
#     print(f"Site {site.id} (ancestral state '{site.ancestral_state}')", alleles[v.genotypes])
#     if site.id >= 4:
#         print("...")
#         break
#
# count = 0
# for v in ts_sweep.variants():
#     site = v.site
#     individual = ts_sweep.node(count).individual
#     alleles = np.array(v.alleles)
#     print(site.id, " : ", individual, " : ", alleles[v.genotypes])
#     count+=1
#     if site.id >= 4:
#         print("...")
#         break

'''
export this data
'''
print("Genotypes")
samp_ids = ts_sweep.samples()
individuals = [ts_sweep.node(i).individual for i in samp_ids]  # diploids - so each individual takes two slots
geno_array = np.empty((0, ts_sweep.sample_size))  # length is all SNPs and width is all individuals (diploids)
for v in ts_sweep.variants():
    # print(f"Site {v.site.id}: {v.genotypes}")
    vars = v.genotypes.reshape((1, ts_sweep.sample_size))
    geno_array = np.append(geno_array, vars, axis=0)  # transpose this

    # if v.site.id >= 4:  # only print up to site ID 4
    #     print("...")
    #     break

print('done')