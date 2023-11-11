"""
Review stdpopsim tutorial - part 3 - SLiM model
Forward-time simulator. More time consuming, but more flexibility. The default engine is msprime, but SLiM
engine is an alternate option. Popsim provides the ability to run a simulation using msprime as the engine
and then turnaround and use SLiM to run the exact same simulation, but with a different engine. This is
pretty cool bc you can compare the backward and forward modeling results, just swap out the engine.

https://popsim-consortium.github.io/stdpopsim-docs/stable/tutorial.html
"""
# sloppy, but the 'warnings' module ignores a couple future warnings in stdpopsim
import warnings
warnings.filterwarnings("ignore")
import stdpopsim
""" 
------------------------------------------------------------------------------------------------------------------------
SPECIES, CONTIG, RECOMBINATION MAP, AND SIM ENGINE
Select the species (humans)
------------------------------------------------------------------------------------------------------------------------
"""
print("-----------------------------------------------------------------------------")
print("SPECIES, CONTIG, RECOMBINATION MAP, AND SIM ENGINE")
print("-----------------------------------------------------------------------------")
species = stdpopsim.get_species("HomSap") # same process as 01_ tutorial
model = species.get_demographic_model("Africa_1T12")
contig = species.get_contig(
    "chr22", left=10e6, right=20e6, mutation_rate=model.mutation_rate
)
# default is a flat genetic map with avg rate across chr22
samples = {"AFR":100}
# getting an error at the line below that appears to be w/in the SLiM code
engine = stdpopsim.get_engine("slim")
ts = engine.simulate(model, contig, samples, slim_scaling_factor=10)