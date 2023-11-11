'''
Tutorial code from getting started with tskit.
Note that msprime CAN simulate a selective sweep as shown below
'''

import msprime
import matplotlib.pyplot as plt
import numpy as np

pop_size=10_000
seq_length=10_000_000

sweep_model = msprime.SweepGenicSelection(
    position=seq_length/2, start_frequency=0.0001, end_frequency=0.9999, s=0.25, dt=1e-6)

ts = msprime.sim_ancestry(
    20,
    model=[sweep_model, msprime.StandardCoalescent()],
    population_size=pop_size,
    sequence_length=seq_length,
    recombination_rate=1e-8,
    random_seed=1234,  # only needed for repeatabilty
    )
# Optionally add finite-site mutations to the ts using the Jukes & Cantor model, creating SNPs
ts = msprime.sim_mutations(ts, rate=1e-8, random_seed=4321)
print(ts)

for tree in ts.trees():
    print(f"Tree {tree.index} covers {tree.interval}")
    if tree.index >= 4:
        print("...")
        break
print(f"Tree {ts.last().index} covers {ts.last().interval}")

kb = [0]  # Starting genomic position
mrca_time = []
for tree in ts.trees():
    kb.append(tree.interval.right/1000)  # convert to kb
    mrca = ts.node(tree.root)  # For msprime tree sequences, the root node is the MRCA
    mrca_time.append(mrca.time)
plt.stairs(mrca_time, kb, baseline=None)
plt.xlabel("Genome position (kb)")
plt.ylabel("Time of root (or MRCA) in generations")
plt.yscale("log")
plt.show()

np.set_printoptions(linewidth=200)  # print genotypes on a single line

print("Genotypes")
for v in ts.variants():
    print(f"Site {v.site.id}: {v.genotypes}")
    if v.site.id >= 4:  # only print up to site ID 4
        print("...")
        break

'''
probably the most helpful code for what I need is copied below
https://tskit.dev/tutorials/getting_started.html#sec-tskit-getting-started-exporting-data
'''

samp_ids = ts.samples()
print(" ID of diploid individual: ", " ".join([f"{ts.node(s).individual:3}" for s in samp_ids]))
print("      ID of (sample) node: ", " ".join([f"{s:3}" for s in samp_ids]))
for v in ts.variants():
    site = v.site
    alleles = np.array(v.alleles)
    print(f"Site {site.id} (ancestral state '{site.ancestral_state}')", alleles[v.genotypes])
    if site.id >= 4:
        print("...")
        break