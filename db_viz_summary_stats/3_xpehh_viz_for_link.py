"""
Script to analyze the sweep footprint based on the xp-ehh statistic. This statistic was used
due to the low occurrence of null values as wells as the distinct signature of the sweep and its footprint.

Steps applied:
1. Query the db to extract all xpehh data greater than 0 over a span of sweep_position +/- 500K bp
2. A 500K span was assumed to be far more than a typical human sweep footprint, which is approx. 100K (lauren s.)
3. Plot the distribution of snp positions from the query.
4. Fit a normal distribution to the observed snp position
5. Use the 95% area (~1.96 * sd) to define the +/- range to use for defining link events.
"""


import sqlite3
import time
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm

'''
-----------------------------------------------------------------------------------------------------------------------
Query the simulation db to determine the unique number of simulations
-----------------------------------------------------------------------------------------------------------------------
'''
start = time.time()
path = 'C:/Users/scott/PycharmProjects/population_simulations/db_build/'
conn = sqlite3.connect(path + 'population_simulation_v2.db')
# enter the stats table name (neutral or sweep)
table_name = 'sweep_simulations_stats'
sweep_pos = 2.5e6

# the sql below extracts the xpehh data between +/-200K bp away from the sweep point of 2.5e6
sql1 = (f"""
       SELECT
        vcf_name,
        uniq_id,
        snp_position,
        xpehh
        
        FROM {table_name}
        WHERE (snp_position >= 2e6) AND
            (snp_position <=  3e6) AND
            (xpehh > 0) --This restrain ensures that only xpehh values typical of a sweep-like signature are included
       """)

# collect a list of the unique simulations
sim = pd.read_sql(sql1, conn)
conn.close()

'''
-----------------------------------------------------------------------------------------------------------------------
Using the distribution of xpehh around the sweep event for all simulations
Histogram showing the distribution of XP-EHH values > 0 around the sweep location
-----------------------------------------------------------------------------------------------------------------------
'''
n_bins = 500
snp_mean = sim['snp_position'].mean()  # mean of snps with xp-ehh values greater than zero
snp_sd = sim['snp_position'].std()  # sd of snps with xp-ehh values greater than zero
x_range = np.linspace(2e6, 3e6, n_bins)  # range to use for normal distribution
norm_snp = norm.pdf(x_range, snp_mean, snp_sd)  # normal dist with mean and sd of snp_mean and snp_sd respectively
norm_ppf = norm.ppf(0.975, loc=0, scale=1)  # z-score to use for defining the 95% range (basically just 1.96)
sweep_foot = round(round(norm_ppf * snp_sd, -4),)

print("-------------------------------------------------------------------------------------------------------------")
print("The +/- range for the sweep footprint is: ", sweep_foot)
print("The above value (rounded to nearest 10k) will be used to define the linked SNPs by bracketing the sweep location with this value")
print("-------------------------------------------------------------------------------------------------------------")

plt.figure(figsize=(12, 5))
plt.hist(sim['snp_position'], bins=n_bins, density=False, label='Sweep Footprint', alpha=0.7, color='blue')
plt.plot(x_range, norm_snp/norm_snp.sum() * len(sim['snp_position']), color='red')
plt.axvline(x=sweep_pos, color='black', linestyle='--', label='Sweep Location')
plt.axvline(x=sweep_pos - norm_ppf*snp_sd, color='red', linestyle='--', label='lower 2.5%')
plt.axvline(x=sweep_pos + norm_ppf*snp_sd, color='red', linestyle='--', label='upper 97.5%')
plt.annotate('Sweep Footprint +/- ' + str(sweep_foot), xy=(2.8e6, 5000))
plt.xlabel('SNP Position')
plt.ylabel('Count of XP-EHH > 0')
plt.title('Sweep Footprint Distribution')
plt.legend()
plt.show()






