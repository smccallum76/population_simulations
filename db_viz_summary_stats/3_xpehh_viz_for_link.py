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
            (snp_position <=  3e6)
       """)

# collect a list of the unique simulations
sim = pd.read_sql(sql1, conn)
conn.close()

sim_pos = sim[sim['xpehh'] > 0]  # greater than zero inferred to be  sweep-like
sim_neg = sim[sim['xpehh'] <= 0]  # less than zero inferred to be neutral-like

'''
-----------------------------------------------------------------------------------------------------------------------
Analysis #1 --> Review the distribution of xpehh > zero +/-500K bp from sweep event

Using the distribution of xpehh around the sweep event for all simulations
Histogram showing the distribution of XP-EHH values > 0 around the sweep location
-----------------------------------------------------------------------------------------------------------------------
'''
n_bins = 500
snp_mean = sim_pos['snp_position'].mean()  # mean of snps with xp-ehh values greater than zero
snp_sd = sim_pos['snp_position'].std()  # sd of snps with xp-ehh values greater than zero
x_range = np.linspace(2e6, 3e6, n_bins)  # range to use for normal distribution
norm_snp = norm.pdf(x_range, snp_mean, snp_sd)  # normal dist with mean and sd of snp_mean and snp_sd respectively
norm_ppf = norm.ppf(0.975, loc=0, scale=1)  # z-score to use for defining the 95% range (basically just 1.96)
sweep_foot = round(round(norm_ppf * snp_sd, -4),)

print("-------------------------------------------------------------------------------------------------------------")
print("The +/- range for the sweep footprint is: ", sweep_foot)
print("The above value (rounded to nearest 10k) will be used to define the linked SNPs by bracketing the sweep location with this value")
print("-------------------------------------------------------------------------------------------------------------")

plt.figure(figsize=(12, 5))
plt.hist(sim_pos['snp_position'], bins=n_bins, density=False, label='Sweep Footprint', alpha=0.7, color='blue')
plt.plot(x_range, norm_snp/norm_snp.sum() * len(sim_pos['snp_position']), color='red')
plt.axvline(x=sweep_pos, color='black', linestyle='--', label='Sweep Location')
plt.axvline(x=sweep_pos - norm_ppf*snp_sd, color='red', linestyle='--', label='lower 2.5%')
plt.axvline(x=sweep_pos + norm_ppf*snp_sd, color='red', linestyle='--', label='upper 97.5%')
plt.annotate('Sweep Footprint +/- ' + str(sweep_foot), xy=(2.8e6, 5000))
plt.xlabel('SNP Position')
plt.ylabel('Count of XP-EHH > 0')
plt.title('Sweep Footprint Distribution, Method 2')
plt.legend()
plt.show()

'''
-----------------------------------------------------------------------------------------------------------------------
Analysis #2 --> Review the fraction of xpehh>0 when including xpehh<0 on a per bin basis
-----------------------------------------------------------------------------------------------------------------------
'''

hist_pos, bin_pos = np.histogram(sim_pos['snp_position'], bins=n_bins, density=False)  # hist and bins for xpehh>0
hist_neg, bin_neg = np.histogram(sim_neg['snp_position'], bins=bin_pos, density=False)  # hist and bins for xpehh<=0
xp_frc = hist_pos / (hist_pos + hist_neg)  # fraction of xpehh > 0

# xp_frc_scale = xp_frc * hist_pos  # scaling to match sweep footprint method 1
xp_frc_sd = sim_pos['snp_position'].std()   # assumed to be the same as method 1
xp_frc_mean = np.mean(bin_pos)  # this should be the same as method 1

x_range2 = np.linspace(2e6, 3e6, n_bins)  # range to use for normal distribution
norm_snp2 = norm.pdf(x_range2, xp_frc_mean, xp_frc_sd)  # normal dist with mean and sd of snp_mean and snp_sd respectively
norm_ppf2 = norm.ppf(0.975, loc=0, scale=1)  # z-score to use for defining the 95% range (basically just 1.96)
sweep_foot2 = round(round(norm_ppf2 * xp_frc_sd, -4),)

plt.figure(figsize=(12, 5))
plt.plot(bin_pos[1::], xp_frc, label='Sweep Footprint Method 2', color='dodgerblue')
# plt.plot(x_range, norm_snp2/np.sum(norm_snp2) * 125, color='red')  # scaling of 125 not robust, was done by trial and error
plt.plot(x_range, norm_snp2 * np.max(xp_frc)/np.max(norm_snp2), color='red')  # alternate scale based on max xp_frc
plt.axvline(x=sweep_pos, color='black', linestyle='--', label='Sweep Location')
plt.axvline(x=sweep_pos - norm_ppf2*xp_frc_sd, color='red', linestyle='--', label='lower 2.5%')
plt.axvline(x=sweep_pos + norm_ppf2*xp_frc_sd, color='red', linestyle='--', label='upper 97.5%')
plt.annotate('Sweep Footprint +/- ' + str(sweep_foot2), xy=(2.8e6, 0.5))
plt.xlabel('SNP Position')
plt.ylabel('Fraction of XP-EHH > 0')
plt.title('Sweep Footprint Distribution, Method 2')
plt.legend()
plt.show()

print('done')



