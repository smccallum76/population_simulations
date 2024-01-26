"""
Code to query the population simulation db and visualize the various statistics. This code is intended as  a QC
step only. It can be used to verify that the sweep location is in fact showing a sweep signature.
"""

import sqlite3
import time
import pandas as pd
import matplotlib.pyplot as plt

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
# enter the simulation to evaluate, must be neutral if table is neutral or sweep if the table is sweep. Changing
# the number at the end of the string extracts a new simulation.
sim_select = 'ts_sweep_14.vcf'  # replace number with value from 0 to 999 inclusive

sql = (f"""
       SELECT *
        FROM {table_name}
        WHERE vcf_name IS {"'" + sim_select + "'"}
       """)
# collect a list of the unique simulations
sim = pd.read_sql(sql, conn)

end = time.time()
mins = round((end - start) / 60, 2)
print("Time to load unique sims [mins] = ", mins)

sweep_pos = 2.5e6
x_min = 1.5e6
x_max = 3.5e6
x_range = sim.index[(sim['snp_position'] >= x_min) &
                    (sim['snp_position'] <= x_max)].tolist()

'''
-----------------------------------------------------------------------------------------------------------------------
PLOTS - Allele Frequency
-----------------------------------------------------------------------------------------------------------------------
'''
allele_pos = str(1)  # 0, 1, or 2

fig, ax = plt.subplots(nrows=1, ncols=1, sharey=False, sharex=True, figsize=(18, 8))
ax.plot(sim.loc[x_range, 'snp_position'], sim.loc[x_range, 'afr_ac_'+ allele_pos], label="AFR ac Sweep", alpha=0.7)
ax.plot(sim.loc[x_range, 'snp_position'], sim.loc[x_range, 'eur_ac_'+ allele_pos]*-1, label="EUR ac Sweep", alpha=0.7, color='red')
ax.axvline(x=sweep_pos, color='black', linestyle='--', label='Sweep Location')
ax.set_xlabel("SNP Number")
ax.set_title("Sweep Sims -- Allele Frequency by SNP")
ax.legend(loc="upper left")
plt.show()

'''
-----------------------------------------------------------------------------------------------------------------------
PLOTS - XP-EHH
-----------------------------------------------------------------------------------------------------------------------
'''

plt.figure(figsize=(15, 3))
plt.plot(sim['snp_position'], sim['xpehh'], label="AFR/EUR Sweep", alpha=0.7)
plt.axhline(y=0, color='black')
plt.axvline(x=sweep_pos, color='black', linestyle='--', label='Sweep Location')
plt.xlabel("Genomic Position")
plt.ylabel("XP-EHH")
plt.title("XP-EHH using AFR and EUR populations")
plt.legend(loc="upper left")
plt.show()


'''
-----------------------------------------------------------------------------------------------------------------------
PLOTS - iHS
-----------------------------------------------------------------------------------------------------------------------
'''
# check to see if all the values are 'None' and if so, do not plot
check_none_afr = sim['ihs_afr_std'].isnull().all()
check_none_eur = sim['ihs_eur_std'].isnull().all()

fig, ax = plt.subplots(nrows=1, ncols=1, sharey=True, figsize=(15, 6))
if check_none_afr == False:
       ax.scatter(sim['snp_position'], abs(sim['ihs_afr_std']), label='AFR Sweep', s=3, color='blue')
if check_none_eur == False:
       ax.scatter(sim['snp_position'], abs(sim['ihs_eur_std'])*-1, label='EUR Sweep', s=3, color='red')
ax.axvline(x=sweep_pos, color='black', linestyle='--', label="Sweep Location")
ax.set_ylabel("|iHS| or |iHS|*-1")
ax.set_title("Sweep Sims -- iHS using AFR and EUR populations")
ax.legend(loc="upper left")
plt.show()

'''
-----------------------------------------------------------------------------------------------------------------------
PLOTS - Fst
-----------------------------------------------------------------------------------------------------------------------
'''

plt.figure(figsize=(15, 5))
plt.plot(sim['snp_position'], sim['fst'], label="Fst")
plt.axvline(x=sweep_pos, color='black', linestyle='--', label='Sweep Location')
plt.xlabel("Genomic Position")
plt.ylabel("Fst")
plt.title("Fst using AFR and EUR populations")
plt.legend(loc="upper left")
plt.show()

print('done')