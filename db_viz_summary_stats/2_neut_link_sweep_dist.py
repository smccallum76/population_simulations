"""
Code to query the population simulation db and visualize the neutral, sweep, and link_left distributions
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
conn = sqlite3.connect(path + 'population_simulation_v3.db')
# enter the stats table name (neutral or sweep)
table_name = 'sweep_simulations_stats'

sql1 = (f"""
       SELECT *
        FROM {table_name}
        WHERE label IS 'sweep'
       """)

sql2 = (f"""
       SELECT *
        FROM {table_name}
        WHERE label IS 'link_left'
       """)

sql3 = (f"""
       SELECT *
        FROM {table_name}
        WHERE label IS 'link_right'
       """)

sql4 = (f"""
       SELECT *
        FROM {table_name}
        WHERE label IS 'neutral'
       """)

# collect a list of the unique simulations
sweep = pd.read_sql(sql1, conn)
link_left = pd.read_sql(sql2, conn)
link_right = pd.read_sql(sql3, conn)
neutral = pd.read_sql(sql4, conn)

end = time.time()
mins = round((end - start) / 60, 2)
print("Time to load unique sims [mins] = ", mins)

'''
-----------------------------------------------------------------------------------------------------------------------
HISTOGRAM - XP-EHH Distribution
-----------------------------------------------------------------------------------------------------------------------
'''
n_bins = 100
plt.figure(figsize=(12, 5))
plt.hist(neutral['xpehh'], density=True, label='neutral', bins=n_bins+500, alpha=0.7, color='green')
plt.hist(link_left['xpehh'], density=True, label='link_left', bins=n_bins+300, alpha=0.7, color='red')
plt.hist(link_right['xpehh'], density=True, label='link_right', bins=n_bins+300, alpha=0.7, color='tomato')
plt.hist(sweep['xpehh'], density=True, label='sweep', bins=n_bins, alpha=0.7, color='blue')
plt.xlabel('XP-EHH Value')
plt.ylabel('Frequency')
plt.title('Histogram of XP-EHH')
plt.legend()
plt.show()

'''
-----------------------------------------------------------------------------------------------------------------------
HISTOGRAM - iHS Distribution
-----------------------------------------------------------------------------------------------------------------------
'''
n_bins = 100

# african population
plt.figure(figsize=(12, 5))
plt.hist(neutral['ihs_afr_std'], density=True, label='neutral', bins=n_bins+500, alpha=0.7, color='green')
plt.hist(link_left['ihs_afr_std'], density=True, label='link_left', bins=n_bins+300, alpha=0.7, color='red')
plt.hist(link_right['ihs_afr_std'], density=True, label='link_right', bins=n_bins+300, alpha=0.7, color='tomato')
plt.hist(sweep['ihs_afr_std'], density=True, label='sweep', bins=n_bins-50, alpha=0.7, color='blue')
plt.xlabel('iHS Value')
plt.ylabel('Frequency')
plt.title('Histogram of African Standard iHS')
plt.legend()
plt.show()

# european population
plt.figure(figsize=(12, 5))
plt.hist(neutral['ihs_afr_std'], density=True, label='neutral', bins=n_bins+500, alpha=0.7, color='green')
plt.hist(link_left['ihs_eur_std'], density=True, label='link_left', bins=n_bins+300, alpha=0.7, color='red')
plt.hist(link_right['ihs_eur_std'], density=True, label='link_right', bins=n_bins+300, alpha=0.7, color='tomato')
plt.hist(sweep['ihs_eur_std'], density=True, label='sweep', bins=n_bins-50, alpha=0.7, color='blue')
plt.xlabel('iHS Value')
plt.ylabel('Frequency')
plt.title('Histogram of European Standard iHS')
plt.legend()
plt.show()

'''
-----------------------------------------------------------------------------------------------------------------------
HISTOGRAM - Fst Distribution
-----------------------------------------------------------------------------------------------------------------------
'''
n_bins = 100

plt.figure(figsize=(12, 5))
plt.hist(neutral['fst'], density=True, label='neutral', bins=n_bins+500, alpha=0.7, color='green')
plt.hist(link_left['fst'], density=True, label='link_left', bins=n_bins+300, alpha=0.7, color='red')
plt.hist(link_right['fst'], density=True, label='link_right', bins=n_bins+300, alpha=0.7, color='tomato')
plt.hist(sweep['fst'], density=True, label='sweep', bins=n_bins-50, alpha=0.7, color='blue')
plt.xlabel('Fst Value')
plt.ylabel('Frequency')
plt.title('Histogram of Fst')
plt.legend()
plt.show()

conn.close()