"""
Code to query the population simulation db and visualize the neutral, sweep, and link distributions
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
table_name = 'neutral_simulations_stats'

sql = (f"""
       SELECT *
        FROM {table_name}
        WHERE label IS 'neutral'
       """)

# collect a list of the unique simulations
neutral = pd.read_sql(sql, conn)

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
plt.xlabel('iHS Value')
plt.ylabel('Frequency')
plt.title('Histogram of African Standard iHS')
plt.legend()
plt.show()

# european population
plt.figure(figsize=(12, 5))
plt.hist(neutral['ihs_afr_std'], density=True, label='neutral', bins=n_bins+500, alpha=0.7, color='green')
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
plt.xlabel('Fst Value')
plt.ylabel('Frequency')
plt.title('Histogram of Fst')
plt.legend()
plt.show()

conn.close()