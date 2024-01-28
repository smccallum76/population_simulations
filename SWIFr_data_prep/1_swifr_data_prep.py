"""
Script pulls data from the population_simulations db and formats it to be used by the SWIFr script. Main operations
that are performed include:
- Replace nulls with -998
- Query the nulls, links, and sweeps and save them as separate files (neutral, link, sweep)
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

sql1 = (f"""
       SELECT
        snp_position,
        vcf_name,
        xpehh,
        fst,
        ihs_afr_std
        
        FROM {table_name}
        WHERE label IS 'sweep'
       """)

sql2 = (f"""
       SELECT
        snp_position,
        vcf_name,
        xpehh,
        fst,
        ihs_afr_std
        
        FROM {table_name}
        WHERE (label IS 'link_left') OR (label IS 'link_right')
       """)

sql3 = (f"""
       SELECT
        snp_position,
        vcf_name,
        xpehh,
        fst,
        ihs_afr_std
        
        FROM {table_name}
        WHERE label IS 'neutral'
       """)

# collect a list of the unique simulations
sweep = pd.read_sql(sql1, conn)
link = pd.read_sql(sql2, conn)
neutral = pd.read_sql(sql3, conn)

''' Print the pct of nans in ihs data'''
nans_sweep = (sweep['ihs_afr_std'].isna().sum()) / len(sweep['ihs_afr_std'])
nans_link = (link['ihs_afr_std'].isna().sum()) / len(link['ihs_afr_std'])
nans_neutral = (neutral['ihs_afr_std'].isna().sum()) / len(neutral['ihs_afr_std'])
print("The percent of the ihs sweep data that are nans is: ", str(round(nans_sweep * 100, 2)) + "%")
print("The percent of the ihs link data that are nans is: ", str(round(nans_link * 100, 2)) + "%")
print("The percent of the ihs neutral data that are nans is: ", str(round(nans_neutral * 100, 2)) + "%")

''' Replace nans with -998 to match SWIFr null format'''
link = link.fillna(-998)
sweep = sweep.fillna(-998)
neutral = neutral.fillna(-998)

''' Save sim stats in their respective folders '''
link.to_csv('simulations_4_swifr/link/link', sep='\t', index=False)
neutral.to_csv('simulations_4_swifr/neutral/neutral', sep='\t', index=False)
sweep.to_csv('simulations_4_swifr/sweep/sweep', sep='\t', index=False)

print('done')
