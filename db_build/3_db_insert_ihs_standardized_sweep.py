'''
Script to run through all the neutral simulations, compute statistics, and load the pop sim db
'''

import sqlite3
import time
import pandas as pd
import allel as sa
import numpy as np
from sqlalchemy import create_engine

# update_db = input("Would you like to update the DB, 'yes' or 'no': ")
db_name = 'population_simulation.db'
table_name = 'sweep_simulations_stats'
db_path = 'C:/Users/scott/PycharmProjects/population_simulations/db_build/'

# if update_db == 'yes':
engine = create_engine('sqlite:///' + db_path + db_name, echo=False)

'''
-----------------------------------------------------------------------------------------------------------------------
Query the simulation db and extract unstandardized iHS and allele counts
-----------------------------------------------------------------------------------------------------------------------
'''

start = time.time()
conn = sqlite3.connect(db_path + 'population_simulation.db')

# define sql query based on the unique simulation name
sql2 = """
SELECT
    uniq_id,
    afr_ac_1,
    eur_ac_1,
    ihs_afr_unstd,
    ihs_eur_unstd

FROM sweep_simulations_stats
"""
# create df that contains the query return
df = pd.read_sql(sql2, conn)

'''
-------------------------------------------------------------------------------------------------------------------
Calculate iHS standardized statistics
-------------------------------------------------------------------------------------------------------------------
'''

''' --- iHS --- '''

# Standardize the iHS scores (this has to be completed after all the unstd scores are calculated
try:
    df['ihs_afr_std'], _ = sa.standardize_by_allele_count(df['ihs_afr_unstd'], df['afr_ac_1'], diagnostics=False)
except TypeError:
    df['ihs_afr_std'] = np.nan
try:
    df['ihs_eur_std'], _ = sa.standardize_by_allele_count(df['ihs_eur_unstd'], df['eur_ac_1'], diagnostics=False)
except TypeError:
    df['ihs_eur_std'] = np.nan


'''
-------------------------------------------------------------------------------------------------------------------
ADD STATS TABLE TO THE SIMULATION DB - not yet complete
-------------------------------------------------------------------------------------------------------------------
'''
# if update_db == 'yes':
#     pass
#     # stat_tbl_batch.to_sql(table_name, con= engine, if_exists='append')

engine.execute('alter table sweep_simulations_stats add column column_name String')

end = time.time()
mins = round((end - start) / 60, 2)
print("Time to load one batch, run all stats, and update the db = ", mins)

conn.close()