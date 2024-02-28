'''
Script to run through all the neutral simulations, compute statistics, and load the pop sim db
'''

import sqlite3
import time
import pandas as pd
import allel as sa
import numpy as np
from sqlalchemy import create_engine

update_db = input("Would you like to update the DB, 'yes' or 'no': ")
db_name = 'population_simulation_v3.db'
table_name = 'sweep_simulations_stats'
db_path = 'C:/Users/scott/PycharmProjects/population_simulations/db_build/'
print('----------------------------------------------------------------------------------')
print('BE SURE THAT YOU SET uniq_id AS THE PRIMARY KEY PRIOR TO RUNNING THIS SCRIPT')
print('----------------------------------------------------------------------------------')
if update_db == 'yes':
    engine = create_engine('sqlite:///' + db_path + db_name, echo=False)

'''
-----------------------------------------------------------------------------------------------------------------------
Query the simulation db and extract unstandardized iHS and allele counts
-----------------------------------------------------------------------------------------------------------------------
'''

start = time.time()
conn = sqlite3.connect(db_path + db_name)

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
    ihs_afr_std, _ = sa.standardize_by_allele_count(df['ihs_afr_unstd'], df['afr_ac_1'], diagnostics=False)
except TypeError:
    ihs_afr_std = np.nan
try:
    ihs_eur_std, _ = sa.standardize_by_allele_count(df['ihs_eur_unstd'], df['eur_ac_1'], diagnostics=False)
except TypeError:
    ihs_eur_std = np.nan

'''
-------------------------------------------------------------------------------------------------------------------
Add the new columns to the existing table
-------------------------------------------------------------------------------------------------------------------
'''

if update_db == 'yes':
    cursor = conn.cursor()
    cursor.execute(
        """
        ALTER TABLE sweep_simulations_stats 
        ADD COLUMN 'ihs_afr_std' 'float'
        """)

    cursor.execute(
        """
        ALTER TABLE sweep_simulations_stats 
        ADD COLUMN 'ihs_eur_std' 'float'
        """)
    print('converting ihs values to tuples')
    ihs_join = np.vstack((ihs_afr_std, ihs_eur_std, df['uniq_id'])).transpose()  # join the two columns
    sqltuples = tuple(map(tuple, ihs_join)) # convert each row of the two columns to a tuple for sql
    # Insert a row into the table
    # IMPORTANT --> BE SURE TO MAKE uniq_id the primary key. This can be done in DB Browser or python
    update_sql = """
        UPDATE sweep_simulations_stats SET ihs_afr_std=?, ihs_eur_std=? WHERE uniq_id=?
    """
    print('Updating the database.')
    cursor.executemany(update_sql, sqltuples)
    print('Committing db changes')
    # Commit the changes
    conn.commit()
    print('Closing cursor')
    # Close the cursor and the connection
    cursor.close()

end = time.time()
mins = round((end - start) / 60, 2)
print("Time to load one batch, run all stats, and update the db = ", mins)

conn.close()