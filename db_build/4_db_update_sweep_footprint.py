"""
This script updates the sweep_simulations and sweep_simulations_stats db to honor the sweep footprint as determined
by the analysis of xp-ehh values > 0 (see script db_viz_summary_stats/3_xpehh_viz_for_link.py).

Note --> this is not the ideal workflow. I should have not labeled the SNPs prior to reviewing the xp-ehh data. Future
updated may include eliminating the labels when the sims are initially pushed into the db.
"""

import sqlite3
import time
import pandas as pd
import allel as sa
import numpy as np
from sqlalchemy import create_engine

update_db = input("Would you like to update the DB, 'yes' or 'no': ")
# update_db = 'yes'
db_name = 'population_simulation_v3.db'
table_name = 'sweep_simulations_stats'  # select sweep_simulations_stats OR sweep_simulations
db_path = 'C:/Users/scott/PycharmProjects/population_simulations/db_build/'
sweep_footprint = 250000  # this was derived using the b_viz_summary_stats/3_xpehh_viz_for_link.py script
sweep_position = 2.5e6
link_left = sweep_position - sweep_footprint
link_right = sweep_position + sweep_footprint

if update_db == 'yes':
    engine = create_engine('sqlite:///' + db_path + db_name, echo=False)

'''
-------------------------------------------------------------------------------------------------------------------
Pull the necessary data from the db
-------------------------------------------------------------------------------------------------------------------
'''
start = time.time()
if update_db == 'yes':
    conn = sqlite3.connect(db_path + db_name)

    # define sql query based on the unique simulation name
    sql = f"""
    SELECT
        vcf_name,
        uniq_id,
        label,
        snp_position
        
    FROM {table_name}
    """
    # create df that contains the query return
    df = pd.read_sql(sql, conn)
    # create the new labels that will update the old labels
    df['new_labels'] = np.select(
                            [
                                df['snp_position'].between(link_left, sweep_position, inclusive='left'),
                                df['snp_position'] == sweep_position,
                                df['snp_position'].between(sweep_position, link_right, inclusive='right' )
                            ],
                            [
                                'link_left',
                                'sweep',
                                'link_right'
                            ],
                            default='neutral'
                        )

    '''
    -------------------------------------------------------------------------------------------------------------------
    Add the new columns to the existing table
    -------------------------------------------------------------------------------------------------------------------
    '''
    cursor = conn.cursor()
    update_sql = f"""
        UPDATE {table_name} SET label=? WHERE uniq_id=?
    """
    print('converting values to tuples')
    new_labels = np.vstack((df['new_labels'], df['uniq_id'])).transpose()  # join the two columns
    sqltuples = tuple(map(tuple, new_labels))  # convert each row of the two columns to a tuple for sql
    print('updating db: ', table_name)
    cursor.executemany(update_sql, sqltuples)
    print('commiting changes')
    # Commit the changes
    conn.commit()
    # Close the cursor and the connection
    cursor.close()
    conn.close()
    # post timing
    end = time.time()
    mins = round((end - start) / 60, 2)
    print("Time to load one batch, run all stats, and update the db = ", mins)

