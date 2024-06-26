"""
Script pulls data from the population_simulations db and formats it to be used by the SWIFr script. Main operations
that are performed include:
- Replace nulls with -998
- Query the nulls, links, and sweeps and save them as separate files (neutral, link_left, sweep)
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
test_count = 20  # number of simulations to use in the test data for swifr
vcf_names = tuple(['ts_sweep_' + str(i) + '.vcf' for i in range(test_count)])  # list of vcf files to include in test

sql1 = (f"""
       SELECT
        snp_position,
        vcf_name,
        xpehh,
        fst,
        ihs_afr_std

        FROM {table_name}
        WHERE (label IS 'sweep') AND (vcf_name NOT IN {vcf_names})
       """)

sql2 = (f"""
       SELECT
        snp_position,
        vcf_name,
        xpehh,
        fst,
        ihs_afr_std
        
        FROM {table_name}
        WHERE (label IS 'neutral' OR label IS 'link_left' OR label IS 'link_right') AND (vcf_name NOT IN {vcf_names})
       """)

sql4 = (f"""
       SELECT
        label,
        vcf_name,
        snp_position,
        xpehh,
        fst,
        ihs_afr_std

        FROM {table_name}
        WHERE vcf_name in {vcf_names}
       """)

# collect a list of the unique simulations
print("Collecting Sweeps")
sweep = pd.read_sql(sql1, conn)
print("Collecting Neutral")
neutral = pd.read_sql(sql2, conn)
neutral = neutral.replace({'label': {'link_left': 'neutral', 'link_right': 'neutral'}})  # rename linked data
print("Collecting Test")
test = pd.read_sql(sql4, conn)

# there are ~22 million rows of neutral data, for now I'm knocking this back to 3 million b/c that's just excessive
neutral = neutral.iloc[0:3000000, :]

''' Print the pct of nans in ihs data'''
nans_sweep = (sweep['ihs_afr_std'].isna().sum()) / len(sweep['ihs_afr_std'])
nans_neutral = (neutral['ihs_afr_std'].isna().sum()) / len(neutral['ihs_afr_std'])
nans_test = (test['ihs_afr_std'].isna().sum()) / len(test['ihs_afr_std'])
print("The percent of the ihs sweep data that are nans is: ", str(round(nans_sweep * 100, 2)) + "%")
print("The percent of the ihs neutral data that are nans is: ", str(round(nans_neutral * 100, 2)) + "%")
print("The percent of the ihs test data that are nans is: ", str(round(nans_test * 100, 2)) + "%")

''' Replace nans with -998 to match SWIFr null format'''
print("Replacing nan with -998")
sweep = sweep.fillna(-998)
neutral = neutral.fillna(-998)
test = test.fillna(-998)

''' Save sim stats in their respective folders '''
print("Saving Data")
neutral.to_csv('simulations_4_swifr/2_class/neutral/neutral', sep='\t', index=False)
sweep.to_csv('simulations_4_swifr/2_class/sweep/sweep', sep='\t', index=False)
test.to_csv('simulations_4_swifr/2_class/test/test', sep='\t', index=False)

# Create a training set with neutral, ll, sweep, lr in the same file (like test, but only with training data)
# this was only needed for the hmm code to use for determining pi and A_trans
neutral['label'] = 'neutral'
sweep['label'] = 'sweep'

train = pd.concat([neutral, sweep], ignore_index=True)
train = train.sort_values(by=['vcf_name', 'snp_position'])
# save the train data
train.to_csv('simulations_4_swifr/2_class/train/train', sep='\t', index=False)

print('done')
