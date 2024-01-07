'''
Script to eval and test sqlite as a db repository for vcf data
'''

import sqlite3
from sqlalchemy import create_engine
import pandas as pd
import allel as sa
import numpy as np
import os
import time

'''
To Do:
- seems to be working for a single file and db update -- use os to loop over all files --> done
- need to add labels --> done
- need to ensure that the file name is modified to honor the file being imported (not hard coded as currently set)
- need to see how duplicate rows are handled in the db...don't want to add data that already exists
- I think keeping the vcf test clean for now is the right idea, but stats will need to be added.  
- looks like a file size reduction from about 18Mb to 10 Mb. 
'''

# update_db = 'no'
update_db = input("Would you like to update the DB, 'yes' or 'no': ")
db_name = 'population_simulation.db'
table_name = 'sweep_simulations'
file_path = 'C:/Users/scott/OneDrive/pop_sim_output/sweep/'
db_path = 'C:/Users/scott/PycharmProjects/population_simulations/db_build/'

afr_samples = 100  # african samples
eur_samples = 100  # european samples
'''
-----------------------------------------------------------------------------------------------------------------------
Import data and reshape
-----------------------------------------------------------------------------------------------------------------------
'''
# import
files = os.listdir(file_path)
total_files = len(files)

if update_db == 'yes':
    engine = create_engine('sqlite:///' + db_path + db_name, echo=False)

count = 0
start = time.time()
for f in files:
    stop = time.time()
    cum_min = round((stop - start) / 60, 2)
    print("Percent Complete: ", round(count/total_files*100, 2), " | File Name: ", f, " | Total Time: ", cum_min, " minutes")
    vcf = sa.read_vcf(file_path + f)
    gt_shape = np.shape(vcf['calldata/GT'])
    gt_len = gt_shape[0]
    gt_ind = gt_shape[1]
    gt_dip = gt_shape[2]

    # the 2d shape below has the following structure:
    # -it has length equal to the genotype array length (~23K rows)
    # - the columns are structured such that we have afr_ind0_diploid0, afr_ind0_dip1, afr_ind1_dip0, afr_ind1_dip1...eur_ind99_dip0, eur_ind99_dip1
    gt_2d = np.array(vcf['calldata/GT']).reshape(gt_len, -1)

    '''
    -----------------------------------------------------------------------------------------------------------------------
    Build column names based on the number of afr and eur samples (100 of each)
    -----------------------------------------------------------------------------------------------------------------------
    '''
    # loop to construct the column names that are consistent with gt_2d. Structure = Demographic_Individual_Diploid
    all_cols = []
    for i in range(gt_ind):
        for d in range(gt_dip):
            if i < afr_samples:
                all_cols.append('afr_' + str(i) + '_' + str(d))
            else:
                all_cols.append('eur_' + str(i-afr_samples) + '_' + str(d))

    # add the genotype data to a dataframe with the new column labels
    gt_df = pd.DataFrame(gt_2d, columns=all_cols)  # diploid 0

    # add the additional columnns that are of interest
    gt_df['snp_position'] = vcf['variants/POS']
    gt_df['vcf_row_id'] = vcf['variants/ID']
    gt_df['reference'] = vcf['variants/REF']
    gt_df['variant_alt0'] = vcf['variants/ALT'][:, 0]
    gt_df['variant_alt1'] = vcf['variants/ALT'][:, 1]
    gt_df['variant_alt2'] = vcf['variants/ALT'][:, 2]
    gt_df['vcf_name'] = f
    gt_df['uniq_id'] = gt_df['vcf_name'] + "_" + gt_df['snp_position'].astype(str)

    '''
    -----------------------------------------------------------------------------------------------------------------------
    Add neutral, link, and sweep labels
    
    For this set of sims the length is 5e6 and the sweep is hard set at 2.5e6. 
    The linkage is assumed to be a fixed buffer of 500K  +/- the sweep location. This linkage range was determined based
    on the XP-EHH plot that showed the increase in the sweep covering this approximate range. 
    -----------------------------------------------------------------------------------------------------------------------
    '''

    gt_df['label'] = 'neutral'
    gt_df.loc[(gt_df['snp_position'] > 2.5e6 - 5e5) & (gt_df['snp_position'] < 2.5e6 + 5e5), 'label'] = 'link'
    gt_df.loc[gt_df['snp_position'] == 2.5e6, 'label'] = 'sweep'

    if update_db == 'yes':
        gt_df.to_sql(table_name, con= engine, if_exists='append')

    count += 1
print('done')