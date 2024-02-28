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

update_db = 'no'
# update_db = input("Would you like to update the DB, 'yes' or 'no': ")
db_name = 'population_simulation_v3.db'
table_name = 'neutral_simulations'
file_path = 'C:/Users/scott/OneDrive/pop_sim_output/neutral_01222024_active/'
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

count=0
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

    # script below is to get rid of monomorphic SNP sites
    gt_2d_width = np.shape(gt_2d)[1]  # get the width of gt_2d, but will be 400 (200 ind with 2 dips)
    # mask the data that does not sum to 0 or 400 (width of gt_2d)
    mask = (np.sum(gt_2d, axis=1) != 0) & (np.sum(gt_2d, axis=1) != gt_2d_width)
    gt_df = gt_df.iloc[mask, :].reset_index(drop=True)  # apply mask

    '''
    -----------------------------------------------------------------------------------------------------------------------
    Add neutral, link, and sweep labels
    
    For the neutral simulations there was not a sweep introduced, therefore all rows are labeled as neutral. 
    -----------------------------------------------------------------------------------------------------------------------
    '''
    gt_df['label'] = 'neutral'

    if update_db == 'yes':
        gt_df.to_sql(table_name, con=engine, if_exists='append')

    count += 1
print('done')