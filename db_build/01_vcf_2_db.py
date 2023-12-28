'''
Script to eval and test sqlite as a db repository for vcf data
'''

import sqlite3
from sqlalchemy import create_engine
import pandas as pd
import allel as sa
import numpy as np

'''
To Do:
- seems to be working for a single file and db update
- need to ensure that the file name is modified to honor the file being imported (not hard coded as currently set)
- need to see how duplicate rows are handled in the db...don't want to add data that already exists
- I think keeping the vcf test clean for now is the right idea, but stats will need to be added.  
- looks like a file size reduction from about 18Mb to 10 Mb. 
'''

afr_samples = 100  # african samples
eur_samples = 100  # european samples

# import
vcf = sa.read_vcf("vcf_files_temp/sweep/ts_sweep_5.vcf")
gt_shape = np.shape(vcf['calldata/GT'])
gt_len = gt_shape[0]
gt_ind = gt_shape[1]
gt_dip = gt_shape[2]

# the 2d shape below has the following structure:
# -it has length equal to the genotype array length (~23K rows)
# - the columns are structured such that we have afr_ind0_diploid0, afr_ind0_dip1, afr_ind1_dip0, afr_ind1_dip1...eur_ind99_dip0, eur_ind99_dip1
gt_2d = np.array(vcf['calldata/GT']).reshape(gt_len, -1)

# loop to construct the column names that are consistent with gt_2d. Structure = Demographic_Individual_Diploid
all_cols = []
for i in range(gt_ind):
    for d in range(gt_dip):
        if i < afr_samples:
            all_cols.append('afr_' + str(i) + '_' + str(d))
        else:
            all_cols.append('eur_' + str(i-afr_samples) + '_' + str(d))

# convert genotype to a dataframe for each diploid
gt_df = pd.DataFrame(gt_2d, columns=all_cols)  # diploid 0

# add the additional columnns that are of interest
gt_df['snp_position'] = vcf['variants/POS']
gt_df['vcf_row_id'] = vcf['variants/ID']
gt_df['reference'] = vcf['variants/REF']
gt_df['variant_alt0'] = vcf['variants/ALT'][:, 0]
gt_df['variant_alt1'] = vcf['variants/ALT'][:, 1]
gt_df['variant_alt2'] = vcf['variants/ALT'][:, 2]
gt_df['vcf_name'] = 'ts_sweep_5.vcf'


engine = create_engine('sqlite:///test.db', echo=False)
gt_df.to_sql('sweep_test', con=engine, if_exists='append')
# conn = sqlite3.connect('test.db')
print('done')