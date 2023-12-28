'''
Script to eval and test sqlite as a db repository for vcf data
'''

import sqlite3
from sqlalchemy import create_engine
import pandas as pd
import allel as sa
import numpy as np

'''
Note about location of the hard sweep:
-The hard sweep is located at the round(len(contig) /2)
- The contig length is 5,000,000, therefore the hard sweep is located
- at position 2,500,000 (SNP #, not row number). 
- This location should exist at every single Sweep simulation, but will not necessarily exist at every
neutral location (in most cases it will not). 
- Therefore,  the stats for the  sweep should correspond to variants/POS  == 2,500,000 in all cases)
- The linkage should be bumped x-base pairs left and right from this position
- Everything else is neutral. 
'''

sweep = sa.read_vcf("vcf_files_temp/sweep/ts_sweep_5.vcf")
gt_sweep = sa.GenotypeArray(sweep['calldata/GT'])

# if data is brought in as a numpy array and then converted to a GenotypeArray then it appears that the functionality
# of scikit-allel is available. Therefore, if the data is sent to a DB as a 2D matrix then when it is brought back into
# a script for stat analysis THEN be sure to honor the expected shape of the needed array. Specifically, keep the
# Genotype array as a 3D matrix (i.e., reshape). Copied below is a basic example of bringing in numpy data and then
# flipping it to a Genotypearray.
gt_test = np.array(gt_sweep)
gt_test2 = sa.GenotypeArray(gt_test)

# allele counts at each position
ac_afr = gt_sweep.count_alleles()
ac_test = gt_test2.count_alleles() # same answer as ac_afr




# engine = create_engine('sqlite:///test.db', echo=False)
# sweep.to_sql('sweep_test', con=engine, if_exists='append')
# conn = sqlite3.connect('test.db')
print('done')