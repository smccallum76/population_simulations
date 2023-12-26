'''
Script to eval and test sqlite as a db repository for vcf data
'''

import sqlite3
from sqlalchemy import create_engine
import pandas as pd
import allel as sa

'''
Note about location of the hard sweep:
-The hard sweep is located at the round(len(contig) /2)
- The contig length is 5,000,000, therefore the hard sweep is located
- at position 2,500,000 (SNP #, not row number). 
- This location should exist at every single Sweep simulation, but will not necessarily exist at every
neutral location (in most cases it will not). 
- Therefore,  the stats for the  sweep should correspond to variants/POS  == 2,500,000 in all cases)
- The linkage should be bumpded x-base pairs left and right from this position
- Everything else is neutral. 
'''

sweep = sa.read_vcf("vcf_files_temp/sweep/ts_sweep_5.vcf")

# engine = sqlalchemy.create_engine('sqlite:///test.db', echo=False)
# df = pd.DataFrame([[1,2],[1,2]], columns=['a', 'b'])
#
# df.to_sql('mytable', con=engine, if_exists='append')
#
# conn = sqlite3.connect('test.db')
print('done')