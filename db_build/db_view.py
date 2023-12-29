'''
Script to eval and test sqlite as a db repository for vcf data
'''

import sqlite3
from sqlalchemy import create_engine
import pandas as pd
import allel as sa
import numpy as np


# engine = create_engine('sqlite:///test.db', echo=False)
conn = sqlite3.connect('test.db')

sql = ("SELECT * FROM sweep_test")

df = pd.read_sql(sql, conn)

print(df.head())


print('done')