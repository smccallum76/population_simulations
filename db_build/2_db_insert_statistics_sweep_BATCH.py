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
db_name = 'population_simulation.db'
table_name = 'sweep_simulations_stats'
db_path = 'C:/Users/scott/PycharmProjects/population_simulations/db_build/'

if update_db == 'yes':
    engine = create_engine('sqlite:///' + db_path + db_name, echo=False)


'''
-----------------------------------------------------------------------------------------------------------------------
Query the simulation db to determine the unique number of simulations
-----------------------------------------------------------------------------------------------------------------------
'''

start = time.time()
conn = sqlite3.connect(db_path + 'population_simulation.db')

sql = ("""
       SELECT DISTINCT
        vcf_name
        FROM sweep_simulations
       """)

# collect a list of the unique simulations
uniq_sims = tuple(pd.read_sql(sql, conn)['vcf_name'])
# uniq_sims = ['ts_sweep_120.vcf']

end = time.time()
mins = round((end - start) / 60, 2)
print("Time to load unique sims [mins] = ", mins)

'''
-----------------------------------------------------------------------------------------------------------------------
Set up batch index list
-----------------------------------------------------------------------------------------------------------------------
'''
# define a batch
batch_size = 50  # number of simulations to use in each batch
batch_rmdr = len(uniq_sims) % batch_size  # account for the leftovers when total sims is not evenly divided by batch
batch_count = int(len(uniq_sims) / batch_size)  # number of batches based on the batch size
if batch_rmdr > 0:  # add on more to the batch to account for the leftovers
    batch_count += 1

idx_batch_list = []  # list of tuples containing the starting and ending idx of each batch
for i in range(batch_count):
    start_idx = i * batch_size
    end_idx = (i+1) * batch_size
    if (i == batch_count - 1) and (batch_rmdr > 0):
        start_idx = i * batch_size
        end_idx = i * batch_size + batch_rmdr
    idx_batch_list.append((start_idx, end_idx))

'''
-----------------------------------------------------------------------------------------------------------------------
Iterate over each unique simulation and calculate the stats
-----------------------------------------------------------------------------------------------------------------------
'''
start = time.time()
for batch in idx_batch_list:
    print('Batch: ', batch)
    sims = uniq_sims[batch[0]:batch[1]]

    # define sql query based on the unique simulation name
    sql2 = f"""
    SELECT *
    FROM sweep_simulations
    WHERE vcf_name IN {sims}
    """
    # create df that contains the query return
    df1 = pd.read_sql(sql2, conn)
    stat_tbl_batch = pd.DataFrame()  # df to hold all the stats from a batch set, gets wiped at each batch reset
    for sim in sims:
        crnt = df1[df1['vcf_name'] == sim]  # current df containing one simulation from the batch
        crnt = crnt.reset_index(drop=True)

        # split the df between the african and european samples using a regex
        afr = crnt.filter(regex="afr_*")
        eur = crnt.filter(regex="eur_*")
        # define dimensions of afr and eur with the assumption that the last dim is 2 since diploid
        d1_afr = len(afr)
        d2_afr = int(len(afr.columns)/2)  # always even b/c diploid
        d3_afr = 2

        d1_eur = len(eur)
        d2_eur = int(len(eur.columns)/2)  # always even b/c diploid
        d3_eur = 2

        # convert the data to a genotype array (first convert to numpy, then reshape, then genotype array)
        afr_ga = sa.GenotypeArray(np.array(afr).reshape((d1_afr, d2_afr, d3_afr)))
        eur_ga = sa.GenotypeArray(np.array(eur).reshape((d1_eur, d2_eur, d3_eur)))

        '''
        -------------------------------------------------------------------------------------------------------------------
        Calculate statistics
        -------------------------------------------------------------------------------------------------------------------
        '''
        ''' --- ALLELE COUNTS --- '''
        afr_ac = afr_ga.count_alleles()
        eur_ac = eur_ga.count_alleles()

        ''' --- XP-EHH --- '''
        # xp-ehh (expects a 2d array, so the genotype array is NOT used here, use the afr and eur df's instead)
        xpehh = sa.xpehh(afr, eur, pos=crnt['snp_position'], include_edges=False, use_threads=True)

        ''' --- Fst --- '''
        num, den = sa.hudson_fst(afr_ac, eur_ac)
        fst = num / den  # fst for each variant individually
        fst_overall = np.sum(num) / np.sum(den)  # fst averaging over all variants

        ''' --- iHS --- '''
        # iHS --> error checking required b/c this stat can crash (ts_sweep_120.vcf will break) ZeroDivision Error
        # standardized first
        try:
            ihs_afr = sa.ihs(afr, pos=crnt['snp_position'], include_edges=False, use_threads=True)
        except ZeroDivisionError:
            ihs_afr = np.nan

        try:
            ihs_eur = sa.ihs(eur, pos=crnt['snp_position'], include_edges=False, use_threads=True)
        except ZeroDivisionError:
            ihs_eur = np.nan

        # # Standardize the iHS scores (this has to be completed after all the unstd scores are calculated
        # try:
        #     ihs_afr_std, _ = sa.standardize_by_allele_count(ihs_afr, afr_ac[:, 1], diagnostics=False)
        # except TypeError:
        #     ihs_afr_std = np.nan
        # try:
        #     ihs_eur_std, _ = sa.standardize_by_allele_count(ihs_eur, eur_ac[:, 1], diagnostics=False)
        # except TypeError:
        #     ihs_eur_std = np.nan

        '''
        -------------------------------------------------------------------------------------------------------------------
        BUILD TABLE OF STATS DATA
        The stat_tbl df will need to be appended for each sim ... not yet done
        -------------------------------------------------------------------------------------------------------------------
        '''
        stat_tbl = pd.DataFrame()

        stat_tbl['snp_position'] = crnt['snp_position']
        stat_tbl['label'] = crnt['label']
        stat_tbl['simulation_type'] = 'sweep_type'
        stat_tbl['vcf_row_id'] = crnt['vcf_row_id']
        stat_tbl['reference'] = crnt['reference']
        stat_tbl['vcf_name'] = crnt['vcf_name']
        stat_tbl['uniq_id'] = crnt['uniq_id']
        stat_tbl['afr_ac_0'] = afr_ac[:, 0]
        stat_tbl['afr_ac_1'] = afr_ac[:, 1]
        stat_tbl['afr_ac_2'] = afr_ac[:, 2]
        stat_tbl['eur_ac_0'] = eur_ac[:, 0]
        stat_tbl['eur_ac_1'] = eur_ac[:, 1]
        stat_tbl['eur_ac_2'] = eur_ac[:, 2]
        stat_tbl['xpehh'] = xpehh
        stat_tbl['fst'] = fst
        stat_tbl['ihs_afr_unstd'] = ihs_afr
        stat_tbl['ihs_eur_unstd'] = ihs_eur

        stat_tbl_batch = pd.concat([stat_tbl_batch, stat_tbl], ignore_index=True)
        '''
        -------------------------------------------------------------------------------------------------------------------
        ADD STATS TABLE TO THE SIMULATION DB - not yet complete
        -------------------------------------------------------------------------------------------------------------------
        '''
    if update_db == 'yes':
        stat_tbl_batch.to_sql(table_name, con= engine, if_exists='append')

    end = time.time()
    mins = round((end - start) / 60, 2)
    print("Time to load one batch, run all stats, and update the db = ", mins)

conn.close()