'''
Script to eval and test sqlite as a db repository for vcf data
'''

"""
To Do:
- revise this code such that (say) 100 sims are queried at a time (the multiple query method is slow b/c each one takes
about 20 seconds)
- still need to add the code to update the db with the sweep_stats
- need to visualize the sweep stats data after the db update is complete
"""


import sqlite3
import time
import pandas as pd
import allel as sa
import numpy as np

start = time.time()
path = 'C:/Users/scott/PycharmProjects/population_simulations/db_build/'
conn = sqlite3.connect(path + 'population_simulation.db')

sql = ("""
       SELECT DISTINCT
        vcf_name
        FROM sweep_simulations
       """)

# collect a list of the unique simulations
uniq_sims = list(pd.read_sql(sql, conn)['vcf_name'])
# uniq_sims = ['ts_sweep_120.vcf']

end = time.time()
mins = round((end - start) / 60, 2)
print("Time to load unique sims [mins] = ", mins)

# iterate over each unique simulation and calculate the stats
start = time.time()
for sim in uniq_sims:
    print('Simulation File: ', sim)
    # define sql query based on the unique simulation name
    sql2 = f"""
    SELECT * 
    FROM sweep_simulations
    WHERE vcf_name LIKE {"'" + sim + "'"}
    """
    # create df that contains the query return
    df1 = pd.read_sql(sql2, conn)
    # split the df between the african and european samples using a regex
    afr = df1.filter(regex="afr_*")
    eur = df1.filter(regex="eur_*")
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

    """ Calculate Statistics """
    # allele counts at each position
    afr_ac = afr_ga.count_alleles()
    eur_ac = eur_ga.count_alleles()

    # xp-ehh (expects a 2d array, so the genotype array is NOT used here, use the afr and eur df's instead)
    xpehh = sa.xpehh(afr, eur, pos=df1['snp_position'], include_edges=False, use_threads=True)

    # Fst
    num, den = sa.hudson_fst(afr_ac, eur_ac)
    fst = num / den  # fst for each variant individually
    fst_overall = np.sum(num) / np.sum(den)  # fst averaging over all variants

    # iHS --> error checking required b/c this stat can crash (ts_sweep_120.vcf will break) ZeroDivisioon Error
    # unstandardized first
    try:
        ihs_afr = sa.ihs(afr, pos=df1['snp_position'], include_edges=False, use_threads=True)
    except ZeroDivisionError:
        ihs_afr = np.nan

    try:
        ihs_eur = sa.ihs(eur, pos=df1['snp_position'], include_edges=False, use_threads=True)
    except ZeroDivisionError:
        ihs_eur = np.nan

    # Standardize the iHS scores
    try:
        ihs_afr_std, _ = sa.standardize_by_allele_count(ihs_afr, afr_ac[:, 1], diagnostics=False)
    except TypeError:
        ihs_afr_std = np.nan
    try:
        ihs_eur_std, _ = sa.standardize_by_allele_count(ihs_eur, eur_ac[:, 1], diagnostics=False)
    except TypeError:
        ihs_eur_std = np.nan

    """ Build table of stats data """
    stat_tbl = pd.DataFrame()

    stat_tbl['snp_position'] = df1['snp_position']
    stat_tbl['label'] = df1['label']
    stat_tbl['simulation_type'] = 'sweep_type'
    stat_tbl['vcf_row_id'] = df1['vcf_row_id']
    stat_tbl['reference'] = df1['reference']
    stat_tbl['vcf_name'] = df1['vcf_name']
    stat_tbl['uniq_id'] = df1['uniq_id']
    stat_tbl['afr_ac_0'] = afr_ac[:, 0]
    stat_tbl['afr_ac_1'] = afr_ac[:, 1]
    stat_tbl['afr_ac_2'] = afr_ac[:, 2]
    stat_tbl['eur_ac_0'] = eur_ac[:, 0]
    stat_tbl['eur_ac_1'] = eur_ac[:, 1]
    stat_tbl['eur_ac_2'] = eur_ac[:, 2]
    stat_tbl['xpehh'] = xpehh
    stat_tbl['fst'] = fst
    stat_tbl['ihs_afr_std'] = ihs_afr_std
    stat_tbl['ihs_eur_std'] = ihs_eur_std

end = time.time()
mins = round((end - start) / 60, 2)
print("Time to load one simulation = ", mins)


print('done')

conn.close()