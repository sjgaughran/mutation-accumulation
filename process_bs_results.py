#!/usr/bin/python

import os
import sys
import fnmatch
import re
import pandas as pd 
import numpy as np 
from collections import defaultdict

def make_dd():
	return defaultdict(make_dd)

rxy_dict = make_dd()
rxy_hom_dict = make_dd()
rprime_dict = make_dd()
rprime_hom_dict = make_dd()
lstat_dict = make_dd()
lstat_hom_dict = make_dd()
count_ns_dict = make_dd()
count_snpeff_dict = make_dd()
count_gran_dict = make_dd()

granth_f_no = 0

# Read in the bs files

for filename in os.listdir(os.getcwd()):
	if fnmatch.fnmatch(filename, 'Lstats_results_*'):
		number = re.search('results_(.+?).txt', filename).group(1)
		df = pd.read_csv(filename, index_col='Species_pair', sep='\t')
		df_dict = df.to_dict()
		for key1 in df_dict:
			lstat_dict[key1]
			for key2 in df_dict[key1]:
				lstat_dict[key1][key2]
		for key1 in lstat_dict:
			for key2 in lstat_dict[key1]:
				lstat_dict[key1][key2][number] = df_dict[key1][key2]
	if fnmatch.fnmatch(filename, 'Lstats_Hom_results_*'):
		number = re.search('results_(.+?).txt', filename).group(1)
		df = pd.read_csv(filename, index_col='Species_pair', sep='\t')
		df_dict = df.to_dict()
		for key1 in df_dict:
			lstat_hom_dict[key1]
			for key2 in df_dict[key1]:
				lstat_hom_dict[key1][key2]
		for key1 in lstat_hom_dict:
			for key2 in lstat_hom_dict[key1]:
				lstat_hom_dict[key1][key2][number] = df_dict[key1][key2]
	if fnmatch.fnmatch(filename, 'Count_sum_N_S_results_*'):
		number = re.search('results_(.+?).txt', filename).group(1)
		df = pd.read_csv(filename, index_col='SPECIES', sep='\t')
		df_dict = df.to_dict()
		for key1 in df_dict:
			count_ns_dict[key1]
			for key2 in df_dict[key1]:
				count_ns_dict[key1][key2]
		for key1 in count_ns_dict:
			for key2 in count_ns_dict[key1]:
				count_ns_dict[key1][key2][number] = df_dict[key1][key2]
	if fnmatch.fnmatch(filename, 'Count_sum_Grantham_results_*'):
		granth_f_no += 1
		print("This is Grantham file number: " + str(granth_f_no))		
		number = re.search('results_(.+?).txt', filename).group(1)
		df = pd.read_csv(filename, index_col='SPECIES', sep='\t')
		df_dict = df.to_dict()
		for key1 in df_dict:
			count_gran_dict[key1]
			for key2 in df_dict[key1]:
				count_gran_dict[key1][key2]
		for key1 in count_gran_dict:
			for key2 in count_gran_dict[key1]:
				count_gran_dict[key1][key2][number] = df_dict[key1][key2]

# Perform bootstrapping

# number of bs files (usually 1000)
bs_f_number = granth_f_no

# number of iterations (1000)
bs_iter = range(1,1001)

# keep track of the species pairs
sp_pairs_l = []

# regular L/R statistics

bs_lstat_dict = make_dd()
bs_lstat_dict_out = make_dd()
with open("Lstats_bs.txt", "w+") as f:
	f.write("Species_pair\tL_xnoty_SYN_lower\tL_xnoty_SYN_avg\tL_xnoty_SYN_upper\tL_xnoty_NS_lower\tL_xnoty_NS_avg\tL_xnoty_NS_upper\tL_xnoty_STOP_lower\tL_xnoty_STOP_avg\tL_xnoty_STOP_upper\tL_xnoty_CON_lower\tL_xnoty_CONS_avg\tL_xnoty_CONS_upper\tL_xnoty_MODC_lower\tL_xnoty_MODC_avg\tL_xnoty_MODC_upper\tL_xnoty_MODR_lower\tL_xnoty_MODR_avg\tL_xnoty_MODR_upper\tL_xnoty_RAD_lower\tL_xnoty_RAD_avg\tL_xnoty_RAD_upper\tL_ynotx_SYN_lower\tL_ynotx_SYN_avg\tL_ynotx_SYN_upper\tL_ynotx_NS_lower\tL_ynotx_NS_avg\tL_ynotx_NS_upper\tL_ynotx_STOP_lower\tL_ynotx_STOP_avg\tL_ynotx_STOP_upper\tL_ynotx_CONS_lower\tL_ynotx_CONS_avg\tL_ynotx_CONS_upper\tL_ynotx_MODC_lower\tL_ynotx_MODC_avg\tL_ynotx_MODC_upper\tL_ynotx_MODR_lower\tL_ynotx_MODR_avg\tL_ynotx_MODR_upper\tL_ynotx_RAD_lower\tL_ynotx_RAD_avg\tL_ynotx_RAD_upper")
	for b in bs_iter:
		bs_array = np.random.choice(bs_iter, bs_f_number)
		for stat_key in lstat_dict:
			for sp_key in lstat_dict[stat_key]:
				count = 0
				for i in bs_array:
					count = count + lstat_dict[stat_key][sp_key][str(i)]
				bs_lstat_dict[stat_key][sp_key][b] = count
	for stat_key in bs_lstat_dict:
		for sp_key in bs_lstat_dict[stat_key]:
			sort_values = sorted(bs_lstat_dict[stat_key][sp_key].items(), key=lambda x: x[1])
			lower_ci = int(len(sort_values)*0.025)
			upper_ci = int(len(sort_values)*0.975)
			count = 0
			for i in range(len(sort_values)):
				count = count + sort_values[i][1]
			average = count/(len(sort_values))
			bs_lstat_dict_out[sp_key][stat_key]['lower_ci'] = sort_values[lower_ci][1]
			bs_lstat_dict_out[sp_key][stat_key]['upper_ci'] = sort_values[upper_ci][1]
			bs_lstat_dict_out[sp_key][stat_key]['avg'] = average	
	for sp_key in bs_lstat_dict_out:
		sp_pairs_l.append(sp_key)
		f.write('\n' + str(sp_key) + '\t')
		for stat_key in bs_lstat_dict_out[sp_key]:
			f.write(str(bs_lstat_dict_out[sp_key][stat_key]['lower_ci'])+'\t'+str(bs_lstat_dict_out[sp_key][stat_key]['avg'])+'\t'+str(bs_lstat_dict_out[sp_key][stat_key]['upper_ci'])+'\t')
	f.close()

bs_rxy_dict = make_dd()
bs_rxy_dict_out = make_dd()
with open("Rxy_bs.txt", "w+") as f:
	f.write("Species_pair\tRxy_SYN_lower\tRxy_SYN_avg\tRxy_SYN_upper\tRxy_NS_lower\tRxy_NS_avg\tRxy_NS_upper\tRxy_STOP_lower\tRxy_STOP_avg\tRxy_STOP_upper\tRxy_CONS_lower\tRxy_CONS_avg\tRxy_CONS_upper\tRxy_MODC_lower\tRxy_MODC_avg\tRxy_MODC_upper\tRxy_MODR_lower\tRxy_MODR_avg\tRxy_MODR_upper\tRxy_RAD_lower\tRxy_RAD_avg\tRxy_RAD_upper")
	for b in bs_iter:
		for sp_key in sp_pairs_l:
			bs_rxy_dict['Rxy_SYN'][sp_key][b] = bs_lstat_dict['L_xnoty_SYN'][sp_key][b]/bs_lstat_dict['L_ynotx_SYN'][sp_key][b]
			bs_rxy_dict['Rxy_NS'][sp_key][b] = bs_lstat_dict['L_xnoty_NS'][sp_key][b]/bs_lstat_dict['L_ynotx_NS'][sp_key][b]
			bs_rxy_dict['Rxy_STOP'][sp_key][b] = bs_lstat_dict['L_xnoty_STOP'][sp_key][b]/bs_lstat_dict['L_ynotx_STOP'][sp_key][b]
			bs_rxy_dict['Rxy_CONS'][sp_key][b] = bs_lstat_dict['L_xnoty_CONS'][sp_key][b]/bs_lstat_dict['L_ynotx_CONS'][sp_key][b]
			bs_rxy_dict['Rxy_MODC'][sp_key][b] = bs_lstat_dict['L_xnoty_MODC'][sp_key][b]/bs_lstat_dict['L_ynotx_MODC'][sp_key][b]
			bs_rxy_dict['Rxy_MODR'][sp_key][b] = bs_lstat_dict['L_xnoty_MODR'][sp_key][b]/bs_lstat_dict['L_ynotx_MODR'][sp_key][b]
			bs_rxy_dict['Rxy_RAD'][sp_key][b] = bs_lstat_dict['L_xnoty_RAD'][sp_key][b]/bs_lstat_dict['L_ynotx_RAD'][sp_key][b]
	for stat_key in bs_rxy_dict:
		for sp_key in bs_rxy_dict[stat_key]:
			sort_values = sorted(bs_rxy_dict[stat_key][sp_key].items(), key=lambda x: x[1])
			lower_ci = int(len(sort_values)*0.025)
			upper_ci = int(len(sort_values)*0.975)
			count = 0
			for i in range(len(sort_values)):
				count = count + sort_values[i][1]
			average = count/(len(sort_values))
			bs_rxy_dict_out[sp_key][stat_key]['lower_ci'] = sort_values[lower_ci][1]
			bs_rxy_dict_out[sp_key][stat_key]['upper_ci'] = sort_values[upper_ci][1]
			bs_rxy_dict_out[sp_key][stat_key]['avg'] = average	
	for sp_key in bs_rxy_dict_out:
		f.write('\n' + str(sp_key) + '\t')
		for stat_key in bs_rxy_dict_out[sp_key]:
			f.write(str(bs_rxy_dict_out[sp_key][stat_key]['lower_ci'])+'\t'+str(bs_rxy_dict_out[sp_key][stat_key]['avg'])+'\t'+str(bs_rxy_dict_out[sp_key][stat_key]['upper_ci'])+'\t')
	f.close()

bs_rprime_dict = make_dd()
bs_rprime_dict_out = make_dd()
with open("Rprime_bs.txt", "w+") as f:
	f.write("Species_pair\tRprime_N_S_lower\tRprime_N_S_avg\tRprime_N_S_upper\tRprime_STOP_S_lower\tRprime_STOP_S_avg\tRprime_STOP_S_upper\tRprime_RAD_S_lower\tRprime_RAD_S_avg\tRprime_RAD_S_upper\tRprime_RAD_CON_lower\tRprime_RAD_CON_avg\tRprime_RAD_CON_upper")
	for b in bs_iter:
		for sp_key in sp_pairs_l:
			bs_rprime_dict['Rprime_N_S'][sp_key][b] = bs_rxy_dict['Rxy_NS'][sp_key][b]/bs_rxy_dict['Rxy_SYN'][sp_key][b]
			bs_rprime_dict['Rprime_STOP_S'][sp_key][b] = bs_rxy_dict['Rxy_STOP'][sp_key][b]/bs_rxy_dict['Rxy_SYN'][sp_key][b]
			bs_rprime_dict['Rprime_RAD_S'][sp_key][b] = bs_rxy_dict['Rxy_RAD'][sp_key][b]/bs_rxy_dict['Rxy_SYN'][sp_key][b]
			bs_rprime_dict['Rprime_RAD_CON'][sp_key][b] = bs_rxy_dict['Rxy_RAD'][sp_key][b]/bs_rxy_dict['Rxy_CONS'][sp_key][b]
	for stat_key in bs_rprime_dict:
		for sp_key in bs_rprime_dict[stat_key]:
			sort_values = sorted(bs_rprime_dict[stat_key][sp_key].items(), key=lambda x: x[1])
			lower_ci = int(len(sort_values)*0.025)
			upper_ci = int(len(sort_values)*0.975)
			count = 0
			for i in range(len(sort_values)):
				count = count + sort_values[i][1]
			average = count/(len(sort_values))
			bs_rprime_dict_out[sp_key][stat_key]['lower_ci'] = sort_values[lower_ci][1]
			bs_rprime_dict_out[sp_key][stat_key]['upper_ci'] = sort_values[upper_ci][1]
			bs_rprime_dict_out[sp_key][stat_key]['avg'] = average	
	for sp_key in bs_rprime_dict_out:
		f.write('\n' + str(sp_key) + '\t')
		for stat_key in bs_rprime_dict_out[sp_key]:
			f.write(str(bs_rprime_dict_out[sp_key][stat_key]['lower_ci'])+'\t'+str(bs_rprime_dict_out[sp_key][stat_key]['avg'])+'\t'+str(bs_rprime_dict_out[sp_key][stat_key]['upper_ci'])+'\t')
	f.close()

# Homozygous state-based statistics

bs_lstat_hom_dict = make_dd()
bs_lstat_hom_dict_out = make_dd()
with open("Lstats_Hom_bs.txt", "w+") as f:
	f.write("Species_pair\tL_xnoty_Hom_SYN_lower\tL_xnoty_Hom_SYN_avg\tL_xnoty_Hom_SYN_upper\tL_xnoty_Hom_NS_lower\tL_xnoty_Hom_NS_avg\tL_xnoty_Hom_NS_upper\tL_xnoty_Hom_STOP_lower\tL_xnoty_Hom_STOP_avg\tL_xnoty_Hom_STOP_upper\tL_xnoty_Hom_CON_lower\tL_xnoty_Hom_CONS_avg\tL_xnoty_Hom_CONS_upper\tL_xnoty_Hom_MODC_lower\tL_xnoty_Hom_MODC_avg\tL_xnoty_Hom_MODC_upper\tL_xnoty_Hom_MODR_lower\tL_xnoty_Hom_MODR_avg\tL_xnoty_Hom_MODR_upper\tL_xnoty_Hom_RAD_lower\tL_xnoty_Hom_RAD_avg\tL_xnoty_Hom_RAD_upper\tL_xnoty_Hom_SYN_lower\tL_xnoty_Hom_SYN_avg\tL_xnoty_Hom_SYN_upper\tL_xnoty_Hom_NS_lower\tL_xnoty_Hom_NS_avg\tL_xnoty_Hom_NS_upper\tL_xnoty_Hom_STOP_lower\tL_xnoty_Hom_STOP_avg\tL_xnoty_Hom_STOP_upper\tL_xnoty_Hom_CONS_lower\tL_xnoty_Hom_CONS_avg\tL_xnoty_Hom_CONS_upper\tL_xnoty_Hom_MODC_lower\tL_xnoty_Hom_MODC_avg\tL_xnoty_Hom_MODC_upper\tL_xnoty_Hom_MODR_lower\tL_xnoty_Hom_MODR_avg\tL_xnoty_Hom_MODR_upper\tL_xnoty_Hom_RAD_lower\tL_xnoty_Hom_RAD_avg\tL_xnoty_Hom_RAD_upper")
	for b in bs_iter:
		bs_array = np.random.choice(bs_iter, bs_f_number)
		for stat_key in lstat_hom_dict:
			for sp_key in lstat_hom_dict[stat_key]:
				count = 0
				for i in bs_array:
					count = count + lstat_hom_dict[stat_key][sp_key][str(i)]
				bs_lstat_hom_dict[stat_key][sp_key][b] = count
	for stat_key in bs_lstat_hom_dict:
		for sp_key in bs_lstat_hom_dict[stat_key]:
			sort_values = sorted(bs_lstat_hom_dict[stat_key][sp_key].items(), key=lambda x: x[1])
			lower_ci = int(len(sort_values)*0.025)
			upper_ci = int(len(sort_values)*0.975)
			count = 0
			for i in range(len(sort_values)):
				count = count + sort_values[i][1]
			average = count/(len(sort_values))
			bs_lstat_hom_dict_out[sp_key][stat_key]['lower_ci'] = sort_values[lower_ci][1]
			bs_lstat_hom_dict_out[sp_key][stat_key]['upper_ci'] = sort_values[upper_ci][1]
			bs_lstat_hom_dict_out[sp_key][stat_key]['avg'] = average	
	for sp_key in bs_lstat_hom_dict_out:
		f.write('\n' + str(sp_key) + '\t')
		for stat_key in bs_lstat_hom_dict_out[sp_key]:
			f.write(str(bs_lstat_hom_dict_out[sp_key][stat_key]['lower_ci'])+'\t'+str(bs_lstat_hom_dict_out[sp_key][stat_key]['avg'])+'\t'+str(bs_lstat_hom_dict_out[sp_key][stat_key]['upper_ci'])+'\t')
	f.close()


bs_rxy_hom_dict = make_dd()
bs_rxy_hom_dict_out = make_dd()
with open("Rxy_Hom_bs.txt", "w+") as f:
	f.write("Species_pair\tRxy_Hom_SYN_lower\tRxy_Hom_SYN_avg\tRxy_Hom_SYN_upper\tRxy_Hom_NS_lower\tRxy_Hom_NS_avg\tRxy_Hom_NS_upper\tRxy_Hom_STOP_lower\tRxy_Hom_STOP_avg\tRxy_Hom_STOP_upper\tRxy_Hom_CONS_lower\tRxy_Hom_CONS_avg\tRxy_Hom_CONS_upper\tRxy_Hom_MODC_lower\tRxy_Hom_MODC_avg\tRxy_Hom_MODC_upper\tRxy_Hom_MODR_lower\tRxy_Hom_MODR_avg\tRxy_Hom_MODR_upper\tRxy_Hom_RAD_lower\tRxy_Hom_RAD_avg\tRxy_Hom_RAD_upper")
	for b in bs_iter:
		for sp_key in sp_pairs_l:
			bs_rxy_hom_dict['Rxy_Hom_SYN'][sp_key][b] = bs_lstat_hom_dict['L_xnoty_Hom_SYN'][sp_key][b]/bs_lstat_hom_dict['L_ynotx_Hom_SYN'][sp_key][b]
			bs_rxy_hom_dict['Rxy_Hom_NS'][sp_key][b] = bs_lstat_hom_dict['L_xnoty_Hom_NS'][sp_key][b]/bs_lstat_hom_dict['L_ynotx_Hom_NS'][sp_key][b]
			bs_rxy_hom_dict['Rxy_Hom_STOP'][sp_key][b] = bs_lstat_hom_dict['L_xnoty_Hom_STOP'][sp_key][b]/bs_lstat_hom_dict['L_ynotx_Hom_STOP'][sp_key][b]
			bs_rxy_hom_dict['Rxy_Hom_CONS'][sp_key][b] = bs_lstat_hom_dict['L_xnoty_Hom_CONS'][sp_key][b]/bs_lstat_hom_dict['L_ynotx_Hom_CONS'][sp_key][b]
			bs_rxy_hom_dict['Rxy_Hom_MODC'][sp_key][b] = bs_lstat_hom_dict['L_xnoty_Hom_MODC'][sp_key][b]/bs_lstat_hom_dict['L_ynotx_Hom_MODC'][sp_key][b]
			bs_rxy_hom_dict['Rxy_Hom_MODR'][sp_key][b] = bs_lstat_hom_dict['L_xnoty_Hom_MODR'][sp_key][b]/bs_lstat_hom_dict['L_ynotx_Hom_MODR'][sp_key][b]
			bs_rxy_hom_dict['Rxy_Hom_RAD'][sp_key][b] = bs_lstat_hom_dict['L_xnoty_Hom_RAD'][sp_key][b]/bs_lstat_hom_dict['L_ynotx_Hom_RAD'][sp_key][b]
	for stat_key in bs_rxy_hom_dict:
		for sp_key in bs_rxy_hom_dict[stat_key]:
			sort_values = sorted(bs_rxy_hom_dict[stat_key][sp_key].items(), key=lambda x: x[1])
			lower_ci = int(len(sort_values)*0.025)
			upper_ci = int(len(sort_values)*0.975)
			count = 0
			for i in range(len(sort_values)):
				count = count + sort_values[i][1]
			average = count/(len(sort_values))
			bs_rxy_hom_dict_out[sp_key][stat_key]['lower_ci'] = sort_values[lower_ci][1]
			bs_rxy_hom_dict_out[sp_key][stat_key]['upper_ci'] = sort_values[upper_ci][1]
			bs_rxy_hom_dict_out[sp_key][stat_key]['avg'] = average	
	for sp_key in bs_rxy_hom_dict_out:
		f.write('\n' + str(sp_key) + '\t')
		for stat_key in bs_rxy_hom_dict_out[sp_key]:
			f.write(str(bs_rxy_hom_dict_out[sp_key][stat_key]['lower_ci'])+'\t'+str(bs_rxy_hom_dict_out[sp_key][stat_key]['avg'])+'\t'+str(bs_rxy_hom_dict_out[sp_key][stat_key]['upper_ci'])+'\t')
	f.close()




bs_rprime_hom_dict = make_dd()
bs_rprime_hom_dict_out = make_dd()
with open("Rprime_Hom_bs.txt", "w+") as f:
	f.write("Species_pair\tRprime_Hom_N_S_lower\tRprime_Hom_N_S_avg\tRprime_Hom_N_S_upper\tRprime_Hom_STOP_S_lower\tRprime_Hom_STOP_S_avg\tRprime_Hom_STOP_S_upper\tRprime_Hom_RAD_S_lower\tRprime_Hom_RAD_S_avg\tRprime_Hom_RAD_S_upper\tRprime_Hom_RAD_CON_lower\tRprime_Hom_RAD_CON_avg\tRprime_Hom_RAD_CON_upper")
	for b in bs_iter:
		for sp_key in sp_pairs_l:
			bs_rprime_hom_dict['Rprime_Hom_N_S'][sp_key][b] = bs_rxy_hom_dict['Rxy_Hom_NS'][sp_key][b]/bs_rxy_hom_dict['Rxy_Hom_SYN'][sp_key][b]
			bs_rprime_hom_dict['Rprime_Hom_STOP_S'][sp_key][b] = bs_rxy_hom_dict['Rxy_Hom_STOP'][sp_key][b]/bs_rxy_hom_dict['Rxy_Hom_SYN'][sp_key][b]
			bs_rprime_hom_dict['Rprime_Hom_RAD_S'][sp_key][b] = bs_rxy_hom_dict['Rxy_Hom_RAD'][sp_key][b]/bs_rxy_hom_dict['Rxy_Hom_SYN'][sp_key][b]
			bs_rprime_hom_dict['Rprime_Hom_RAD_CON'][sp_key][b] = bs_rxy_hom_dict['Rxy_Hom_RAD'][sp_key][b]/bs_rxy_hom_dict['Rxy_Hom_CONS'][sp_key][b]
	for stat_key in bs_rprime_hom_dict:
		for sp_key in bs_rprime_hom_dict[stat_key]:
			sort_values = sorted(bs_rprime_hom_dict[stat_key][sp_key].items(), key=lambda x: x[1])
			lower_ci = int(len(sort_values)*0.025)
			upper_ci = int(len(sort_values)*0.975)
			count = 0
			for i in range(len(sort_values)):
				count = count + sort_values[i][1]
			average = count/(len(sort_values))
			bs_rprime_hom_dict_out[sp_key][stat_key]['lower_ci'] = sort_values[lower_ci][1]
			bs_rprime_hom_dict_out[sp_key][stat_key]['upper_ci'] = sort_values[upper_ci][1]
			bs_rprime_hom_dict_out[sp_key][stat_key]['avg'] = average	
	for sp_key in bs_rprime_hom_dict_out:
		f.write('\n' + str(sp_key) + '\t')
		for stat_key in bs_rprime_hom_dict_out[sp_key]:
			f.write(str(bs_rprime_hom_dict_out[sp_key][stat_key]['lower_ci'])+'\t'+str(bs_rprime_hom_dict_out[sp_key][stat_key]['avg'])+'\t'+str(bs_rprime_hom_dict_out[sp_key][stat_key]['upper_ci'])+'\t')
	f.close()


# Count-based statistics 

bs_count_ns_dict = make_dd()
bs_count_ns_dict_out = make_dd()
with open("Count_sum_N_S_bs.txt", "w+") as f:
	f.write("SPECIES\tGRANTHAM_SUM_LOWER\tGRANTHAM_SUM_AVG\tGRANTHAM_SUM_UPPER\tCOUNT_HOMO_DER_LOWER\tCOUNT_HOMO_DER_AVG\tCOUNT_HOMO_DER_UPPER\tCOUNT_HET_DER_LOWER\tCOUNT_HET_DER_AVG\tCOUNT_HET_DER_UPPER\tCOUNT_ALL_DER_LOWER\tCOUNT_ALL_DER_AVG\tCOUNT_ALL_DER_UPPER\tCOUNT_HOMO_SYN_LOWER\tCOUNT_HOMO_SYN_AVG\tCOUNT_HOMO_SYN_UPPER\tCOUNT_HET_SYN_LOWER\tCOUNT_HET_SYN_AVG\tCOUNT_HET_SYN_UPPER\tCOUNT_ALL_SYN_LOWER\tCOUNT_ALL_SYN_AVG\tCOUNT_ALL_SYN_UPPER\tCOUNT_HOMO_NS_LOWER\tCOUNT_HOMO_NS_AVG\tCOUNT_HOMO_NS_LOWER\tCOUNT_HET_NS_LOWER\tCOUNT_HET_NS_AVG\tCOUNT_HET_NS_UPPER\tCOUNT_ALL_NS_LOWER\tCOUNT_ALL_NS_AVG\tCOUNT_ALL_NS_UPPER")
	for b in bs_iter:
		bs_array = np.random.choice(bs_iter, bs_f_number)
		for stat_key in count_ns_dict:
			for sp_key in count_ns_dict[stat_key]:
				count = 0
				for i in bs_array:
					count = count + count_ns_dict[stat_key][sp_key][str(i)]
				bs_count_ns_dict[stat_key][sp_key][b] = count
	for stat_key in bs_count_ns_dict:
		for sp_key in bs_count_ns_dict[stat_key]:
			sort_values = sorted(bs_count_ns_dict[stat_key][sp_key].items(), key=lambda x: x[1])
			lower_ci = int(len(sort_values)*0.025)
			upper_ci = int(len(sort_values)*0.975)
			count = 0
			for i in range(len(sort_values)):
				count = count + sort_values[i][1]
			average = count/(len(sort_values))
			bs_count_ns_dict_out[sp_key][stat_key]['lower_ci'] = sort_values[lower_ci][1]
			bs_count_ns_dict_out[sp_key][stat_key]['upper_ci'] = sort_values[upper_ci][1]
			bs_count_ns_dict_out[sp_key][stat_key]['avg'] = average	
	for sp_key in bs_count_ns_dict_out:
		f.write('\n' + str(sp_key) + '\t')
		for stat_key in bs_count_ns_dict_out[sp_key]:
			f.write(str(bs_count_ns_dict_out[sp_key][stat_key]['lower_ci'])+'\t'+str(bs_count_ns_dict_out[sp_key][stat_key]['avg'])+'\t'+str(bs_count_ns_dict_out[sp_key][stat_key]['upper_ci'])+'\t')
	f.close()

bs_count_snpeff_dict = make_dd()
bs_count_snpeff_dict_out = make_dd()
with open("Count_sum_SNPeff_bs.txt", "w+") as f:
	f.write("SPECIES\tGRANTHAM_SUM_LOWER\tGRANTHAM_SUM_AVG\tGRANTHAM_SUM_UPPER\tCOUNT_HOMO_STOP_LOWER\tCOUNT_HOMO_STOP_AVG\tCOUNT_HOMO_STOP_UPPER\tCOUNT_HET_STOP_LOWER\tCOUNT_HET_STOP_AVG\tCOUNT_HET_STOP_UPPER\tCOUNT_ALL_STOP_LOWER\tCOUNT_ALL_STOP_AVG\tCOUNT_ALL_STOP_UPPER")
	for b in bs_iter:
		bs_array = np.random.choice(bs_iter, bs_f_number)
		for stat_key in count_snpeff_dict:
			for sp_key in count_snpeff_dict[stat_key]:
				count = 0
				for i in bs_array:
					count = count + count_snpeff_dict[stat_key][sp_key][str(i)]
				bs_count_snpeff_dict[stat_key][sp_key][b] = count
	for stat_key in bs_count_snpeff_dict:
		for sp_key in bs_count_snpeff_dict[stat_key]:
			sort_values = sorted(bs_count_snpeff_dict[stat_key][sp_key].items(), key=lambda x: x[1])
			lower_ci = int(len(sort_values)*0.025)
			upper_ci = int(len(sort_values)*0.975)
			count = 0
			for i in range(len(sort_values)):
				count = count + sort_values[i][1]
			average = count/(len(sort_values))
			bs_count_snpeff_dict_out[sp_key][stat_key]['lower_ci'] = sort_values[lower_ci][1]
			bs_count_snpeff_dict_out[sp_key][stat_key]['upper_ci'] = sort_values[upper_ci][1]
			bs_count_snpeff_dict_out[sp_key][stat_key]['avg'] = average	
	for sp_key in bs_count_snpeff_dict_out:
		f.write('\n' + str(sp_key) + '\t')
		for stat_key in bs_count_snpeff_dict_out[sp_key]:
			f.write(str(bs_count_snpeff_dict_out[sp_key][stat_key]['lower_ci'])+'\t'+str(bs_count_snpeff_dict_out[sp_key][stat_key]['avg'])+'\t'+str(bs_count_snpeff_dict_out[sp_key][stat_key]['upper_ci'])+'\t')
	f.close()

bs_count_gran_dict = make_dd()
bs_count_gran_dict_out = make_dd()
with open("Count_sum_Grantham_bs.txt", "w+") as f:
	f.write("SPECIES\tGRANTHAM_SUM_LOWER\tGRANTHAM_SUM_AVG\tGRANTHAM_SUM_UPPER\tCOUNT_HOMO_CONS_LOWER\tCOUNT_HOMO_CONS_AVG\tCOUNT_HOMO_CONS_UPPER\tCOUNT_HET_CONS_LOWER\tCOUNT_HET_CONS_AVG\tCOUNT_HET_CONS_UPPER\tCOUNT_ALL_CONS_LOWER\tCOUNT_ALL_CONS_AVG\tCOUNT_ALL_CONS_UPPER\tCOUNT_HOMO_MODC_LOWER\tCOUNT_HOMO_MODC_AVG\tCOUNT_HOMO_MODC_UPPER\tCOUNT_HET_MODC_LOWER\tCOUNT_HET_MODC_AVG\tCOUNT_HET_MODC_UPPER\tCOUNT_ALL_MODC_LOWER\tCOUNT_ALL_MODC_AVG\tCOUNT_ALL_MODC_UPPER\tCOUNT_HOMO_MODR_LOWER\tCOUNT_HOMO_MODR_AVG\tCOUNT_HOMO_MODR_UPPER\tCOUNT_HET_MODR_LOWER\tCOUNT_HET_MODR_AVG\tCOUNT_HET_MODR_UPPER\tCOUNT_ALL_MODR_LOWER\tCOUNT_ALL_MODR_AVG\tCOUNT_ALL_MOD_UPPER\tCOUNT_HOMO_RAD_LOWER\tCOUNT_HOMO_RAD_AVG\tCOUNT_HOMO_RAD_UPPER\tCOUNT_HET_RAD_LOWER\tCOUNT_HET_RAD_AVG\tCOUNT_HET_RAD_UPPER\tCOUNT_ALL_RAD_LOWER\tCOUNT_ALL_RAD_AVG\tCOUNT_ALL_RAD_UPPER")
	for b in bs_iter:
		bs_array = np.random.choice(bs_iter, bs_f_number)
		for stat_key in count_gran_dict:
			for sp_key in count_gran_dict[stat_key]:
				count = 0
				for i in bs_array:
					count = count + count_gran_dict[stat_key][sp_key][str(i)]
				bs_count_gran_dict[stat_key][sp_key][b] = count
	for stat_key in bs_count_gran_dict:
		for sp_key in bs_count_gran_dict[stat_key]:
			sort_values = sorted(bs_count_gran_dict[stat_key][sp_key].items(), key=lambda x: x[1])
			lower_ci = int(len(sort_values)*0.025)
			upper_ci = int(len(sort_values)*0.975)
			count = 0
			for i in range(len(sort_values)):
				count = count + sort_values[i][1]
			average = count/(len(sort_values))
			bs_count_gran_dict_out[sp_key][stat_key]['lower_ci'] = sort_values[lower_ci][1]
			bs_count_gran_dict_out[sp_key][stat_key]['upper_ci'] = sort_values[upper_ci][1]
			bs_count_gran_dict_out[sp_key][stat_key]['avg'] = average	
	for sp_key in bs_count_gran_dict_out:
		f.write('\n' + str(sp_key) + '\t')
		for stat_key in bs_count_gran_dict_out[sp_key]:
			f.write(str(bs_count_gran_dict_out[sp_key][stat_key]['lower_ci'])+'\t'+str(bs_count_gran_dict_out[sp_key][stat_key]['avg'])+'\t'+str(bs_count_gran_dict_out[sp_key][stat_key]['upper_ci'])+'\t')
	f.close()
