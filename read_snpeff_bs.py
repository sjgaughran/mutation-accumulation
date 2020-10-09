#!/usr/bin/python

import sys
import os
import re
import traceback
import itertools
import csv
import math
import vcf
from operator import itemgetter
from collections import defaultdict

if len(sys.argv) < 2:
	print("\nUSE: python program.py snpeff.vcf pop_file.txt outgroup_sample_name \nExiting--please try again! \n")
	sys.exit()

vcf_name = sys.argv[1]
vcf_reader = vcf.Reader(open(vcf_name, 'r'))
bs_number = re.search('bs_(.+?).vcf', vcf_name).group(1)

pop_f_name = sys.argv[2]
pop_dict = defaultdict(list)
pop_l_dict = defaultdict(list)
pop_names = []

with open(pop_f_name, 'r') as pops:
	for line in pops:
		pop_names.append(line.strip().split(':')[0])
		pop_l_dict[line.strip().split(':')[0]] = line.strip().split(':')[1].replace(" ","").split(',')

outgroup_name = sys.argv[3]

sample_names = []
sample_dict = defaultdict(list)
pairwise_dict = defaultdict(list)
rstats_dict = defaultdict(list)

number_snps = 0
err_num = 0

def grantham(aa1, aa2):
	g_distance = {
	'A': {
			'A': '0',
			'C': '195',
			'D': '126',
			'E': '107',
			'F': '113',
			'G': '60',
			'H': '86',
			'I': '94',
			'K': '106',
			'L': '96',
			'M': '84',
			'N': '111',
			'P': '27',
			'Q': '91',
			'R': '112',
			'S': '99',	
			'T': '58',
			'V': '64',
			'W': '148',
			'Y': '112',
			},
	'C': {
			'A': '195',
			'C': '0',
			'D': '154',
			'E': '170',
			'F': '205',
			'G': '159',
			'H': '174',
			'I': '198',
			'K': '202',
			'L': '198',
			'M': '196',
			'N': '139',
			'P': '169',
			'Q': '154',
			'R': '180',
			'S': '112',
			'T': '149',
			'V': '192',
			'W': '215',
			'Y': '194',
			},
	'D': {
			'A': '126',
			'C': '154',
			'D': '0',
			'E': '45',
			'F': '177',
			'G': '94',
			'H': '81',
			'I': '168',
			'K': '101',
			'L': '172',
			'M': '160',
			'N': '23',
			'P': '108',
			'Q': '61',
			'R': '96',
			'S': '65',
			'T': '85',
			'V': '152',
			'W': '181',
			'Y': '160',
			},
	'E': {
			'A': '107',
			'C': '170',
			'D': '45',
			'E': '0',
			'F': '140',
			'G': '98',
			'H': '40',
			'I': '134',
			'K': '56',
			'L': '138',
			'M': '126',
			'N': '42',
			'P': '93',
			'Q': '29',
			'R': '54',
			'S': '80',
			'T': '65',
			'V': '121',
			'W': '152',
			'Y': '122',
			},
	'F': {
			'A': '113',
			'C': '205',
			'D': '177',
			'E': '140',
			'F': '0',
			'G': '153',
			'H': '100',
			'I': '21',
			'K': '102',
			'L': '22',
			'M': '28',
			'N': '158',
			'P': '114',
			'Q': '116',
			'R': '97',
			'S': '155',
			'T': '103',
			'V': '50',
			'W': '40',
			'Y': '22',
			},
	'G': {
			'A': '60',
			'C': '159',
			'D': '94',
			'E': '98',
			'F': '153',
			'G': '0',
			'H': '98',
			'I': '135',
			'K': '127',
			'L': '138',
			'M': '127',
			'N': '80',
			'P': '42',
			'Q': '87',
			'R': '125',
			'S': '56',
			'T': '59',
			'V': '109',
			'W': '184',
			'Y': '147',
			},
	'H': {		
			'A': '86',
			'C': '174',
			'D': '81',
			'E': '40',
			'F': '100',
			'G': '98',
			'H': '0',
			'I': '94',
			'K': '32',
			'L': '99',
			'M': '87',
			'N': '68',
			'P': '77',
			'Q': '24',
			'R': '29',
			'S': '89',
			'T': '47',
			'V': '84',
			'W': '115',
			'Y': '83',
			},
	'I': {
			'A': '94',
			'C': '198',
			'D': '168',
			'E': '134',
			'F': '21',
			'G': '135',
			'H': '94',
			'I': '0',
			'K': '102',
			'L': '5',
			'M': '10',
			'N': '149',
			'P': '95',
			'Q': '109',
			'R': '97',
			'S': '142',
			'T': '89',
			'V': '29',
			'W': '61',
			'Y': '33',
			},
	'K': {
			'A': '106',
			'C': '202',
			'D': '101',
			'E': '56',
			'F': '102',
			'G': '127',
			'H': '32',
			'I': '102',
			'K': '0',
			'L': '107',
			'M': '95',
			'N': '94',
			'P': '103',
			'Q': '53',
			'R': '26',
			'S': '121',
			'T': '78',
			'V': '97',
			'W': '110',
			'Y': '85',
			},
	'L': {
			'A': '96',
			'C': '198',
			'D': '172',
			'E': '138',
			'F': '22',
			'G': '138',
			'H': '99',
			'I': '5',
			'K': '107',
			'L': '0',
			'M': '15',
			'N': '153',
			'P': '98',
			'Q': '113',
			'R': '102',
			'S': '145',
			'T': '92',
			'V': '32',
			'W': '61',
			'Y': '36',
			},
	'M':{
			'A': '84',
			'C': '196',
			'D': '160',
			'E': '126',
			'F': '28',
			'G': '127',
			'H': '87',
			'I': '10',
			'K': '95',
			'L': '15',
			'M': '0',
			'N': '142',
			'P': '87',
			'Q': '101',
			'R': '91',
			'S': '135',
			'T': '81',
			'V': '21',
			'W': '67',
			'Y': '36',
			},
	'N': {
			'A': '111',
			'C': '139',
			'D': '23',
			'E': '42',
			'F': '158',
			'G': '80',
			'H': '68',
			'I': '149',
			'K': '94',
			'L': '153',
			'M': '142',
			'N': '0',
			'P': '91',
			'Q': '46',
			'R': '86',
			'S': '46',
			'T': '65',
			'V': '133',
			'W': '174',
			'Y': '143',
			},
	'P': {
			'A': '27',
			'C': '169',
			'D': '108',
			'E': '93',
			'F': '114',
			'G': '42',
			'H': '77',
			'I': '95',
			'K': '103',
			'L': '98',
			'M': '87',
			'N': '91',
			'P': '0',
			'Q': '76',
			'R': '103',
			'S': '74',
			'T': '38',
			'V': '68',
			'W': '147',
			'Y': '110',
			},
	'Q': {
			'A': '91',
			'C': '154',
			'D': '61',
			'E': '29',
			'F': '116',
			'G': '87',
			'H': '24',
			'I': '109',
			'K': '53',
			'L': '113',
			'M': '101',
			'N': '46',
			'P': '76',
			'Q': '0',
			'R': '43',
			'S': '68',
			'T': '42',
			'V': '96',
			'W': '130',
			'Y': '99',
			},
	'R': {		
			'A': '112',
			'C': '180',
			'D': '96',
			'E': '54',
			'F': '97',
			'G': '125',
			'H': '29',
			'I': '97',
			'K': '26',
			'L': '102',
			'M': '91',
			'N': '86',
			'P': '103',
			'Q': '43',
			'R': '0',
			'S': '110',
			'T': '71',
			'V': '96',
			'W': '101',
			'Y': '77',
			},
	'S': {
			'A': '99',
			'C': '112',
			'D': '65',
			'E': '80',
			'F': '155',
			'G': '56',
			'H': '89',
			'I': '142',
			'K': '121',
			'L': '145',
			'M': '135',
			'N': '46',
			'P': '74',
			'Q': '68',
			'R': '110',
			'S': '0',
			'T': '58',
			'V': '124',
			'W': '177',
			'Y': '144',
			},
	'T': {
			'A': '58',
			'C': '149',
			'D': '85',
			'E': '65',
			'F': '103',
			'G': '59',
			'H': '47',
			'I': '89',
			'K': '78',
			'L': '92',
			'M': '81',
			'N': '65',
			'P': '38',
			'Q': '42',
			'R': '71',
			'S': '58',
			'T': '0',
			'V': '69',
			'W': '128',
			'Y': '92',
			},
	'V': {
			'A': '64',
			'C': '192',
			'D': '152',
			'E': '121',
			'F': '50',
			'G': '109',
			'H': '84',
			'I': '29',
			'K': '97',
			'L': '32',
			'M': '21',
			'N': '133',
			'P': '68',
			'Q': '96',
			'R': '96',
			'S': '124',
			'T': '69',
			'V': '0',
			'W': '88',
			'Y': '55',
			},
	'W': {
			'A': '148',
			'C': '215',
			'D': '181',
			'E': '152',
			'F': '40',
			'G': '184',
			'H': '115',
			'I': '61',
			'K': '110',
			'L': '61',
			'M': '67',
			'N': '174',
			'P': '147',
			'Q': '130',
			'R': '101',
			'S': '177',
			'T': '128',
			'V': '88',
			'W': '0',
			'Y': '37',
			},
	'Y': {
			'A': '112',
			'C': '194',
			'D': '160',
			'E': '122',
			'F': '22',
			'G': '147',
			'H': '83',
			'I': '33',
			'K': '85',
			'L': '36',
			'M': '36',
			'N': '143',
			'P': '110',
			'Q': '99',
			'R': '77',
			'S': '144',
			'T': '92',
			'V': '55',
			'W': '37',
			'Y': '0'
			}
	}
	return g_distance[aa1][aa2];

#Amino acid names

aa_dict = {
	"Ala": "A",
	"Arg": "R",
	"Asn": "N",
	"Asp": "D",
	"Asx": "B",
	"Cys": "C",
	"Glu": "E",
	"Gln": "Q",
	"Glx": "Z",
	"Gly": "G",
	"His": "H",
	"Ile": "I",
	"Leu": "L",
	"Lys": "K",
	"Met": "M",
	"Phe": "F",
	"Pro": "P",
	"Ser": "S",
	"Thr": "T",
	"Trp": "W",
	"Tyr": "Y",
	"Val": "V"
}

with open(vcf_name) as vcf_f:
	for line in vcf_f:
		if line.startswith("##"):
			continue
		elif line.startswith("#"):
			line = line.strip().split()
			len_l = len(line)
			for i in range(9, len_l):
				s_name = line[i]
				if outgroup_name in s_name:
					outgroup_number = (i-9)
				sample_names.append(s_name)
				sample_dict[s_name] = defaultdict(list)
				sample_dict[s_name]["ANNOTATION"]
				sample_dict[s_name]["ANN_IMPACT"]
				sample_dict[s_name]["GRANTHAM_SCORE"]
				sample_dict[s_name]["COUNT_HOMO_DER"]
				sample_dict[s_name]["COUNT_HET_DER"]
				sample_dict[s_name]["COUNT_ALL_DER"]
				sample_dict[s_name]["COUNT_HOMO_SYN"]
				sample_dict[s_name]["COUNT_HET_SYN"]
				sample_dict[s_name]["COUNT_ALL_SYN"]
				sample_dict[s_name]["COUNT_HOMO_NS"]
				sample_dict[s_name]["COUNT_HET_NS"]
				sample_dict[s_name]["COUNT_ALL_NS"]
				sample_dict[s_name]["COUNT_HOMO_STOP"]
				sample_dict[s_name]["COUNT_HET_STOP"]
				sample_dict[s_name]["COUNT_ALL_STOP"]
				sample_dict[s_name]["COUNT_HOMO_CONS"]
				sample_dict[s_name]["COUNT_HET_CONS"]
				sample_dict[s_name]["COUNT_ALL_CONS"]
				sample_dict[s_name]["COUNT_HOMO_MODC"]
				sample_dict[s_name]["COUNT_HET_MODC"]
				sample_dict[s_name]["COUNT_ALL_MODC"]
				sample_dict[s_name]["COUNT_HOMO_MODR"]
				sample_dict[s_name]["COUNT_HET_MODR"]
				sample_dict[s_name]["COUNT_ALL_MODR"]
				sample_dict[s_name]["COUNT_HOMO_RAD"]
				sample_dict[s_name]["COUNT_HET_RAD"]
				sample_dict[s_name]["COUNT_ALL_RAD"]
			for p_name in pop_names:
				pop_dict[p_name] = defaultdict(list)
				pop_dict[p_name]["ANNOTATION"]
				pop_dict[p_name]["ANN_IMPACT"]
				pop_dict[p_name]["GRANTHAM_SCORE"]
				pop_dict[p_name]["COUNT_HOMO_DER"]
				pop_dict[p_name]["COUNT_HET_DER"]
				pop_dict[p_name]["COUNT_ALL_DER"]
				pop_dict[p_name]["COUNT_HOMO_SYN"]
				pop_dict[p_name]["COUNT_HET_SYN"]
				pop_dict[p_name]["COUNT_ALL_SYN"]
				pop_dict[p_name]["COUNT_HOMO_NS"]
				pop_dict[p_name]["COUNT_HET_NS"]
				pop_dict[p_name]["COUNT_ALL_NS"]
				pop_dict[p_name]["COUNT_HOMO_STOP"]
				pop_dict[p_name]["COUNT_HET_STOP"]
				pop_dict[p_name]["COUNT_ALL_STOP"]
				pop_dict[p_name]["COUNT_HOMO_CONS"]
				pop_dict[p_name]["COUNT_HET_CONS"]
				pop_dict[p_name]["COUNT_ALL_CONS"]
				pop_dict[p_name]["COUNT_HOMO_MODC"]
				pop_dict[p_name]["COUNT_HET_MODC"]
				pop_dict[p_name]["COUNT_ALL_MODC"]
				pop_dict[p_name]["COUNT_HOMO_MODR"]
				pop_dict[p_name]["COUNT_HET_MODR"]
				pop_dict[p_name]["COUNT_ALL_MODR"]
				pop_dict[p_name]["COUNT_HOMO_RAD"]
				pop_dict[p_name]["COUNT_HET_RAD"]
				pop_dict[p_name]["COUNT_ALL_RAD"]
			pairs = itertools.combinations(pop_names, 2)
			for pair in pairs:
				pairwise_dict[pair] = defaultdict(list)
				pairwise_dict[pair]["L_xnoty"] = defaultdict(list)
				pairwise_dict[pair]["L_ynotx"] = defaultdict(list)
				pairwise_dict[pair]["L_xnoty"]["SYN"]
				pairwise_dict[pair]["L_xnoty"]["NS"]
				pairwise_dict[pair]["L_xnoty"]["STOP"]
				pairwise_dict[pair]["L_ynotx"]["SYN"]
				pairwise_dict[pair]["L_ynotx"]["NS"]
				pairwise_dict[pair]["L_ynotx"]["STOP"]
				pairwise_dict[pair]["L_xnoty"]["CONS"]
				pairwise_dict[pair]["L_xnoty"]["MODC"]
				pairwise_dict[pair]["L_xnoty"]["MODR"]
				pairwise_dict[pair]["L_xnoty"]["RAD"]
				pairwise_dict[pair]["L_ynotx"]["CONS"]
				pairwise_dict[pair]["L_ynotx"]["MODC"]
				pairwise_dict[pair]["L_ynotx"]["MODR"]
				pairwise_dict[pair]["L_ynotx"]["RAD"]
				pairwise_dict[pair]["L_xnoty_Hom"] = defaultdict(list)
				pairwise_dict[pair]["L_ynotx_Hom"] = defaultdict(list)
				pairwise_dict[pair]["L_xnoty_Hom"]["SYN"]
				pairwise_dict[pair]["L_xnoty_Hom"]["NS"]
				pairwise_dict[pair]["L_xnoty_Hom"]["STOP"]
				pairwise_dict[pair]["L_ynotx_Hom"]["SYN"]
				pairwise_dict[pair]["L_ynotx_Hom"]["NS"]
				pairwise_dict[pair]["L_ynotx_Hom"]["STOP"]
				pairwise_dict[pair]["L_xnoty_Hom"]["CONS"]
				pairwise_dict[pair]["L_xnoty_Hom"]["MODC"]
				pairwise_dict[pair]["L_xnoty_Hom"]["MODR"]
				pairwise_dict[pair]["L_xnoty_Hom"]["RAD"]
				pairwise_dict[pair]["L_ynotx_Hom"]["CONS"]
				pairwise_dict[pair]["L_ynotx_Hom"]["MODC"]
				pairwise_dict[pair]["L_ynotx_Hom"]["MODR"]
				pairwise_dict[pair]["L_ynotx_Hom"]["RAD"]
		else:
			try:
				record = next(vcf_reader)
				line = line.strip().split()
				if "WARNING" in line:
					pass
				elif record.genotype(outgroup_name)['GT'] == "0/1":
					pass
				else:
					number_snps += 1
					parts = line[7].strip().split('ANN=')
					info = parts[1].strip().split(',')
					ann_info = info[0].strip().split('|')
					ANNOTATION = ann_info[1]
					ANN_IMPACT = ann_info[2]
					HGVS_P = ann_info[10].strip().split('.')
					if ANNOTATION == "missense_variant":
						three_aa1 = HGVS_P[1][0:3]
						three_aa2 = HGVS_P[1][-3:]
						one_aa1 = aa_dict[three_aa1]
						one_aa2 = aa_dict[three_aa2]
						GRANTHAM_SCORE = grantham(one_aa1, one_aa2)
					# target_index = line.index('GT:PL')
					# genotypes = line[target_index+1:]
					# outgroup_gt = genotypes[outgroup_number]
					# outgroup_state = outgroup_gt.strip().split(':')
					ratio_i = defaultdict(list)
					hom_ratio_i = defaultdict(list)
					for p_name in pop_names:
						s_n = 0
						d_i = 0
						n_i = 0
						het_i = 0
						hom_i = 0
						der_i = 0
						hom_syn_i = 0
						het_syn_i = 0
						der_syn_i = 0
						hom_mis_i = 0
						het_mis_i = 0
						der_mis_i = 0
						hom_con_i = 0
						het_con_i = 0
						der_con_i = 0
						hom_mdc_i = 0
						het_mdc_i = 0
						der_mdc_i = 0
						hom_mdr_i = 0
						het_mdr_i = 0
						der_mdr_i = 0
						hom_rad_i = 0
						het_rad_i = 0
						der_rad_i = 0
						hom_stop_i = 0
						het_stop_i = 0
						der_stop_i = 0

						pop_list = pop_l_dict[p_name]

						if record.genotype(outgroup_name)['GT'] == "0/0":
							for sample in pop_list:
								s_n += 1
								n_i += 2
	
								if record.genotype(sample)['GT'] == "0/1":
									d_i += 1
									het_i += 1
									der_i += 1
									sample_dict[sample]["ANNOTATION"].append(ANNOTATION)
									sample_dict[sample]["COUNT_HOMO_DER"].append(0)
									sample_dict[sample]["COUNT_HET_DER"].append(1)
									sample_dict[sample]["COUNT_ALL_DER"].append(1)
									if ANNOTATION == "synonymous_variant":
										het_syn_i += 1
										der_syn_i += 1
										sample_dict[sample]["COUNT_HOMO_SYN"].append(0)
										sample_dict[sample]["COUNT_HET_SYN"].append(1)
										sample_dict[sample]["COUNT_ALL_SYN"].append(1)
									elif ANNOTATION == "missense_variant":
										het_mis_i += 1
										der_mis_i += 1
										sample_dict[sample]["COUNT_HOMO_NS"].append(0)
										sample_dict[sample]["COUNT_HET_NS"].append(1)
										sample_dict[sample]["COUNT_ALL_NS"].append(1)
										if float(GRANTHAM_SCORE) < 51:
											het_con_i += 1
											der_con_i += 1
											sample_dict[sample]["COUNT_HOMO_CONS"].append(0)
											sample_dict[sample]["COUNT_HET_CONS"].append(1)
											sample_dict[sample]["COUNT_ALL_CONS"].append(1)
										elif 50 < float(GRANTHAM_SCORE) < 101:
											het_mdc_i += 1
											der_mdc_i += 1
											sample_dict[sample]["COUNT_HOMO_MODC"].append(0)
											sample_dict[sample]["COUNT_HET_MODC"].append(1)
											sample_dict[sample]["COUNT_ALL_MODC"].append(1)
										elif 100 < float(GRANTHAM_SCORE) < 151:
											het_mdr_i += 1
											der_mdr_i += 1
											sample_dict[sample]["COUNT_HOMO_MODR"].append(0)
											sample_dict[sample]["COUNT_HET_MODR"].append(1)
											sample_dict[sample]["COUNT_ALL_MODR"].append(1)
										elif float(GRANTHAM_SCORE) > 150:
											het_rad_i += 1
											der_rad_i += 1
											sample_dict[sample]["COUNT_HOMO_RAD"].append(0)
											sample_dict[sample]["COUNT_HET_RAD"].append(1)
											sample_dict[sample]["COUNT_ALL_RAD"].append(1)
									elif ANNOTATION == "stop_gained":
										het_stop_i += 1
										der_stop_i += 1
										sample_dict[sample]["COUNT_HOMO_STOP"].append(0)
										sample_dict[sample]["COUNT_HET_STOP"].append(1)
										sample_dict[sample]["COUNT_ALL_STOP"].append(1)
								
								elif record.genotype(sample)['GT'] == "1/1":
									d_i += 2
									hom_i += 1
									der_i += 2
									sample_dict[sample]["ANNOTATION"].append(ANNOTATION)
									sample_dict[sample]["COUNT_HOMO_DER"].append(1)
									sample_dict[sample]["COUNT_HET_DER"].append(0)
									sample_dict[sample]["COUNT_ALL_DER"].append(2)
									if ANNOTATION == "synonymous_variant":
										hom_syn_i += 1
										der_syn_i += 2
										sample_dict[sample]["COUNT_HOMO_SYN"].append(1)
										sample_dict[sample]["COUNT_HET_SYN"].append(0)
										sample_dict[sample]["COUNT_ALL_SYN"].append(2)
									elif ANNOTATION == "missense_variant":
										hom_mis_i += 1
										der_mis_i += 2
										sample_dict[sample]["COUNT_HOMO_NS"].append(1)
										sample_dict[sample]["COUNT_HET_NS"].append(0)
										sample_dict[sample]["COUNT_ALL_NS"].append(2)
										if float(GRANTHAM_SCORE) < 51:
											hom_con_i += 1
											der_con_i += 2
											sample_dict[sample]["COUNT_HOMO_CONS"].append(1)
											sample_dict[sample]["COUNT_HET_CONS"].append(0)
											sample_dict[sample]["COUNT_ALL_CONS"].append(2)
										elif 50 < float(GRANTHAM_SCORE) < 101:
											hom_mdc_i += 1
											der_mdc_i += 2
											sample_dict[sample]["COUNT_HOMO_MODC"].append(1)
											sample_dict[sample]["COUNT_HET_MODC"].append(0)
											sample_dict[sample]["COUNT_ALL_MODC"].append(2)
										elif 100 < float(GRANTHAM_SCORE) < 151:
											hom_mdr_i += 1
											der_mdr_i += 2
											sample_dict[sample]["COUNT_HOMO_MODR"].append(1)
											sample_dict[sample]["COUNT_HET_MODR"].append(0)
											sample_dict[sample]["COUNT_ALL_MODR"].append(2)
										elif float(GRANTHAM_SCORE) > 150:
											hom_rad_i += 1
											hom_rad_i += 2
											sample_dict[sample]["COUNT_HOMO_RAD"].append(1)
											sample_dict[sample]["COUNT_HET_RAD"].append(0)
											sample_dict[sample]["COUNT_ALL_RAD"].append(2)
									elif ANNOTATION == "stop_gained":
										hom_stop_i += 1
										der_stop_i += 2
										sample_dict[sample]["COUNT_HOMO_STOP"].append(1)
										sample_dict[sample]["COUNT_HET_STOP"].append(0)
										sample_dict[sample]["COUNT_ALL_STOP"].append(2)


						elif record.genotype(outgroup_name)['GT'] == "1/1":
							for sample in pop_list:
								s_n += 1
								n_i += 2
	
								if record.genotype(sample)['GT'] == "0/1":
									d_i += 1
									het_i += 1
									der_i += 1
									sample_dict[sample]["ANNOTATION"].append(ANNOTATION)
									sample_dict[sample]["COUNT_HOMO_DER"].append(0)
									sample_dict[sample]["COUNT_HET_DER"].append(1)
									sample_dict[sample]["COUNT_ALL_DER"].append(1)
									if ANNOTATION == "synonymous_variant":
										het_syn_i += 1
										der_syn_i += 1
										sample_dict[sample]["COUNT_HOMO_SYN"].append(0)
										sample_dict[sample]["COUNT_HET_SYN"].append(1)
										sample_dict[sample]["COUNT_ALL_SYN"].append(1)
									elif ANNOTATION == "missense_variant":
										het_mis_i += 1
										der_mis_i += 1
										sample_dict[sample]["COUNT_HOMO_NS"].append(0)
										sample_dict[sample]["COUNT_HET_NS"].append(1)
										sample_dict[sample]["COUNT_ALL_NS"].append(1)
										if float(GRANTHAM_SCORE) < 51:
											het_con_i += 1
											der_con_i += 1
											sample_dict[sample]["COUNT_HOMO_CONS"].append(0)
											sample_dict[sample]["COUNT_HET_CONS"].append(1)
											sample_dict[sample]["COUNT_ALL_CONS"].append(1)
										elif 50 < float(GRANTHAM_SCORE) < 101:
											het_mdc_i += 1
											der_mdc_i += 1
											sample_dict[sample]["COUNT_HOMO_MODC"].append(0)
											sample_dict[sample]["COUNT_HET_MODC"].append(1)
											sample_dict[sample]["COUNT_ALL_MODC"].append(1)
										elif 100 < float(GRANTHAM_SCORE) < 151:
											het_mdr_i += 1
											der_mdr_i += 1
											sample_dict[sample]["COUNT_HOMO_MODR"].append(0)
											sample_dict[sample]["COUNT_HET_MODR"].append(1)
											sample_dict[sample]["COUNT_ALL_MODR"].append(1)
										elif float(GRANTHAM_SCORE) > 150:
											het_rad_i += 1
											der_rad_i += 1
											sample_dict[sample]["COUNT_HOMO_RAD"].append(0)
											sample_dict[sample]["COUNT_HET_RAD"].append(1)
											sample_dict[sample]["COUNT_ALL_RAD"].append(1)
									elif ANNOTATION == "stop_gained":
										het_stop_i += 1
										der_stop_i += 1
										sample_dict[sample]["COUNT_HOMO_STOP"].append(0)
										sample_dict[sample]["COUNT_HET_STOP"].append(1)
										sample_dict[sample]["COUNT_ALL_STOP"].append(1)
								
								elif record.genotype(sample)['GT'] == "0/0":
									d_i += 2
									hom_i += 1
									der_i += 2
									sample_dict[sample]["ANNOTATION"].append(ANNOTATION)
									sample_dict[sample]["COUNT_HOMO_DER"].append(1)
									sample_dict[sample]["COUNT_HET_DER"].append(0)
									sample_dict[sample]["COUNT_ALL_DER"].append(2)
									if ANNOTATION == "synonymous_variant":
										hom_syn_i += 1
										der_syn_i += 2
										sample_dict[sample]["COUNT_HOMO_SYN"].append(1)
										sample_dict[sample]["COUNT_HET_SYN"].append(0)
										sample_dict[sample]["COUNT_ALL_SYN"].append(2)
									elif ANNOTATION == "missense_variant":
										hom_mis_i += 1
										der_mis_i += 2
										sample_dict[sample]["COUNT_HOMO_NS"].append(1)
										sample_dict[sample]["COUNT_HET_NS"].append(0)
										sample_dict[sample]["COUNT_ALL_NS"].append(2)
										if float(GRANTHAM_SCORE) < 51:
											hom_con_i += 1
											der_con_i += 2
											sample_dict[sample]["COUNT_HOMO_CONS"].append(1)
											sample_dict[sample]["COUNT_HET_CONS"].append(0)
											sample_dict[sample]["COUNT_ALL_CONS"].append(2)
										elif 50 < float(GRANTHAM_SCORE) < 101:
											hom_mdc_i += 1
											der_mdc_i += 2
											sample_dict[sample]["COUNT_HOMO_MODC"].append(1)
											sample_dict[sample]["COUNT_HET_MODC"].append(0)
											sample_dict[sample]["COUNT_ALL_MODC"].append(2)
										elif 100 < float(GRANTHAM_SCORE) < 151:
											hom_mdr_i += 1
											der_mdr_i += 2
											sample_dict[sample]["COUNT_HOMO_MODR"].append(1)
											sample_dict[sample]["COUNT_HET_MODR"].append(0)
											sample_dict[sample]["COUNT_ALL_MODR"].append(2)
										elif float(GRANTHAM_SCORE) > 150:
											hom_rad_i += 1
											hom_rad_i += 2
											sample_dict[sample]["COUNT_HOMO_RAD"].append(1)
											sample_dict[sample]["COUNT_HET_RAD"].append(0)
											sample_dict[sample]["COUNT_ALL_RAD"].append(2)
									elif ANNOTATION == "stop_gained":
										hom_stop_i += 1
										der_stop_i += 2
										sample_dict[sample]["COUNT_HOMO_STOP"].append(1)
										sample_dict[sample]["COUNT_HET_STOP"].append(0)
										sample_dict[sample]["COUNT_ALL_STOP"].append(2)

						ratio_i[p_name] = (d_i/n_i)
						hom_ratio_i[p_name] = (2*d_i*(n_i-d_i))/(n_i*(n_i-1))
						pop_dict[p_name]["COUNT_HOMO_DER"].append(hom_i/s_n)
						pop_dict[p_name]["COUNT_HET_DER"].append(het_i/s_n)
						pop_dict[p_name]["COUNT_ALL_DER"].append(der_i/s_n)

					## calc l stats
					for pair in pairwise_dict:
						x = pair[0]
						y = pair[1]
	

						lxnoty = ratio_i[x] * (1 - ratio_i[y])
						lynotx = ratio_i[y] * (1 - ratio_i[x])
						lxnoty_hom = hom_ratio_i[x] * (1 - hom_ratio_i[y])
						lynotx_hom = hom_ratio_i[y] * (1 - hom_ratio_i[x])
	
						if ANNOTATION == "synonymous_variant":
							pairwise_dict[pair]["L_xnoty"]["SYN"].append(lxnoty)
							pairwise_dict[pair]["L_ynotx"]["SYN"].append(lynotx)
							pairwise_dict[pair]["L_xnoty_Hom"]["SYN"].append(lxnoty_hom)
							pairwise_dict[pair]["L_ynotx_Hom"]["SYN"].append(lynotx_hom)
						elif ANNOTATION == "missense_variant":
							pairwise_dict[pair]["L_xnoty"]["NS"].append(lxnoty)
							pairwise_dict[pair]["L_ynotx"]["NS"].append(lynotx)
							pairwise_dict[pair]["L_xnoty_Hom"]["NS"].append(lxnoty_hom)
							pairwise_dict[pair]["L_ynotx_Hom"]["NS"].append(lynotx_hom)
							if float(GRANTHAM_SCORE) < 51:
								pairwise_dict[pair]["L_xnoty"]["CONS"].append(lxnoty)
								pairwise_dict[pair]["L_ynotx"]["CONS"].append(lynotx)
								pairwise_dict[pair]["L_xnoty_Hom"]["CONS"].append(lxnoty_hom)
								pairwise_dict[pair]["L_ynotx_Hom"]["CONS"].append(lynotx_hom)
							elif 50 < float(GRANTHAM_SCORE) < 101:
								pairwise_dict[pair]["L_xnoty"]["MODC"].append(lxnoty)
								pairwise_dict[pair]["L_ynotx"]["MODC"].append(lynotx)
								pairwise_dict[pair]["L_xnoty_Hom"]["MODC"].append(lxnoty_hom)
								pairwise_dict[pair]["L_ynotx_Hom"]["MODC"].append(lynotx_hom)
							elif 100 < float(GRANTHAM_SCORE) < 151:
								pairwise_dict[pair]["L_xnoty"]["MODR"].append(lxnoty)
								pairwise_dict[pair]["L_ynotx"]["MODR"].append(lynotx)
								pairwise_dict[pair]["L_xnoty_Hom"]["MODR"].append(lxnoty_hom)
								pairwise_dict[pair]["L_ynotx_Hom"]["MODR"].append(lynotx_hom)
							elif float(GRANTHAM_SCORE) > 150:
								pairwise_dict[pair]["L_xnoty"]["RAD"].append(lxnoty)
								pairwise_dict[pair]["L_ynotx"]["RAD"].append(lynotx)
								pairwise_dict[pair]["L_xnoty_Hom"]["RAD"].append(lxnoty_hom)
								pairwise_dict[pair]["L_ynotx_Hom"]["RAD"].append(lynotx_hom)
						elif ANNOTATION == "stop_gained":
							pairwise_dict[pair]["L_xnoty"]["STOP"].append(lxnoty)
							pairwise_dict[pair]["L_ynotx"]["STOP"].append(lynotx)
							pairwise_dict[pair]["L_xnoty_Hom"]["STOP"].append(lxnoty_hom)
							pairwise_dict[pair]["L_ynotx_Hom"]["STOP"].append(lynotx_hom)
			except IndexError:
				err_num += 1
				print("error!")
				pass
			
## Perform bootstrapping

print("Number of sites is: " + str(number_snps))
print("Number of errors was: " + str(err_num))


## make R stat output files

# out_file_Rxy_name = str("Rxy_results_{}.txt").format(bs_number)
# out_file_Rxy = open(out_file_Rxy_name, "w+")
# out_file_Rxy.write("Species_pair\tRxy_SYN\tRxy_NS\tRxy_STOP\tRxy_CONS\tRxy_MODC\tRxy_MODR\tRxy_RAD\n")

# out_file_RHom_name = str("Rxy_Hom_results_{}.txt").format(bs_number)
# out_file_RHom = open(out_file_RHom_name, "w+")
# out_file_RHom.write("Species_pair\tRxy_Hom_SYN\tRxy_Hom_NS\tRxy_Hom_STOP\tRxy_Hom_CONS\tRxy_Hom_MODC\tRxy_Hom_MODR\tRxy_Hom_RAD\n")

# out_file_Rprime_name = str("Rprime_results_{}.txt").format(bs_number)
# out_file_Rprime = open(out_file_Rprime_name, "w+")
# out_file_Rprime.write("Species_pair\tRprime_N_S\tRprime_STOP_S\tRprime_RAD_S\tRprime_RAD_CON\n")

# out_file_RprimeHom_name = str("Rprime_Hom_results_{}.txt").format(bs_number)
# out_file_RprimeHom = open(out_file_RprimeHom_name, "w+")
# out_file_RprimeHom.write("Species_pair\tRprime_Hom_N_S\tRprime_Hom_STOP_S\tRprime_Hom_RAD_S\tRprime_Hom_RAD_CON\n")


## make L stat output files

out_file_L_name = str("Lstats_results_{}.txt").format(bs_number)
out_file_L = open(out_file_L_name, "w+")
out_file_L.write("Species_pair\tL_xnoty_SYN\tL_xnoty_NS\tL_xnoty_STOP\tL_xnoty_CONS\tL_xnoty_MODC\tL_xnoty_MODR\tL_xnoty_RAD\tL_ynotx_SYN\tL_ynotx_NS\tL_ynotx_STOP\tL_ynotx_CONS\tL_ynotx_MODC\tL_ynotx_MODR\tL_ynotx_RAD\n")

out_file_LHom_name = str("Lstats_Hom_results_{}.txt").format(bs_number)
out_file_LHom = open(out_file_LHom_name, "w+")
out_file_LHom.write("Species_pair\tL_xnoty_Hom_SYN\tL_xnoty_Hom_NS\tL_xnoty_Hom_STOP\tL_xnoty_Hom_CONS\tL_xnoty_Hom_MODC\tL_xnoty_Hom_MODR\tL_xnoty_Hom_RAD\tL_ynotx_Hom_SYN\tL_ynotx_Hom_NS\tL_ynotx_Hom_STOP\tL_ynotx_Hom_CONS\tL_ynotx_Hom_MODC\tL_ynotx_Hom_MODR\tL_ynotx_Hom_RAD\n")


## make count output files

out_file_countNS_name = str("Count_sum_N_S_results_{}.txt").format(bs_number)
out_file_countNS = open(out_file_countNS_name, "w+")
out_file_countNS.write("SPECIES\tGRANTHAM_SUM\tCOUNT_HOMO_DER\tCOUNT_HET_DER\tCOUNT_ALL_DER\tCOUNT_HOMO_SYN\tCOUNT_HET_SYN\tCOUNT_ALL_SYN\tCOUNT_HOMO_NS\tCOUNT_HET_NS\tCOUNT_ALL_NS\n")

out_file_counteff_name = str("Count_sum_SNPeff_results_{}.txt").format(bs_number)
out_file_counteff = open(out_file_counteff_name, "w+")
out_file_counteff.write("SPECIES\tGRANTHAM_SUM\tCOUNT_HOMO_STOP\tCOUNT_HET_STOP\tCOUNT_ALL_STOP\n")

out_file_countg_name = str("Count_sum_Grantham_results_{}.txt").format(bs_number)
out_file_countg = open(out_file_countg_name, "w+")
out_file_countg.write("SPECIES\tGRANTHAM_SUM\tCOUNT_HOMO_CONS\tCOUNT_HET_CONS\tCOUNT_ALL_CONS\tCOUNT_HOMO_MODC\tCOUNT_HET_MODC\tCOUNT_ALL_MODC\tCOUNT_HOMO_MODR\tCOUNT_HET_MODR\tCOUNT_ALL_MODR\tCOUNT_HOMO_RAD\tCOUNT_HET_RAD\tCOUNT_ALL_RAD\n")

for pair in pairwise_dict:
	if outgroup_name not in pair[0] and outgroup_name not in pair[1]:
		sp_pair = str(pair[0] + ", " + pair[1])
		# rstats_dict[pair] = defaultdict(list)
		# if sum(pairwise_dict[pair]["L_ynotx"]["SYN"]) > 0:
		# 	rstats_dict[pair]["Rxy_SYN"].append(sum(pairwise_dict[pair]["L_xnoty"]["SYN"]) / sum(pairwise_dict[pair]["L_ynotx"]["SYN"]))
		# else:
		# 	rstats_dict[pair]["Rxy_SYN"].append(1)
		# if sum(pairwise_dict[pair]["L_ynotx"]["NS"]) > 0:
		# 	rstats_dict[pair]["Rxy_NS"].append(sum(pairwise_dict[pair]["L_xnoty"]["NS"]) / sum(pairwise_dict[pair]["L_ynotx"]["NS"]))
		# else:
		# 	rstats_dict[pair]["Rxy_NS"].append(1)
		# if sum(pairwise_dict[pair]["L_ynotx"]["STOP"]) > 0:
		# 	rstats_dict[pair]["Rxy_STOP"].append(sum(pairwise_dict[pair]["L_xnoty"]["STOP"]) / sum(pairwise_dict[pair]["L_ynotx"]["STOP"]))
		# else: 
		# 	rstats_dict[pair]["Rxy_STOP"].append(1)
		# if sum(pairwise_dict[pair]["L_ynotx"]["CONS"]) > 0:
		# 	rstats_dict[pair]["Rxy_CONS"].append(sum(pairwise_dict[pair]["L_xnoty"]["CONS"]) / sum(pairwise_dict[pair]["L_ynotx"]["CONS"]))
		# else:
		# 	rstats_dict[pair]["Rxy_CONS"].append(1)
		# if sum(pairwise_dict[pair]["L_ynotx"]["MODC"]) > 0:
		# 	rstats_dict[pair]["Rxy_MODC"].append(sum(pairwise_dict[pair]["L_xnoty"]["MODC"]) / sum(pairwise_dict[pair]["L_ynotx"]["MODC"]))
		# else:
		# 	rstats_dict[pair]["Rxy_MODC"].append(1)
		# if sum(pairwise_dict[pair]["L_ynotx"]["MODR"]) > 0:
		# 	rstats_dict[pair]["Rxy_MODR"].append(sum(pairwise_dict[pair]["L_xnoty"]["MODR"]) / sum(pairwise_dict[pair]["L_ynotx"]["MODR"]))
		# else:
		# 	rstats_dict[pair]["Rxy_MODR"].append(1)
		# if sum(pairwise_dict[pair]["L_ynotx"]["RAD"]) > 0:
		# 	rstats_dict[pair]["Rxy_RAD"].append(sum(pairwise_dict[pair]["L_xnoty"]["RAD"]) / sum(pairwise_dict[pair]["L_ynotx"]["RAD"]))
		# else: 
		# 	rstats_dict[pair]["Rxy_RAD"].append(1)
		# if sum(pairwise_dict[pair]["L_ynotx_Hom"]["SYN"]) > 0:
		# 	rstats_dict[pair]["Rxy_Hom_SYN"].append(sum(pairwise_dict[pair]["L_xnoty_Hom"]["SYN"]) / sum(pairwise_dict[pair]["L_ynotx_Hom"]["SYN"]))
		# else:
		# 	rstats_dict[pair]["Rxy_Hom_SYN"].append(1)
		# if sum(pairwise_dict[pair]["L_ynotx_Hom"]["NS"]) > 0:
		# 	rstats_dict[pair]["Rxy_Hom_NS"].append(sum(pairwise_dict[pair]["L_xnoty_Hom"]["NS"]) / sum(pairwise_dict[pair]["L_ynotx_Hom"]["NS"]))
		# else:
		# 	rstats_dict[pair]["Rxy_Hom_NS"].append(1)
		# if sum(pairwise_dict[pair]["L_ynotx_Hom"]["STOP"]) > 0:
		# 	rstats_dict[pair]["Rxy_Hom_STOP"].append(sum(pairwise_dict[pair]["L_xnoty_Hom"]["STOP"]) / sum(pairwise_dict[pair]["L_ynotx_Hom"]["STOP"]))
		# else: 
		# 	rstats_dict[pair]["Rxy_Hom_STOP"].append(1)
		# if sum(pairwise_dict[pair]["L_ynotx_Hom"]["CONS"]) > 0:
		# 	rstats_dict[pair]["Rxy_Hom_CONS"].append(sum(pairwise_dict[pair]["L_xnoty_Hom"]["CONS"]) / sum(pairwise_dict[pair]["L_ynotx_Hom"]["CONS"]))
		# else:
		# 	rstats_dict[pair]["Rxy_Hom_CONS"].append(1)
		# if sum(pairwise_dict[pair]["L_ynotx_Hom"]["MODC"]) > 0:
		# 	rstats_dict[pair]["Rxy_Hom_MODC"].append(sum(pairwise_dict[pair]["L_xnoty_Hom"]["MODC"]) / sum(pairwise_dict[pair]["L_ynotx_Hom"]["MODC"]))
		# else:
		# 	rstats_dict[pair]["Rxy_Hom_MODC"].append(1)
		# if sum(pairwise_dict[pair]["L_ynotx_Hom"]["MODR"]) > 0:
		# 	rstats_dict[pair]["Rxy_Hom_MODR"].append(sum(pairwise_dict[pair]["L_xnoty_Hom"]["MODR"]) / sum(pairwise_dict[pair]["L_ynotx_Hom"]["MODR"]))
		# else:
		# 	rstats_dict[pair]["Rxy_Hom_MODR"].append(1)
		# if sum(pairwise_dict[pair]["L_ynotx_Hom"]["RAD"]) > 0:
		# 	rstats_dict[pair]["Rxy_Hom_RAD"].append(sum(pairwise_dict[pair]["L_xnoty_Hom"]["RAD"]) / sum(pairwise_dict[pair]["L_ynotx_Hom"]["RAD"]))
		# else:
		# 	rstats_dict[pair]["Rxy_Hom_RAD"].append(1) 
		# if rstats_dict[pair]["Rxy_SYN"][0] > 0:
		# 	rstats_dict[pair]["Rprime_N_S"].append(rstats_dict[pair]["Rxy_NS"][0] / rstats_dict[pair]["Rxy_SYN"][0])
		# else:
		# 	rstats_dict[pair]["Rprime_N_S"].append(1)
		# if rstats_dict[pair]["Rxy_SYN"][0] > 0:
		# 	rstats_dict[pair]["Rprime_STOP_S"].append(rstats_dict[pair]["Rxy_STOP"][0] / rstats_dict[pair]["Rxy_SYN"][0])
		# else:
		# 	rstats_dict[pair]["Rprime_STOP_S"].append(1)
		# if rstats_dict[pair]["Rxy_SYN"][0] > 0:
		# 	rstats_dict[pair]["Rprime_RAD_S"].append(rstats_dict[pair]["Rxy_RAD"][0] / rstats_dict[pair]["Rxy_SYN"][0])
		# else:
		# 	rstats_dict[pair]["Rprime_RAD_S"].append(1)
		# if rstats_dict[pair]["Rxy_CONS"][0] > 0:	
		# 	rstats_dict[pair]["Rprime_RAD_CON"].append(rstats_dict[pair]["Rxy_RAD"][0] / rstats_dict[pair]["Rxy_CONS"][0])
		# else:
		# 	rstats_dict[pair]["Rprime_RAD_CON"].append(1)
		# if rstats_dict[pair]["Rxy_Hom_SYN"][0] > 0: 	
		# 	rstats_dict[pair]["Rprime_Hom_N_S"].append(rstats_dict[pair]["Rxy_Hom_NS"][0] / rstats_dict[pair]["Rxy_Hom_SYN"][0])
		# else:
		# 	rstats_dict[pair]["Rprime_Hom_N_S"].append(1)
		# if rstats_dict[pair]["Rxy_Hom_SYN"][0] > 0:
		# 	rstats_dict[pair]["Rprime_Hom_HIGH_S"].append(rstats_dict[pair]["Rxy_Hom_HIGH"][0] / rstats_dict[pair]["Rxy_Hom_SYN"][0])
		# else:
		# 	rstats_dict[pair]["Rprime_Hom_HIGH_S"].append(1)
		# if rstats_dict[pair]["Rxy_Hom_LOW"][0] > 0: 
		# 	rstats_dict[pair]["Rprime_Hom_HIGH_LOW"].append(rstats_dict[pair]["Rxy_Hom_HIGH"][0] / rstats_dict[pair]["Rxy_Hom_LOW"][0])
		# else:
		# 	rstats_dict[pair]["Rprime_Hom_HIGH_LOW"].append(1)
		# if rstats_dict[pair]["Rxy_Hom_SYN"][0] > 0: 
		# 	rstats_dict[pair]["Rprime_Hom_STOP_S"].append(rstats_dict[pair]["Rxy_Hom_STOP"][0] / rstats_dict[pair]["Rxy_Hom_SYN"][0])
		# else: 
		# 	rstats_dict[pair]["Rprime_Hom_STOP_S"].append(1)
		# if rstats_dict[pair]["Rxy_Hom_SYN"][0] > 0:
		# 	rstats_dict[pair]["Rprime_Hom_RAD_S"].append(rstats_dict[pair]["Rxy_Hom_RAD"][0] / rstats_dict[pair]["Rxy_Hom_SYN"][0])
		# else: 
		# 	rstats_dict[pair]["Rprime_Hom_RAD_S"].append(1)
		# if rstats_dict[pair]["Rxy_Hom_CONS"][0] > 0: 
		# 	rstats_dict[pair]["Rprime_Hom_RAD_CON"].append(rstats_dict[pair]["Rxy_Hom_RAD"][0] / rstats_dict[pair]["Rxy_Hom_CONS"][0])
		# else:
		# 	rstats_dict[pair]["Rprime_Hom_RAD_CON"].append(1)

		# out_file_Rxy = open(out_file_Rxy_name, "a")
		# out_file_Rxy.write(str(sp_pair) + "\t" + str(rstats_dict[pair]["Rxy_SYN"][0]) + "\t" + str(rstats_dict[pair]["Rxy_NS"][0]) + "\t" + str(rstats_dict[pair]["Rxy_LOW"][0]) + "\t" +
		#  str(rstats_dict[pair]["Rxy_MOD"][0]) + "\t" + str(rstats_dict[pair]["Rxy_HIGH"][0]) + "\t" + str(rstats_dict[pair]["Rxy_STOP"][0]) + "\t" + 
		#  str(rstats_dict[pair]["Rxy_CONS"][0]) + "\t" + str(rstats_dict[pair]["Rxy_MODC"][0]) + "\t" + str(rstats_dict[pair]["Rxy_MODR"][0]) + "\t" + str(rstats_dict[pair]["Rxy_RAD"][0]) + "\n")
		# out_file_Rxy.close()

		# out_file_RHom = open(out_file_RHom_name, "a")
		# out_file_RHom.write(str(sp_pair) + "\t" + 
		#  str(rstats_dict[pair]["Rxy_Hom_SYN"][0]) + "\t" + str(rstats_dict[pair]["Rxy_Hom_NS"][0]) + "\t" + str(rstats_dict[pair]["Rxy_Hom_LOW"][0]) + "\t" + 
		#  str(rstats_dict[pair]["Rxy_Hom_MOD"][0]) + "\t" + str(rstats_dict[pair]["Rxy_Hom_HIGH"][0]) + "\t" + str(rstats_dict[pair]["Rxy_Hom_STOP"][0]) + "\t" + 
		#  str(rstats_dict[pair]["Rxy_Hom_CONS"][0]) + "\t" + str(rstats_dict[pair]["Rxy_Hom_MODC"][0]) + "\t" + str(rstats_dict[pair]["Rxy_Hom_MODR"][0]) + "\t" + 
		#  str(rstats_dict[pair]["Rxy_Hom_RAD"][0]) + "\n")
		# out_file_RHom.close()

		# out_file_Rprime = open(out_file_Rprime_name, "a")
		# out_file_Rprime.write(str(sp_pair) + "\t" + str(rstats_dict[pair]["Rprime_N_S"][0]) + "\t" + str(rstats_dict[pair]["Rprime_HIGH_S"][0]) + "\t" +
		#  str(rstats_dict[pair]["Rprime_HIGH_LOW"][0]) + "\t" + str(rstats_dict[pair]["Rprime_STOP_S"][0]) + "\t" + str(rstats_dict[pair]["Rprime_RAD_S"][0]) + "\t" +
		#  str(rstats_dict[pair]["Rprime_RAD_CON"][0]) + "\n")
		# out_file_Rprime.close()

		# out_file_RprimeHom = open(out_file_RprimeHom_name, "a")
		# out_file_RprimeHom.write(str(sp_pair) + "\t" + str(rstats_dict[pair]["Rprime_Hom_N_S"][0]) + "\t" + str(rstats_dict[pair]["Rprime_Hom_HIGH_S"][0]) + "\t" +
		#  str(rstats_dict[pair]["Rprime_Hom_HIGH_LOW"][0]) + "\t" + str(rstats_dict[pair]["Rprime_Hom_STOP_S"][0]) + "\t" + str(rstats_dict[pair]["Rprime_Hom_RAD_S"][0]) + "\t" +
		#  str(rstats_dict[pair]["Rprime_Hom_RAD_CON"][0]) + "\n")
		# out_file_RprimeHom.close()

		out_file_L = open(out_file_L_name, "a")
		out_file_L.write(str(sp_pair) + "\t" + str(sum(pairwise_dict[pair]["L_xnoty"]["SYN"])) +  "\t" + str(sum(pairwise_dict[pair]["L_xnoty"]["NS"])) + "\t" + 
		 str(sum(pairwise_dict[pair]["L_xnoty"]["STOP"])) + "\t" + str(sum(pairwise_dict[pair]["L_xnoty"]["CONS"])) + "\t" + str(sum(pairwise_dict[pair]["L_xnoty"]["MODC"])) + "\t" + 
		 str(sum(pairwise_dict[pair]["L_xnoty"]["MODR"])) + "\t" + str(sum(pairwise_dict[pair]["L_xnoty"]["RAD"])) + "\t" + str(sum(pairwise_dict[pair]["L_ynotx"]["SYN"])) +  "\t" + 
		 str(sum(pairwise_dict[pair]["L_ynotx"]["NS"])) + "\t" + str(sum(pairwise_dict[pair]["L_ynotx"]["STOP"])) + "\t" + str(sum(pairwise_dict[pair]["L_ynotx"]["CONS"])) + "\t" + 
		 str(sum(pairwise_dict[pair]["L_ynotx"]["MODC"])) + "\t" + str(sum(pairwise_dict[pair]["L_ynotx"]["MODR"])) + "\t" + str(sum(pairwise_dict[pair]["L_ynotx"]["RAD"])) + "\n")
		out_file_L.close()

		out_file_LHom = open(out_file_LHom_name, "a")
		out_file_LHom.write(str(sp_pair) + "\t" + str(sum(pairwise_dict[pair]["L_xnoty_Hom"]["SYN"])) + "\t" + str(sum(pairwise_dict[pair]["L_xnoty_Hom"]["NS"])) + "\t" + 
		 str(sum(pairwise_dict[pair]["L_xnoty_Hom"]["STOP"])) + "\t" + str(sum(pairwise_dict[pair]["L_xnoty_Hom"]["CONS"])) + "\t" + 
		 str(sum(pairwise_dict[pair]["L_xnoty_Hom"]["MODC"])) + "\t" + str(sum(pairwise_dict[pair]["L_xnoty_Hom"]["MODR"])) + "\t" + str(sum(pairwise_dict[pair]["L_xnoty_Hom"]["RAD"])) + "\t" + 
		 str(sum(pairwise_dict[pair]["L_ynotx_Hom"]["SYN"])) + "\t" + str(sum(pairwise_dict[pair]["L_ynotx_Hom"]["NS"])) + "\t" + str(sum(pairwise_dict[pair]["L_ynotx_Hom"]["STOP"])) + "\t" + 
		 str(sum(pairwise_dict[pair]["L_ynotx_Hom"]["CONS"])) + "\t" + str(sum(pairwise_dict[pair]["L_ynotx_Hom"]["MODC"])) + "\t" + str(sum(pairwise_dict[pair]["L_ynotx_Hom"]["MODR"])) + "\t" + 
		 str(sum(pairwise_dict[pair]["L_ynotx_Hom"]["RAD"])) + "\n")
		out_file_LHom.close()

	else:
		pass

def calsum(l):
	return  str(sum([float(i) for i in l if type(i) != str]))

for name in sample_names:
	out_file_countNS = open(out_file_countNS_name, "a")
	out_file_countNS.write(name + "\t" + calsum(sample_dict[name]["GRANTHAM"])+ "\t" + calsum(sample_dict[name]["COUNT_HOMO_DER"]) + "\t" + calsum(sample_dict[name]["COUNT_HET_DER"]) + "\t" + 
	calsum(sample_dict[name]["COUNT_ALL_DER"]) + "\t" + calsum(sample_dict[name]["COUNT_HOMO_SYN"]) + "\t" + calsum(sample_dict[name]["COUNT_HET_SYN"]) + "\t" + 
	calsum(sample_dict[name]["COUNT_ALL_SYN"]) + "\t" + calsum(sample_dict[name]["COUNT_HOMO_NS"]) + "\t" + calsum(sample_dict[name]["COUNT_HET_NS"]) + "\t" + 
	calsum(sample_dict[name]["COUNT_ALL_NS"]) + "\n")
	out_file_countNS.close()

	out_file_counteff = open(out_file_counteff_name, "a")
	out_file_counteff.write(name + "\t" + calsum(sample_dict[name]["GRANTHAM"]) + "\t" + calsum(sample_dict[name]["COUNT_HOMO_STOP"]) + "\t" + calsum(sample_dict[name]["COUNT_HET_STOP"]) + "\t" + 
	calsum(sample_dict[name]["COUNT_ALL_STOP"]) + "\n")
	out_file_counteff.close()

	out_file_countg = open(out_file_countg_name, "a")
	out_file_countg.write(name + "\t" + calsum(sample_dict[name]["GRANTHAM"])+ "\t" + calsum(sample_dict[name]["COUNT_HOMO_CONS"]) + "\t" + calsum(sample_dict[name]["COUNT_HET_CONS"]) + "\t" + 
	calsum(sample_dict[name]["COUNT_ALL_CONS"]) + "\t" + calsum(sample_dict[name]["COUNT_HOMO_MODC"]) + "\t" + calsum(sample_dict[name]["COUNT_HET_MODC"]) + "\t" + 
	calsum(sample_dict[name]["COUNT_ALL_MODC"]) + "\t" + calsum(sample_dict[name]["COUNT_HOMO_MODR"]) + "\t" + calsum(sample_dict[name]["COUNT_HET_MODR"]) + "\t" + 
	calsum(sample_dict[name]["COUNT_ALL_MODR"]) + "\t" + calsum(sample_dict[name]["COUNT_HOMO_RAD"]) + "\t" + calsum(sample_dict[name]["COUNT_HET_RAD"]) + "\t" + 
	calsum(sample_dict[name]["COUNT_ALL_RAD"]) + "\n")
	out_file_countg.close()

