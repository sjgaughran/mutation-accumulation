#!/usr/bin/python

import sys
import os
import re
import traceback
import itertools
import csv
import math
from operator import itemgetter
from collections import defaultdict

if len(sys.argv) < 3:
	print("\nUSE: python program.py snpeff.vcf file_length \nExiting--please try again! \n")
	sys.exit()

vcf_name = sys.argv[1]

number_lines = int(sys.argv[2])

number_head_lines = 0
err_num = 0

header_lines = []
vcf_lines = 0 
file_number = 0
small_vcf = None


with open(vcf_name) as vcf_f:
	for line in vcf_f:
		if line.startswith("##"):
			number_head_lines += 1
			header_lines.append(line)
		elif line.startswith("#"):
			number_head_lines += 1
			header_lines.append(line)
		else:
			block_size = (number_lines-number_head_lines)//1000
			if vcf_lines % block_size == 0:
				if small_vcf:
					small_vcf.close()
				file_number += 1
				small_vcf_name = 'snpeff_bs_{}.vcf'.format(file_number)
				small_vcf = open(small_vcf_name, 'w')
				for item in header_lines:
					small_vcf.write(item)
			vcf_lines += 1
			small_vcf.write(line)

if small_vcf:
	small_vcf.close()

