from __future__ import division
import sys
import argparse
import time
start_time = time.time()

parser=argparse.ArgumentParser()
parser.add_argument("-i","--INPUT",help="The input file containing the stats for the gene types you wish to compare to the background.",type=str,required=True)
parser.add_argument("-o","--OUTPUT",help="The output generated containing resampled differences for each stat.",type=str,required=True)
parser.add_argument("-l","--LABEL",help="The label of the ortholog column, to pair with melanogaster ID.",type=str)
parser.add_argument("-d","--DATABASE",help="path to a database containing orthologs relative to a focal species",type=str)
parser.add_argument("-s","--SPECIES",help="species examined for the ortholog database, usually the species of the reference genome mapped to.",type=str)
args=parser.parse_args()

"""

"""
import os
import pandas as pd
import numpy as np
from collections import defaultdict

if args.DATABASE is not None:
	db = pd.read_table(args.DATABASE)
	db2 = db.loc[db['Species'] == args.SPECIES]
	id_db = dict(db2[["ID", args.LABEL]].copy().set_index(args.LABEL).to_dict())["ID"]
	name_db =db2[["Name", args.LABEL]].copy().set_index(args.LABEL).to_dict()["Name"]
	out = open(args.OUTPUT,"w")
	for line in open(args.INPUT):
		line = line.rstrip()
		for j in line.rstrip().split("\t")[3].split(","):
			if j in id_db:
				groups = [line.split("\t")[0],line.split("\t")[1],line.split("\t")[2],j,id_db[j],name_db[j]]
				groups.extend(line.split("\t")[4:])
				out.write("\t".join(map(str,groups)) + "\n")
			elif j not in id_db:
				groups = [line.split("\t")[0],line.split("\t")[1],line.split("\t")[2],j,"NA","NA"]
				groups.extend(line.split("\t")[4:])
				out.write("\t".join(map(str,groups)) + "\n")
	out.close()

elif args.DATABASE is None:
	out = open(args.OUTPUT,"w")
	for line in open(args.INPUT):
		line = line.rstrip()
		for j in line.rstrip().split("\t")[3].split(","):
			groups = [line.split("\t")[0],line.split("\t")[1],line.split("\t")[2],j,"NA","NA"]
			groups.extend(line.split("\t")[4:])
			out.write("\t".join(groups) + "\n")
	out.close()
