#!/usr/bin/python3

"""
Main secondary metabolism pipeline. Imports modules and processes data to write a dataframe which can later be processed with R.
"""

import os
import csv
import sys
import argparse
import pandas as pd
import shutil

from smModule.smServerSide import tmpSmBiTable, createBidirSmurf, mysqlSmChecker
from smModule.aspSMDl import dlSMdata
import smModule.bioSlim3 as bio # tupleToFasta
# from smModule.processMibig3 import processMibig, dlSmurfProteins, writeMibigFormatted, processBlastResult, dlMibig

from smModule.misc import readConfig



parser=argparse.ArgumentParser(description='''
Script to execute the seconary metabolite analysis pipeline. In case you want to leave out some analysis you have to modify the script.\n

Example:

python3 MAIN.py -o nigri_orgs.txt -bibase biblast_table -biFinal smurf_bidir_hits_test -t species_tree.nwk -l run.log -od out_dir
''')
parser.add_argument("--orgs", "-o",
					dest="filename",
					required=True,
					help="Input file with jgi names (one per row) of organisms", metavar="FILE")
parser.add_argument("--treeFile", "-t",
					dest="tree",
					required=True,
					help="Tree file in newick format",
					metavar="FILE")
parser.add_argument("--biblastTable", "-bibase",
					dest="bibase",
					required=True,
					help="Specify the original biblast table to use as basis for the smurf bidirectional hits table", metavar="CHAR")
parser.add_argument("--clusterBiblast", "-biFinal",
					dest="biFinal",
					required=True,
					help="Specify name for smurf bidirectional hits table", metavar="CHAR")
parser.add_argument("--log", "-l",
					dest="logFile",
					required=True,
					help="Specify the name for log file", metavar="CHAR")
parser.add_argument("--outdir", "-od",
					dest="sn",
					required=True,
					help="Output directory, preferably the name of your set", metavar="CHAR")
parser.add_argument("--otherBlast", "-ob",
					dest="otherBlast",
					action="store_true",
					default=False,
					help="If problems are encountered with blast use this command for another query")
# parser.add_argument("--fetch_raw_data", "-fr",
# 					dest="fetch_raw_data",
# 					action="store_true",
# 					default=False,
# 					help="This command uses non subsetted tables")
parser.add_argument("--cluster_only", "-co",
					dest="cluster_only",
					action="store_true",
					default=False,
					help="If you only want to rerun clustering, e.g. after updating blast data, choose this option")
args=parser.parse_args()




with open("configNew.txt") as c:
	config = readConfig(c.readlines())


########
# LOADING ORGS

filename = args.filename 

biblastBaseTable = args.bibase
smurfBidirHitsName = args.biFinal 
testLogName = args.logFile 
treeFile = args.tree 
setName = args.sn 

with open(filename, "r") as tmp:
	orgSet = [item.strip() for item in tmp.readlines()]

if setName not in os.listdir():
	os.mkdir(setName)

else:
	input("A folder for the specified set is already available. If you want to rerun the analysis press Enter, else ctrl+c/d\n")

shutil.copy2(treeFile, setName)

def download_and_processing():
	os.chdir(setName)
	# Checking data
	mysqlSmChecker(orgSet, testLogName)

	# Creating smurf bidir hits table for dataset.
	tmpSmBiTable(smtable = biblastBaseTable) # Creating a temporary table


	createBidirSmurf(smtable = smurfBidirHitsName) # Creating a sm cluster table
	# 

	# DOWNLOADING SM DATA
	input("Did you create a smurfbiblast table? If yes, press Enter, if not, ctrl+c/d")


	smIpGf = dlSMdata(orgSet)

	smIpGf.to_csv("sm_data_"+setName+".tsv", sep = '\t', encoding = "UTF-8", index = False, na_rep='none')



def clustering_and_output():
	if os.getcwd().split("/")[-1] != setName:
		os.chdir(setName) 

	cluster_blast_all =bio.dbFetch("""SELECT *, ROUND(
                 ( COALESCE( ( pident_tailoring /(clust_size-q_max_bb)  )*0.35,0)+
                 COALESCE( (pident_bb/q_max_bb)*0.65, 0)
                 ),2) AS pident_score FROM (

                 SELECT bidir.q_org, bidir.q_clust_id, bidir.h_org, bidir.h_clust_id,
                 SUM(CASE WHEN sm_short != 'none' then pident else 0 end) AS pident_bb,
                 smurf.clust_size,
                 COALESCE(SUM(sm_short != 'none'),0) AS count_bb,
                 tqmax.q_max_bb,
                 SUM(CASE WHEN sm_short = 'none' then pident else 0 end) AS pident_tailoring,
                 COALESCE(SUM(sm_short = 'none'),0) AS count_tailoring

                 FROM (
                 SELECT bia.* FROM (
                 SELECT * FROM %s) bia
                 JOIN (
                 SELECT * FROM %s ) bib
                 ON bia.h_clust_id = bib.q_clust_id
                 WHERE bia.q_clust_id = bib.h_clust_id
                 GROUP BY bia.q_clust_id, bia.q_protein_id, bia.h_clust_id, bia.h_protein_id
                 ) AS bidir

                 JOIN smurf
                 ON bidir.q_org = smurf.org_id AND bidir.q_protein_id = smurf.sm_protein_id AND bidir.q_clust_id != bidir.h_clust_id

                 JOIN (SELECT CONCAT(org_id, '_' , clust_backbone,'_', clust_size) AS q_clust_id, SUM(sm_short != 'none') AS q_max_bb FROM smurf
                 GROUP BY q_clust_id) tqmax
                 ON bidir.q_clust_id = tqmax.q_clust_id

                 GROUP BY q_clust_id, h_clust_id ) ta;""" %  ("publication_smurf_bidir_hits_nigri", "publication_smurf_bidir_hits_nigri"))
	with open("clusterBlastAll.csv", "w") as handle:
		cluster_writer = csv.writer(handle, delimiter = ",")
		cluster_writer.writerow(["q_org", "q_clust_id", "h_org", "h_clust_id", "pident_bb", "clust_size", "count_bb", "q_max_bb", "pident_tailoring", "count_tailoring", "pident_score"])
		cluster_writer.writerows(cluster_blast_all)
	
	smFile = "sm_data_"+setName+".tsv"


	os.chdir("..")

	cmdString = "Rscript clusterData.R %s %s" % (smFile, setName) # File will be renamed with a _c extension
	print(cmdString)
	print("Calculating cluster families")

	os.system(cmdString)


	print("Searching for unique secondary metabolic gene clusters in tree")
	cmdString = "Rscript uniquesAtNodes.R %s %s %s" % (smFile.replace(".tsv","")+"_c.tsv", setName, treeFile)

	os.system(cmdString)


if __name__ == '__main__' and not args.cluster_only:
	download_and_processing()
	clustering_and_output()

if __name__ == '__main__' and args.cluster_only:
	clustering_and_output()
