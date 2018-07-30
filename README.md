# Creating a virtual environment and installing python packages

```bash
virtualenv /home/seth/virtualenvs/secMetPipeline -p /usr/bin/python3.5

pip install biopython pandas ete3 mysqlclient
```

# R packages

The following packages are required by R (R version 3.4.1) and need to be installed:
```
optparse
igraph
parallel
ape
ggtree
ggplot2
```

# Query example:

Before running the query, add credentials in configNew.txt.

```bash
python3 MAIN.py -o nigri_set.txt -bibase publication_nigri_biblast -biFinal publication_smurf_bidir_hits_nigri -t smBasicTree.nwk -l run.log -od nigri_test
```
