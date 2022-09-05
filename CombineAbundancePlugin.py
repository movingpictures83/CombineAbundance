# Objective:
# Combine abundance, PTR, metadata, and AMR

import pandas as pd
import numpy as np
import os
import math
import pdb
import PyPluMA
# Step 1 - Combine abundance

class CombineAbundancePlugin:
 def input(self, infile):
        inputfile = open(infile, 'r')
        self.parameters = dict()
        for line in inputfile:
            contents = line.strip().split('\t')
            self.parameters[contents[0]] = contents[1]

 def run(self):
     pass

 def output(self, outputfile):
  amr_abundance_dir = PyPluMA.prefix()+"/"+self.parameters["abundancedir"]
  out_amr = outputfile
  samples_list_file = PyPluMA.prefix()+"/"+self.parameters["sampleslist"]

  samples_list = [x.strip('\n') for x in open(samples_list_file).readlines()]
  samples_list = sorted(samples_list)

  amr_dict={"sample":[]}
  all_genes=[]

  # Get list of all genes:
  for sample in samples_list:
    # Read abundance file

    f = os.path.join(amr_abundance_dir,sample+"-abundance.txt")
    if os.path.exists(f):
        with open(f, 'r') as re:
            re.readline()
            for line in re.readlines():
                line = line.strip("\n")
                row = line.split("\t")
                gene, abundance = row[0], row[3]
                if gene not in all_genes:
                    all_genes.append(gene)

  # assign dictionary:
  for gene in all_genes:
    amr_dict[gene] = []

  # Get all values
  for sample in samples_list:
    # Read abundance file
    amr_dict["sample"].append(sample)
    f = os.path.join(amr_abundance_dir,sample+"-abundance.txt")
    if os.path.exists(f):
        genes_found = []
        with open(f, 'r') as re:
            re.readline()
            for line in re.readlines():
                line = line.strip("\n")
                row = line.split("\t")
                gene, abundance = row[0], float(row[3])
                genes_found.append(gene)
                amr_dict[gene].append(abundance)
        for gene in all_genes:
            if gene not in genes_found:
                amr_dict[gene].append(0)
    else:
        for gene in all_genes:
            amr_dict[gene].append(np.nan)

  amr_df = pd.DataFrame(amr_dict)
  amr_df.index = amr_df["sample"]
  amr_df = amr_df.loc[:, (amr_df != 0).any(axis=0)]

  amr_df = amr_df.drop("sample", axis=1)

  amr_df = amr_df.T

  amr_df.to_csv(out_amr)
