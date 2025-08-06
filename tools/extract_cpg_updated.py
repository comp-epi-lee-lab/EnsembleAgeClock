#!/usr/bin/env python


import pandas as pd
import time
import threading
import sys

#the following code does the following:

#Imports file you want to have cpg extracted
#Find overlapping with cpg with union_out_corrected
#Makes a new csv file with just the overlapping cpgs for clock analysis. 

# Define file paths (recommended that index contains some sort of sample identification)
csv_file = sys.argv[1]

#this is the overlap files
union_file = r"union_out_corrected.csv"

# Load the union and csv set of CpGs
union_cpgs = pd.read_csv(union_file, header=None)[0].tolist()
csv_analysis = pd.read_csv(csv_file)
valid_cpgs = list(set(union_cpgs) & set(csv_analysis.iloc[:,0].tolist()))

#create filtered csv with overlapping cpgs
csv_analysis_T = csv_analysis.T
csv_analysis_T.columns = csv_analysis_T.iloc[0,:]
filtered_csv = csv_analysis_T[valid_cpgs].T
filtered_csv = filtered_csv.iloc[: , 1:]

# Write the output csv file, it will be transposed for the EnsembleMeAgingClock
output_file = sys.argv[2]
filtered_csv.to_csv(output_file)

