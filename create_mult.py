#!/usr/bin/python
import sys
import os
from os import listdir
import shutil
import ROOT
import numpy as np
from ROOT import TLorentzVector
import csv
import pandas as pd
from ROOT import TFile, TTree
from rootpy.io import root_open
from rootpy.tree import Tree, TreeChain
from rootpy.plotting import Hist, Hist2D
from rootpy.extern.six.moves import range
from root_numpy import hist2array, root2array
from functions import  Delphes_analysis

if len(sys.argv) < 2:
  print " Usage: Example1.py input_Dir"
  sys.exit(1)

ROOT.gSystem.Load("/home/felipe/madanalysis5_1_5/tools/delphes/libDelphes")

inputDir = sys.argv[1]

file_names = listdir(inputDir)

result_table = os.path.join(inputDir, 'result_table')

inputs = []
for i in file_names:
	if i != 'result_table' and i != 'run_card.dat':
		inputs.append(os.path.join(inputDir, i))

res_T = {}
with open(result_table, 'r') as res_tab:
	for line in res_tab:
		k,cHW,tcHW,xsec = line.split()
		res_T[k +'.root'] = [cHW,tcHW,xsec]

######## loop through all inputs###############

for i in inputs:
	Delphes_analysis().load_analysis(i,res_T)

