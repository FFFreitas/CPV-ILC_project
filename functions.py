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

class Delphes_analysis():
	def load_analysis(self,inputs,res_T):
		cHW = res_T[inputs.split('/')[-1]][0]
		tcHW = res_T[inputs.split('/')[-1]][1]
		xsec = res_T[inputs.split('/')[-1]][2]
		# Create chain of root trees
		chain1 = ROOT.TChain("Delphes")
		chain1.Add(inputs)

		# Create object of class ExRootTreeReader
		treeReader = ROOT.ExRootTreeReader(chain1)
		numberOfEntries = treeReader.GetEntries()


		# create new root file
		root_name = 'cHW_{}_tcHW_{}.root'.format(cHW,tcHW)
		csv_name = 'cHW_{}_tcHW_{}.csv'.format(cHW,tcHW)
		f = root_open(root_name, "recreate")
		tree = Tree("cHW_{}_tcHW_{}".format(cHW,tcHW))
		tree.create_branches({'PT_l1': 'F',
			'PT_l2': 'F',
			'PT_ll': 'F',
			'Cos_lZ': 'F',
			'DPHI_ll': 'F',
			'PT_j1': 'F',
			'PT_j2': 'F',
			'PT_b1': 'F',
			'PT_b2': 'F',
			'Eta_H': 'F',
			'phi_H': 'F',
			'M_H': 'F',
			'Cos_Hb1': 'F',
			'PT_H': 'F',
			'PT_ZH': 'F',
			'M_Z': 'F',
			'M_ZH': 'F',
			'cHW': 'F',
			'tcHW': 'F',
			'xsec': 'F'
			})

		# Get pointers to branches used in this analysis
		branchJet = treeReader.UseBranch("Jet")
		branchElectron = treeReader.UseBranch("Electron")
		branchMuon = treeReader.UseBranch("Muon")
		branchPhoton = treeReader.UseBranch("Photon")
		branchMET = treeReader.UseBranch("MissingET")
		# Loop over all events
		for entry in range(0, numberOfEntries):
		  # Load selected branches with data from specified event
			treeReader.ReadEntry(entry)
			muons = []
			for n in xrange(branchMuon.GetEntries()):
				muons.append(branchMuon.At(n))

			if len(muons) >= 2:
				muons = sorted(branchMuon, key=lambda Muon: Muon.P4().Pt(), reverse=True)
			else:
				continue

			missing = sorted(branchMET, key=lambda MisingET: MisingET.MET, reverse=True)
			muon1 = muons[0]
			muon2 = muons[1]
			Muon1 = ROOT.TLorentzVector()
			Muon2 = ROOT.TLorentzVector()
			Muon1.SetPtEtaPhiE(muon1.P4().Pt(),muon1.P4().Eta(),muon1.P4().Phi(),muon1.P4().E())
			Muon2.SetPtEtaPhiE(muon2.P4().Pt(),muon2.P4().Eta(),muon2.P4().Phi(),muon2.P4().E())
			met = ROOT.TLorentzVector()
			met.SetPtEtaPhiE(missing[0].P4().Pt(),missing[0].P4().Eta(),missing[0].P4().Phi(),missing[0].P4().E())
			bjato1 = ROOT.TLorentzVector()
			bjato2 = ROOT.TLorentzVector()
			jato1 = ROOT.TLorentzVector()
			jato2 = ROOT.TLorentzVector()
		####################################################################################
			bjets, ljets = [], []
			for n in xrange(branchJet.GetEntries()):
				if branchJet.At(n).BTag == 1:
				        bjets.append(branchJet.At(n))
				else:
				        ljets.append(branchJet.At(n))

			if len(bjets) >= 2:
				bjets = sorted(bjets, key=lambda BJet:  BJet.P4().Pt(), reverse=True)
			else:
				continue

			ljets = sorted(ljets, key=lambda Jet:  Jet.P4().Pt(), reverse=True)

			try:
				jato1.SetPtEtaPhiE(ljets[0].P4().Pt(),ljets[0].P4().Eta(),ljets[0].P4().Phi(),ljets[0].P4().E())
			except IndexError:
				tree.PT_j1 = -999

			try:
				jato2.SetPtEtaPhiE(ljets[1].P4().Pt(),ljets[1].P4().Eta(),ljets[1].P4().Phi(),ljets[1].P4().E())
			except IndexError:
				tree.PT_j2 = -999

		####################################################################################
			bjato1.SetPtEtaPhiE(bjets[0].P4().Pt(),bjets[0].P4().Eta(),bjets[0].P4().Phi(),bjets[0].P4().E())
			bjato2.SetPtEtaPhiE(bjets[1].P4().Pt(),bjets[1].P4().Eta(),bjets[1].P4().Phi(),bjets[1].P4().E())

		###################################################################################################
			if 95 < (bjato1 + bjato2).M() < 135:
				tree.PT_l1 = Muon1.Pt()
				tree.PT_l2 = Muon2.Pt()
				tree.PT_ll = (Muon1 + Muon2).Pt()
				tree.PT_b1 = bjato1.Pt()
				tree.PT_b2 = bjato2.Pt()
				tree.PT_j1 = jato1.Pt()
				tree.PT_j2 = jato2.Pt()
				Z = ROOT.TLorentzVector()
				H = ROOT.TLorentzVector()
				ZH = ROOT.TLorentzVector()
				Z = (Muon1 + Muon2)
				H = (bjato1 + bjato2)
				ZH = Z + H
				tree.phi_H = H.Phi()
				tree.PT_ZH = ZH.Pt()
				tree.M_ZH = ZH.M()
				tree.PT_H = H.Pt()
				tree.Eta_H = H.Eta()
				tree.M_H = H.M()
				tree.M_Z = Z.M()
				tree.DPHI_ll = np.abs(Muon1.DeltaPhi(Muon2))
		########################## boosted objects  ############################################
				Ztob = ROOT.TLorentzVector()
				Ztob.SetPxPyPzE(Z.Px(),Z.Py(),Z.Pz(),Z.E())
				Zboost = ROOT.TVector3()
				Zboost = Ztob.BoostVector()
				v = Zboost.Unit()
				Muon1.Boost(-Zboost)
				Htob = ROOT.TLorentzVector()
				Htob.SetPxPyPzE(H.Px(),H.Py(),H.Pz(),H.E())
				Hboost = ROOT.TVector3()
				Hboost = Htob.BoostVector()
				ang = Hboost.Unit()
				bjato1.Boost(-Hboost)
				tree.Cos_Hb1 = np.cos(bjato1.Angle(ang))
				tree.Cos_lZ = np.cos(Muon1.Angle(v))
		##########################################################################################
				tree.cHW = cHW
				tree.tcHW = tcHW
				tree.xsec = xsec
				tree.Fill()


		tree.write()
		f.close()

		#create the csv output

		to_convert = root2array(root_name,"cHW_{}_tcHW_{}".format(cHW,tcHW))

		df_conv = pd.DataFrame(to_convert)

		df_conv.to_csv( csv_name, index=False, header= df_conv.keys(), mode='w', sep=' ')


		### move everything
		if not os.path.exists('500GeV_res'):
			os.makedirs('500GeV_res')
			os.makedirs('500GeV_res/roots')
			os.makedirs('500GeV_res/csv')

		shutil.move(root_name, '500GeV_res/roots')
		shutil.move(csv_name, '500GeV_res/csv')

