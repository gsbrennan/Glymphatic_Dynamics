#=============================================================================
#
#	 _____
#	|     \
#	| | ) |   A Python 3  Proteopathy-on-Graphs Solver Software Suite
#	|  ___/
#	| |
#	| |rYon
#
#	----------------------------------------------------
#	PrYon is distributed under the GNU GPL V3 License
#	----------------------------------------------------
#	Mathematical Institute, Oxford University
#	Oxford, United Kingdom
#
#	Authors: Pavanjit Chaggar, T.B. Thompson
#	Copyright (c) 2020 A. Goriely, and T.B. Thompson
#
#	GNU GPL V3: https://www.gnu.org/licenses/gpl-3.0.html
#	-----------------------------------------------
import matplotlib.pyplot as plt
from pryon import *
import pandas as pd
#=============================================================================
# PrYon Example Script: Solve a Brennan-Goriely (BG) model problem on a
#			Brain Connectome Graph with PrYon
#=============================================================================





################# CHOOSE A SCALE CONNECTOME, SET UP SOLVER ##################


bg = solverBG('master-std500.graphml')


# --set error / status reporting on (True: on, False: off)
bg.setVerbose(True)

# --set up the solver--
bg.setup(0.0,300.0,1,1e-8)

# --enable VTK (paraview visualization) output
bg.enableVisualization('mybg','./output/')

# -- You can call these functions to get solver-specific information
# bg.discoverSolutionFields()	# what solution fields are supported by the solver
params = bg.discoverParameters()	# what parameter names does the solver use
labels = bg.discoverRegionLabels()	# what are the region labels available





############# MAKE DICTIONARY TO  CHANGE REGIONAL PARAMS EASILY ###################

d = dict.fromkeys(labels)
lambda_0 = d


# load the data
col_names = ['name','BNdata']

# loop through all patients - here is example for patient 228 only
data1 = pd.read_csv('228.csv', names=col_names, header=None)
voxels = pd.read_csv('voxels_allpatients_ascending.csv', names=col_names, header=None)
numregions = pd.read_csv('nodes_rois_ascending.csv',names=col_names, header=None)

datamatrix = data1.to_numpy()
weightmatrix = voxels.to_numpy()
roinodesmatrix = numregions.to_numpy()

# assign lambda_0 values to dict
for i in range(int(len(datamatrix))):
	lambda_0[datamatrix[i,0]] = datamatrix[i,1]


############################ GENERATE STAT FILES ########################################

# -- Add output for the global (full connectome) misfolded protein
#	concentration solver results
bg.addGlobalResult('Misfolded-Protein-Concentration')
bg.addGlobalResult('Clearance')
bg.addGlobalResult('Damage-Percentage')

# -- regional stat files --
# here is code for misflded protein concentration - can also do 'Clearance' and 'Damage-Percentage' regionally by changing the option
# ADNI BRAAK REGIONS

# Braak 1
bg.addRegionalResult(['cortical.entorhinal.right','cortical.entorhinal.left'],'Misfolded-Protein-Concentration', optionalIdTag='braak1')

# Braak 2
bg.addRegionalResult(['subcortical.Right-Hippocampus.right','subcortical.Left-Hippocampus.left'],'Misfolded-Protein-Concentration', optionalIdTag='braak2')

# Braak 3
bg.addRegionalResult(['cortical.parahippocampal.right','cortical.parahippocampal.left','cortical.fusiform.right','cortical.fusiform.left','cortical.lingual.right','cortical.lingual.left','subcortical.Left-Amygdala.left','subcortical.Right-Amygdala.right'],'Misfolded-Protein-Concentration', optionalIdTag='braak3')

# Braak 4
bg.addRegionalResult(['cortical.rostralanteriorcingulate.right','cortical.caudalanteriorcingulate.right','cortical.rostralanteriorcingulate.left','cortical.caudalanteriorcingulate.left','cortical.middletemporal.left','cortical.middletemporal.right','cortical.posteriorcingulate.left','cortical.posteriorcingulate.right','cortical.isthmuscingulate.right','cortical.isthmuscingulate.left','cortical.insula.right','cortical.insula.left','cortical.inferiortemporal.right','cortical.inferiortemporal.left','cortical.temporalpole.right','cortical.temporalpole.left'],'Misfolded-Protein-Concentration', optionalIdTag='braak4')

# Braak 5 and 6
bg.addRegionalResult(['cortical.lateraloccipital.right','cortical.lateraloccipital.left','cortical.superiorfrontal.left','cortical.superiorfrontal.right','cortical.lateralorbitofrontal.left','cortical.lateralorbitofrontal.right','cortical.medialorbitofrontal.left','cortical.medialorbitofrontal.right','cortical.frontalpole.left','cortical.frontalpole.right','cortical.caudalmiddlefrontal.left','cortical.caudalmiddlefrontal.right','cortical.rostralmiddlefrontal.right','cortical.rostralmiddlefrontal.left','cortical.parsopercularis.right','cortical.parsopercularis.left','cortical.parsorbitalis.right','cortical.parsorbitalis.left','cortical.parstriangularis.left','cortical.parstriangularis.right','cortical.supramarginal.right','cortical.supramarginal.left','cortical.inferiorparietal.right','cortical.inferiorparietal.left','cortical.superiortemporal.right','cortical.superiortemporal.left','cortical.superiorparietal.right','cortical.superiorparietal.left','cortical.precuneus.right','cortical.precuneus.left','cortical.bankssts.right','cortical.bankssts.left','cortical.transversetemporal.right','cortical.transversetemporal.left'],'Misfolded-Protein-Concentration', optionalIdTag='braak5')

bg.addRegionalResult(['cortical.cuneus.right','cortical.pericalcarine.right','cortical.cuneus.left','cortical.pericalcarine.left','cortical.postcentral.left','cortical.postcentral.right','cortical.precentral.left','cortical.precentral.right','cortical.paracentral.left','cortical.paracentral.right'],'Misfolded-Protein-Concentration', optionalIdTag='braak6')


# LOBES

# Frontal
bg.addRegionalResult(['cortical.lateralorbitofrontal.right','cortical.parsorbitalis.right','cortical.frontalpole.right','cortical.medialorbitofrontal.right','cortical.parstriangularis.right','cortical.parsopercularis.right','cortical.rostralmiddlefrontal.right','cortical.superiorfrontal.right','cortical.caudalmiddlefrontal.right','cortical.precentral.right'],'Misfolded-Protein-Concentration', optionalIdTag='frontal')

# Parietal
bg.addRegionalResult(['cortical.paracentral.right','cortical.postcentral.right','cortical.supramarginal.right','cortical.superiorparietal.right','cortical.inferiorparietal.right','cortical.precuneus.right'],'Misfolded-Protein-Concentration', optionalIdTag='parietal')

# Limbic
bg.addRegionalResult(['cortical.rostralanteriorcingulate.right','cortical.caudalanteriorcingulate.right','cortical.posteriorcingulate.right','cortical.isthmuscingulate.right','cortical.parahippocampal.right','cortical.entorhinal.right'],'Misfolded-Protein-Concentration', optionalIdTag='limbic')

# Occipital
bg.addRegionalResult(['cortical.precuneus.right','cortical.pericalcarine.right','cortical.lateraloccipital.right','cortical.lingual.right'],'Misfolded-Protein-Concentration', optionalIdTag='occipital')

# Temporal
bg.addRegionalResult(['cortical.fusiform.right','cortical.temporalpole.right','cortical.inferiortemporal.right','cortical.middletemporal.right','cortical.bankssts.right','cortical.superiortemporal.right','cortical.transversetemporal.right','cortical.insula.right','subcortical.Right-Hippocampus.right'],'Misfolded-Protein-Concentration', optionalIdTag='temporal')

# Basal Ganglia
bg.addRegionalResult(['subcortical.Right-Thalamus-Proper.right','subcortical.Right-Caudate.right','subcortical.Right-Putamen.right','subcortical.Right-Pallidum.right','subcortical.Right-Accumbens-area.right','subcortical.Right-Amygdala.right'],'Misfolded-Protein-Concentration', optionalIdTag='basal')


# Generate stat files for each patient for all 83 ROIs:

bg.addRegionalResult(['subcortical.Brain-Stem.left'],'Misfolded-Protein-Concentration', optionalIdTag='subcortical.Brain-Stem.left')
bg.addRegionalResult(['subcortical.Left-Accumbens-area.left'],'Misfolded-Protein-Concentration', optionalIdTag='subcortical.Left-Accumbens-area.left')
bg.addRegionalResult(['subcortical.Left-Caudate.left'],'Misfolded-Protein-Concentration', optionalIdTag='subcortical.Left-Caudate.left')
bg.addRegionalResult(['subcortical.Left-Pallidum.left'],'Misfolded-Protein-Concentration', optionalIdTag='subcortical.Left-Pallidum.left')
bg.addRegionalResult(['subcortical.Left-Putamen.left'],'Misfolded-Protein-Concentration', optionalIdTag='subcortical.Left-Putamen.left')
bg.addRegionalResult(['subcortical.Left-Thalamus-Proper.left'],'Misfolded-Protein-Concentration', optionalIdTag='subcortical.Left-Thalamus-Proper.left')
bg.addRegionalResult(['subcortical.Right-Accumbens-area.right'],'Misfolded-Protein-Concentration', optionalIdTag='subcortical.Right-Accumbens-area.right')
bg.addRegionalResult(['subcortical.Right-Caudate.right'],'Misfolded-Protein-Concentration', optionalIdTag='subcortical.Right-Caudate.right')
bg.addRegionalResult(['subcortical.Right-Pallidum.right'],'Misfolded-Protein-Concentration', optionalIdTag='subcortical.Right-Pallidum.right')
bg.addRegionalResult(['subcortical.Right-Putamen.right'],'Misfolded-Protein-Concentration', optionalIdTag='subcortical.Right-Putamen.right')
bg.addRegionalResult(['subcortical.Right-Thalamus-Proper.right'],'Misfolded-Protein-Concentration', optionalIdTag='subcortical.Right-Thalamus-Proper.right')
bg.addRegionalResult(['cortical.entorhinal.right'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.entorhinal.right')
bg.addRegionalResult(['cortical.entorhinal.left'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.entorhinal.left')
bg.addRegionalResult(['subcortical.Right-Hippocampus.right'],'Misfolded-Protein-Concentration', optionalIdTag='subcortical.Right-Hippocampus.right')
bg.addRegionalResult(['subcortical.Left-Hippocampus.left'],'Misfolded-Protein-Concentration', optionalIdTag='subcortical.Left-Hippocampus.left')
bg.addRegionalResult(['cortical.parahippocampal.right'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.parahippocampal.right')
bg.addRegionalResult(['cortical.parahippocampal.left'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.parahippocampal.left')
bg.addRegionalResult(['cortical.fusiform.right'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.fusiform.right')
bg.addRegionalResult(['cortical.fusiform.left'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.fusiform.left')
bg.addRegionalResult(['cortical.lingual.right'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.lingual.right')
bg.addRegionalResult(['cortical.lingual.left'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.lingual.left')
bg.addRegionalResult(['subcortical.Left-Amygdala.left'],'Misfolded-Protein-Concentration', optionalIdTag='subcortical.Left-Amygdala.left')
bg.addRegionalResult(['subcortical.Right-Amygdala.right'],'Misfolded-Protein-Concentration', optionalIdTag='subcortical.Right-Amygdala.right')
bg.addRegionalResult(['cortical.rostralanteriorcingulate.right'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.rostralanteriorcingulate.right')
bg.addRegionalResult(['cortical.caudalanteriorcingulate.right'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.caudalanteriorcingulate.right')
bg.addRegionalResult(['cortical.rostralanteriorcingulate.left'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.rostralanteriorcingulate.left')
bg.addRegionalResult(['cortical.caudalanteriorcingulate.left'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.caudalanteriorcingulate.left')
bg.addRegionalResult(['cortical.middletemporal.left'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.middletemporal.left')
bg.addRegionalResult(['cortical.middletemporal.right'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.middletemporal.right')
bg.addRegionalResult(['cortical.posteriorcingulate.left'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.posteriorcingulate.left')
bg.addRegionalResult(['cortical.posteriorcingulate.right'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.posteriorcingulate.right')
bg.addRegionalResult(['cortical.isthmuscingulate.right'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.isthmuscingulate.right')
bg.addRegionalResult(['cortical.isthmuscingulate.left'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.isthmuscingulate.left')
bg.addRegionalResult(['cortical.insula.right'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.insula.right')
bg.addRegionalResult(['cortical.insula.left'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.insula.left')
bg.addRegionalResult(['cortical.inferiortemporal.right'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.inferiortemporal.right')
bg.addRegionalResult(['cortical.inferiortemporal.left'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.inferiortemporal.left')
bg.addRegionalResult(['cortical.temporalpole.right'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.temporalpole.right')
bg.addRegionalResult(['cortical.temporalpole.left'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.temporalpole.left')
bg.addRegionalResult(['cortical.cuneus.right'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.cuneus.right')
bg.addRegionalResult(['cortical.pericalcarine.right'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.pericalcarine.right')
bg.addRegionalResult(['cortical.lateraloccipital.right'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.lateraloccipital.right')
bg.addRegionalResult(['cortical.cuneus.left'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.cuneus.left')
bg.addRegionalResult(['cortical.pericalcarine.left'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.pericalcarine.left')
bg.addRegionalResult(['cortical.lateraloccipital.left'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.lateraloccipital.left')
bg.addRegionalResult(['cortical.postcentral.left'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.postcentral.left')
bg.addRegionalResult(['cortical.postcentral.right'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.postcentral.right')
bg.addRegionalResult(['cortical.precentral.left'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.precentral.left')
bg.addRegionalResult(['cortical.precentral.right'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.precentral.right')
bg.addRegionalResult(['cortical.paracentral.left'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.paracentral.left')
bg.addRegionalResult(['cortical.paracentral.right'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.paracentral.right')
bg.addRegionalResult(['cortical.superiorfrontal.left'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.superiorfrontal.left')
bg.addRegionalResult(['cortical.superiorfrontal.right'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.superiorfrontal.right')
bg.addRegionalResult(['cortical.lateralorbitofrontal.left'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.lateralorbitofrontal.left')
bg.addRegionalResult(['cortical.lateralorbitofrontal.right'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.lateralorbitofrontal.right')
bg.addRegionalResult(['cortical.medialorbitofrontal.left'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.medialorbitofrontal.left')
bg.addRegionalResult(['cortical.medialorbitofrontal.right'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.medialorbitofrontal.right')
bg.addRegionalResult(['cortical.frontalpole.left'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.frontalpole.left')
bg.addRegionalResult(['cortical.frontalpole.right'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.frontalpole.right')
bg.addRegionalResult(['cortical.caudalmiddlefrontal.left'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.caudalmiddlefrontal.left')
bg.addRegionalResult(['cortical.caudalmiddlefrontal.right'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.caudalmiddlefrontal.right')
bg.addRegionalResult(['cortical.rostralmiddlefrontal.right'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.rostralmiddlefrontal.right')
bg.addRegionalResult(['cortical.rostralmiddlefrontal.left'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.rostralmiddlefrontal.left')
bg.addRegionalResult(['cortical.parsopercularis.right'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.parsopercularis.right')
bg.addRegionalResult(['cortical.parsopercularis.left'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.parsopercularis.left')
bg.addRegionalResult(['cortical.parsorbitalis.right'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.parsorbitalis.right')
bg.addRegionalResult(['cortical.parsorbitalis.left'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.parsorbitalis.left')
bg.addRegionalResult(['cortical.parstriangularis.left'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.parstriangularis.left')
bg.addRegionalResult(['cortical.parstriangularis.right'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.parstriangularis.right')
bg.addRegionalResult(['cortical.supramarginal.right'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.supramarginal.right')
bg.addRegionalResult(['cortical.supramarginal.left'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.supramarginal.left')
bg.addRegionalResult(['cortical.inferiorparietal.right'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.inferiorparietal.right')
bg.addRegionalResult(['cortical.inferiorparietal.left'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.inferiorparietal.left')
bg.addRegionalResult(['cortical.superiortemporal.right'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.superiortemporal.right')
bg.addRegionalResult(['cortical.superiortemporal.left'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.superiortemporal.left')
bg.addRegionalResult(['cortical.superiorparietal.right'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.superiorparietal.right')
bg.addRegionalResult(['cortical.superiorparietal.left'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.superiorparietal.left')
bg.addRegionalResult(['cortical.precuneus.right'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.precuneus.right')
bg.addRegionalResult(['cortical.precuneus.left'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.precuneus.left')
bg.addRegionalResult(['cortical.bankssts.right'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.bankssts.right')
bg.addRegionalResult(['cortical.bankssts.left'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.bankssts.left')
bg.addRegionalResult(['cortical.transversetemporal.right'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.transversetemporal.right')
bg.addRegionalResult(['cortical.transversetemporal.left'],'Misfolded-Protein-Concentration', optionalIdTag='cortical.transversetemporal.left')






################## ASSIGN UNIFORM PARAMETER VALUES ####################################

	#set the asymptotic minimal clearance (l_{\infinity} in the equations) to 1e-6 in all regions
bg.setUniformParameter('Asymptotic-Minimal-Clearance',1e-6)
 	#set the critical clearance (mu in the equations) to 0.72 in all regions
bg.setUniformParameter('Critical-Clearance',0.72)
	# set the Diffusion coefficient (rho in the equations) to 1e-2 in all regions
#bg.setUniformParameter('Diffusion-Coefficient',1e-2)
	# set the linear growth coefficient (the coefficient G in G*(mu-lambda) in the equations) to 2.5 in all regions
bg.setUniformParameter('Linear-Growth-Coefficient',1)
	# set the non-local (deafferentation) degradation coefficient (tau in the damage and clearance equations) to 0.0 in all regions
bg.setUniformParameter('Nonlocal-Degradation-Rate',0.0)
	# set the Saturation growth coefficient (alpha in the equations) to 2.1 in all regions
bg.setUniformParameter('Saturation-Growth-Coefficient',2.1)
	# set the Toxic Degradiation Rate (the coefficient B in B*p_j in the damage and clearance equations) to 1.0
bg.setUniformParameter('Toxic-Degradation-Rate',1.0)






################## ASSIGN REGIONAL PARAMETER VALUES (NOT CLEARANCE RELATED) #########

# -- set all soln fields to zero and then seed individual fields as
#	needed. We will use the discoverSolutionFields() to first return a list
#	of all solution fields and then loop over the items in the list to set
# 	each of them individually
fields = bg.discoverSolutionFields()

for x in fields:
	bg.setUniformInitialValue(x,0.0,False)



################## ASSIGN REGIONAL INITIAL CLEARANCE AND DIFFUSION CO-EFF ##############################

for i in range(int(len(datamatrix))):
	bg.setInitialValue('Clearance',datamatrix[i,0],datamatrix[i,1], False)

	bg.setParameter('Diffusion-Coefficient', weightmatrix[i,0], 1e-2*(1.0/float(weightmatrix[i,1]))*(1.0/float(roinodesmatrix[i,1])))





################################## SEEDING ###################################

# -- To seed the model, set the initial value for
#	the misfolded protein concentration in the entorhinal cortex to a higher
#	value.  We set `True' to distribute this initial value equally to all
#	nodes in the region so that the total initial value in the region
#	coincides with 0.075.

# For Tau seeding:
bg.setInitialValue('Misfolded-Protein-Concentration','cortical.entorhinal.right',0.1,False)  #75e-3,True)
bg.setInitialValue('Misfolded-Protein-Concentration','cortical.entorhinal.left',0.1,False) #75e-3,True)

# solve
bg.solve()




########### PLOT RESULTS ###################################################

plt.figure(1)

rglob = resultsReader()
rglob.loadResults('Brennan-Goriely-Model-Solver-Global-Damage-Percentage.stat')
globTimes = rglob.readTime()
globVals  = rglob.readMeanValue()
plt.plot(globTimes,globVals, label ='Damage')

rglob2 = resultsReader()
rglob2.loadResults('Brennan-Goriely-Model-Solver-Global-Clearance.stat')
globTimes2 = rglob2.readTime()
globVals2  = rglob2.readMeanValue()
plt.plot(globTimes2,globVals2, label = 'Clearance')

rglob3 = resultsReader()
rglob3.loadResults('Brennan-Goriely-Model-Solver-Global-Misfolded-Protein-Concentration.stat')
globTimes3 = rglob3.readTime()
globVals3  = rglob3.readMeanValue()

plt.plot(globTimes3,globVals3, label = 'Toxic Mass')
plt.xlabel("Time")
plt.legend(loc="upper right")
plt.show()
plt.savefig("brennan-goriely-global-clearance.png")
