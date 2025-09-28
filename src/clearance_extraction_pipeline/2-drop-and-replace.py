# --------------------------------------------------------
#
#  ***Oxford Mathematical Brain Modeling Group***
#   This script is the second in the clearance pipeline
#   and does two things
#
#   1. drops undesired fields from patient data CSV files.
#
#   2. [optional] replaces IDs with different IDs according to
#       a renaming dictionary
#
#   The format of the processed file is assumed to be
#   the following:
#
#   Row 1: Header (kept)
#   Rows 2 to N: [Region name] data, data, data
#
#   A row with region [Region name] is kept if [Region Name]
#   is not indicated as a dropped field
#
#   [Region Name] can also be replaced with a new name if
#       it exists in the renaming dictionary
#
#  Authors:
#  ================================================
#       Georgia S. Brennan      brennan@maths.ox.ac.uk
#                   ----
#       Travis B. Thompson      thompsont@maths.ox.ac.uk
#                   ----
#       Marie E. Rognes         meg@simula.no
#                   ----
#       Vegard Vinje            vegard@simula.no
#                   ----
#       Alain Goriely           goriely@maths.ox.ac.uk
# ---------------------------------------------------------

import os
import csv
import sys
import shutil


def writeRenamedOnly(input, output, renamingList):
    line = 0

    with open(input) as incsv:
        csv_reader = csv.reader(incsv, delimiter=',')

        with open(output, mode='w') as outcsv:
            csv_writer = csv.writer(outcsv, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)

            for row in csv_reader:
                if line == 0:
                    # write the header directly
                    csv_writer.writerow(row)
                else:
                    # The region should be the first field of the row
                    region = row[0]
                    if region in renamingList:
                        newname = renamingList[region]
                        row[0] = newname
                        csv_writer.writerow(row)
                line += 1


# ------------------------------------------------------------------------------------------------------------
#                                                Configuration
# ------------------------------------------------------------------------------------------------------------

# --..--..--..--.. Input / Output ..--..--..--..--
# relative directory where the original patient files reside
inputdirectory = "./reformatted-data/1-extracted-mean/"

# relative directory where you want the reformatted patient files to go
#   Note: You should create this directory if it does not already exist
outputdirectory = "./reformatted-data/2-dropped-fields/"

# --..--..--..--.. Dropped Fields ..--..--..--..--
droplist = ['Left-Cerebral-White-Matter',
            'Left-Lateral-Ventricle',
            'Left-Inf-Lat-Vent',
            'Left-Cerebellum-White-Matter',
            'Left-Cerebellum-Cortex',
            '3rd-Ventricle',
            '4th-Ventricle',
            'CSF',
            'Left-VentralDC',
            'Left-vessel',
            'Left-choroid-plexus',
            'Right-Cerebral-White-Matter',
            'Right-Lateral-Ventricle',
            'Right-Inf-Lat-Vent',
            'Right-Cerebellum-White-Matter',
            'Right-Cerebellum-Cortex',
            'Right-VentralDC',
            'Right-vessel',
            'Right-choroid-plexus',
            'WM-hypointensities',
            'Optic-Chiasm',
            'CC_Posterior',
            'CC_Mid_Posterior',
            'CC_Central',
            'CC_Mid_Anterior',
            'CC_Anterior',
            'ctx-lh-unknown',
            'ctx-rh-unknown',
            'wm-lh-bankssts',
            'wm-lh-caudalanteriorcingulate',
            'wm-lh-caudalmiddlefrontal',
            'wm-lh-cuneus',
            'wm-lh-entorhinal',
            'wm-lh-fusiform',
            'wm-lh-inferiorparietal',
            'wm-lh-inferiortemporal',
            'wm-lh-isthmuscingulate',
            'wm-lh-lateraloccipital',
            'wm-lh-lateralorbitofrontal',
            'wm-lh-lingual',
            'wm-lh-medialorbitofrontal',
            'wm-lh-middletemporal',
            'wm-lh-parahippocampal',
            'wm-lh-paracentral',
            'wm-lh-parsopercularis',
            'wm-lh-parsorbitalis',
            'wm-lh-parstriangularis',
            'wm-lh-pericalcarine',
            'wm-lh-postcentral',
            'wm-lh-posteriorcingulate',
            'wm-lh-precentral',
            'wm-lh-precuneus',
            'wm-lh-rostralanteriorcingulate',
            'wm-lh-rostralmiddlefrontal',
            'wm-lh-superiorfrontal',
            'wm-lh-superiorparietal',
            'wm-lh-superiortemporal',
            'wm-lh-supramarginal',
            'wm-lh-frontalpole',
            'wm-lh-temporalpole',
            'wm-lh-transversetemporal',
            'wm-lh-insula',
            'wm-rh-bankssts',
            'wm-rh-caudalanteriorcingulate',
            'wm-rh-caudalmiddlefrontal',
            'wm-rh-cuneus',
            'wm-rh-entorhinal',
            'wm-rh-fusiform',
            'wm-rh-inferiorparietal',
            'wm-rh-inferiortemporal',
            'wm-rh-isthmuscingulate',
            'wm-rh-lateraloccipital',
            'wm-rh-lateralorbitofrontal',
            'wm-rh-lingual',
            'wm-rh-medialorbitofrontal',
            'wm-rh-middletemporal',
            'wm-rh-parahippocampal',
            'wm-rh-paracentral',
            'wm-rh-parsopercularis',
            'wm-rh-parsorbitalis',
            'wm-rh-parstriangularis',
            'wm-rh-pericalcarine',
            'wm-rh-postcentral',
            'wm-rh-posteriorcingulate',
            'wm-rh-precentral',
            'wm-rh-precuneus',
            'wm-rh-rostralanteriorcingulate',
            'wm-rh-rostralmiddlefrontal',
            'wm-rh-superiorfrontal',
            'wm-rh-superiorparietal',
            'wm-rh-superiortemporal',
            'wm-rh-supramarginal',
            'wm-rh-frontalpole',
            'wm-rh-temporalpole',
            'wm-rh-transversetemporal',
            'wm-rh-insula',
            'Left-UnsegmentedWhiteMatter',
            'Right-UnsegmentedWhiteMatter']

# --..--..--..--.. Replaced Field Names ..--..--..--..--
rename = {'Left-Thalamus-Proper': 'subcortical.Left-Thalamus-Proper.left',
          'Left-Caudate': 'subcortical.Left-Caudate.left',
          'Left-Putamen': 'subcortical.Left-Putamen.left',
          'Left-Pallidum': 'subcortical.Left-Pallidum.left',
          'Brain-Stem': 'subcortical.brainstem.right',
          'Left-Hippocampus': 'subcortical.Left-Hippocampus.left',
          'Left-Amygdala': 'subcortical.Left-Amygdala.left',
          'Left-Accumbens-area': 'subcortical.Left-Accumbens-area.left',
          'Right-Thalamus-Proper': 'subcortical.Right-Thalamus-Proper.right',
          'Right-Caudate': 'subcortical.Right-Caudate.right',
          'Right-Putamen': 'subcortical.Right-Putamen.right',
          'Right-Pallidum': 'subcortical.Right-Pallidum.right',
          'Right-Hippocampus': 'subcortical.Right-Hippocampus.right',
          'Right-Amygdala': 'subcortical.Right-Amygdala.right',
          'Right-Accumbens-area': 'subcortical.Right-Accumbens-area.right',
          'ctx-lh-bankssts': 'cortical.bankssts.left',
          'ctx-lh-caudalanteriorcingulate': 'cortical.caudalanteriorcingulate.left',
          'ctx-lh-caudalmiddlefrontal': 'cortical.caudalmiddlefrontal.left',
          'ctx-lh-cuneus': 'cortical.cuneus.left',
          'ctx-lh-entorhinal': 'cortical.entorhinal.left',
          'ctx-lh-fusiform': 'cortical.fusiform.left',
          'ctx-lh-inferiorparietal': 'cortical.inferiorparietal.left',
          'ctx-lh-inferiortemporal': 'cortical.inferiortemporal.left',
          'ctx-lh-isthmuscingulate': 'cortical.isthmuscingulate.left',
          'ctx-lh-lateraloccipital': 'cortical.lateraloccipital.left',
          'ctx-lh-lateralorbitofrontal': 'cortical.lateralorbitofrontal.left',
          'ctx-lh-lingual': 'cortical.lingual.left',
          'ctx-lh-medialorbitofrontal': 'cortical.medialorbitofrontal.left',
          'ctx-lh-middletemporal': 'cortical.middletemporal.left',
          'ctx-lh-parahippocampal': 'cortical.parahippocampal.left',
          'ctx-lh-paracentral': 'cortical.paracentral.left',
          'ctx-lh-parsopercularis': 'cortical.parsopercularis.left',
          'ctx-lh-parsorbitalis': 'cortical.parsorbitalis.left',
          'ctx-lh-parstriangularis': 'cortical.parstriangularis.left',
          'ctx-lh-pericalcarine': 'cortical.pericalcarine.left',
          'ctx-lh-postcentral': 'cortical.postcentral.left',
          'ctx-lh-posteriorcingulate': 'cortical.posteriorcingulate.left',
          'ctx-lh-precentral': 'cortical.precentral.left',
          'ctx-lh-precuneus': 'cortical.precuneus.left',
          'ctx-lh-rostralanteriorcingulate': 'cortical.rostralanteriorcingulate.left',
          'ctx-lh-rostralmiddlefrontal': 'cortical.rostralmiddlefrontal.left',
          'ctx-lh-superiorfrontal': 'cortical.superiorfrontal.left',
          'ctx-lh-superiorparietal': 'cortical.superiorparietal.left',
          'ctx-lh-superiortemporal': 'cortical.superiortemporal.left',
          'ctx-lh-supramarginal': 'cortical.supramarginal.left',
          'ctx-lh-frontalpole': 'cortical.frontalpole.left',
          'ctx-lh-temporalpole': 'cortical.temporalpole.left',
          'ctx-lh-transversetemporal': 'cortical.transversetemporal.left',
          'ctx-lh-insula': 'cortical.insula.left',
          'ctx-rh-bankssts': 'cortical.bankssts.right',
          'ctx-rh-caudalanteriorcingulate': 'cortical.caudalanteriorcingulate.right',
          'ctx-rh-caudalmiddlefrontal': 'cortical.caudalmiddlefrontal.right',
          'ctx-rh-cuneus': 'cortical.cuneus.right',
          'ctx-rh-entorhinal': 'cortical.entorhinal.right',
          'ctx-rh-fusiform': 'cortical.fusiform.right',
          'ctx-rh-inferiorparietal': 'cortical.inferiorparietal.right',
          'ctx-rh-inferiortemporal': 'cortical.inferiortemporal.right',
          'ctx-rh-isthmuscingulate': 'cortical.isthmuscingulate.right',
          'ctx-rh-lateraloccipital': 'cortical.lateraloccipital.right',
          'ctx-rh-lateralorbitofrontal': 'cortical.lateralorbitofrontal.right',
          'ctx-rh-lingual': 'cortical.lingual.right',
          'ctx-rh-medialorbitofrontal': 'cortical.medialorbitofrontal.right',
          'ctx-rh-middletemporal': 'cortical.middletemporal.right',
          'ctx-rh-parahippocampal': 'cortical.parahippocampal.right',
          'ctx-rh-paracentral': 'cortical.paracentral.right',
          'ctx-rh-parsopercularis': 'cortical.parsopercularis.right',
          'ctx-rh-parsorbitalis': 'cortical.parsorbitalis.right',
          'ctx-rh-parstriangularis': 'cortical.parstriangularis.right',
          'ctx-rh-pericalcarine': 'cortical.pericalcarine.right',
          'ctx-rh-postcentral': 'cortical.postcentral.right',
          'ctx-rh-posteriorcingulate': 'cortical.posteriorcingulate.right',
          'ctx-rh-precentral': 'cortical.precentral.right',
          'ctx-rh-precuneus': 'cortical.precuneus.right',
          'ctx-rh-rostralanteriorcingulate': 'cortical.rostralanteriorcingulate.right',
          'ctx-rh-rostralmiddlefrontal': 'cortical.rostralmiddlefrontal.right',
          'ctx-rh-superiorfrontal': 'cortical.superiorfrontal.right',
          'ctx-rh-superiorparietal': 'cortical.superiorparietal.right',
          'ctx-rh-superiortemporal': 'cortical.superiortemporal.right',
          'ctx-rh-supramarginal': 'cortical.supramarginal.right',
          'ctx-rh-frontalpole': 'cortical.frontalpole.right',
          'ctx-rh-temporalpole': 'cortical.temporalpole.right',
          'ctx-rh-transversetemporal': 'cortical.transversetemporal.right',
          'ctx-rh-insula': 'cortical.insula.right'
          }

# Execution starts here
if __name__ == "__main__":

    if not os.path.exists(inputdirectory):
        print(f"The relative (to this script) input directory {inputdirectory} does not exist")
        sys.exit()

    # clear any existing files and remake the target output directory
    if os.path.exists(outputdirectory):
        shutil.rmtree(outputdirectory)
    os.mkdir(outputdirectory)

    dirlevel = 0

    # ----------- Read in and parse all subject files --------------------
    for rootdir, subjectdirs, files in os.walk(inputdirectory):

        totalsubj = len(files)
        thissubj = 0
        dirlevel += 1

        # only process the top level subdirectories
        if dirlevel == 1:
            for subj in files:
                thissubj += 1
                print(f"Processing file {subj}")
                infile = inputdirectory + subj
                outfile = outputdirectory + subj
                writeRenamedOnly(infile, outfile, rename)

    if not os.path.exists(outputdirectory):
        print(f"The relative (to this script) output directory {outputdirectory} does not exist")
        sys.exit()


