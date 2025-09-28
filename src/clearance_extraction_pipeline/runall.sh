#!/bin/bash

#------------------------------------------------------------
# Oxford Mathematical Brain Modeling Group
#
#    Clearance Pipeline: Version 2.0
#-----------------------------------------------------------
#
#
# Authors
# Georgia S. Brennan  	(brennan@maths.ox.ac.uk)
# Travis B. Thompson 	(thompsont@maths.ox.ac.uk)
# Vegard Vinje		(vegard@simula.no)
# Marie E. Rognes	(meg@simula.no)
# Alain Goriely		(goriely@maths.ox.ac.uk)
#----------------------------------------------------------

python3 ./1-extract-field.py
python3 ./2-drop-and-replace.py
python3 ./3-cull-data.py
python3 ./4-compute-clearance.py
python3 ./4b-average-computed-clearance.py

