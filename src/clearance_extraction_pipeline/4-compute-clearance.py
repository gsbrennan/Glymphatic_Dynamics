# --------------------------------------------------------
#
#  ***Oxford Mathematical Brain Modeling Group***
#
#   It is assumed that clearance in each region should fit
#   an exponential decay model although we do not strictly
#   enforce this point of view.
#
#   The general model is I(t) = a*exp(-kt) + b
#
#   We are interested in computing k.  The data for this
#   version of the pipeline is assumed to be three values
#   (~24 hours, ~48 hours, ~30 day baseline).  Data are
#   normalized to the 30 day baseline and fitted to the
#   model I(t) above using the first two (normalized)
#   points and taking b=1.
#
#   (This is the fourth script in the pipeline)
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
import shutil
import math

import numpy as np
from pylab import *
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


# returns the following sum of squares errors:
#   sum of squares error
#   mean square error
#   relative sum of squares error, i.e. ((actual-data)/actual)^2
#   relative mean square error
def sserrs(actual, fitted):

    differ = [math.pow(actual[i] - fitted[i], 2) for i in range(len(actual))]
    reldiffer = [math.pow( (actual[i] - fitted[i])/actual[i], 2) for i in range(len(actual))]

    sse = np.sum(differ)
    mse = sse / len(fitted)

    relsse = np.sum(reldiffer)
    relmse = relsse / len(fitted)

    return sse, mse, relsse, relmse

# plot the data and the (lambda)
# function f with options opts
def plotit(xdat, ydat, f, opts):
    xmin = np.min(xdat)
    xmax = np.max(xdat)

    xsynth = np.linspace(xmin, xmax)
    ysynth = f(xsynth, *opts)

    plt.plot(xsynth, ysynth, color='red', linestyle='dashed', label='model fit')
    plt.scatter(xdat,ydat, s=25, c='blue', alpha=0.6, label='data')
    plt.legend()
    #plt.show()


# this function assumes that the xvalues are three times (in days)
#   at ~24 hours, ~48 hours and ~30 days.  The y-values correspond to
#   the measurements at those times.
def fitted(xv, yv, plotres=False):

    #------------------------------------------
    def lambdalin(X, Y, plotfit=False):
        # For a linear fit, we don't fix any points
        # as this would fully determine the behaviour
        f = lambda x, a, b: np.asarray(a) * x + b
        popt, pcov = curve_fit(f, X, Y)

        if plotfit:
            plotit(X, Y, f, popt)

        # get the errors for the prediction
        ypred = f(X, *popt)
        sse, mse, relsse, relmse = sserrs(Y, ypred)

        # clearance should be the negative of the slope
        clearance = -1.0 * popt[0]

        return relmse, clearance

    # ------------------------------------------
    def lambdaexp(X, Y, increasing=False, plotfit=False):

        # fit a lambda function of the form
        # y(x) = a exp^( -k(x-x0) ) + 1 where
        # sign(k) = 1 if increasing == True
        # sign(k) =-1 if increasing == False
        # it is assumed that the data in Y is
        # decreasing in X and that it has been
        # normalized (by the last point) so
        # that the asymptotic value is 1
        f = None

        # x0 is determined by the initial point
        # of the X vector
        xz = X[0]

        # We assume that the data has been normalized
        # by a baseline value at larget time so that
        # the asymptotic value expected is 1.
        B = 1.0

        # It is assumed that the Y vector contains two
        # values: Y[0] = y(t0) and Y[1] = y(t1) where
        # t1 is not far from t0 (e.g. t1 is not an
        # asymptotic return-to-baseline time).

        # We enforce y(x0) = y0 by setting a = Y[0]-B = Y[0]-1
        A = Y[0] - B

        # Now create the lambda function to fit against
        # we only have one parameter to determine at this
        # point and that is k
        f = lambda x, k: np.asarray(A) * np.exp(np.asarray(-1.0 * k) * (x - np.asarray(xz))) + B

        # initial guess for k0
        k0 = (np.log(np.max(Y)) - np.log(np.min(Y))) / (np.max(X) - np.min(X))

        popt, pcov = curve_fit(f, X, Y, p0=(k0))

        if plotfit:
            # for debugging visualization
            # set renorm=True (in debugger)
            # to see all three points
            renorm = False
            if renorm:
                x = list(X)
                y = list(Y)
                x.append(30.0)
                y.append(1.0)
                plotit(x, y, f, popt)
            else:
                plotit(X, Y, f, popt)

        # get the errors for the prediction
        ypred = f(X, *popt)
        sse, mse, relsse, relmse = sserrs(Y, ypred)

        clearance = popt[0]

        return relmse, clearance
    #------------------------------------------

    clearance = 0.00
    modeltype = 'Exponential'

    # normalize the y value vector
    ynorm = np.asarray(yv) / np.asarray(yv[2])

    # make sure that the maximum is at the first point
    ymax = np.max(yv)
    if ymax != yv[0]:
        # the maximum occurs at 48 hours or at "baseline"
        # we fit a linear model to this type of unexpected
        # data and use the slope as the clearance.  Note,
        # this occurs mostly in the white matter or subcortical
        # areas where the signal changes very little and noise
        # tends to dominate
        error, clearance = lambdalin(xv, ynorm, plotfit=False)
        #error, clearance = lambdalin(xv[:2], ynorm[:2], plotfit=False)
        modeltype = 'Linear'
        print('degenerate type 1')
        
    elif yv[2]>yv[1]:
        error, clearance = lambdalin(xv, ynorm, plotfit=False)
        #error, clearance = lambdalin(xv[:2], ynorm[:2], plotfit=False)
        modeltype = 'Linear'
        print('degenerate type 2')

    else:
        error, clearance = lambdaexp(xv[:2], ynorm[:2], plotfit=False)

    return clearance, modeltype
#------------------------------------------
#------------------------------------------


#---------------------------------------------------------
# In version 2 of the pipeline, the data has been culled
# (by 3-cull-data.py) such that every patient dataset,
# which has progress to this point, has three values for
# every anatomical region.
#   1: Time at ~24 hours
#   2: Time at ~48 hours
#   3: Time at ~30 days (baseline)
#
#   Thus, version 2 of the clearance pipeline no longer
#   searches for maxima, etc, in order to fit based on
#   peak location.  We now normalize the first two values
#   by the baseline value and then ask two questions
#   1: Does the data decrease or increase?
#       * Determines exponential or linear model
#   2: What is the fit?
#---------------------------------------------------------

def writeClearance(input, output):
    line = 0

    headertimes = []
    filedata = {}
    clearancedata = {}
    clearancemodel = {}

    with open(input) as incsv:
        csv_reader = csv.reader(incsv, delimiter=',')

        for row in csv_reader:
            if line == 0:
                # save the header time data
                headertimes = row[1:]
            else:
                # bin the data
                filedata[row[0]] = row[1:]
            line += 1

    if len(headertimes) < 2:
        print(f"patient file {input} cannot be processed due to data paucity (at least 3 data points are needed)")
    else:
        print(headertimes)
        xvals = [float(t) for t in headertimes]

        # fit each anatomical field in the file
        for field in filedata:
            yvals = [float(d) for d in filedata[field]]
            lyv = len(yvals)

            if lyv != 3:
                # this is an error and should never happen
                print(f"V2 of the clearance pipeline requires exactly 3 data points.  Please double check patient data {input}")
                clearancedata[field] = 0.00
                clearancemodel[field] = 'Irregular Data'
            else:
                # fit the data
                clearance, modeltype = fitted(xvals, yvals, plotres=False)

                clearancedata[field] = clearance
                clearancemodel[field] = modeltype


        with open(output, mode='w') as outcsv:
            csv_writer = csv.writer(outcsv, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)

            # write the header
            row = ['StructName'] + ['Clearance', 'Model Type']
            csv_writer.writerow(row)

            # write the fit parameters
            for field in clearancedata:
                row = [field] + [clearancedata[field], clearancemodel[field]]
                csv_writer.writerow(row)

# ------------------------------------------------------------------------------------------------------------
#                                                Configuration
# ------------------------------------------------------------------------------------------------------------

# --..--..--..--.. Input / Output ..--..--..--..--
# relative directory where the original patient files reside
inputdirectory = "./reformatted-data/3-culled-data/"

# relative directory where you want the reformatted patient files to go
#   Note: You should create this directory if it does not already exist
outputdirectory = "./reformatted-data/4-clearance-initial/"



# Execution starts here
if __name__ == "__main__":
    # ---------------------------------------------------------------------
    if not os.path.exists(inputdirectory):
        print(f"The relative (to this script) input directory {inputdirectory} does not exist")
        sys.exit()

    # clear any existing files and remake the target output directory
    if os.path.exists(outputdirectory):
        shutil.rmtree(outputdirectory)
    os.mkdir(outputdirectory)

    #---------------------------------------------------------------------


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
                writeClearance(infile, outfile)

    if not os.path.exists(outputdirectory):
        print(f"The relative (to this script) output directory {outputdirectory} does not exist (please create it first)")
        sys.exit()
