# --------------------------------------------------------
#
#  ***Oxford Mathematical Brain Modeling Group***
#
#  This (optional) script can be used following
#  4-compute-clearance.py.  This script replaces
#  the "linear model" results with averaged results
#  from nearby regions where the fit was exponential.
#
#  This script implements two views of "nearby"
#    1. Nearby via location in space: based on (x,y,z)
#       coordinates of the ROIs
#
#    2. Nearby via structural graph connectivity: uses
#       the edge weights of the graph laplacian with
#       diffusive weighting
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
import math
import shutil
import xml.etree.ElementTree as ET

import numpy as np
from pylab import *
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

class node:

    def __init__(self, nid):
        self.x = 0.0
        self.y = 0.0
        self.z = 0.0

        self.clearance = 0.0

        self.nodeid = int(nid)
        self.edges = []

        # node related structural information
        self.region = ""                #node tag d4 in connectome file
        self.freesurfername = ""        #node tag d5 in connectome file
        self.hemisphere = ""            #node tag d7 in connectome file
        self.__nodename = ""

        self.bValidClearance = False

    def setRegion(self, regiond4):
        self.region = regiond4

    def setFreesurferName(self, fsn):
        self.freesurfername = fsn

    def setHemisphere(self, hs):
        self.hemisphere = hs

    def makeNodeString(self):
        # there is a slight naming difference between the way the
        # connectome names the brain stem naming and the way
        # the clearance files name the brain stem.
        if self.freesurfername == "Brain-Stem":
            self.__nodename = "subcortical.brainstem.right"
        else:
            self.__nodename = self.region + "." + self.freesurfername + "." + self.hemisphere


    def setClearance(self, clr):
        self.clearance = clr

    def setClearanceValid(self, bValidity):
        self.bValidClearance = bValidity

    def associateEdge(self, e):
        self.edges.append(e)

    def setxcoord(self, x):
        self.x = float(x)

    def setycoord(self, y):
        self.y = float(y)

    def setzcoord(self, z):
        self.z = float(z)

    # --------------------------
    def getNodeString(self):
        return self.__nodename

    def getID(self):
        return self.nodeid

    def getxcoord(self):
        return self.x

    def getycoord(self):
        return self.y

    def getzcoord(self):
        return self.z

    def getClearance(self):
        return self.clearance

    def getIsClearanceValid(self):
        return self.bValidClearance
    # ---------------------------

    # compare this node to another one
    def compareNID(self, n):
        return (self.nid == n.getID())

class edge:

    def __init__(self, srcid, trgid):
        self.n = 0
        self.l = 0

        self.sourceid = srcid
        self.targetid = trgid

    def setfibercount(self, n):
        self.n = float(n)

    def setfiberlength(self, l):
        self.l = float(l)

    def getSrcId(self):
        return self.source

    def getTrgId(self):
        return self.target

    # return the edge weight.
    # Default: n/l^2
    def getWeight(self, lpow=2):
        return float(self.n / math.pow(self.l, lpow))


class connectome:

    def __init__(self):
        self.__reset()

    def __reset(self):
        self.edgeset = []

        self.edgesByNodeID = {}
        self.nodeneighbors = {}
        self.nodesbyID = {}
        self.edgemap = {}
        self.nodeStridtoNid = {}


        self.nodeRadialProximity = -1.0
        self.nodesByProximity = {}


    def __getNSKeyedStr(self, st):
        ns = 'http://graphml.graphdrawing.org/xmlns'
        ky = '{' + ns + '}' + st
        return ky

    def __addGraphmlNode(self, n):
        # get the node ID attribute
        iid = int(n.get('id'))

        # the node we will add to the collection
        newn = node(iid)

        # extract the data keys
        datky = self.__getNSKeyedStr('data')

        for dat in n.findall(datky):
            # get the code in the file
            code = dat.get('key')
            code = code.strip()

            if code == "d0":
                # x coordinate
                newn.setxcoord(dat.text)
            if code == "d1":
                # y coordinate
                newn.setycoord(dat.text)
            if code == "d2":
                # z coordinate
                newn.setzcoord(dat.text)
            if code == "d4":
                # cortical / subcortical
                newn.setRegion(dat.text)
            if code == "d5":
                # freesurfer parcellation id
                newn.setFreesurferName(dat.text)
            if code == "d7":
                # hemisphere (left / right)
                newn.setHemisphere(dat.text)

        # make the coded node string for this node
        # ( this string can be matched to the freesurfer
        #   naming convention in the clearance data )
        newn.makeNodeString()
        nodestrid = newn.getNodeString()

        self.nodeStridtoNid[nodestrid] = newn.getID()

        # add to the internal list of nodes
        if iid not in self.nodesbyID:
            self.nodesbyID[iid] = newn
        else:
            print(f"[Parsing Error] Node with id {iid} already exists in the connectome.")


    def __addGraphmlEdge(self, e):
        srcid = int(e.get('source'))
        trgid = int(e.get('target'))

        newe = edge(srcid, trgid)

        # extract the data keys
        datky = self.__getNSKeyedStr('data')

        # ---
        for dat in e.findall(datky):
            # get the code in the file
            code = dat.get('key')
            code = code.strip()

            if code == "d9":
                # number of fibers (n)
                newe.setfibercount(dat.text)

            if code == "d12":
                # average fiber length
                newe.setfiberlength(dat.text)
        # ---

        if srcid not in self.nodeneighbors:
            self.nodeneighbors[srcid] = [trgid]
        else:
            self.nodeneighbors[srcid].append(trgid)

        if trgid not in self.nodeneighbors:
            self.nodeneighbors[trgid] = [srcid]
        else:
            self.nodeneighbors[trgid].append(srcid)

        # We assume that all of the nodes have already
        # been added to the connectome.  If not, parse
        # the nodes first before parsing the edges.
        srn = self.nodesbyID[srcid]
        tgn = self.nodesbyID[trgid]

        srn.associateEdge(newe)
        tgn.associateEdge(newe)

        # list building
        self.edgeset.append(newe)

        if srcid not in self.edgesByNodeID:
            self.edgesByNodeID[srcid] = [newe]
        else:
            self.edgesByNodeID[srcid].append(newe)

        if trgid not in self.edgesByNodeID:
            self.edgesByNodeID[trgid] = [newe]
        else:
            self.edgesByNodeID[trgid].append(newe)

        # add an entry to the tuple-based edge map
        self.edgemap[(srcid, trgid)] = newe

    def __getNearestProximalNeighbor(self, nid):
        rProx = 1.0e6
        nNearest = None

        if nid in self.nodesbyID:
            n = self.nodesbyID[nid]

            # get the key node's x,y,z coordinates
            nx = n.getxcoord()
            ny = n.getycoord()
            nz = n.getzcoord()

            for nextid in self.nodesbyID:
                if nextid != nid:
                    nodenxt = self.nodesbyID[nextid]
                    nxx = nodenxt.getxcoord()
                    nxy = nodenxt.getycoord()
                    nxz = nodenxt.getzcoord()

                    rn = math.pow(nx-nxx, 2) + math.pow(ny-nxy, 2) + math.pow(nz-nxz, 2)
                    rn = math.sqrt(rn)

                    if rn < rProx:
                        rProx = rn
                        nNearest = nodenxt
        else:
            print(f"[ERROR] node id {nid} not found in node list")

        return rProx, nNearest

    # --------------------------------------------
    # Call this function to replace the node with
    # string id `strid' with an average based on
    # its neighbors determined by spatial proximity
    #
    # The average is taken over neighboring nodes with
    # a non-zero clearance value
    # ---------------------------------------------
    def __getAverageClearanceByProximity(self, strid):
        success = False
        clearavg = 0.0

        if strid in self.nodeStridtoNid:
            nid = self.nodeStridtoNid[strid]
            ngbrlist = self.nodesByProximity[nid]

            clearsum = 0.0
            clearcount = 0

            for ngbrid in ngbrlist:
                ngbrnd = self.nodesbyID[ngbrid]

                bValidNode = ngbrnd.getIsClearanceValid()
                clearval = ngbrnd.getClearance()

                # We only want to average over nodes are / were originally
                # nodes with valid clearance values.  We don't want to
                # include any nodes which started off invalid in the
                # averaging process.
                if bValidNode and clearval > 0.0:
                    clearsum += clearval
                    clearcount += 1

            if clearcount > 0:
                clearavg = float(clearsum / clearcount)
                success = True
            else:
                print(f"[WARNING] Node {strid} has no valid proximal neighbors with nonzero clearance.")
        else:
            print(f"[ERROR] The string ID {strid} does not match any known nodes")

        return success, clearavg


    def __getAverageClearanceByConnectivity(self, strid, bWeighted):
        success = False
        clearavg = 0.0

        if strid in self.nodeStridtoNid:
            nid = self.nodeStridtoNid[strid]

            # use the neighbor list built when reading in
            # the connectome (i.e. connectivity)
            ngbrlist = self.nodeneighbors[nid]

            clearsum = 0.0
            weightsum = 0
            edgeErrors = 0

            for ngbrid in ngbrlist:
                ngbrnd = self.nodesbyID[ngbrid]

                bValidNode = ngbrnd.getIsClearanceValid()
                clearval = ngbrnd.getClearance()

                # We only want to average over nodes are / were originally
                # nodes with valid clearance values.  We don't want to
                # include any nodes which started off invalid in the
                # averaging process.
                if bValidNode and clearval > 0.0:
                    weight = 1.0

                    if bWeighted == True:
                        # get the edge corresponding to the (nid, ngbrid) pair.
                        # this edge could be (nid, ngbrid) or (ngbrid, nid)
                        edg = None
                        if (nid, ngbrid) in self.edgemap:
                            edg = self.edgemap[(nid, ngbrid)]
                        elif (ngbrid, nid) in self.edgemap:
                            edg = self.edgemap[(ngbrid, nid)]
                        else:
                            print(f"[ERROR] Cannot find edge ({nid}, {ngbrid}) or ({ngbrid}, {nid}) in the edge map")
                            # set the clearance to 0.0 so that it is not averaged in
                            clearval = 0.0
                            edgeErrors += 1
                        if edg != None:
                            # get the diffusive weight
                            weight = edg.getWeight(lpow=2)

                    clearsum += weight*clearval
                    weightsum += weight

            if weightsum > 0 and edgeErrors == 0:
                success = True
                clearavg = float(clearsum / weightsum)

        return success, clearavg

    # ---------------------------------------
    # Find all connectome neighbors for a
    # anatomically labeled region
    # The string id should be that of the clearance
    # files.  E.g. "cortical.medialtemporal.right"
    # ---------------------------------------
    def getGraphNeighbors(self, strid):
        bRet = None
        # get the numerical ID of the node
        if strid in self.nodeStridtoNid:
            nid = self.nodeStridtoNid[strid]
            bRet = self.nodeneighbors[nid]

        return bRet

    # ------------------------------------
    # Gets the average radial distance between
    # all nodes in the connectome
    # ------------------------------------
    def getAverageNodeRadialProximity(self):
        rAccum = 0.0

        for nV in self.nodesbyID:
            rProx, nnode = self.__getNearestProximalNeighbor(nV)
            if rProx > 0.0:
                rAccum += rProx
            else:
                print(f"[ERROR] Node ID with no valid nearest neighbor encountered")

        nnodes = len(self.nodesbyID)
        avgProx = float(rAccum / nnodes)
        return avgProx

    # parse a scale-33 connectome.  Note: it is
    # important that this connectome is scale-33
    # and not scale-500, etc.
    def parseConnectome(self, xmlfile):
        self.__reset()

        tree = ET.parse(xmlfile)
        root = tree.getroot()

        # get the graph main object
        gkey = self.__getNSKeyedStr('graph')
        graph = root.find(gkey)

        # the graph consists of node and
        #   edge children.  We process the
        # nodes first
        nkey = self.__getNSKeyedStr('node')
        allnodes = graph.findall(nkey)

        for node in allnodes:
            self.__addGraphmlNode(node)

        # now get the edges from the xml tree
        # edge parsing depends on nodes existing
        # so it needs to go second
        nkey = self.__getNSKeyedStr('edge')
        alledges = graph.findall(nkey)

        for edge in alledges:
            self.__addGraphmlEdge(edge)

    # -------------------------------------
    # build a proximity neighbor list by specifying
    # a distance (in mm) between nodes to
    # consider them as "proximal".  This list can
    # be used to average values by proximity
    #
    # This function returns True if every node is
    # proximal to at least one other node and false
    # otherwise.
    #
    # The second value that this function
    # returns is the minimal size of the proximal
    # neighbor list over all nodes.
    #
    # The third return value is the maximal
    # size of the proximal neighbor list
    #
    # The fourth value is the average size
    # of the proximal neighbor list
    # --------------------------------------
    def groupNodesByProximity(self, r=10.0):
        bAllNonempty = True
        iMinSize = len(self.nodesbyID)
        iMaxSize = -1

        self.nodeRadialProximity = r
        self.nodesByProximity = {}

        for na in self.nodesbyID:
            nda = self.nodesbyID[na]

            self.nodesByProximity[na] = []

            nax = nda.getxcoord()
            nay = nda.getycoord()
            naz = nda.getzcoord()

            for nb in self.nodesbyID:
                ndb = self.nodesbyID[nb]

                if na != nb:
                    nbx = ndb.getxcoord()
                    nby = ndb.getycoord()
                    nbz = ndb.getzcoord()

                    rval = math.sqrt( math.pow(nax-nbx, 2) + math.pow(nay - nby, 2) + math.pow(naz - nbz, 2) )

                    if rval <= self.nodeRadialProximity:
                        self.nodesByProximity[na].append(nb)

        proxCount = len(self.nodesByProximity)
        proxSum = 0
        for chk in self.nodesByProximity:
            proxSz = len(self.nodesByProximity[chk])
            if proxSz == 0:
                bAllNonempty = False
                iMinSize = 0
            if proxSz < iMinSize:
                iMinSize = proxSz
            if proxSz > iMaxSize:
                iMaxSize = proxSz

            proxSum += proxSz

        return bAllNonempty, iMinSize, iMaxSize, float(proxSum/proxCount)


    # ---------------------------------------
    # Set the clearance value of a node
    # based on anatomical string id.
    # ---------------------------------------
    def setNodeClearance(self, strid, clearanceVal, bClearanceValid = True):
        bSuccess = False

        if strid in self.nodeStridtoNid:
            nid = self.nodeStridtoNid[strid]
            n = self.nodesbyID[nid]
            n.setClearance(clearanceVal)
            n.setClearanceValid(bClearanceValid)
            bSuccess = True

        return bSuccess



    #-------------------------------------------
    # Call this function to reset all clearance
    # values to zero.
    #-------------------------------------------
    def resetNodalClearanceValues(self):
        for n in self.nodesbyID:
            nd = self.nodesbyID[n]
            nd.setClearance(0.00)


    # ------------------------------------------
    #  Call this function to average any invalid
    #  clearance values by the values of their
    #  proximity neighbors
    # ------------------------------------------
    def averageInvalidClearanceByProximity(self):

        for n in self.nodesbyID:
            nd = self.nodesbyID[n]
            if nd.getIsClearanceValid() == False:
                bSuccess, avgClearance = self.__getAverageClearanceByProximity(nd.getNodeString())

                if bSuccess:
                    nd.setClearance(avgClearance)
                else:
                    print(f"[ERROR] Could not set average proximity clearance for node {nd.getNodeString()}")


    # -------------------------------------------------
    # Call this function to average any invalid
    # clearance values by the values of their
    # axonally connected (graph) neighbors.
    #
    # Optional:
    # bWeighted = True (Default): computes an average
    #   weighted by diffusive edge weightings
    # bWeighted = False: computes an unweighted average
    #---------------------------------------------------
    def averageInvalidClearanceByConnectivity(self, bWeighted=True):

        for n in self.nodesbyID:
            nd = self.nodesbyID[n]
            if nd.getIsClearanceValid() == False:
                bSuccess, avgClearance = self.__getAverageClearanceByConnectivity(nd.getNodeString(), bWeighted)

                if bSuccess:
                    nd.setClearance(avgClearance)
                else:
                    print(f"[ERROR] Could not set average proximity clearance for node {nd.getNodeString()}")


    #--------------------------------------------
    # Get a count of the number of invalid clearance
    # nodes currently in the connectome (call after
    # loading a patient CSV into the connectome object)
    #--------------------------------------------
    def getInvalidClearanceCount(self):
        invalid = 0
        for n in self.nodesbyID:
            nd = self.nodesbyID[n]
            if nd.getIsClearanceValid() == False:
                invalid += 1

        return invalid


    #------------------------------------------------
    # Write the current content of the clearance values
    # to a desired output CSV file
    #------------------------------------------------
    def writeClearanceToCSV(self, outputcsv):

        with open(outputcsv, mode='w') as outcsv:
            csv_writer = csv.writer(outcsv, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)

            # write the header
            row = ['StructName'] + ['Clearance', 'Model Type']
            csv_writer.writerow(row)

            for n in self.nodesbyID:
                node = self.nodesbyID[n]
                structname = node.getNodeString()
                clearance  = node.getClearance()
                model = "Exponential"

                bValid = node.getIsClearanceValid()
                if bValid == False:
                    model = "Averaged"

                row = [structname, str(clearance), model]
                csv_writer.writerow(row)



# ------------------------------------------------------------------------------------------------------------
#                                                Configuration
# ------------------------------------------------------------------------------------------------------------

# --..--..--..--.. Input / Output ..--..--..--..--
# relative directory where the original patient files reside
inputdirectory = "./reformatted-data/4-clearance-initial/"

# relative directory where you want the reformatted patient files to go
#   Note: You should create this directory if it does not already exist
outputdirectoryroot = "./reformatted-data/4-clearance-initial/averaged/"

# Path to the scale-33 connectome graph
scale33Connectome = "./master-std33.graphml"
# -------------------------------------------------------------------------------------------------------------



# Load clearance values from a patient CSV into a connectome object
def loadClearanceCSV(connectome, inputcsv):

    # reset any currently stored clearance values to zero
    connectome.resetNodalClearanceValues()

    line = 0
    header = None

    with open(inputcsv) as incsv:
        csv_reader = csv.reader(incsv, delimiter=',')

        for row in csv_reader:
            if line == 0:
                # save the header time data
                header = row
                if header != ['StructName', 'Clearance', 'Model Type']:
                    print(f"[ERROR] Unexpected header format in Patient file {incsv}")
                    sys.exit()
            else:
                # bin the clearance data
                # The format for this file is expected to be the output format of
                # the script 4-compute-clearance.py.  This should be
                # ['StructName', 'Clearance', 'Model Type']
                strid = row[0].strip()
                clearance = float(row[1])
                ModelType = row[2].strip()

                if ModelType == "Exponential":
                    if connectome.setNodeClearance(strid, clearance, bClearanceValid=True) == False:
                        print(f"[ERROR] Could not set the clearance status for {strid}")
                else:
                    #if we could not fit an exponential model, this clearance value is invalid and
                    #needs to be averaged out somehow
                    if connectome.setNodeClearance(strid, clearance, bClearanceValid=False) == False:
                        print(f"[ERROR] Could not set the clearance status for {strid}")
            line += 1



# Execution starts here
if __name__ == "__main__":

    if not os.path.exists(inputdirectory):
        print(f"The relative (to this script) input directory {inputdirectory} does not exist")
        sys.exit()

    outdirProximity = outputdirectoryroot + "proximity-averaged"
    outdirConnectivity = outputdirectoryroot + "connectivity-averaged"

    # clear any existing files and remake the target output directory
    if os.path.exists(outputdirectoryroot):
        shutil.rmtree(outputdirectoryroot)
    os.mkdir(outputdirectoryroot)

    os.mkdir(outdirProximity)
    os.mkdir(outdirConnectivity)

    objConnectome = connectome()
    objConnectome.parseConnectome(scale33Connectome)

    # The graph neighbors are built automatically by the connectome.
    # We now establish the radial proximity neighbor list
    avgProx = objConnectome.getAverageNodeRadialProximity()

    # Now we build a proximal neighbor list by defining a radius
    # around each node equal to some percentage of the average
    # nearest neighbor radius.  We consider two nodes to be proximal
    # if their distance is less than or equal to twice the average
    # of the distance between each node and its nearest spatial neighbor
    groupval = 2.2
    allGrouped, minProximal, maxProximal, avgProximal = objConnectome.groupNodesByProximity(r=groupval*avgProx)

    print("")
    print(f"*** Proximity set to {groupval} times average nearest neighbor distance")
    print(f"*** Minimum proximity connectome grouping contains {minProximal} neighbors")
    print(f"*** Maximal proximity connectome grouping contains {maxProximal} neighbors")
    print(f"*** The average proximity connectome grouping contains {avgProximal} neighbors")
    print("")

    if allGrouped == False:
        print(f"[ERROR] The nodal proximity threshold is too low to facilitate clearance averaging by region")
        sys.exit()



    dirlevel = 0

    # ----------- create normalized files --------------------
    for rootdir, subjectdirs, files in os.walk(inputdirectory):

        totalsubj = len(files)
        thissubj = 0
        dirlevel += 1

        # only process the top level subdirectories
        if dirlevel == 1:
            for subj in files:
                thissubj += 1
                infile = inputdirectory + subj
                #outfile = outputdirectory + subj
                print(f"Processing file {subj}")

                loadClearanceCSV(objConnectome, infile)

                iInvalid = objConnectome.getInvalidClearanceCount()

                # if iInvalid == 0:
                #     print(f"All regions in patient file {subj} are valid (Exponential) model clearances")
                #     continue
                # else:
                #     print(f"Patient file {subj} contains {iInvalid} invalid clearance regions (i.e. linear model fitted).")
                #     print(f"Repairing subject file {subj} by averaging valid (Exponential) neighbors using two different methods")

                print(f"Patient file {subj} contains {iInvalid} invalid clearance regions (i.e. linear model fitted).")
                print(f"Repairing subject file {subj} by averaging valid (Exponential) neighbors using two different methods")

                # Now we average invalid values (these correspond to the linear model)
                # by proximity and output the result
                objConnectome.averageInvalidClearanceByProximity()

                # Write the averaged normalization
                proximityOutput = outdirProximity + f"/{subj}"
                objConnectome.writeClearanceToCSV(proximityOutput)

                # Now average invalid values (these correspond to the linear model)
                # by graph connectivity and output the result
                objConnectome.averageInvalidClearanceByConnectivity()

                connectivityOutput = outdirConnectivity + f"/{subj}"
                objConnectome.writeClearanceToCSV(connectivityOutput)
