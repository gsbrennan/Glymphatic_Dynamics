# --------------------------------------------------------
#
#  ***Oxford Mathematical Brain Modeling Group***
#
#   This script extracts a data field from the structured
#   clearance data provided by Simula Simula Research
#   Laboratory and Oslo Research Hospital Partners
#
#   This script writes a set of files for each raw-data
#   patient.  Each file contains columns whose numeric
#   value indicates days post-injection and whose rows
#   are the regions and associated (extracted) data field
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
from datetime import datetime


#-------------------------------------------------------
# Class for manipulating patient data
#-------------------------------------------------------



class patient:

    # construct a patient using an (integer) ID
    # and a path to a csv file for this patient.
    # [optional] datafield: the data field you want to extract for this patient into the
    #                       reformatted file.  (Default = 'Mean')
    # [optional] timepointstr: this field indicates the search
    #               string corresponding to a data point.  Extraction
    #               occurs based on finding fields that contain this
    #               substring.
    def __init__(self,id, csvfile, datafield='Median', timepointstr='[date, time]'):
        self.patid = int(id)
        self.csvin = csvfile

        self.records = {}
        self.timepointid = timepointstr

        self.extractval = datafield

        self.normalized = False

        self.valndxs = []
        self.timevals = []

        self.regionndx = -1

        self.data = {}

        # read in the patient data from the indicated CSV file
        self.__readdata()

    def __populateTimesFromHeader(self,header):
        def converttime(h):
            # [[ extract the date / time substring and return a datetime object]]
            # -- find the date / time substring index
            dtndx = h.find(self.timepointid)
            # -- extract everything before this point
            dtstr = h[:dtndx]
            # -- trim any whitespace from beginning and end of the string
            dtstr = dtstr.strip()

            # -- create the datetime object from the record
            dt = datetime.strptime(dtstr, "%Y%m%d %H:%M:%S")

            return dt

        timen = 0
        initDate = None

        for h in header:
            if self.timepointid in h:
                if (timen == 0):
                    initDate = converttime(h)

                    # save the numeric index (position) for time zero
                    self.timevals.append(float(0.00))
                else:
                    timep = converttime(h)

                    # get time in days
                    dt = timep - initDate

                    # convert and save into times dictionary
                    dayval = dt.days + (dt.seconds) / (60 * 60 * 24)

                    self.timevals.append(float(dayval))

                timen += 1

    # this function extracts the row indices for the datafield specified
    #   by self.extractval.  It is assumed that each patient file format
    #   has the same set of value indices for each time (for instance,
    #   each time has an associated data field called 'Mean', 'Median', 'StdDev', etc)
    def __extractValueIndices(self, header, regionName='StructName'):

        # we expect to extract one field (with value matching self.extractval) for
        # each time in the array self.timevals)
        nfields = len(self.timevals)

        searchndx = 0

        for h in header:
            if(h == regionName):
                self.regionndx = header.index(regionName)
            else:
                if (h == self.extractval):
                    #thisndx = header.index(self.extractval, searchndx)
                    self.valndxs.append(searchndx)

            searchndx += 1

    # (Private) parse a data record according to the stored indices
    def __parsedata(self, record):
        # This line is assumed to have the following format
        #  ... [Field name] .... [data value] ... [data value] ... [data value]
        # where
        # [Field name] is the anatomical region identifier (e.g. Left-Hippocampus, etc)
        # ... represent intermediate data values (which we are not interested in)
        # [data value] is the data corresponding to the index entries for the heading
        #   of interest whose identifier is stored in self.extractval.  Examples include
        #   'Mean', 'StdDev', 'Min', etc

        if self.regionndx == -1:
            print(f"[Error] Incorrect region index for patient {self.patid}.  Something went wrong")
        else:
            # create the data dictionary entry for this region
            region = record[self.regionndx]

            # extract the data field and the values
            self.data[region] = []

            for i in range(len(self.valndxs)):
                idxat = self.valndxs[i]
                self.data[region].append(record[idxat])

    # (Private) read in data from the raw csv file whose path is
    #   contained in the internal field self.csvin
    def __readdata(self):
        with open(self.csvin) as csvfile:
            csv_reader = csv.reader(csvfile, delimiter=',')
            line=0

            for row in csv_reader:
                if line==0:
                    # If this is the header (first line) we populate the times
                    # dictionary with the numeric time (in days since the first
                    # measurement) as the keys and the indices of the column as
                    # the values
                    self.__populateTimesFromHeader(row)
                elif line == 1:
                    # This is the second line of the header.  This is where we
                    # extract the field name index
                    self.__extractValueIndices(row)
                else:
                    # These are all data records.  We read the values for this row
                    # and populate the patient structure according to the desired ID.
                    self.__parsedata(row)
                line +=1

    # writes the extracted patient data to a CSV file.
    # dirout: the output directory (must end with a '/' such as '/home/user/output/')
    def writedata(self, dirout):

        csvout = ""

        if self.normalized:
            #csvout = dirout + str(self.patid) + "-" + self.extractval + "-normalized.csv"
            csvout = dirout + str(self.patid) +".csv"
        else:
            csvout = dirout + str(self.patid) + "-" + self.extractval + ".csv"

        linecount = 0
        with open(csvout, mode='w') as outcsv:
            data_writer = csv.writer(outcsv, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)

            # write the header first
            # The header should be: 'StructName' time, time, time, time
            header = ['StructName'] + self.timevals
            data_writer.writerow(header)

            # write the data
            for region in self.data:
                # the data row should be: Region, data, data, data, .. , data
                towrite = [region] + self.data[region]
                data_writer.writerow(towrite)

    def getID(self):
        return self.patid

    # Normalize the data for this patient according to the values in `normalizedValues'
    # Note: a normalization object for this patient can be retrieved from a normalizer
    #   object using this patient's ID
    # [Optional] valstr is the prefix of the expected keys for the normalizedValues object
    #       they should be 'TP1', 'TP2', etc where 'TP' is the default valstr
    def normalizePatientData(self, normalizedValues, valstr='TP'):
        expectedValues = len(self.timevals)
        nvalues = len(normalizedValues)
        allfound = True

        if expectedValues != nvalues:
            print(f"[ERROR] Expected {expectedValues} for normalization of patient {self.patid} but recieved {nvalues}")
        else:
            for i in range(expectedValues):
                keystr = valstr + str(i+1)
                if keystr in normalizedValues:
                    nval =  float(normalizedValues[keystr])

                    # Normalize the data in every region by the normalization value
                    for region in self.data:
                        rdat = self.data[region][i]
                        self.data[region][i] = str(float(rdat) / 1)#nval)
                else:
                    print(f"Missing expected key {keystr} in normalized values object")
                    allfound = False
            if allfound:
                self.normalized = True
            else:
                print(f"Failed to normalize all values for patient {self.patid}")

#-------------------------------------------------------
# Class for normalization
#-------------------------------------------------------
class normalization:

    def __init__(self, csvfile):

        # Each patient
        self.patients = {}
        self.indices = {}
        self.__readpatientdata(csvfile)



    def __readpatientdata(self, csvfile, patid = 'PatID', ignore=['COMMENT']):

        with open(csvfile) as csvfile:
            csv_reader = csv.reader(csvfile, delimiter=',')

            line = 0
            for row in csv_reader:

                if line == 0:
                    # -- process the header
                    for h in row:
                        if h not in ignore:
                            if patid in h:
                                self.indices[patid] = row.index(h)
                            elif 'TP' in h:
                                # we want the ref index corresponding to this 'TP' entry
                                strh = h.strip()
                                self.indices[strh] = row.index(h) + 1

                else:
                    # -- get the patient ID and create the patient
                    pidx = int(self.indices[patid])
                    pid = int(row[pidx])
                    self.patients[pid] = {}

                    for sidx in self.indices:
                        # if sidx is a label containing 'TP' we
                        # want to associate this label to the 'ref'
                        # value that comes after it
                        if 'TP' in sidx:
                            tmval = row[self.indices[sidx]]
                            if tmval != '':
                                self.patients[pid][sidx] = float(tmval)
                line += 1

    # Retrieve a list of the patient IDs in the normalization object
    def getPatientIDList(self):
        return list(self.patients.keys())

    # Get a specific patient's normalization value
    def getPatientNormalization(self, id):
        p = None
        if id in self.patients:
            p = self.patients[id]
        else:
            print(f"a record of patient {id} does not exist in this normalization object")

        return p
#-------------------------------------------------------


# ------------------------------------------------------------------------------------------------------------
#                                                Configuration
# ------------------------------------------------------------------------------------------------------------

# --..--..--..--.. Input / Output ..--..--..--..--
# relative directory where the original patient files reside
patientinputdirectory = "./raw-data/clearance_data_excel/"

# relative directory where you want the reformatted patient files to go
#   Note: You should create this directory if it does not already exist
outputroot = "./reformatted-data/"
outputdirectory = outputroot + "1-extracted-mean/"

# --..--..--..--.. Normalization Options ..--..--..--..--
normalize = False
normalizationcsv = "./raw-data/ref-ROI-values.csv"

#--------------------------------------------------------------------------------------------------------------

# Execution starts here
if __name__ == "__main__":

    if not os.path.exists(patientinputdirectory):
        print(f"The relative (to this script) input directory {patientinputdirectory} does not exist")
        sys.exit()

    # clear any existing files and remake the target output directory
    # you will need to re-run the full pipeline
    if os.path.exists(outputroot):
        shutil.rmtree(outputroot)
    os.mkdir(outputroot)
    os.mkdir(outputdirectory)

    dirlevel = 0
    patientList = []

    # ----------- Read in and parse all subject files --------------------
    for rootdir, subjectdirs, files in os.walk(patientinputdirectory):

        totalsubj = len(files)
        thissubj = 0
        dirlevel += 1

        # only process the top level subdirectories
        if dirlevel == 1:
            for subj in files:
                thissubj += 1
                print(f"Processing file {subj}")
                sidx = subj.find("-")
                sid = subj[:sidx]
                id = int(sid)

                filepath = patientinputdirectory + subj

                p = patient(id, filepath)
                patientList.append(p)

    if not os.path.exists(outputdirectory):
        print(f"The relative (to this script) output directory {outputdirectory} does not exist")
        sys.exit()


    # --------------- Normalize patient data -------------- #

    normalizer = normalization(normalizationcsv)
    normalizedIDs = normalizer.getPatientIDList()





    # --------------- output all patients -------------------------------------
    for p in patientList:
        pid = p.getID()
        if pid in normalizedIDs:
            normvals = normalizer.getPatientNormalization(pid)
            p.normalizePatientData(normvals)
        else:
            print(f"No normalization data is available for patient {pid}")
        p.writedata(outputdirectory)
