# --------------------------------------------------------
#
#  ***Oxford Mathematical Brain Modeling Group***
#   This script is the third script in the clearance pipeline.
#   The purpose of this script is to capitalize on the following
#   observations
#   1. Most of the patient data measurements occur in the following
#       way:
#       1a. Several measurements within the first 10-15 hours
#       1b. One followup around 24 hours ( 20h <= t <= 35h )
#       1c. One followup around 48 hours ( 40h <= t <= 60h )
#       1d. A baseline measurement circa 30 days
#       1e. Possibly an additional measurement farther out
#
#   2. This script aims to keep the three measurements from
#       1b, 1c and 1d and drop the rest.  If a patient does
#       not have all three, they are dropped.
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


class patient:
    def __init__(self, pid):
        self.pid = pid
        self.isValid = False
        self.times = {'first': 0.0, 'second': 0.0, 'third': 0.0}
        self.btimes = {'first': False, 'second': False, 'third': False}
        self.itimes = {'first': 0, 'second': 0, 'third': 0}
        self.data = {}

    def __checkTimes(self, csvheader):
        # the header is assumed to have the format
        # [Struct Name] [time] [time] ... [time]
        # These times are in terms of *days* and we
        # need to determine the first and second bounds
        # in terms of *hours*
        ttimes = csvheader[1:]

        for st in ttimes:
            timed = float(st)   # the time in days
            timeh = timed*24    # the time in hours

            # the first time slot is the timepoint between
            #   20 and 35 hours
            if 20.0 <= timeh and timeh <= 35.0:
                if self.btimes['first'] == True:
                    print(f"Found duplicate first time at {timed} ({timeh} hours)")
                else:
                    print(f"First time for patient {self.pid} is {timed} days ({timeh} hours)")
                    self.times['first'] = timed
                    self.btimes['first'] = True

            if 35.0 < timeh and timeh <= 60.0:
                if self.btimes['second'] == True:
                    print(f"Found duplicate second time at {timed} ({timeh} hours)")
                else:
                    print(f"Second time for patient {self.pid} is {timed} days ({timeh} hours)")
                    self.times['second'] = timed
                    self.btimes['second'] = True

            # now we look for the ~1 month checkup
            if 20.0 < timed and timed <= 50.0:
                if self.btimes['third'] == True:
                    print(f"Found duplicate third time at {timed}")
                else:
                    print(f"Third time for patient {self.pid} is {timed} days")
                    self.times['third'] = timed
                    self.btimes['third'] = True

        self.isValid = (self.btimes['first'] and self.btimes['second'] and self.btimes['third'])

        if self.isValid:
            print(f"Patient {self.pid} time extraction successful")
            for tkey in self.times:
                # store the indices
                self.itimes[tkey] = csvheader.index(str(self.times[tkey]))

        return self.isValid


    def importdata(self,input):
        imported = True
        line = 0

        with open(input) as incsv:
            csv_reader = csv.reader(incsv, delimiter=',')

            for row in csv_reader:
                if line == 0:
                    bstat = self.__checkTimes(row)
                    if bstat == False:
                        imported = False
                        break
                else:
                    idxf = self.itimes['first']
                    idxs = self.itimes['second']
                    idxt = self.itimes['third']

                    self.data[row[0]] = {'first': row[idxf], 'second': row[idxs], 'third': row[idxt]}

                line +=1

        return imported


    def writepatient(self,outdir):
        if self.isValid == False:
            return

        flnm = str(self.pid) + ".csv"
        csvout = outdir + flnm


        with open(csvout, mode='w') as outcsv:
            csv_writer = csv.writer(outcsv, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)

            #write the header
            header = ['StructName'] + [self.times['first'], self.times['second'], self.times['third']]
            csv_writer.writerow(header)

            for fld in self.data:
                row = [fld] + [self.data[fld]['first'], self.data[fld]['second'], self.data[fld]['third']]
                csv_writer.writerow(row)



# --..--..--..--.. Input / Output ..--..--..--..--
# relative directory where the original patient files reside
inputdirectory = "./reformatted-data/2-dropped-fields/"

# relative directory where you want the reformatted patient files to go
#   Note: You should create this directory if it does not already exist
outputdirectory = "./reformatted-data/3-culled-data/"


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
    patientdata = []

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

                # we assume that the filenames are XXX.csv
                #   where XXX is the patient ID (e.g. 7.csv
                #   or 124.csv etc)
                ipid = subj.find('.csv')
                pid = int(subj[:ipid])

                p = patient(pid)
                bValid = p.importdata(infile)

                if bValid:
                    patientdata.append(p)

    for pd in patientdata:
        pd.writepatient(outputdirectory)






