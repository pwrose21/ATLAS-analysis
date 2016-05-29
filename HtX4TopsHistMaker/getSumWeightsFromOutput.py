import ROOT, sys, commands
from ROOT import *

class Sample:
    def __init__(self, dsid, tag, weight):
        self.dsid = int(dsid)
        self.tag = tag
        self.weight = float(weight)

    """
    def __add__(self, other):
        total_weight = self.weight + other.weight
        return Sample(self.dsid, self.tag, total_weight)
    """
    def add(self, other):
        if isinstance(other, self.__class__):
            self.weight = self.weight + other.weight
            return True
        return False

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return (self.dsid == other.dsid and self.tag == other.tag)
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)
    
    def __str__(self):
        return ("DSID : " + str(self.dsid) + "\nTag : " + 
                str(self.tag) + "\nWeight : " + str(self.weight) + "\n")


dsDir = sys.argv[1] + '/'
status, fileDirs = commands.getstatusoutput('ls ' + dsDir)
fileDirs = fileDirs.split('\n')

samples = []

for iDir in fileDirs:

    ## skip log files
    if not 'output.root' in iDir:
        continue

    ## skip data files
    if 'physics_Main' in iDir:
        continue
    
    ## debug
    #print "Starting dir", iDir

    ## get info about the files within this directory
    dsid = iDir.split('.')[2]
    tag = iDir.split('.')[-3]
    
    ## iterate through all output files for this sample
    status, files = commands.getstatusoutput("ls " + dsDir + iDir)
    files = files.split('\n')
    for iFile in files:
        if not iFile:
            continue
        if 'part' in iFile:
            continue
        newSample = True
        fileSumWeights = 0
        f = ROOT.TFile(dsDir + iDir + '/' + iFile)
        t = f.Get('sumWeights')
        if not t:
            continue
        for iEvt in t:
            if int(iEvt.dsid) != int(dsid):
                print "ERROR : Mismatch in dsid between sample name and tree. Investigate!"
                print "Tree dsid:", iEvt.dsid
                print "File name dsid:", dsid
                sys.exit()
            fileSumWeights = fileSumWeights + float(iEvt.totalEventsWeighted)

        print dsid, tag, fileSumWeights
            
        thisSample = Sample(dsid, tag, fileSumWeights)

        print thisSample

        for iSample in samples:
            if thisSample == iSample:
                iSample.add(thisSample)
                newSample = False
        if newSample:
            samples.append(thisSample)


f = open("sumWeights.txt", 'w')
f.write("# dsid ## tag ## weight #\n")
for i in samples:
    f.write(str(i.dsid) + " " + str(i.tag) + " " + str(i.weight) + "\n")
f.close()

"""

fileDir = sys.argv[1]
status, files = commands.getstatusoutput('ls ' + fileDir)
files = files.split('\n')
dsidWeights = [['410000', 0.]]

for iFile in files:
    if not '.root' in iFile:
        continue
    print 'Working on file:', iFile
    f = ROOT.TFile(fileDir+'/'+iFile, 'read')
    t = f.Get('sumWeights')
    print t
    if not t:
        print 'Could not get sumWeights tree -- is this a data file?'
        continue
    if t.GetEntries() > 1:
        print 'WARNING::FILE', iFile, 'HAS MORE THAN 1 ENTRY IN SUMWEIGHTS TREE. INVESTIGATE!!!!'
    for iEvt in t:
        dsid = str(iEvt.dsid)
        if dsid == '0'and not 'physics_Main' in iFile:
            print 'Found dsid==0 in file:', iFile
            sys.exit('failed')
        if dsid == '0'and 'physics_Main' in iFile:
            continue
        sumWeight = float(iEvt.totalEventsWeighted)
        addNew = True
        for iSample in dsidWeights:
            if iSample[0] == dsid:
                print 'DSID already exists in list of DSID weights -- adding to current value!'
                iSample[1] = iSample[1] + sumWeight
                addNew = False
        if addNew:
            print 'Adding new DSID to list of DSID weights'
            dsidWeights.append([dsid, sumWeight])


y = open('sampleYields.txt', 'w')
for iSample in dsidWeights:
    print str(iSample[0]), str(iSample[1]), '\n'
    y.write(str(iSample[0]) + ' ' + str(iSample[1]) + '\n')

"""

"""
#input file
f = ROOT.TFile(sys.argv[1], 'read')
#sum weights
t = f.Get('sumWeights')

nEntries = t.GetEntries()

print 'There are', str(nEntries), 'in the sumWeights tree'

dsid = 0
nEvt = 0
for iEvt in t:
    dsid = iEvt.dsid
    print dsid
    nEvt = iEvt.totalEventsWeighted
    print nEvt
f.Close()
"""
