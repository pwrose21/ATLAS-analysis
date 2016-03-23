

## /////// ##
## imports ##
## /////// ##

import sys, os, subprocess
import re, time, copy, math, array
import ROOT


## /////// ##
## classes ##
## /////// ##
class plotVar:
    def __init__(self, branch_name, title, nbinsx, xlow, xup):
        self.branch_name = branch_name
        self.title = title
        self.nbinsx = nbinsx
        self.xlow = xlow
        self.xup = xup


## /////// ##
## globals ##
## /////// ##

## common to all files ----------------------------------
## ------------------------------------------------------

## trees to loop over -----------------------------------
tree_nominal = ['nominal']
tree_sys     = ['JET_19NP_JET_Flavor_Response__1down', 'JET_19NP_JET_Flavor_Composition__1up', 'JET_19NP_JET_EtaIntercalibration_Modelling__1up', 'JET_19NP_JET_Flavor_Composition__1down', 'JET_19NP_JET_EtaIntercalibration_Modelling__1down', 'JET_19NP_JET_EffectiveNP_3__1down', 'JET_19NP_JET_EffectiveNP_6restTerm__1up', 'JET_19NP_JET_EffectiveNP_6restTerm__1down', 'JET_19NP_JET_EffectiveNP_4__1up', 'JET_19NP_JET_BJES_Response__1up', 'JET_19NP_JET_Pileup_RhoTopology__1up', 'JET_19NP_JET_EtaIntercalibration_TotalStat__1up', 'JET_19NP_JET_BJES_Response__1down', 'JET_19NP_JET_Flavor_Response__1up', 'JET_19NP_JET_EffectiveNP_4__1down', 'JET_19NP_JET_Pileup_OffsetNPV__1down', 'JET_19NP_JET_Pileup_OffsetMu__1up', 'JET_19NP_JET_Pileup_PtTerm__1down', 'JET_19NP_JET_PunchThrough_MC15__1down', 'JET_19NP_JET_Pileup_RhoTopology__1down', 'JET_19NP_JET_PunchThrough_MC15__1up', 'JET_19NP_JET_Pileup_OffsetNPV__1up', 'JET_19NP_JET_EffectiveNP_1__1down', 'JET_19NP_JET_EffectiveNP_5__1down', 'JET_19NP_JET_Pileup_OffsetMu__1down', 'JET_19NP_JET_SingleParticle_HighPt__1down', 'JET_19NP_JET_EffectiveNP_2__1up']

tree_all     = tree_nominal + tree_sys

## weight systematics -----------------------------------
weight_sys   = ['weight_pileup_UP']

## variables to plot ------------------------------------
plot_vars = [plotVar('HT_all', 'HT_all', 23, 0, 4600*1000.)]

## branches needed --------------------------------------
branches = set(['lep_n', 'el_n', 'mu_n',
            'jet_n', 'jet_isB_77_n', 'ljet_isHOT_100_n',
            'Mbb_MindR', 'HT_all', 'weight_mc',
            'weight_leptonSF', 'weight_bTagSF_77', 'weight_pileup'])
            

## normalize to this luminosity -------------------------
lumi = 0.001

## file containing list of cross sections by dsid -------
xsec_file = "/export/home/prose/2t2bAnalysis/HQTFourTopTools-2.3.37/HQTFourTopTools/scripts/histMaker/PMGXSections.txt"

## file containing sum of event weights by dsid ---------
sumweights_file = "/export/home/prose/2t2bAnalysis/HQTFourTopTools-2.3.37/HQTFourTopTools/scripts/histMaker/sampleYields.txt"


## file specific ----------------------------------------
## ------------------------------------------------------

## file names -------------------------------------------
input_file  = ""
output_file = "outputfile.root"
output_directory = "output-" + time.strftime("%d_%m_%Y")

## the tfile
tFile = ""

## dsid -------------------------------------------------
dsid = 0

## data / mc --------------------------------------------
isMC = True

## file xsec --------------------------------------------
xsec = 1.

## file sum weights -------------------------------------
sumweights = 1.

## sample name ------------------------------------------
sample_name = ""

## list of histograms -----------------------------------
hist_dict = {}


def main(argv):

    global tFile
    global hist_dict

    tFile = ROOT.TFile(argv[1], "READ")
    f2 = ROOT.TFile(output_file, 'RECREATE')
    tFile.cd()
    if not tFile:
        sys.exit("ERROR : Input file could not be read. Exiting.")

    GetFileMetaData(tFile)
    print "\nINFO : Summary of MetaData --"
    print "       DSID       :", dsid
    print "       isMC       :", isMC
    print "       XSec       :", xsec
    print "       SumWeights :", sumweights

    ## analyze all trees
    for iTree in tree_all:
        ## run the analysis on each tree
        AnalyzeTree(iTree)

        #print "Number of hists:", len(hist_dict)
        ## save to file and clear hist_dict
        f2.cd()
        for iH in hist_dict:
            hist_dict[iH].Write()
            hist_dict[iH].Delete()
        #f2.Close()
        tFile.cd()
        hist_dict = {}
    
    tFile.Close()
    f2.Close()

def AnalyzeTree(tree_name):
    global hist_dict

    print "Analyzing tree:", tree_name
    t = tFile.Get(tree_name)
    if not t:
        print "\nWARNING : The current tree is not available. Skipping"
    t.SetBranchStatus('*', 0)
    [t.SetBranchStatus(branchname, 1) for branchname in branches]
    if tree_name == 'nominal':
        [t.SetBranchStatus(branchname, 1) for branchname in weight_sys]
 
    nPreSel = 0
    for iEvt in range(t.GetEntries()):
        t.GetEntry(iEvt)
        
        if iEvt % 1000 == 0:
            print "INFO :", iEvt, "events have been analyzed"

        ## ========================================================== ##
        ## ================== your analysis here ==================== ##
        ## ========================================================== ##

        ## preselection -----------------------------------------------
        if not DoPreSelection(t):
            continue
        nPreSel = nPreSel + 1
        #print "Passed pre selection"
        ## ------------------------------------------------------------

        ## categorize events -------------------------------------------
        ## -------------------------------------------------------------
        event_cat_list = GetEventCategories(t)

        ## --------------------------------------------------------------

        ## dictionary of event weights by name
        ## --------------------------------------------------------------
        if tree_name == 'nominal':
            weight_dict = FillEventWeightDict(t)
        else:
            weight_dict = {}
            weight_dict[tree_name] = ( xsec * lumi * t.weight_mc *
                                       t.weight_leptonSF *
                                       t.weight_bTagSF_77 / sumweights)

        ## create and/or fill histograms
        ## -------------------------------------------------------------
        for iVar in plot_vars:
            FillHists(t, iVar, event_cat_list, weight_dict, hist_dict)
        ## ========================================================== ##
        ## ========================== end =========================== ##
        ## ========================================================== ##
    print "nPreSel", nPreSel

def FillHists(tree, plot_var, event_cat_list, weight_dict, hist_dict):
    val = getattr(tree, plot_var.branch_name)
    if type(val) in [float, int]:
        val = [val]
    for iCat in event_cat_list:
        for iSys in weight_dict:
            h_name = ('h_' + sample_name + '_'+ 
                      plot_var.title + '_' + iCat + '_' + iSys)
            h = hist_dict.get(h_name, '')
            if not h:
                h = ROOT.TH1F(h_name, h_name, plot_var.nbinsx, plot_var.xlow, plot_var.xup)
                h.GetXaxis().SetTitle(plot_var.title + ' ' + GetUnits(plot_var.title))
                h.GetYaxis().SetTitle("Entries / " + "{0:.1f}".format((plot_var.xup - plot_var.xlow) / plot_var.nbinsx ) 
                                      + " " + GetUnits(plot_var.title))
                h.SetDirectory(0)
                hist_dict[h_name] = h
            # now that it is made, fill
            for iVal in val:
                #print "Filling hist with value:", iVal, "weight:", weight_dict[iSys]
                #print "Weight sys / value:", iSys, weight_dict[iSys]
                h.Fill(iVal, weight_dict[iSys])

def GetUnits(str):
    return "[]"

def FillEventWeightDict(tree):
    weight_dict = {}
    wPU = tree.weight_pileup
    wLSF = tree.weight_leptonSF
    wBTag = tree.weight_bTagSF_77
    wBase = xsec * lumi * tree.weight_mc / sumweights
    weight_dict['nominal'] = wBase * wPU * wLSF * wBTag
    for iW in weight_sys:
        #tree.SetBranchStatus(iW, 1)
        if 'weight_pileup' in iW:
            weight_dict[iW] = wBase * getattr(tree, iW) * wLSF * wBTag
        if 'weight_leptonSF' in iW:
            weight_dict[iW] = wBase * wPU * getattr(tree, iW) * wBTag
        if 'weight_bTagSF' in iW:
            if 'eigenvars' in iW:
                eigenvars = getattr(tree, iW)
                for i,iVar in enumerate(eigenvars):
                    weight_dict[iW+'_NP'+str(i)] = wBase * wPU * wLSF * iVar
            else:
                weight_dict[iW] = wBase * wPU * wLSF * getattr(tree, iW)
    return weight_dict
        

def GetEventCategories(tree):
    event_cat_list = [] ## given event might be in more than one category
        ##                     i.e. inclusive/exclusive ntag
    ## lepton category
    lep_cat = '' + '1el'*(tree.el_n==1) + '1mu'*(tree.mu_n==1)
    ## jet category
    jet_cat = '' + '5j'*(tree.jet_n==5) + 'ge6j'*(tree.jet_n>=6) 
    ## b category
    b_cat = '' + '2b'*(tree.jet_isB_77_n==2) + '3b'*(tree.jet_isB_77_n==3) + 'ge4b'*(tree.jet_isB_77_n>=4)
    ## HOT category
    ljet_cat = '' + '0lj'*(tree.ljet_isHOT_100_n==0) + '1lj'*(tree.ljet_isHOT_100_n==1) + 'ge2lj'*(tree.ljet_isHOT_100_n>=2)
    ## mbb cat
    mbb_cat = '' + 'LM'*(tree.Mbb_MindR<100000.) + 'HM'*(tree.Mbb_MindR>=100000.)
    ## event cat
    event_cat = lep_cat + jet_cat + b_cat + ljet_cat
    if all(a in event_cat for a in ['ge6j', '1lj']):
        if any(a in event_cat for a in ['3b', 'ge4b']):
            event_cat = event_cat + mbb_cat
    event_cat_list.append(event_cat)

    return event_cat_list

def DoPreSelection(tree):
    #print "Number of jets", tree.jet_n
    #print "Number of b jets", tree.jet_isB_77_n
    if (tree.lep_n != 1) or (tree.jet_n < 5) or (tree.jet_isB_77_n < 2):
        return False
    return True


def GetFileMetaData(root_file):
    global isMC
    GetDSIDFromSumWeightsTree(root_file)
    isMC = bool(dsid)
    if isMC:
        GetXSecFromFile(dsid, xsec_file)
        GetSumWeightsFromFile(dsid, sumweights_file)
    GetSampleNameFromDSID(dsid)

def GetSampleNameFromDSID(dsid):
    global sample_name
    # data
    if dsid == 0:
        print "This sample is data. Did you expect this for file", input_file, " ?\n"
        sample_name = 'data'

    # backgrounds
    elif dsid == 410000:
        sample_name = 'ttbar'
    elif dsid in range(407009, 407012):
        sample_name = 'ttbar'
    elif dsid in range(410001, 410004+1):
        sample_name = 'ttbarSys'
    elif dsid in range(361300, 361371+1):
        sample_name = 'wjets'
    elif dsid in range(361372, 361467+1):
        sample_name = 'zjets'
    elif dsid in range(361520, 361534+1):
        sample_name = 'wjetsSys'
    elif dsid in range(361500, 361519+1):
        sample_name = 'zjetsSys'
    elif dsid in [410013, 410014, 410011, 410012, 410025, 410026]:
        sample_name = 'singletop'
    elif dsid in [410017, 410018, 410019, 410020, 410099, 410100, 410101, 410102]:
        sample_name = 'singletopSys'
    elif dsid in range(361081, 361087+1):
        sample_name = 'diboson'
    elif dsid in [410066, 410067, 410068, 410069, 410070, 410073, 410074, 410075, 410081]:
        sample_name = 'ttV'
    elif dsid in [341177, 341270, 341271]:
        sample_name = 'ttH'

    # signals
    
    elif dsid == 410080:
        sample_name = 'fourtopSM'
    elif dsid == 302777:
        sample_name = 'fourtopCI'
    elif dsid in range(302055, 302059+1):
        sample_name = '2uedmKK' + str(1000 + 200*(dsid-302055))
    elif dsid in range(302470, 302480+1):
        sample_name = 'TTS' + str(700 + 50*(dsid-302470))
    elif dsid == 302469:
        sample_name = 'TTS600'
    elif dsid == 302481:
        sample_name = 'TTS1300'
    elif dsid == 302482:
        sample_name = 'TTS1400'
    elif dsid == 302483:
        sample_name = 'TTD700'
    elif dsid == 302484:
        sample_name = 'TTD950'
    elif dsid == 302485:
        sample_name = 'TTD1200'
    else:
        print "Could not determine the sample for input file", input_file, "\n"
        sample_name = 'unknown'

def GetSumWeightsFromFile(dsid, sumweights_file):
    global sumweights
    f = open(sumweights_file, 'r')
    for line in f:
        if line.startswith('#'):
            continue
        if int(line.split()[0]) == dsid:
            sumweights = float(line.split()[1])
            break
    f.close()

def GetXSecFromFile(dsid, xsec_file):
    global xsec
    f = open(xsec_file, 'r')
    for line in f:
        if line.startswith('#'):
            continue
        if int(line.split()[0]) == dsid:
            xsec = float(line.split()[6])
            break
    f.close()

def GetDSIDFromSumWeightsTree(root_file):
    global dsid
    sumweights = 0
    t = root_file.Get("sumWeights")
    if not t:
        sys.exit("ERROR : sumWeights tree not present.\n        Could not determine DSID. Exiting.")
    t.SetBranchStatus('dsid', 1)
    for i in range(t.GetEntries()):
        t.GetEntry(i)
        dsid = t.dsid
        break


def GetXSecFromDSID(xsec_file, dsid):
    pass

if __name__ == "__main__" : main(sys.argv)
