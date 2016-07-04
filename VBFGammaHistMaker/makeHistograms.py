#!/usr/bin/env python

## /////// ##
## imports ##
## /////// ##

import sys, os, subprocess, getopt
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

# TODO -- put some of this in a config file!!!

## ------------------------------------------------------
## ------------------------------------------------------
## common to all files ----------------------------------
## ------------------------------------------------------
## ------------------------------------------------------

## Convention : globals start with "m_"

## trees to loop over -----------------------------------
m_tree_nominal = ['Nominal']
m_tree_sys     = ['EG_RESOLUTION_ALL__1down', 'EG_RESOLUTION_ALL__1up',
                  'EG_SCALE_ALL__1down', 'EG_SCALE_ALL__1up', 
                  'JET_19NP_JET_BJES_Response__1down', 'JET_19NP_JET_BJES_Response__1up',
                  'JET_19NP_JET_EffectiveNP_1__1down', 'JET_19NP_JET_EffectiveNP_1__1up', 
                  'JET_19NP_JET_EffectiveNP_2__1down', 'JET_19NP_JET_EffectiveNP_2__1up', 
                  'JET_19NP_JET_EffectiveNP_3__1down', 'JET_19NP_JET_EffectiveNP_3__1up',
                  'JET_19NP_JET_EffectiveNP_4__1down', 'JET_19NP_JET_EffectiveNP_4__1up', 
                  'JET_19NP_JET_EffectiveNP_5__1down', 'JET_19NP_JET_EffectiveNP_5__1up', 
                  'JET_19NP_JET_EffectiveNP_6restTerm__1down', 'JET_19NP_JET_EffectiveNP_6restTerm__1up',
                  'JET_19NP_JET_EtaIntercalibration_Modelling__1down', 'JET_19NP_JET_EtaIntercalibration_Modelling__1up',
                  'JET_19NP_JET_EtaIntercalibration_TotalStat__1down', 'JET_19NP_JET_EtaIntercalibration_TotalStat__1up',
                  'JET_19NP_JET_Flavor_Composition__1down', 'JET_19NP_JET_Flavor_Composition__1up',
                  'JET_19NP_JET_Flavor_Response__1down', 'JET_19NP_JET_Flavor_Response__1up', 'JET_19NP_JET_GroupedNP_1__1down',
                  'JET_19NP_JET_GroupedNP_1__1up', 'JET_19NP_JET_Pileup_OffsetMu__1down', 'JET_19NP_JET_Pileup_OffsetMu__1up',
                  'JET_19NP_JET_Pileup_OffsetNPV__1down', 'JET_19NP_JET_Pileup_OffsetNPV__1up', 'JET_19NP_JET_Pileup_PtTerm__1down',
                  'JET_19NP_JET_Pileup_PtTerm__1up', 'JET_19NP_JET_Pileup_RhoTopology__1down', 'JET_19NP_JET_Pileup_RhoTopology__1up',
                  'JET_19NP_JET_PunchThrough_MC15__1down', 'JET_19NP_JET_PunchThrough_MC15__1up', 
                  'JET_19NP_JET_SingleParticle_HighPt__1down', 'JET_19NP_JET_SingleParticle_HighPt__1up', 'JET_JER_SINGLE_NP__1up',
                  'MUONS_ID__1down', 'MUONS_ID__1up', 'MUONS_MS__1down', 'MUONS_MS__1up', 'MUONS_SCALE__1down', 'MUONS_SCALE__1up',
                  ]
## weight systematics -----------------------------------
#m_weight_sys   = ['weight_pileup_UP']

#testing
m_tree_sys = ['JET_JER_SINGLE_NP__1up','EG_RESOLUTION_ALL__1down','EG_RESOLUTION_ALL__1up']
m_tree_sys = []
m_weight_sys = []

## variables to plot ------------------------------------
#m_plot_vars = [plotVar('mBB_Regression', 'mBBReg', 22, 40*1000., 260*1000.)]
m_plot_vars = [plotVar('mBB_PtRecollbbOneMuPartonBukinNew', 'mBBReg', 35, 50., 400.),
               plotVar('pTB1', 'pTB1', 36, 40., 400.),
               plotVar('pTB2', 'pTB2', 36, 40., 400.)
               ]

#plotVar('mBB_PtRecollbbOneMuPartonBukinNew', 'mBBReg', 22, 40*1000., 260*1000.)]
## branches needed --------------------------------------
m_branches = set(['mBB_Regression', 'MV2c20B1', 'MV2c20B2', 'nJ', 'pTBB', 'MCWeight', 
                  'mJJ', 'pTJ1', 'pTB1', 'dEtaJJ', 
                  'dRJJ','dRBB', 'dRB1J1', 'dRB1J2', 'dRB2J1', 'dRB2J2',
                  'dRJ1Ph','dRJ2Ph', 'dRB1Ph', 'dRB2Ph'])
#m_branchesMC = set([ 'weight_mc', 'weight_leptonSF', 'weight_bTagSF_77', 'weight_pileup'])
m_branchesMC = set([])
   
## normalize to this luminosity -------------------------
m_lumi = 1000. # convert x-sec from pb to fb
#m_lumi = m_lumi * 100. # 100 fb-1

## file containing list of cross sections by dsid -------
m_xsec_file_str = "/export/home/prose/ATLAS-analysis/VBFGammaHistMaker/crossx.txt"

## file containing sum of event weights by dsid ---------
m_sumweights_file_str = "/export/home/prose/ATLAS-analysis/VBFGammaHistMaker/sumweights.txt"

##
m_runDir = "/export/home/prose/ATLAS-analysis/VBFGammaHistMaker/"

## 
#m_bcut = -0.4434 #77%
m_bcut = -0.0436
## ------------------------------------------------------
## file specific ----------------------------------------
## ------------------------------------------------------
m_input_file_str  = "" # file name / location string
m_file_tag_str = '' # tag information extracted from file name
m_dsid = -2 # dataset ID
m_isMC = True # data/MC
m_xsec = 1. # cross section UNITS ?
m_sumweights = 1. # sample sum of event weights
m_sample_name = "" # sample name
m_sample_cat = "" # sample category -- signal, background, data
m_sample_sys = m_tree_nominal[0] # sample systematic

def main(argv):

    # get input_file, sumweights_file, xsec_file, lumi
    ParseCommandLineArguments(argv)

    # open the TFile, check that it can be read
    in_tFile = ROOT.TFile(m_input_file_str, "READ")
    if not in_tFile:
        sys.exit("ERROR : Input file could not be read. Exiting.")

    # get metadata frm the TFile
    GetFileMetaData(in_tFile)
    
    # make output file
    n_out = 0
    output_file_str = 'outputHists_' + str(n_out) + '_' + m_input_file_str.split('/')[-1]
    while os.path.exists(output_file_str):
        output_file_str = output_file_str.replace('_' + str(n_out) + '_', '_' +str(n_out+1) + '_')
        n_out = n_out + 1
        
    print "Writing histograms to output file:", output_file_str
    out_tFile = ROOT.TFile(output_file_str, 'RECREATE')
    out_tFile.cd()
    d = out_tFile.mkdir('Systematics')
    #cdtof = top->mkdir("tof")
    #d = ROOT.TDirectory('Systematics', 'Systematics')
    #d.SetMother(out_tFile)
    # cd back to input file
    in_tFile.cd()

    # define which trees to analyze
    # -----------------------------
    # systematic trees for MC
    if m_isMC:
        tree_all = m_tree_nominal + m_tree_sys
    # only nominal trees for Data
    else:
        tree_all = m_tree_nominal

    # analyze trees
    # -------------
    for iTree in tree_all:
        # do not analyze systematic trees for systematic samples
        if iTree in m_tree_sys and m_sample_sys != m_tree_nominal[0]:
            continue

        # get the TTree
        this_tTree = in_tFile.Get(iTree)
        if not this_tTree:
            continue

        ## run the analysis on each tree
        hist_dict = AnalyzeTree(this_tTree)

        ## save to file and clear hist_dict
        out_tFile.cd()
        for iH in hist_dict:
            if '_Sys' in hist_dict[iH].GetName():
                d.cd()
            #SetBinErrors(hist_dict[iH],0.1)
            hist_dict[iH].Write()
            hist_dict[iH].Delete()
            out_tFile.cd()
        in_tFile.cd()

    # close the input and output files
    in_tFile.Close()
    out_tFile.Close()


def SetBinErrors(h, frac):
    for i in range(1, h.GetNbinsX()+1):
        h.SetBinError(i, h.GetBinContent(i)*frac)


def AnalyzeTree(t):
    hist_dict = {}
    tree_name = t.GetName()

    print "Analyzing tree:", tree_name

    # activate branches
    # -----------------
    "Basically need all"
    """
    t.SetBranchStatus('*', 0)
    branches = copy.deepcopy(m_branches)
    if m_isMC:
        branches.update(m_branchesMC)
    [t.SetBranchStatus(branchname, 1) for branchname in branches]
    if tree_name == m_tree_nominal[0] and m_isMC:
        [t.SetBranchStatus(branchname, 1) for branchname in m_weight_sys]
    """
    n_selected = 0
    # loop over all events in the file
    # --------------------------------
    for iEvt in range(t.GetEntries()):

        # get the ith entry
        t.GetEntry(iEvt)

        # information on progress
        if iEvt % 1000 == 0:
            print "INFO :", iEvt, "events have been analyzed"

        ## ========================================================== ##
        ## ================== your analysis here ==================== ##
        ## ========================================================== ##

        ## preselection ------------------------------------------------
        ## -------------------------------------------------------------
        if not DoPreSelection(t):
            continue
        n_selected = n_selected + 1
        ## categorize events -------------------------------------------
        ## -------------------------------------------------------------
        event_cat_list = GetEventCategories(t)
        if not event_cat_list:
            continue

        ## dictionary of event weights by name
        ## -------------------------------------------------------------
        weight_dict = {}
        if tree_name == m_tree_nominal[0]:
            if m_sample_sys == m_tree_nominal[0]:
                if m_isMC:
                    # for nominal tree, nominal sample_sys, MC, fill all event weights
                    weight_dict = FillEventWeightDict(t)
                else:
                    # data should have nominal tree, nominal sample_sys -- all weights = 1
                    weight_dict[m_tree_nominal[0]] = 1
            else:
                # for sys samples, use nominal weights, but label hists by the sample_sys
                weight_dict[sample_sys] = ( m_xsec * m_lumi * t.weight_mc *
                                          t.weight_leptonSF *
                                          t.weight_bTagSF_77 / m_sumweights)
    
        else:
            # for systematic tree, use nominal weights, but label hists by the tree_sys
            weight_dict[tree_name] =  m_xsec * m_lumi * t.MCWeight / m_sumweights
            #weight_dict[tree_name] = ( m_xsec * m_lumi * t.weight_mc *
            #                           t.weight_leptonSF *
            #                           t.weight_bTagSF_77 / m_sumweights)

        ## create and/or fill histograms
        ## -------------------------------------------------------------
        for iVar in m_plot_vars:
            FillHists(t, iVar, event_cat_list, weight_dict, hist_dict)
        
        ## ========================================================== ##
        ## ========================== end =========================== ##
        ## ========================================================== ##
    print "Number of events selected:", n_selected
    return hist_dict

def DoPreSelection(t):
    # cut based from arXiv
    #if not t.passTrig:
    #    return False

    if not t.mJJ > 800000.:
        return False

    if not t.MV2c20B1 > m_bcut:
        return False
    if not t.MV2c20B2 > m_bcut:
        return False
    
    if not t.nJ >= 4:
        return False

    #if not t.pTPh > 26000.:
    #    return False

    #if not t.BDT > 0:
    #    return False

    if not t.pTBB > 80000.:
        return False

    #if not t.BDT > 0.1:
    #    return False

    #if not t.pTJ2 > 40000.:
    #    return False
    #if not t.pTB2 > 40000.:
    #    return False

    """
    if not t.dEtaJJ > 4.:
        return False

    if any(getattr(t,dR) < 0.7 for dR in ['dRJJ','dRBB', 'dRB1J1', 'dRB1J2', 'dRB2J1', 'dRB2J2']):
        return False
    if any(getattr(t,dR) < 1.4 for dR in ['dRJ1Ph','dRJ2Ph', 'dRB1Ph', 'dRB2Ph']):
        return False

    """

    return True

def GetEventCategories(tree):
    event_cat_list = [] ## given event might be in more than one category
                        ## i.e. inclusive/exclusive ntag

    n_tag = 0
    n_tag = n_tag + int(tree.MV2c20B1>m_bcut) + int(tree.MV2c20B2>m_bcut)
    if not n_tag == 2:
        return []

    n_jet = tree.nJ
    if n_jet >=4:
        n_jet = '4p'

    """
    pt_bb = tree.pTBB
    if pt_bb > 80000.:
        pt_bb = '100'
    else:
        pt_bb = '0'
    """
    pt_bb = tree.BDT
    if pt_bb > 0.1:
        pt_bb = '1'
    else:
        pt_bb = '0'
    
    reg = ''
    
    reg = "mBBcr"
    if n_tag == 2 and tree.mBB_PtRecollbbOneMuPartonBukinNew > 110000. and tree.mBB_PtRecollbbOneMuPartonBukinNew < 140000.:
        #round to binning#tree.mBB_Regression > 112500. and tree.mBB_Regression < 137500.:
        reg = "SR"
    
    #reg = 'SR'

    event_cat = ( str(n_tag) + 'tag' + str(n_jet) + 'jet_' +
                  str(pt_bb) + 'BDT_' + reg)

    event_cat_list.append(event_cat)

    event_cat_list.append(event_cat.replace('BDT', 'ptv'))


    ## Z CR
    reg = "mBBcrZ"
    if n_tag ==2 and tree.mBB_PtRecollbbOneMuPartonBukinNew > 80000. and tree.mBB_PtRecollbbOneMuPartonBukinNew < 100000.:
        # round to binning#tree.mBB_Regression > 81000. and tree.mBB_Regression < 99000.
        reg = "SRZ"

    event_cat = ( str(n_tag) + 'tag' + str(n_jet) + 'jet_' +
                  str(pt_bb) + 'BDT_' + reg)

    event_cat_list.append(event_cat)
    event_cat_list.append(event_cat.replace('BDT', 'ptv'))
    #event_cat_list.append(event_cat.replace(reg, 'allmBB'))

    """
    if n_jet == 4:
        event_cat_list.append(event_cat)
    
        
    event_cat_list.append(event_cat.replace(str(n_jet) + 'jet_', '4pjet_'))
    """

    return event_cat_list


def FillHists(tree, plot_var, event_cat_list, weight_dict, hist_dict):
    val = getattr(tree, plot_var.branch_name)
    if type(val) in [float, int]:
        val = [val]
    for iCat in event_cat_list:
        for iSys in weight_dict:
            h_name = m_sample_name + '_' + iCat + '_' + plot_var.title + (iSys != m_tree_nominal[0]) * ('_Sys' + iSys)
            #if (iSys != m_tree_nominal[0]):
            #    h_name = 'Systematics/' + h_name
            h = hist_dict.get(h_name, '')
            if not h:
                h = ROOT.TH1F(h_name, h_name, plot_var.nbinsx, plot_var.xlow, plot_var.xup)
                h.Sumw2(ROOT.kFALSE)
                h.GetXaxis().SetTitle(plot_var.title + ' ' + GetUnits(plot_var.title))
                h.GetYaxis().SetTitle("Entries / " + "{0:.1f}".format((plot_var.xup - plot_var.xlow) / plot_var.nbinsx ) 
                                      + " " + GetUnits(plot_var.title))
                h.SetDirectory(0)
                hist_dict[h_name] = h
            h = hist_dict[h_name]
            # now that it is made, fill
            for iVal in val:
                #print "Filling hist with value:", iVal, "weight:", weight_dict[iSys], 'for sample', m_sample_name
                #print "Weight sys / value:", iSys, weight_dict[iSys]
                h.Fill(iVal/1000., weight_dict[iSys])

def GetUnits(str):
    return "[]"

def FillEventWeightDict(tree):
    # testing
    return {m_tree_nominal[0] : m_xsec * m_lumi * tree.MCWeight / m_sumweights}

    weight_dict = {}
    wPU = tree.weight_pileup
    wLSF = tree.weight_leptonSF
    wBTag = tree.weight_bTagSF_77
    wBase = m_xsec * m_lumi * tree.weight_mc / m_sumweights
    weight_dict[m_tree_nominal[0]] = wBase * wPU * wLSF * wBTag
    for iW in m_weight_sys:
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

def ParseCommandLineArguments(argv):
    global m_input_file_str
    global m_xsec_file_str
    global m_sumweights_file_str
    global m_lumi

    ## parse options 
    _short_options = 'hf:x:s:l:'
    _long_options = ['help', 'inputFile=', 'xSecFile=', 'sumWeightsFile=','lumi=']
    try:
        opts, args = getopt.gnu_getopt(argv, _short_options, _long_options)
    except getopt.GetoptError:
        print 'getopt.GetoptError\n'
        print __doc__
        sys.exit(2)
    for opt, val in opts:
        if opt in ('-h', '--help'):
            print __doc__
            sys.exit()
        if opt in ('-f', '--inputFile'):
            m_input_file_str = val
        if opt in ('-x', '--xSecFile'):
            m_xsec_file_str = val
        if opt in ('-s', '--sumWeightsFile'):
            m_sumweights_file_str = val
        if opt in ('-l', '--lumi'):
            m_lumi = float(val)

    # check that we have everything
    if not os.path.isfile(m_input_file_str):
        sys.exit("Error : please enter a valid input file with '-f <input_file>' option")
    if not os.path.isfile(m_xsec_file_str):
        sys.exit("Error : please enter a valid xsec file with '-x <sec_file>' option")
    if not os.path.isfile(m_sumweights_file_str):
        sys.exit("Error : please enter a valid sumweights file with '-s <sumweights_file>' option")

    print "\nINFO : Summary of configuration for this run --"
    print "       Input file      :", m_input_file_str
    print "       X-sec file      :", m_xsec_file_str
    print "       Sumweights file :", m_sumweights_file_str
    print "       Luminosity      :", m_lumi,"\n"

def GetFileMetaData(root_file): # root_file is an opened TFile
    # dsid ----------------------------------------------------------
    global m_dsid
    m_dsid = GetDSIDFromNominalTree(root_file)
    if m_dsid < -1:
        sys.exit("Error : Could not get dsid for this sample")

    # isMC ----------------------------------------------------------
    global m_isMC
    m_isMC = not(m_dsid==-1)
    if not m_isMC:
        print "Data sample detected.  Did you expect this?"

    # tag -----------------------------------------------------------
    """
    global m_file_tag_str
    m_file_tag_str = GetTagFromFileName(m_input_file_str)
    if not m_file_tag_str:
        sys.exit("Error : could not determine file tag from file name: " + m_input_file_str)
    """
    # MC only -------------------------------------------------------
    if m_isMC:
        # x-sec -----------------------------------------------------
        global m_xsec
        m_xsec = GetXSecFromFile(m_dsid, m_xsec_file_str)
        if m_xsec < 0:
            sys.exit("Error : Could not get xsec for this sample")

        # sum of event weights --------------------------------------
        global m_sumweights
        m_sumweights = GetSumWeightsFromFile(m_dsid, m_file_tag_str, m_sumweights_file_str)
        if m_sumweights < 0:
            sys.exit("Error : Could not get sumweights for this sample")

    # sample name, file systematic
    global m_sample_name
    global m_sample_sys
    m_sample_name, m_sample_sys = GetSampleNameAndSysFromDSIDAndTag(m_dsid, m_file_tag_str)
    if not m_sample_name:
        sys.exit("Error : Could not determine the sample name for dsid " + str(dsid)
                 + " and tag " + m_file_tag_str)


    print "\nINFO : Summary of MetaData --"
    print "       DSID        :", m_dsid
    print "       isMC        :", m_isMC
    print "       XSec        :", m_xsec
    print "       SumWeights  :", m_sumweights
    print "       Sample name :", m_sample_name
    print "       Sample sys  :", m_sample_sys, '\n'

def GetDSIDFromNominalTree(root_file):
    dsid = -2
    t = root_file.Get("Nominal")
    if not t:
        sys.exit("ERROR : Nominal tree not present.\n        Could not determine DSID. Exiting.")
    t.SetBranchStatus('MCChannelNumber', 1)
    for i in range(t.GetEntries()):
        t.GetEntry(i)
        dsid = t.MCChannelNumber
        break
    return dsid

def GetTagFromFileName(file_name_str):
    #user.prose.410066.MadGraphPythia8EvtGen.DAOD_TOPQ1.e4111_s2608_s2183_r7326_r6282_p2516.HtX4Tops_00-00-01_rt
    try:
        file_tag_str = file_name_str.split('TOPQ1.')[1].split('.HtX4Tops')[0]
        return file_tag_str
    except:
        return ''

def GetXSecFromFile(dsid, xsec_file):
    xsec = -1.
    foundXsec = False
    f = open(xsec_file, 'r')
    for line in f:
        if line.startswith('#'):
            continue
        if int(line.split()[0]) == dsid:
            xsec = float(line.split()[1])
            foundXsec = True
            break
    f.close()
    if not foundXsec:
        f = open(m_runDir + '/missingXsec.txt', 'a')
        f.write(str(dsid) + '\n')
        f.close()
    return xsec

def GetSumWeightsFromFile(dsid, file_tag, sumweights_file):
    sumweights = -1
    foundWeights = False
    f = open(sumweights_file, 'r')
    for line in f:
        if line.startswith('#'):
            continue
        if int(line.split()[0]) == dsid:
            sumweights = float(line.split()[1])
            foundWeights = True
            break
    f.close()
    if not foundWeights:
        f = open(m_runDir + '/missingWeight.txt', 'a')
        f.write(str(dsid) + ' ' + file_tag + '\n')
        f.close()
    return sumweights

def GetSampleNameAndSysFromDSIDAndTag(dsid, file_tag_str):
    sample_name = ''
    sample_sys  = m_tree_nominal[0]

    # data
    if dsid == -1:
        print "This sample is data. Did you expect this for file", m_input_file_str, " ?\n"
        sample_name = 'data'

    # backgrounds
    elif dsid == 343388:
        sample_name = 'HbbjjaSM125'
    elif dsid == 343390:
        sample_name = 'ZbbjjaEWK'
    elif dsid == 343391:
        sample_name = 'ZbbjjaQCD'
    elif dsid == 343392:
        sample_name = 'NonResbbjja'
    else:
        print "Could not determine the sample for input file", m_input_file_str, "\n"
        
    return sample_name, sample_sys


if __name__ == "__main__" : main(sys.argv)
