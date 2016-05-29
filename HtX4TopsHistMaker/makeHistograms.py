#!/usr/bin/env python

## /////// ##
## imports ##
## /////// ##

import sys, os, subprocess, getopt
import re, time, copy, math, array
import ROOT

# local imports
from sampleCategories import *


## /////// ##
## classes ##
## /////// ##
class plotVar:
    def __init__(self, branch_name, title, nbinsx, xlow, xup, units = ''):
        self.branch_name = branch_name
        self.title = title
        self.nbinsx = nbinsx
        self.xlow = xlow
        self.xup = xup
        self.units = units


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
m_tree_nominal = ['nominal_Loose']
m_tree_sys     = []

## weight systematics -----------------------------------
m_weight_sys   = ['weight_pileup_UP', 'weight_pileup_DOWN',
                  'weight_indiv_SF_EL_Trigger_UP', 'weight_indiv_SF_EL_Trigger_DOWN',
                  'weight_indiv_SF_EL_Reco_UP', 'weight_indiv_SF_EL_Reco_DOWN',
                  'weight_indiv_SF_EL_ID_UP', 'weight_indiv_SF_EL_ID_DOWN',
                  'weight_indiv_SF_EL_Isol_UP', 'weight_indiv_SF_EL_Isol_DOWN',
                  'weight_indiv_SF_MU_Trigger_STAT_UP', 'weight_indiv_SF_MU_Trigger_STAT_DOWN',
                  'weight_indiv_SF_MU_Trigger_SYST_UP', 'weight_indiv_SF_MU_Trigger_SYST_DOWN',
                  'weight_indiv_SF_MU_ID_STAT_UP', 'weight_indiv_SF_MU_ID_STAT_DOWN',
                  'weight_indiv_SF_MU_ID_SYST_UP', 'weight_indiv_SF_MU_ID_SYST_DOWN',
                  'weight_indiv_SF_MU_Isol_STAT_UP', 'weight_indiv_SF_MU_Isol_STAT_DOWN',
                  'weight_indiv_SF_MU_Isol_SYST_UP', 'weight_indiv_SF_MU_Isol_SYST_DOWN',
                  'weight_indiv_SF_MU_TTVA_STAT_UP', 'weight_indiv_SF_MU_TTVA_STAT_DOWN',
                  'weight_indiv_SF_MU_TTVA_SYST_UP', 'weight_indiv_SF_MU_TTVA_SYST_DOWN',
                  'weight_jvt_UP', 'weight_jvt_DOWN',
                  ]

##testing
m_weight_sys = []

m_trf_sys = ['trf_weight_77_2ex_eigenvars_B_up', 'trf_weight_77_2ex_eigenvars_B_down',
             'trf_weight_77_2ex_eigenvars_C_up', 'trf_weight_77_2ex_eigenvars_C_down',
             'trf_weight_77_2ex_eigenvars_Light_up', 'trf_weight_77_2ex_eigenvars_Light_down',
             'trf_weight_77_2ex_extrapolation_up', 'trf_weight_77_2ex_extrapolation_down',
             'trf_weight_77_3ex_eigenvars_B_up', 'trf_weight_77_3ex_eigenvars_B_down',
             'trf_weight_77_3ex_eigenvars_C_up', 'trf_weight_77_3ex_eigenvars_C_down',
             'trf_weight_77_3ex_eigenvars_Light_up', 'trf_weight_77_3ex_eigenvars_Light_down',
             'trf_weight_77_3ex_extrapolation_up', 'trf_weight_77_3ex_extrapolation_down',
             'trf_weight_77_4in_eigenvars_B_up', 'trf_weight_77_4in_eigenvars_B_down',
             'trf_weight_77_4in_eigenvars_C_up', 'trf_weight_77_4in_eigenvars_C_down',
             'trf_weight_77_4in_eigenvars_Light_up', 'trf_weight_77_4in_eigenvars_Light_down',
             'trf_weight_77_4in_extrapolation_up', 'trf_weight_77_4in_extrapolation_down',
             ]

m_trf_sys = []

m_btag_sys = []

## variables to plot ------------------------------------
## branch name can be something already in the tree OR
##   the name of a global variable computed in CalculateEventVariables 
m_plot_vars = [plotVar('Meff', 'meff', 60, 0, 3000., 'GeV'),
               plotVar('m_btagjet_n', 'nBTags', 10, 0, 10),
               #plotVar('HT_all', 'HT', 60, 0, 3000., 'GeV'),
               #plotVar('m_jet_n', 'nJets', 10, 0, 10, '')
               ]

## branches needed --------------------------------------
#m_branches = set(['lep_n', 'el_n', 'mu_n',
#                  'jet_n', 'jet_isB_77_n', 'ljet_isHOT_100_n',
#                  'Mbb_MindR', 'Meff',])
#m_branchesMC = set([ 'weight_mc', 'weight_leptonSF', 'weight_bTagSF_77', 'weight_pileup'])
    
## normalize to this luminosity -------------------------
## xsec in pb
m_lumi = 1.

## file containing list of cross sections by dsid -------
m_xsec_file_str = "/export/home/prose/ATLAS-analysis/HtX4TopsHistMaker/XSection-MC15-13TeV.data"

## file containing sum of event weights by dsid ---------
m_sumweights_file_str = "/export/home/prose/ATLAS-analysis/HtX4TopsHistMaker/sumWeights.txt"

##
m_runDir = "/export/home/prose/ATLAS-analysis/HtX4TopsHistMaker/"

##
m_doTRF = False

## ------------------------------------------------------
## file specific ----------------------------------------
## ------------------------------------------------------
m_input_file_str  = "" # file name / location string
m_file_tag_str = '' # tag information extracted from file name
m_dsid = -1 # dataset ID
m_isMC = True # data/MC
m_isttbar = False # needed to save ttl, ttc, ttb in separate files
m_xsec = 0. # cross section UNITS ?
m_sumweights = 0. # sample sum of event weights
m_sample_name = "" # sample name
m_sample_sys = 'nominal' # sample systematic

## ------------------------------------------------------
## selection cuts ---------------------------------------
## ------------------------------------------------------
GeV             = 1000.
invGeV          = 0.001
m_WP            = '77'
m_cut_btag      = -0.4434 
m_cut_jet_pt    = 30. * GeV
m_cut_jet_eta   = 2.5
m_cut_lep_pt    = 30. * GeV
m_cut_lep_eta   = 2.5
m_cut_rcjet_pt  = 300. * GeV
m_cut_rcjet_eta = 2.0
m_cut_rcjet_m   = 100. * GeV
m_cut_rcjet_nsub = 2

## ------------------------------------------------------
## event variables --------------------------------------
## ------------------------------------------------------
m_lep_n       = -1
m_el_n        = -1
m_mu_n        = -1
m_jet_n       = -1
m_rcjet_n     = -1
m_hotjet_n    = -1
m_toptagjet_n = -1


m_btag_vars = ['m_mbb_mindr', 'm_btagjet_pt', 'm_btagjet_eta', 'm_btagjet_n']
m_mbb_mindr   = []
m_btagjet_pt  = []
m_btagjet_eta = []
m_btagjet_n   = ['','','']


def main(argv):
    # get input_file, sumweights_file, xsec_file, lumi
    ParseCommandLineArguments(argv)

    ## ---------------input / output file stuff------------------------------------
    # open the TFile, check that it can be read
    in_tFile = ROOT.TFile(m_input_file_str, "READ")
    if not in_tFile:
        sys.exit("ERROR : Input file could not be read. Exiting.")

    # get metadata frm the TFile
    GetFileMetaData(in_tFile)

    # ttbar handled special
    global m_isttbar
    m_isttbar = (m_sample_name == 'ttbar')

    # make output file(s)
    #output_file_str = m_sample_name + 'output_' + str(m_dsid) + '_' + m_file_tag_str + '_' + m_input_file_str.split('._')[1]
    output_file_str = m_sample_name + '.root'
    print "Writing histograms to output file:", output_file_str
    if not m_isttbar:
        out_tFile = ROOT.TFile(output_file_str, 'RECREATE')
    else:
        out_tFile_l = ROOT.TFile(output_file_str.replace('ttbar', 'ttl'), 'RECREATE')
        out_tFile_c = ROOT.TFile(output_file_str.replace('ttbar', 'ttc'), 'RECREATE')
        out_tFile_b = ROOT.TFile(output_file_str.replace('ttbar', 'ttb'), 'RECREATE')

    # cd back to input file
    in_tFile.cd()
    ## ----------------------------------------------------------------------------

    # define which trees to analyze
    # -----------------------------
    # systematic trees for MC
    if m_isMC:
        tree_all = m_tree_nominal + m_tree_sys
    # only nominal trees for Data
    else:
        tree_all = m_tree_nominal

    # define which b-tag systematics
    # ------------------------------
    global m_weight_sys
    if m_doTRF:
        m_weight_sys = m_weight_sys + m_trf_sys
    else:
        m_weight_sys = m_weight_sys + m_btag_sys

    # analyze trees
    # -------------
    for iTree in tree_all:
        # do not analyze systematic trees for systematic samples
        if iTree in m_tree_sys and m_sample_sys != 'nominal':
            continue

        # get the TTree
        this_tTree = in_tFile.Get(iTree)
        if not this_tTree:
            print 'WARNING ::', iTree, 'not available in this file!'
            continue

        ## run the analysis on each tree
        hist_dict = AnalyzeTree(this_tTree)

        ## save to file and clear hist_dict
        if not m_isttbar:
            out_tFile.cd()
            for iH in hist_dict:
                hist_dict[iH].Write()
                hist_dict[iH].Delete()
            in_tFile.cd()

        else:
            for iH in hist_dict:
                hist = hist_dict[iH]
                if hist.GetName().startswith('ttl'):
                    out_tFile_l.cd()
                    hist.SetName(hist.GetName().replace('ttl',''))
                    hist.SetTitle(hist.GetName())
                    hist.Write()
                    hist.Delete()
                elif hist.GetName().startswith('ttc'):
                    out_tFile_c.cd()
                    hist.SetName(hist.GetName().replace('ttc',''))
                    hist.SetTitle(hist.GetName())
                    hist.Write()
                    hist.Delete()
                elif hist.GetName().startswith('ttb'):
                    out_tFile_b.cd()
                    hist.SetName(hist.GetName().replace('ttb',''))
                    hist.SetTitle(hist.GetName())
                    hist.Write()
                    hist.Delete()
                else:
                    sys.exit("ERROR : ttbar hist not properly categorized")
            in_tFile.cd

    # close the input and output files
    in_tFile.Close()
    if not m_isttbar:
        out_tFile.Close()
    else:
        out_tFile_l.Close()
        out_tFile_c.Close()
        out_tFile_b.Close()

def AnalyzeTree(t):
    hist_dict = {}
    tree_name = t.GetName()

    print "Analyzing tree:", tree_name

    # activate branches
    # -----------------
    # we need mostly everything
    t.SetBranchStatus('*', 1)

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

        ## sample slice veto
        if not ApplySampleSliceVeto(t):
            continue

        ## preselection ------------------------------------------------
        ## -------------------------------------------------------------
        if not DoPreSelection(t):
            continue

        ## categorize events -------------------------------------------
        ## -------------------------------------------------------------
        event_cat_dict = GetEventCategories(t)
        if not event_cat_dict:
            continue




        ## dictionary of event weights by name
        ## needs to be done after pre-selection to get lepton SF correct
        ## -------------------------------------------------------------
        weight_dict = {}
        weight_lepton_SF = 1

        if m_el_n == 1:
            weight_lepton_SF = (t.weight_indiv_SF_EL_Trigger * t.weight_indiv_SF_EL_Reco *
                                t.weight_indiv_SF_EL_ID * t.weight_indiv_SF_EL_Isol)
        elif m_mu_n == 1:
            weight_lepton_SF = (t.weight_indiv_SF_MU_Trigger * t.weight_indiv_SF_MU_ID *
                                t.weight_indiv_SF_MU_Isol * t.weight_indiv_SF_MU_TTVA)

        else:
            print "WARNING :: did not find 1 el or 1 mu in a 1 lep event!!"

        # reweight VLQ samples according to decay type
        # right now returns 1
        weight_vlq = GetVLQDecayWeight(t)
        print 'weight_vlq', weight_vlq

        weight_nominal = (m_xsec * m_lumi * t.weight_mc * weight_lepton_SF * 
                          t.weight_pileup * weight_vlq / m_sumweights)
        print 'weight_pileup', t.weight_pileup
        print 'weight_nominal', weight_nominal


        if tree_name == m_tree_nominal[0]:
            if m_sample_sys == 'nominal':
                if m_isMC:
                    # for nominal tree, nominal sample_sys, MC, fill all event weights
                    weight_dict = FillEventWeightDict(t, weight_nominal)
                else:
                    # data should have nominal tree, nominal sample_sys -- all weights = 1
                    weight_dict['nominal'] = 1
            else:
                # for sys samples, use nominal weights, but label hists by the sample_sys
                weight_dict[m_sample_sys] = weight_nominal
        else:
            # for systematic tree, use nominal weights, but label hists by the tree_sys
            weight_dict[tree_name.replace('_Loose','')] = weight_nominal

        ## create and/or fill histograms
        ## -------------------------------------------------------------
        for iVar in m_plot_vars:
            FillHists(t, iVar, event_cat_dict, weight_dict, hist_dict)
        
        ## ========================================================== ##
        ## ========================== end =========================== ##
        ## ========================================================== ##
    return hist_dict

def GetVLQDecayWeight(tree):
    return 1.

def CalculateEventVariables(tree):

    ## number of jets / bjets
    global m_jet_n
    #global m_btagjet_n
    m_jet_n  = 0
    #m_btagjet_n = 0

    jet_pt     = tree.jet_pt
    jet_eta    = tree.jet_eta
    jet_mv2c20 = tree.jet_mv2c20
    for i in xrange(len(jet_pt)):
        if ( jet_pt[i] > m_cut_jet_pt and
             abs(jet_eta[i]) < m_cut_jet_eta):
            m_jet_n = m_jet_n + 1
            #if jet_mv2c20[i] > m_cut_btag:
            #    m_btagjet_n = m_btagjet_n + 1

    ## number of leptons
    global m_lep_n
    m_lep_n = 0
    
    global m_el_n
    m_el_n = 0
    el_pt = tree.el_pt
    el_eta = tree.el_eta
    for i in xrange(len(el_pt)):
        if ( el_pt[i] > m_cut_lep_pt and
             abs(el_eta[i]) < m_cut_lep_eta ):
            m_el_n = m_el_n + 1
            m_lep_n = m_lep_n + 1

    global m_mu_n
    m_mu_n = 0
    mu_pt = tree.mu_pt
    mu_eta = tree.mu_eta
    for i in xrange(len(mu_pt)):
        if ( mu_pt[i] > m_cut_lep_pt and
             abs(mu_eta[i]) < m_cut_lep_eta):
            m_mu_n = m_mu_n + 1
            m_lep_n = m_lep_n + 1
    
    ## number of mass-tagged rc jets
    global m_rcjet_n
    m_rcjet_n = 0
    rcjet_pt  = tree.reclustered_jets_pt
    rcjet_eta = tree.reclustered_jets_eta
    rcjet_m   = tree.reclustered_jets_m
    rcjet_nsub = tree.reclustered_jets_constituents_n
    for i in xrange(len(rcjet_pt)):
        if ( rcjet_pt[i] > m_cut_rcjet_pt and
             abs(rcjet_eta[i]) < m_cut_rcjet_eta and
             rcjet_m[i] > m_cut_rcjet_m  and rcjet_nsub[i]>=m_cut_rcjet_nsub):
            m_rcjet_n = m_rcjet_n + 1



def GetListOfBTagJets(tree, trf_cat):
    btag_jets = []

    jet_pt     = tree.jet_pt
    jet_eta    = tree.jet_eta
    jet_phi    = tree.jet_phi
    jet_m      = 0

    if trf_cat == '2ex' or trf_cat == '3ex' or trf_cat == '4in':
        jet_istagged = getattr(tree, 'trf_tagged_' + m_WP + '_' + trf_cat)
        for i in range(len(jet_istagged)):
            if jet_istagged[i]:
                tlv = ROOT.TLorentzVector(0,0,0,0)
                tlv.SetPtEtaPhiM(jet_pt[i], jet_eta[i], jet_phi[i], 0)
                btag_jets.append(tlv)

    else:
        jet_mv2c20 = tree.jet_mv2c20
        for i in range(len(jet_mv2c20)):
            if jet_mv2c20[i] > m_cut_btag:
                tlv = ROOT.TLorentzVector(0,0,0,0)
                tlv.SetPtEtaPhiM(jet_pt[i], jet_eta[i], jet_phi[i], 0)
                btag_jets.append(tlv)

    return btag_jets

def CalculateBTagVars(tree, doTRF):
    global m_btagjet_n

    if not doTRF:
        btag_jets = GetListOfBTagJets(tree, '')
        m_btagjet_n = len(btag_jets)

    else:
        for i,cat in enumerate(['2ex', '3ex', '4in']):
            btag_jets = GetListOfBTagJets(tree, cat)
            m_btagjet_n[i] = len(btag_jets)


def DoPreSelection(tree):

    if not(getattr(tree, '1el5j') or getattr(tree, '1mu5j')):
        return False

    CalculateEventVariables(tree)

    CalculateBTagVars(tree, m_doTRF)

    #print globals()["m_lep_n"]
    #print eval('m_lep_n')

    if (m_lep_n != 1) or (m_jet_n < 5):
        return False

    if (not m_doTRF) and (m_btagjet_n < 2):
        return False

    if (tree.met_met < 20000.) or (tree.met_met + tree.mTW < 60000.):
        return False

    return True

def GetEventCategories(tree):
    global m_isttbar
    event_cat_dict = {} ## given event might be in more than one category
                        ## i.e. inclusive/exclusive ntag
                        ## entries are {category : trf_weight}
    ## lepton category
    lep_cat = '' + 'el'*(m_el_n==1) + 'mu'*(m_mu_n==1)
    ## jet category
    jet_cat = '' + '5j'*(m_jet_n==5) + '6j'*(m_jet_n>=6) 
    ## b category
    if not m_doTRF:
        b_cat = '' + '2b'*(m_btagjet_n==2) + '3b'*(m_btagjet_n==3) + '4b'*(m_btagjet_n>=4)
    else:
        b_cat = ''
    ## HOT category
    ljet_cat = '' + '0TTRCLooser'*(m_rcjet_n==0) + '1TTRCLooser'*(m_rcjet_n==1) + '2TTRCLooser'*(m_rcjet_n>=2)

    event_cat_base = 'c' + str(m_lep_n) + 'l' + ljet_cat + jet_cat

    if m_isttbar:
        hfcat = tree.HF_SimpleClassification
        prefix = 'ttb'*(hfcat==1) + 'ttl'*(hfcat==0) + 'ttc'*(hfcat==-1)
        event_cat_base = prefix + event_cat_base

    # no trf, use btag category
    if not m_doTRF:
        event_cat = event_cat_base + b_cat
        if ( all(a in event_cat for a in ['6j', '1TTRCLooser'])
             and any(b in event_cat for b in ['3b', '4b'])  ):
            #print "Determining mBB cat for event_cat:", event_cat
            mbb_cat = getMbbCat(tree, m_cut_btag)
            event_cat = event_cat + mbb_cat
            #print event_cat
        event_cat_dict[event_cat] = getattr(tree, 'weight_bTagSF_' + m_WP)
        event_cat_dict[event_cat + '_' + lep_cat] = getattr(tree, 'weight_bTagSF_' + m_WP)

    # do trf, use all categories
    else:
        for [i, br] in [['2b','trf_weight_'+m_WP+'_2ex'], ['3b','trf_weight_'+m_WP+'_3ex'], ['4b','trf_weight_'+m_WP+'_4in']]:
            event_cat = event_cat_base + i
            if ( all(a in event_cat for a in ['6j', '1TTRCLooser'])
                 and any(a in event_cat for a in ['3b', '4b'])  ):
                #print "Determining mBB cat for event_cat:", event_cat
                mbb_cat = getMbbCatTRF(tree, i, m_WP)
                event_cat = event_cat + mbb_cat
                #print event_cat
            event_cat_dict[event_cat] = getattr(tree, br)
            event_cat_dict[event_cat + '_' + lep_cat] = getattr(tree, br)
    
    return event_cat_dict

def getMbbCat(tree, tagCut):

    jet_pt       = tree.jet_pt
    jet_eta      = tree.jet_eta
    jet_phi      = tree.jet_phi
    jet_m        = 0 #tree.jet_m
    jet_mv2c20   = tree.jet_mv2c20

    mBB = 0
    mindR = 200

    for i in range(0, len(jet_pt)-1):
        if not jet_mv2c20[i] > tagCut:
            continue
        tlv1 = ROOT.TLorentzVector(0,0,0,0)
        tlv1.SetPtEtaPhiM(jet_pt[i], jet_eta[i], jet_phi[i], 0)
        for j in range(i+1, len(jet_pt)):
            if not jet_mv2c20[j] > tagCut:
                continue
            tlv2 = ROOT.TLorentzVector(0,0,0,0)
            tlv2.SetPtEtaPhiM(jet_pt[j], jet_eta[j], jet_phi[j], 0)

            if tlv1.DeltaR(tlv2) < mindR:
                mindR = tlv1.DeltaR(tlv2)
                mBB = (tlv1 + tlv2).M()

    #print mBB
    if mBB > 100. * GeV:
        return "HighMbb"
    return "LowMbb"


def getMbbCatTRF(tree, ntagstr, WP):
    if ntagstr == '2b':
        tag_br = 'trf_tagged_' + WP + '_2ex'
    if ntagstr == '3b':
        tag_br = 'trf_tagged_' + WP + '_3ex'
    if ntagstr == '4b':
        tag_br = 'trf_tagged_' + WP + '_4in'

    jet_pt       = tree.jet_pt
    jet_eta      = tree.jet_eta
    jet_phi      = tree.jet_phi
    jet_m        = 0 # tree.jet_m
    jet_isTagged = getattr(tree, tag_br)


    mBB = 0
    mindR = 20
    for i in range(0, len(jet_pt)-1):
        if not jet_isTagged[i]:
            continue
        tlv1 = ROOT.TLorentzVector(0,0,0,0)
        tlv1.SetPtEtaPhiM(jet_pt[i], jet_eta[i], jet_phi[i], 0)
        for j in range(i+1, len(jet_pt)):
            if not jet_isTagged[j]:
                continue
            tlv2 = ROOT.TLorentzVector(0,0,0,0)
            tlv2.SetPtEtaPhiM(jet_pt[j], jet_eta[j], jet_phi[j], 0)

            if tlv1.DeltaR(tlv2) < mindR:
                mindR = tlv1.DeltaR(tlv2)
                mBB = (tlv1 + tlv2).M()
    #print mBB
    if mBB > 100. * GeV:
        return "HighMbb"
    return "LowMbb"

def FillHists(tree, plot_var, event_cat_dict, weight_dict, hist_dict):
    for iCat in event_cat_dict:

        # sample values for all systematics
        try:
            val = getattr(tree, plot_var.branch_name)
        except:
            val = globals()[plot_var.branch_name]
        
        if m_doTRF and plot_var.branch_name in m_btag_vars:
            if '2b' in iCat:
                val = val[0]
            if '3b' in iCat:
                val = val[1]
            if '4b' in iCat:
                val = val[2]

        if plot_var.units == 'GeV':
            val = val / GeV
        if type(val) in [float, int]:
            val = [val]

        for iSys in weight_dict:
            if m_doTRF and 'trf_weight' in iSys and not matchCatToSysForTRF(iCat, iSys):
                #print "Skipping trf sys for sys and category:", iSys, iCat
                continue
            #h_name = ('h_' + sample_name + '_'+ plot_var.title + '_' + iCat + '_' + iSys)
            sysName = GetStandardSysName(iSys)
            h_name = iCat + '_' + plot_var.title + (iSys != 'nominal') * ('_' + sysName)
            h = hist_dict.get(h_name, '')
            if not h:
                h = ROOT.TH1F(h_name, h_name, plot_var.nbinsx, plot_var.xlow, plot_var.xup)
                h.GetXaxis().SetTitle(plot_var.title + ' [' + plot_var.units + ']')
                h.GetYaxis().SetTitle("Entries / " + "{0:.1f}".format((plot_var.xup - plot_var.xlow) / plot_var.nbinsx ) 
                                      + " " + GetUnits(plot_var.title))
                h.SetDirectory(0)
                hist_dict[h_name] = h
            h = hist_dict[h_name]
            # now that it is made, fill
            for iVal in val:
                print iSys
                print 'weight_dict_weight:', weight_dict[iSys]
                print 'event_cat_weight', event_cat_dict[iCat]
                h.Fill(iVal, weight_dict[iSys] * event_cat_dict[iCat])

def GetStandardSysName(systematic_name):
    return systematic_name

def matchCatToSysForTRF(cat, sys):
    if '2b' in cat and '2ex' in sys:
        return True
    if '3b' in cat and '3ex' in sys:
        return True
    if '4b' in cat and '4in' in sys:
        return True
    return False

def GetUnits(str):
    return "[]"

def FillEventWeightDict(tree, weight_nominal):
    weight_dict = {}

    wLSF = 1
    if m_el_n == 1:
            wLSF = (tree.weight_indiv_SF_EL_Trigger * tree.weight_indiv_SF_EL_Reco *
                    tree.weight_indiv_SF_EL_ID * tree.weight_indiv_SF_EL_Isol)
    elif m_mu_n == 1:
            wLSF = (tree.weight_indiv_SF_MU_Trigger * tree.weight_indiv_SF_MU_ID *
                    tree.weight_indiv_SF_MU_Isol * tree.weight_indiv_SF_MU_TTVA)
    wPU = tree.weight_pileup

    weight_dict['nominal'] = weight_nominal
    for iW in m_weight_sys:
        if 'weight_pileup' in iW:
            weight_dict[iW] = weight_nominal * getattr(tree, iW) / wPU

        if 'weight_indiv_SF_EL' in iW and m_el_n == 1:
            wSplit = iW.split('_')
            nom_indiv_weight = ( wSplit[0] + '_' + wSplit[1] + '_' + wSplit[2] + 
                                 '_' + wSplit[3] + '_' + wSplit[4] )
            #print "Found nom indiv weight", nom_indiv_weight, " for weight sys", iW
            weight_dict[iW] = weight_nominal * getattr(tree, iW) / getattr(tree, nom_indiv_weight)
        
        if 'weight_indiv_SF_MU' in iW and m_mu_n == 1:
            wSplit = iW.split('_')
            nom_indiv_weight = ( wSplit[0] + '_' + wSplit[1] + '_' + wSplit[2] +
                                 '_' + wSplit[3] + '_' + wSplit[4] )
            #print "Found nom indiv weight", nom_indiv_weight, " for weight sys", iW
            weight_dict[iW] = weight_nominal * getattr(tree, iW) / getattr(tree, nom_indiv_weight)

        
        if 'weight_bTagSF' in iW and not m_doTRF:
            if 'eigenvars' in iW:
                eigenvars = getattr(tree, iW)
                for i,iVar in enumerate(eigenvars):
                    weight_dict[iW+'_NP'+str(i)] = weight_nominal * iVar / getattr(tree, 'weight_bTagSF_' + m_WP)
            else:
                weight_dict[iW] = weight_nominal * getattr(tree, iW) / getattr(tree, 'weight_bTagSF_' + m_WP)
        #TRF handled differently
        if 'trf_weight' in iW and m_doTRF:
            wSplit = iW.split('_')
            nom_trf_weight = (wSplit[0] + '_' + wSplit[1] + '_' + wSplit[2] + '_' +  wSplit[3]) #trf_weight_77_2ex_eigenvars_B_up
            if 'eigenvars' in iW:
                eigenvars = getattr(tree, iW)
                for i, iVar in enumerate(eigenvars):
                    weight_dict[iW+'_NP'+str(i)] = weight_nominal * iVar / getattr(tree, nom_trf_weight)
            else:
                weight_dict[iW] = weight_nominal * iVar / getattr(tree, nom_trf_weight)

    return weight_dict

def ApplySampleSliceVeto(t):

    ## use HT filter for now!
    if (m_dsid == 410000):
        if tree.HT_truth > 600*GeV:
            return False

    ## dont use b-filtered sample
    if(m_dsid == 410120):
        return False

    ## dont use MET filtered sample
    if(m_dsid == 407012):
        return False

    return True


def ParseCommandLineArguments(argv):
    global m_input_file_str
    global m_xsec_file_str
    global m_sumweights_file_str
    global m_lumi
    global m_doTRF

    ## parse options 
    _short_options = 'hf:x:s:l:t'
    _long_options = ['help', 'inputFile=', 'xSecFile=', 'sumWeightsFile=','lumi=', 'trf']
    try:
        opts, args = getopt.gnu_getopt(argv, _short_options, _long_options)
    except getopt.GetoptError:
        print 'ERROR : getopt.GetoptError\n'
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
        if opt in ('-t', '--trf'):
            m_doTRF = True

    print "\nINFO : Summary of configuration for this run --"
    print "       Input file      :", m_input_file_str
    print "       X-sec file      :", m_xsec_file_str
    print "       Sumweights file :", m_sumweights_file_str
    print "       Luminosity      :", m_lumi,"\n"
    print "       Apply TRF       :", boolTOstr(m_doTRF)

    # check that we have everything
    if not os.path.isfile(m_input_file_str):
        sys.exit("ERROR : please enter a valid input file with '-f <input_file>' option")
    if not os.path.isfile(m_xsec_file_str):
        sys.exit("ERROR : please enter a valid xsec file with '-x <sec_file>' option")
    if not os.path.isfile(m_sumweights_file_str):
        sys.exit("ERROR : please enter a valid sumweights file with '-s <sumweights_file>' option")


def boolTOstr(bool):
    if bool:
        return "True"
    return "False"

def GetFileMetaData(root_file): # root_file is an opened TFile
    print "\nINFO : Summary of MetaData --"

    # dsid ----------------------------------------------------------
    global m_dsid
    m_dsid = GetDSIDFromSumWeightsTree(root_file)
    if m_dsid < 0:
        sys.exit("ERROR : Could not get dsid for this sample")

    print "       DSID        :", m_dsid

    # isMC ----------------------------------------------------------
    global m_isMC
    m_isMC = bool(m_dsid)

    print "       isMC        :", m_isMC

    # tag -----------------------------------------------------------
    global m_file_tag_str
    m_file_tag_str = GetTagFromFileName(m_input_file_str)
    if not m_file_tag_str:
        sys.exit("ERROR : could not determine file tag from file name: " + m_input_file_str)

    # MC only -------------------------------------------------------
    if m_isMC:
        # x-sec -----------------------------------------------------
        global m_xsec
        m_xsec = GetXSecFromFile(m_dsid, m_xsec_file_str)
        if m_xsec < 0:
            sys.exit("ERROR : Could not get xsec for this sample")

        print "       XSec        :", m_xsec

        # sum of event weights --------------------------------------
        global m_sumweights
        m_sumweights = GetSumWeightsFromFile(m_dsid, m_file_tag_str, m_sumweights_file_str)
        if m_sumweights < 0:
            sys.exit("ERROR : Could not get sumweights for this sample")

        print "       SumWeights  :", m_sumweights

    # sample name, file systematic
    global m_sample_name
    global m_sample_sys
    m_sample_name, m_sample_sys = GetSampleNameAndSysFromDSIDAndTag(m_dsid, m_file_tag_str)
    if not m_sample_name:
        sys.exit("ERROR : Could not determine the sample name for dsid " + str(m_dsid)
                 + " and tag " + m_file_tag_str)

    print "       Sample name :", m_sample_name
    print "       Sample sys  :", m_sample_sys, '\n'

def GetDSIDFromSumWeightsTree(root_file):
    dsid = -1
    t = root_file.Get("sumWeights")
    if not t:
        sys.exit("ERROR : sumWeights tree not present.\n        Could not determine DSID. Exiting.")
    t.SetBranchStatus('dsid', 1)
    for i in range(t.GetEntries()):
        t.GetEntry(i)
        dsid = t.dsid
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
        if len(line.split()) < 3:
            continue
        if int(line.split()[0]) == dsid:
            xsec = float(line.split()[1]) * float(line.split()[2])
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
        if int(line.split()[0]) == dsid and line.split()[1] == file_tag:
            sumweights = float(line.split()[2])
            foundWeights = True
            break
    f.close()
    if not foundWeights:
        f = open(m_runDir + '/missingWeight.txt', 'a')
        f.write(str(dsid) + ' ' + file_tag + '\n')
        f.close()
    return sumweights




if __name__ == "__main__" : main(sys.argv)





"""
def GetSampleNameAndSysFromDSIDAndTag(dsid, file_tag_str):
    sample_name = ''
    sample_sys  = 'nominal'
    global m_isttbar

    # data
    if dsid == 0:
        print "This sample is data. Did you expect this for file", input_file, " ?\n"
        sample_name = 'data'

    # backgrounds
    elif dsid == 410000:
        sample_name = 'ttbar'
        if file_tag_str == 'e3698_a766_a810_r6282_p2516':
            sample_sys = 'afii'
    elif dsid in range(407009, 407012):
        sample_name = 'ttbar'
    elif dsid in range(410001, 410004+1):
        sample_name = 'ttbar'
        if dsid == 410001:
            sample_sys = 'radHi'
        if dsid == 410002:
            sample_sys = 'radLo'
        if dsid == 410003:
            sample_sys = 'aMCAtNLOHerwigpp'
        if dsid == 410004:
            sample_sys = 'Herwigpp'
    elif dsid in range(361300, 361371+1):
        sample_name = 'wjets'
    elif dsid in range(361372, 361467+1):
        sample_name = 'zjets'
    elif dsid in range(361520, 361534+1):
        sample_name = 'wjets'
        sample_sys = 'MG5Py8'
    elif dsid in range(361500, 361519+1):
        sample_name = 'zjets'
        sample_sys = 'MG5Py8'
    elif dsid in [410013, 410014, 410011, 410012, 410025, 410026]:
        sample_name = 'singletop'
    elif dsid in [410017, 410018, 410019, 410020, 410099, 410100, 410101, 410102]:
        sample_name = 'singletop'
        if dsid == 410017:
            sample_sys = 'radLo'
        if dsid == 410018:
            sample_sys = 'radHi'
        if dsid == 410019:
            sample_sys = 'radLo'
        if dsid == 410020:
            sample_sys = 'radHi'
        if dsid == 410099:
            sample_sys = 'radHi'
        if dsid == 410100:
            sample_sys = 'radLo'
        if dsid == 410101:
            sample_sys = 'radHi'
        if dsid == 410102:
            sample_sys = 'radLo'
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
    
    if sample_name == 'ttbar':
        m_isttbar = True
        
    return sample_name, sample_sys
"""
