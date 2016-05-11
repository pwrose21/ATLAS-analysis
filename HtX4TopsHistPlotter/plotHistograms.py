#!/usr/bin/env python

## /////// ##
## imports ##
## /////// ##
import ROOT
import sys, commands, copy
import math
## /////// ##
## globals ##
## /////// ##

doSysts = False
logPlots = True
histMarginFact = 1.5
signal_names = ['UED1400']
signal_draw_names = ['UED1400']

## /////////// ##
## atlas style ##
## /////////// ##
ROOT.gROOT.LoadMacro("/export/home/prose/myPrograms/scripts/atlasstyle-00-03-05/AtlasUtils.C")
ROOT.gROOT.LoadMacro("/export/home/prose/myPrograms/scripts/atlasstyle-00-03-05/AtlasStyle.C")
ROOT.gROOT.LoadMacro("/export/home/prose/myPrograms/scripts/atlasstyle-00-03-05/AtlasLabels.C")
ROOT.SetAtlasStyle()

def main(argv):
    ROOT.gROOT.SetBatch()

    ## testing
    m_input_file_dir = argv[1]

    # this is a dictionary of dictionaries
    # the keys in m_hists are the histogram regions
    #     i.e. c1lnTTRCblah_meff_systs
    # the values in m_hists are dicts containing 
    #     whose keys are the various samples and
    #     whose values are the actual histograms

    # Get list of input files
    m_files = GetListOfInputFiles(m_input_file_dir)

    # Get list of regions
    m_regions = GetListOfRegions(m_files)

    # canvas containing two pads
    # p1 for the stack histogram and signals and data
    # p2 for data / mc
    c, p1, p2 = SetupCanvasForPlotting()

    m_hists = {}
    # one region at a time
    for iReg in m_regions:
        ##testing
        if not iReg == 'c1l1TTRCLooser6j4b_meff':
            continue
        ##-------

        # m_hists[sample] = TH1F
        m_hists  = GetHistsFromFiles(m_files, iReg)
        m_hstack = MakeStackHist(m_hists, iReg)
        m_hratio = ''
        m_legend = MakeLegend(m_hists)
 
        if m_hists.get('data', '') and m_hists.get('background', ''):
            # testing
            m_hists['data'].Scale(1.5*4)
            m_hists['background'].Scale(4)
            m_hratio = MakeDataOverBkgdHist(m_hists)

        ## -- plotting -- ##
        p2.cd()
        m_hratio.Draw("")

        p1.cd()
        DrawPlotsToUpperPad(m_hists, m_hstack, iReg)

        ## -- legend and labels -- ##
        m_legend.Draw()
        ROOT.myText(       0.41,  0.85, 1, "#sqrt{s}= 13 TeV")
        ROOT.myText(       0.41,  0.80, 1, iReg)
        ROOT.ATLASLabel(0.41,0.75,"Internal")

        
        c.cd()
        c.Print('plots/'+iReg+'.eps')


    sys.exit()


##-------------------------end-of-main----------------------------##
"""
def FindXForLabels(hist):
    nbins = hist.GetNBinsX():
    frac = nbins *
    maxX = nbins / 2
""" 

def MakeLegend(hists):

    ## --need to sort again to get same order as added-- ##
    sortedHistTuple = sorted(hists.items(),
                             key=lambda x: x[1].GetBinContent(x[1].GetMaximumBin()))

    l = ROOT.TLegend(0.2, 0.6, 0.4, 0.8)
    for iHTuple in sortedHistTuple:
        if not SampleExcludedFromStack(iHTuple[0]):
            l.AddEntry(iHTuple[1], iHTuple[0], 'f')
        elif iHTuple[0] in signal_draw_names:
            l.AddEntry(iHTuple[1], iHTuple[0], 'l')
    l.AddEntry(hists['data'], 'data', 'p')

    l.SetBorderSize(0)
    l.SetFillColor(0)
    l.SetTextSize(0.04)

    return l


#------------------------------------------------------------------------
#------------------------------------------------------------------------

def DrawPlotsToUpperPad(hists, stackHist, region):
    global logPlots
    global histMarginFact

    # use the sum of of backgrounds to define the range for the plot
    # todo -- instead, use signal(s), background, data
    bkgdHist = hists.get('background', '')    
    if not logPlots:
        bkgdHist.SetMaximum(histMarginFact * 
                            (bkgdHist.GetMaximum() - bkgdHist.GetMinimum()) +
                            bkgdHist.GetMinimum())
    else:
        minBinVal = min(0.1*GetSmallestBinValAboveZero(bkgdHist), 0.1)
        bkgdHist.SetMinimum(minBinVal)
        bkgdHist.SetMaximum(minBinVal *
                            math.pow(bkgdHist.GetMaximum() / minBinVal, histMarginFact))

    

    PrepHistForLineDraw(bkgdHist)
    bkgdHist.Draw()
    stackHist.Draw("HISTsame")

    hists['data'].Draw("same")
    bkgdHist.Draw("axissame")

def GetSmallestBinValAboveZero(hist):
    minBinVal = 10
    for i in range(1, hist.GetNbinsX()+1):
        binContent = hist.GetBinContent(i)
        if (binContent > 0) and (binContent < minBinVal):
            minBinVal = binContent
    return minBinVal

#------------------------------------------------------------------------
#------------------------------------------------------------------------

def SetupCanvasForPlotting():
    global logPlots

    c = ROOT.TCanvas("c", "c", 0, 0, 600, 600)

    # normal hists
    p1 = ROOT.TPad("p1", "p1", 0, 0.3, 1, 1.0);
    p1.SetBottomMargin(0) # Upper and lower plot are joined
    #p1.SetGridx()         # Vertical grid
    if logPlots:
        p1.SetLogy(ROOT.kTRUE)
    p1.Draw()            # Draw the upper pad: p1

    # data / background hists
    c.cd()          # Go back to the main canvas before defining p2
    #p2 = ROOT.TPad("p2", "p2", 0, 0.05, 1, 0.3)
    p2 = ROOT.TPad("p2", "p2", 0, 0.0, 1, 0.3)
    p2.SetTopMargin(0)
    #p2.SetBottomMargin(0.2)
    p2.SetBottomMargin(0.3)
    #p2.SetGridx() # vertical grid
    p2.Draw()
    c.cd()

    return c, p1, p2


#------------------------------------------------------------------------
#------------------------------------------------------------------------

def MakeDataOverBkgdHist(hists):
    dataHist = hists.get('data', '')
    bkgdHist = hists.get('background', '')
    ratioHist = ''

    if dataHist and bkgdHist:
        ratioHist = copy.deepcopy(dataHist)
        ratioHist.SetDirectory(0)
        ratioHist.Divide(bkgdHist)

    # -- set range -- #
    histMin = GetSmallestBinValAboveZero(ratioHist)
    ratioHist.SetMinimum(histMin*0.8)
    ratioHist.SetMaximum(ratioHist.GetMaximum()*1.2)

    # -- labels and formatting -- #
    ratioHist.GetYaxis().SetTitle("Data / MC")
    ratioHist.GetXaxis().SetTitleSize(ratioHist.GetXaxis().GetTitleSize()*0.7/0.3)
    ratioHist.GetXaxis().SetTitleOffset(ratioHist.GetXaxis().GetTitleOffset()*0.3*1.7/0.7)
    ratioHist.GetXaxis().SetLabelSize(ratioHist.GetXaxis().GetLabelSize()*0.7/0.3)

    ratioHist.GetYaxis().SetLabelSize(ratioHist.GetYaxis().GetLabelSize()*0.7/0.3)
    ratioHist.GetYaxis().SetLabelOffset(ratioHist.GetYaxis().GetLabelOffset()*5)
    ratioHist.GetYaxis().SetTitleSize(ratioHist.GetYaxis().GetTitleSize()*0.7/0.3)
    ratioHist.GetYaxis().SetTitleOffset(ratioHist.GetYaxis().GetTitleOffset()*0.3/0.7)
    ratioHist.GetYaxis().SetNdivisions(506)

    return ratioHist

def MakeStackHist(hists, region):
    hs = ROOT.THStack(region, region)
    
    sortedHistTuple = sorted(hists.items(),
                             key=lambda x: x[1].GetBinContent(x[1].GetMaximumBin()))

    for iHTuple in sortedHistTuple:
        sample = iHTuple[0]
        if SampleExcludedFromStack(sample):
            continue
        print "adding hist to stack for sample", sample
        h = iHTuple[1]
        PrepHistForStack(h, sample)
        hs.Add(h)
    
    return hs
        


def GetHistsFromFiles(file_list, region):
    hists = {}
    for iFile in file_list:
        f = ROOT.TFile(iFile, 'READ')
        # get the sample type
        # a given sample might have 2 types e.g. ttbar is also bkgd
        sample_type_list = GetSampleTypesFromFileName(iFile)
        for iSample in sample_type_list:
            # try to retrieve the hist from the dict
            hist_from_dict = hists.get(iSample, '')
            # if it already exists, add the contribution from this file
            if hist_from_dict:
                hist_from_dict.Add(f.Get(region))
            # otherwise, we need to add it to the dict
            else:
                h = f.Get(region)
                h.SetDirectory(0)
                hists[iSample] = h
        f.Close()
    return hists

#------------------------------------------------------------------------
#------------------------------------------------------------------------                

def GetListOfRegions(file_list):
    regions = []
    for iFile in file_list:
        f = ROOT.TFile(iFile, 'READ')

        for iKey in f.GetListOfKeys():
            # just do TH1 ... not fitting TH2s
            if not 'TH1' in iKey.GetClassName():
                continue
            name = iKey.GetName()
            
            # check that it's not a cutflow or other hist
            if not 'c1l' in name:
                continue

            # might want to only do nominal!
            if doSysts:
                print "Systematics no yet supported. Doing nominal only"
            # if not doSysts
            if True: # for now, no systs
                tmpName = name.replace('_e_', '').replace('_mu_', '').replace('b_', '')
                if '_' in tmpName:
                    # debug
                    #print "Skipping hist", name, "because it was found to be a systematic hist"
                    continue
            #print "Proceeding with hist:", name
            regions.append(name)

        f.Close()

    return regions

########################################################################################################
#########################################new dev #######################################################
########################################################################################################







def PrepHistForStack(hist, sample):
    hist.SetFillColor(GetHistSampleColor(sample))
    hist.SetMarkerStyle(1)
    hist.SetMarkerSize(1.0)

#------------------------------------------------------------------------
#------------------------------------------------------------------------

def PrepHistForLineDraw(hist):
    hist.SetMarkerStyle(1)
    hist.SetMarkerSize(1.0)

#------------------------------------------------------------------------
#------------------------------------------------------------------------

def GetHistSampleColor(sample):
    if sample == 'tt+bb':
        return ROOT.kBlue
    if sample == 'tt+cc':
        return ROOT.kRed
    if sample == 'tt+light':
        return ROOT.kCyan
    if sample == 'others':
        return ROOT.kGreen

    return ROOT.kBlack
#------------------------------------------------------------------------
#------------------------------------------------------------------------

def SampleExcludedFromStack(name):
    global signal_names

    names = ['signal', 'background', 'data']
    ## exclude
    if any(a==name for a in names):
        return True
    if any(a==name for a in signal_names):
        return True

    ## otherwise include
    return False
    
#------------------------------------------------------------------------
#------------------------------------------------------------------------

def UpdateHistDictFromFile(file_name_str, histDict):

    # ttbb, signal, background, data, ttcc, diboson, others
    sample_type_list = GetSampleTypesFromFileName(file_name_str)

    # get the actual TFile
    f = ROOT.TFile(file_name_str, 'READ')
    
    # Find out what's here
    keys = f.GetListOfKeys()

    for iKey in keys:
        # just do TH1 ... not fitting TH2s
        if not 'TH1' in iKey.GetClassName():
            continue
        name = iKey.GetName()
        
        # check that it's not a cutflow or other hist
        if not 'c1l' in name:
            continue

        # might want to only do nominal!
        if doSysts:
            print "Systematics no yet supported. Doing nominal only"
        #if not doSysts:
        if True: # for now, no systs
            tmpName = name.replace('_e_', '').replace('_mu_', '').replace('b_', '')
            if '_' in tmpName:
                # debug
                #print "Skipping hist", name, "because it was found to be a systematic hist"
                continue
        #print "Proceeding with hist:", name

        # Get the histogram
        fileHist = f.Get(name)

        for iT in sample_type_list:
            histRegDict = histDict.get(name, '')
            if not histRegDict:
                histDict[name] = {}
                histRegDict = histDict[name]
            dictHist = histRegDict.get(iT, '')
            if not dictHist:
                h = f.Get(name)
                h.SetDirectory(0)
                #h.SetFillColor(ROOT.kBlue)
                histRegDict[iT] = h
            else:
                h = histRegDict[iT]
                h.Add(f.Get(name))                
    #print histDict

    f.Close()
    
    return True

#------------------------------------------------------------------------
#------------------------------------------------------------------------

def GetSampleTypesFromFileName(file_name_str):
    sample_types = []

    if 'UED' in file_name_str:
        if '1400' in file_name_str:
            sample_types.append('UED1400')

    sample_types.append('tt+bb')
    sample_types.append('tt+cc')
    sample_types.append('tt+light')
    sample_types.append('others')
    sample_types.append('signal')
    sample_types.append('background')
    sample_types.append('data')

    return sample_types

#------------------------------------------------------------------------
#------------------------------------------------------------------------

def GetListOfInputFiles(dir):
    if not dir.endswith('/'):
        dir = dir + '/'
    file_list = []
    status, files = commands.getstatusoutput('ls ' + dir)
    files = files.split('\n')
    for iFile in files:
        if not '.root' in iFile:
            continue
        file_list.append(dir + iFile)

    return file_list

#------------------------------------------------------------------------
#------------------------------------------------------------------------
#------------------------------------------------------------------------
#------------------------------------------------------------------------
    
if __name__ == "__main__" : main(sys.argv)