# ATLAS-analysis

## ---------------------------------------------------------------------------------------------------------------------------------
## HtX4TopsHistMaker ---------------------------------------------------------------------------------------------------------------
## ---------------------------------------------------------------------------------------------------------------------------------

Runs on the output of HtX4TopsNtuple

Usage : ./makeHistograms -f <input_file> -x <xsec_file> -s <sumweights_file> -l <luminosity>

Global variables :

These are set in ParseCommandLineOptions(argv):
-----------------------------------------
m_input_file_str -- location of the input root file
m_xsec_file_str  -- location of the file containing sample cross sections
m_sumweights_file_str -- location of the file containing the sum of event weights for each sample
m_lumi -- the luminosity value to scale to -- default 0.001 


These are set in GetFileMetaData(root_file):
---------------------------------
m_dsid -- the dataset ID for this file.  0 if data
       -- GetDSIDFromSumWeightsTree(root_file)
m_isMC -- 0 if data, 1 if MC
       -- bool(m_dsid)
m_file_tag_str -- the file tag e.g. e4110_r9999_s3245
	       -- GetTagFromFileName(m_input_file_str)
m_xsec -- the cross section for the sample. 0 if data
       -- GetXSecFromFile(m_dsid, m_xsec_file_str)
m_sumweights -- the sum of event weights for this sample. 0 if data
	     -- GetSumWeightsFromFile(m_dsid, m_file_tag_str, m_sumweights_file_str)
m_sample_name -- see below
m_sample_sys -- GetSampleNameAndSysFromDSIDAndTag(m_dsid, m_file_tag_str)


These are set above main:
-----------------------------
m_tree_nominal -- the nominal TTree
m_tree_sys     -- the systematic TTrees
m_weight_sys   -- list of weight systematics
m_plot_vars    -- list of variables to plot
m_branches     -- list of branches needed that are common to data and mc
m_branchesMC   -- list of branches needed only for MC
