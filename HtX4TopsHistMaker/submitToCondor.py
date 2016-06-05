import condor, os, commands, time, sys

# -- the program to run --
exe = '/export/home/prose/ATLAS-analysis/HtX4TopsHistMaker/makeHistograms.py'

# -- master directory and subdirectories --
input_dir = '/export/share/diskvault2/prose/2t2bAnalysis/HtX4TopsNtuples/00-00-04/'
status, subDirs = commands.getstatusoutput('ls ' + input_dir)
subDirs = subDirs.split('\n')


"""
rundirs = [
    'user.farooque.410066.MadGraphPythia8EvtGen.DAOD_TOPQ1.e4111_s2608_s2183_r7326_r6282_p2516.HtX4Tops_00-00-04_V0_output.root',
    'user.farooque.410067.MadGraphPythia8EvtGen.DAOD_TOPQ1.e4111_s2608_s2183_r7326_r6282_p2516.HtX4Tops_00-00-04_V0_output.root',
    'user.farooque.410068.MadGraphPythia8EvtGen.DAOD_TOPQ1.e4111_s2608_s2183_r7326_r6282_p2516.HtX4Tops_00-00-04_V0_output.root',
    'user.farooque.410073.MadGraphPythia8EvtGen.DAOD_TOPQ1.e4631_s2726_r7326_r6282_p2516.HtX4Tops_00-00-04_V0_output.root',
    'user.farooque.410074.MadGraphPythia8EvtGen.DAOD_TOPQ1.e4631_s2726_r7326_r6282_p2516.HtX4Tops_00-00-04_V0_output.root',
    'user.farooque.410080.MadGraphPythia8EvtGen.DAOD_TOPQ1.e4111_s2608_s2183_r7326_r6282_p2516.HtX4Tops_00-00-04_V0_output.root',
    'user.farooque.410081.MadGraphPythia8EvtGen.DAOD_TOPQ1.e4111_s2608_s2183_r7326_r6282_p2516.HtX4Tops_00-00-04_V0_output.root',
    ]

rundirs = [
    'user.prose.410000.PowhegPythiaEvtGen.DAOD_TOPQ1.e3698_s2608_s2183_r7267_r6282_p2516.HtX4Tops_00-00-06_output.root',
    'user.prose.407009.PowhegPythiaEvtGen.DAOD_TOPQ1.e4023_s2608_r7326_r6282_p2516.HtX4Tops_00-00-06_output.root',
    'user.prose.407010.PowhegPythiaEvtGen.DAOD_TOPQ1.e4023_s2608_r7326_r6282_p2516.HtX4Tops_00-00-06_output.root',
    'user.prose.407011.PowhegPythiaEvtGen.DAOD_TOPQ1.e4023_s2608_r7326_r6282_p2516.HtX4Tops_00-00-06_output.root',
    'user.prose.407012.PowhegPythiaEvtGen.DAOD_TOPQ1.e4023_s2608_r7326_r6282_p2516.HtX4Tops_00-00-06_output.root',
    ]
"""

# -- loop subdirectories --
for i,iDir in enumerate(subDirs):

    if not 'output.root' in iDir:
        continue

    if not '410000' in iDir:
        continue
    #if not iDir in rundirs:
    #    continue

    #  -- arg template --
    # for 2ued, no TRF
    arg_template = ' -f %s'
    if (any(a in iDir for a in ['302055', '302056', '302057', '302058', '302059'])
        or 'physics_Main' in iDir):
        pass
    else:
        pass
        #arg_template = arg_template + ' -t'

    if 'TOPQ4' in iDir:
        continue

    #if not 'physics_Main' in iDir:
    #    continue

    # -- dir name --
    # separate condor dir for each sample
    # data example : user.prose.00278968.physics_Main.DAOD_TOPQ4.f628_m1497_p2452.HtX4Tops_00-00-04_output.root/
    # mc example :   user.prose.302057.MadGraphPythia8EvtGen.DAOD_TOPQ1.e4017_s2608_s2183_r6869_r6282_p2516.HtX4Tops_00-00-04_output.root
    dirname = 'condorRun-' + time.strftime("%d_%m_%Y") + '-' + iDir.split('.')[2] + '-' + iDir.split('.')[5] + '-output'
    if '-t' in arg_template:
        dirname = dirname + '_withTRF'

    input_files = []
    status, files = commands.getstatusoutput('ls ' + input_dir + iDir)
    files = files.split('\n')
    
    for iFile in files:
        input_files.append(input_dir + iDir + '/' + iFile)

    try:
        condor.run(exe, arg_template, input_files, dirname = dirname, n_files = 1)
    except:
        print "Problem submitting from directory:", iDir
        print "Need to inspect files:", input_files
    
