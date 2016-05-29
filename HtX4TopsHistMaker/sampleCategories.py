def GetSampleNameAndSysFromDSIDAndTag(dsid, file_tag_str):
    sample_name = ''
    sample_sys  = 'nominal'

    # data
    if dsid == 0:
        print "INFO : This sample is data. Did you expect this?"
        sample_name = 'data'


    # ttbar -- non all had
    elif dsid in [410000, 410001, 410002, 410003,410004, 
                  410120, 407009, 407010, 407011, 407012]:
        sample_name = 'ttbar'
        
        if dsid == 410001: sample_sys = 'radHi'
        if dsid == 410002: sample_sys = 'radLo'
        if dsid == 410003: sample_sys = 'aMcAtNloHerwigpp'
        if dsid == 410004: sample_sys = 'PowhegHerwigpp'


    elif dsid in range(410011, 410016+1) or dsid in [410025, 410026]:
        sample_name = 'singletop'

    elif dsid in range(410066, 410068+1):
        sample_name = 'ttW'

    elif dsid in range(410073, 410075+1):
        sample_name = 'ttZ'

    elif dsid == 410080:
        sample_name = 'fourtopSM'

    elif dsid == 410081:
        sample_name = 'ttWW'

    elif dsid in [410111, 410112]:
        sample_name = 'ttee'

    elif dsid in [410113, 410113]:
        sample_name = 'ttmumu'

    elif dsid in [410115, 410116]:
        sample_name = 'tttautau'

    elif dsid in range(343365,343367+1):
        sample_name = 'ttH'

    elif ( (dsid in range(361331, 361354+1)) or 
           (dsid in range(363436, 363483+1)) ):
        sample_name = 'wjets-sherpa22'

    elif ( (dsid in range(363102, 363122+1)) or
           (disd in range(363361, 363435+1)) ):
        sample_name = 'zjets-sherpa22'

    elif dsid in range(361300, 361371+1):
        sample_name = 'wjets-sherpa21'

    elif dsid in range(361372, 361467+1):
        sample_name = 'zjets-sherpa21'

    elif dsid in range(361091, 361097+1):
        sample_name = 'diboson'

    elif dsid in range(302468, 302470+1): 
        sample_name = 'tts' + str(500 + 100*(dsid - 302468))

    elif dsid in range(302471, 302480+1):
        sample_name = 'tts' + str(750 + 50*(dsid - 302471))

    elif dsid in [302481, 302482]:
        sample_name = 'tts' + str(1300 + 100*(dsid - 302481))

    elif dsid in range(302483, 302485+1):
        sample_name = 'ttd' + str(700 + 250*(dsid - 302483))

    elif dsid in range(302486, 302488+1): 
        sample_name = 'bbs' + str(500 + 100*(dsid - 302486))

    elif dsid in range(302489, 302498+1):
        sample_name = 'bbs' + str(750 + 50*(dsid - 302489))

    elif dsid in [302499, 302500]:
        sample_name = 'bbs' + str(1300 + 100*(dsid - 302499))

    elif dsid in range(302501, 302503+1):
        sample_name = 'bbd' + str(700 + 250*(dsid - 302501))

    elif dsid == 502504:
        sample_name = 'xx500'

    elif dsid in range(302505, 302510+1):
        sample_name = 'xx' + str(700 + 100*(dsid - 302505))

    elif dsid == 302511:
        sample_name = 'xx1400'

    elif dsid == 302512:
        sample_name = 'yy500'

    elif dsid in range(302513, 302518+1):
        sample_name = 'yy' + str(700 + 100*(dsid - 302513))

    elif dsid == 302519:
        sample_name = 'yy1400'

    elif dsid in range(302055, 302059+1):
        sample_name = 'ued' + str(1000 + 200*(dsid - 302055))

    elif dsid == 302777:
        sample_name = 'fourtopCI'

    else:
        print "Could not determine the sample for this file"

    return sample_name, sample_sys



"""
######################################################################










    # ttbar
    elif dsid == 410000 or dsid==407009:
        sample_name = 'ttbar-nah'
        #ttbar afii
        if file_tag_str == 'e3698_a766_a810_r6282_p2516':
            sample_sys = 'afii'

    # ttbar systematics
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
    
    # w/z+jets sherpa 2.1
    elif dsid in range(361300, 361371+1):
        sample_name = 'wjets'
    elif dsid in range(361372, 361467+1):
        sample_name = 'zjets'

    # w/z+jets Madgraph5+Pythia
    elif dsid in range(361520, 361534+1):
        sample_name = 'wjets'
        sample_sys = 'MG5Py8'
    elif dsid in range(361500, 361519+1):
        sample_name = 'zjets'
        sample_sys = 'MG5Py8'

    # single top
    elif dsid in [410013, 410014, 410011, 410012, 410025, 410026]:
        sample_name = 'singletop'

    # single top systematics
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
    # diboson
    elif dsid in range(361081, 361087+1):
        sample_name = 'diboson'

    # ttV
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
        print "Could not determine the sample for this file"

    return sample_name, sample_sys
"""
