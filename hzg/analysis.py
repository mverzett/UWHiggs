#!/usr/bin/env python

#ROOT and low level analysis objects
import ROOT
from ROOT import TTree, TFile, TLorentzVector, TVector3, \
     TRandom3
from MOOSEY.cuts import CompositeCutflow
from MOOSEY.trees import tree_manager

#general cutsets
from hzg_cuts import muon_triggers_data,muon_triggers_mc, \
     electron_triggers_data, electron_triggers_mc, \
     mu_cuts_data,mu_cuts_mc,e_cuts_data,e_cuts_mc, photon_cuts_data, \
     photon_cuts_noifsr, mumu_cuts,ee_cuts,ell_gamma_dr, event_has_ifsr, \
     photon_cuts_mc,photon_cuts_plj

#correction layer
from corrections import setup_corrections

#python standard things
from optparse import OptionParser
from math import fabs
import os,sys
import numpy as np

#definitions of various tree types
import result_trees as outTrees

#processes one input file
#plans are to have concatenation
#options are:
# --runType={data,A,B,AB}
# --leptonType={muon,electron}
#special modes
# --calcCS, implies two more options, --dataInput --signalMC

muonBranches = ['nMu','muPt','muEta','muPhi','muType',
                'muNumberOfValidMuonHits',
                'muNumberOfValidTrkHits','muNumberOfValidPixelHits',
                'muStations','muChi2NDF','muPVD0','muPVDz',
                'muIsoTrk','muIsoEcal','muIsoHcal',
                'nZmumu','ZmumuLeg1Index','ZmumuLeg2Index','ZmumuMass']

electronBranches = ['nEle','elePt','eleEta','elePhi',
                    'eleE3x3','eleSCEta','eleSCEn','eleSCRawEn',
                    'eleConvMissinghit','eleConvDist','eleConvDcot',
                    'eleSigmaIEtaIEta','eledPhiAtVtx','eledEtaAtVtx','elePVD0',
                    'elePVDz','eleIsoTrkDR03','eleIsoEcalDR03',
                    'eleIsoHcalSolidDR03','eleGenIndex','eleGenGMomPID',
                    'nZee','ZeeLeg1Index','ZeeLeg2Index','ZeeMass']

commonBranches = ['run','lumis','event',
                  'HLT','HLTIndex','HLTprescale','nGoodVtx','IsTracksGood',
                  'rho','nPU',
                  'mcIndex','mcPt','mcEta','mcPhi','mcMass',
                  'nMC','mcPID','mcGMomPID',
                  'nPho','phoEt','phoEta','phoPhi','phoSCEta',
                  'phoHoverE','phohasPixelSeed','phoTrkIsoHollowDR04',
                  'phoEcalIsoDR04','phoHcalIsoDR04',
                  'phoSigmaIEtaIEta','phoSigmaIPhiIPhi','phoR9',
                  'phoGenIndex','phoGenGMomPID','phoGenMomPID',
                  'eleTrg','muTrg']

z_infos = {'electron':{'nZ':'nZee',
                       'ell1':'ZeeLeg1Index',
                       'ell2':'ZeeLeg2Index',
                       'mass':'ZeeMass'},
           'muon':{'nZ':'nZmumu',
                   'ell1':'ZmumuLeg1Index',
                   'ell2':'ZmumuLeg2Index',
                   'mass':'ZmumuMass'}}

lepton_masses = {'electron':5.11e-4,'muon':1.06e-1} # in GeV
lepton_infos = {'electron':{'nEll':'nEle',
                            'pt':'elePt',
                            'corpt':'eleCorPt',
                            'eta':'eleEta',
                            'phi':'elePhi'},
                'muon':{'nEll':'nMu',
                        'pt':'muPt',
                        'corpt':'muPt',
                        'eta':'muEta',
                        'phi':'muPhi'}}

run_lumi = {'electron':{'A':2252.0,
                        'B':2709.0,
                        'AB':2252.0+2709.0},
            'muon':{'A':2289.9,
                    'B':2709.0,
                    'AB':2289.9+2709.0}}
rng = TRandom3(0)
             
Z_POLE = 91.188
ell1,ell2,pho = TLorentzVector(),TLorentzVector(),TLorentzVector()
thez, thezg = TLorentzVector(),TLorentzVector()
def run_analysis(options,input_file):
    pwd = ROOT.gDirectory.GetPath()
    in_file = TFile.Open(input_file)
    ROOT.gDirectory.cd(pwd)
    leptonType = options.leptonType
    
    tm = tree_manager()
    tree = None
    nEvents_sample = 0    
    specific = None
    treeName = 'Ntuple'
    if not options.allBranches:    
        if leptonType == 'muon':
            #specific = muonBranches+commonBranches
            mmg = in_file.Get('mmg')
            tree = mmg.Get('final').Get(treeName)
            nEvents_sample = mmg.Get('eventCount').GetBinContent(1)
        elif leptonType == 'electron':
            #specific = electronBranches+commonBranches
            mmg = in_file.Get('eeg')
            tree = in_file.Get('mmg').Get('final').Get(treeName)
            nEvents_sample = eeg.Get('eventCount').GetBinContent(1)
        else:
            raise Exception('invalid lepton type: %s'%options.leptonType)

    total_events = tree.GetEntriesFast()    
    tick_size = int(total_events//100.0)
    selected_events = []

    tm.importTree(treeName,tree,specific)    

    tm.cloneTree(treeName,'%s_zs'%treeName,specific)
    tm.cloneTree(treeName,'%s_zgs'%treeName,specific)
    tm.cloneTree(treeName,'%s_zgs_nosihih'%treeName,specific)

    #setup process dependent stuff
    cuts = setupCuts(options)
    #lepton_mass = lepton_masses[options.leptonType]
    #z_info = z_infos[options.leptonType]
    #lepton_info = lepton_infos[options.leptonType]

    #setup pu-reweighing
    #pu_weight = pu_weight_nominal 
    #if options.isSummer11:
    #    pu_weigt = pu_weight_summer11

    procWeight = 1.0
    if options.datType != 'data':
        procWeight = (run_lumi[options.leptonType][options.datType]*
                      options.crossSection/nEvents_sample)

    # setup the corrector (this links the appropriate four momenta
    # into a common naming scheme
    correct = setup_corrections(options.runYear   , options.runType,
                                options.leptonType, options.datType,
                                options.leptonCor , options.gamCor,
                                options.vanilla                       )

    ievent = long(0)
    for event in tm[treeName]:
        ievent+=1
        #print ievent,total_events,fmod(ievent/total_events,0.01)
        if not (ievent+1)%tick_size or ievent+1 == total_events:  
            sys.stdout.write('\r%3.0f%s complete! (%i/%i)'%(((ievent+1)/
                                                             tick_size),
                                                            '%',
                                                            ievent+1,
                                                            total_events))
            sys.stdout.flush()

        # setup the common event momentum
        # ell1 = lepton1, ell2 = lepton2
        # gam = photon, Z = dilepton, Zg = Z+photon        
        correct(event)
        
        run_idx = getRunIndex(event.run,options.datType,options.leptonType)
        setattr(event,'procWeight',procWeight)
        setattr(event,'puWeight',1.0)
        if options.datType != 'data':            
            setattr(event,'eventFraction',float(ievent+1)/total_events)
            #event.event/nEvents_sample)
            event.puWeight = pu_weight(event.nPU,options.runType)

        selected_z = []
        #selected_pho_nosihih = []
        selected_pho = []        
        bad_leptons = []

        if options.datType == 'data':
            # kill run 170722
            # kill the obvious pile up combinatorial event
            if ( event.run[0] == 170722 or
                 (event.run[0]   == 166512 and
                  event.lumi[0] == 1678 and
                  event.evt[0] == 1822682238) ):
                continue

        bestLLG = None
        bestZdiff = -1
        for i in range(event.N_PATFinalState):            

            cuts.getCutflow('trigger')(event,i)
            if options.exlProc and not cuts.getCutflow('trigger') :
                continue        
            
            cuts.getCutflow('pho')(event,i)            
            cuts.getCutflow('leptons')(event,i) 
            if ( options.exlProc and
                 not cuts.getCutflow('leptons') and
                 not cuts.getCutflow('pho') ):
                continue

            cuts.getCutflow('z')(event,i)
            if options.exlProc and not cuts.getCutflow('z'):
                continue

            if ( cuts.getCutflow('trigger') and
                 cuts.getCutflow('leptons') and
                 cuts.getCutflow('pho') and
                 cuts.getCutflow('z') ):
                thisZdiff = abs(Z_POLE-event.Z[i].M())
                if( cuts.getCut('mindr')[0](event.ell1[i],
                                            event.ell2[i],
                                            event.gam[i]) and
                    (thisZdiff < bestZdiff or bestLLG is None) ):
                    bestZdiff = thisZdiff
                    bestLLG = i
                    
        
        #event object selection done        
        if options.exlProc and bestLLG is None:
            continue
        
        bestZ = bestLLG
        
        if bestLLG is not None:
            setattr(event,'ell1SF',1.0)
            setattr(event,'ell2SF',1.0)
            #if run_idx != -1:
            #    event.ell1SF = leptonSF_nominal(event.nGoodVtx,
            #                          getattr(event,lepton_info['pt'])[idx1],
            #                          getattr(event,lepton_info['eta'])[idx1],
            #                                    run_idx,
            #                                    options.leptonType)
            #    event.ell2SF = leptonSF_nominal(event.nGoodVtx,
            #                          getattr(event,lepton_info['pt'])[idx2],
            #                          getattr(event,lepton_info['eta'])[idx2],
            #                                    run_idx,
            #                                    options.leptonType)
            setattr(event,'bestZ',bestLLG)
            outTrees.bestZTree(event,tm)
            tm.fillTree('%s_zs'%treeName,{})
        
        #bestPhoNoSihih = None
        #bestPhoNoSihihPt = -1
        #for idxph in selected_pho_nosihih:
        #    pho.SetPtEtaPhiM(event.phoCorEt[idxph],
        #                     event.phoEta[idxph],
        #                     event.phoPhi[idxph],
        #                     0.0)
        #    if bestPhoNoSihih is None or pho.Pt() > bestPhoNoSihihPt:
        #        if ( bestZ is not None and
        #             cuts.getCut('mindr')[0](event.bestZLeg1,
        #                                     event.bestZLeg2,
        #                                     pho) ):
        #            bestPhoNoSihih = idxph
        #            bestPhoNoSihihPt = pho.Pt()
        #        elif( bestZ is None ):
        #            bestPhoSihih = idxph
        #            bestPhoSihihPt = pho.Pt()
        
        bestPho = bestLLG

        if  options.exlProc and bestPho is None: #and bestPhoNoSihih is None:
            continue

        if bestPho is not None:            
            setattr(event,'phoSF',1.0)
            #if run_idx != -1:
            #    event.phoSF = phoSF_nominal(event.nGoodVtx,
            #                                event.phoEt[bestPho],
            #                                event.phoEta[bestPho],
            #                                run_idx)            
            setattr(event,'bestPho',bestPho)
            
            #if bestPhoNoSihih is not None:            
            #    setattr(event,'phoNoSihihSF',event.phoSF)                
            #    setattr(event,'bestPhoNoSihihIdx',bestPho)
            #    setattr(event,'bestPhoNoSihih',pho)
        #elif bestPhoNoSihih is not None:
        #    pho.SetPtEtaPhiM(event.phoCorEt[bestPhoNoSihih],
        #                     event.phoEta[bestPhoNoSihih],
        #                     event.phoPhi[bestPhoNoSihih],
        #                     0.0)
        #    setattr(event,'phoNoSihihSF',1.0)
        #    if run_idx != -1:
        #        event.phoNoSihihSF = phoSF_nominal(event.nGoodVtx,
        #                                         event.phoEt[bestPhoNoSihih],
        #                                        event.phoEta[bestPhoNoSihih],
        #                                           run_idx)
        #    setattr(event,'bestPhoNoSihihIdx',bestPhoNoSihih)
        #    setattr(event,'bestPhoNoSihih',pho)
                            
        if bestLLG is not None:
            #if bestPhoNoSihih is not None:
            #    thezg = event.bestZ + event.bestPhoNoSihih
            #    setattr(event,'bestZGNoSihih',thezg)
            #    outTrees.bestZGTreeNoSihih(event,tm)
            #    tm.fillTree('EventTree_zgs_nosihih',{})
            thezg = event.Zg[bestLLG]            
            selected_events.append((event.run[bestLLG],
                                    event.lumi[bestLLG],
                                    event.evt[bestLLG]))
            setattr(event,'bestZG',thezg)
            outTrees.bestZGTree(event,tm)
            tm.fillTree('%s_zgs'%treeName,{})
    
    #make a nice file name
    nameparts = input_file[input_file.rfind('/')+1:]
    nameparts = nameparts.split('.')
    
    #output selected event numbers to file if needed
    print
    print 'Selected %i events after processing!'%(len(selected_events))

    try:
        os.makedirs(options.prefix)
        print 'Created prefix directory: %s'%options.prefix
    except os.error:
        print 'Prefix directory: %s already exists'%options.prefix
    
    if options.dumpSelectedEvents:
        #sort events
        selected_events.sort(key=lambda event: event[2])
        selected_events.sort(key=lambda event: event[1])
        selected_events.sort(key=lambda event: event[0])
        evf = open(options.prefix +
                   '%s_%s_%s_processed.%s'%('_'.join(nameparts[:-1]),
                                            options.runType,
                                            options.leptonType,
                                            'txt'),
                   'w')
        for ev in selected_events:
            evf.write('Run, Lumi, Event #:\t%i\t%i\t%i\n'%ev)
        evf.close()
        
    #push all of our output trees to a file
    outf = TFile.Open(options.prefix +
                      '%s_%s_%s_processed.%s'%('_'.join(nameparts[:-1]),
                                               options.runType,
                                               options.leptonType,
                                               nameparts[-1]),
                      'RECREATE')
    outf.cd()
    tm.write()
    outf.Close()

#determine cutflow given input data type
def setupCuts(options):
    datType = options.datType
    leptonType = options.leptonType
    cuts = CompositeCutflow()
    cuts.addCut('mindr',ell_gamma_dr)
    if leptonType == 'muon':
        if datType == 'data':
            cuts.addCutflow('trigger',muon_triggers_data)
            cuts.addCutflow('leptons',mu_cuts_data)
        else:
            cuts.addCutflow('trigger',muon_triggers_mc)
            cuts.addCutflow('leptons',mu_cuts_mc)
        cuts.addCutflow('z',mumu_cuts)
    elif leptonType == 'electron':
        if datType == 'data':
            cuts.addCutflow('trigger',electron_triggers_data)
            cuts.addCutflow('leptons',e_cuts_data)
        else:
            cuts.addCutflow('trigger',electron_triggers_mc)
            cuts.addCutflow('leptons',e_cuts_mc)
        cuts.addCutflow('z',ee_cuts)
    else:
        raise Exception('Invalid lepton type! {muon,electron} are valid')

    if datType == 'data':
        cuts.addCutflow('pho',photon_cuts_data)    
    elif options.vetoIFSR:
        cuts.addCutflow('pho',photon_cuts_noifsr)
    else:
        cuts.addCutflow('pho',photon_cuts_mc)
    
    return cuts

#this function returns a run period 
run_prob = {'electron':run_lumi['electron']['A']/(run_lumi['electron']['A'] +
                                                  run_lumi['electron']['B']),
            'muon':run_lumi['muon']['A']/(run_lumi['muon']['A'] +
                                          run_lumi['muon']['B'])}
def getRunIndex(run,runType,leptonType):
    if runType == 'data':
        return -1
    elif runType == 'A':
        return 0
    elif runType == 'B':
        return 1
    elif runType == 'AB':
        return int(rng.Rndm() > run_prob[leptonType])


parser = OptionParser(description='%prog : configurable v\gamma analysis',
                      usage='%prog --runType={data,A,B,AB} --leptonType=muon')

parser.add_option('--runYear',dest='runYear',
                  type='int',help='dataset year')
parser.add_option('--runType',dest='runType',                 
                  type='string',help='run era')
parser.add_option('--datType',dest='datType',                 
                  type='string',help='mc or data run')
parser.add_option('--leptonCor',dest='leptonCor',                 
                  type='string',help='lepton correction name')
parser.add_option('--photonCor',dest='gamCor',                 
                  type='string',help='photon correction name',
                  default='PHOSPHOR')
parser.add_option('--leptonType',dest='leptonType',
                  type='string',help='lepton sel.')
parser.add_option('--exclusiveProcessing',dest='exlProc',
                  action='store_true',default=False,
                  help='stop processing if cut failed')
parser.add_option('--allBranches',dest='allBranches',action='store_true',
                  default=False,help='store all input branches.')
parser.add_option('--crossSection',dest='crossSection',
                  type='float',help='MC process cross section, in pb.')
parser.add_option('--calcCS',dest='calcCS',
                  action='store_true',default=False,help='calculate sigma')
parser.add_option('--vanilla',dest='vanilla',
                  action='store_true',default=False,
                  help='true off corrections')
parser.add_option('--dataInput',dest='cs_data_input',
                  type='string',help='the input real data for CS calc')
parser.add_option('--signalMC',dest='cs_mc_input',
                  type='string',help='the input mc data for CS calc')
parser.add_option('--vetoIFSR',dest='vetoIFSR',default=False,
                  action='store_true',help='activate I/FSR vetos.')
parser.add_option('--dumpSelectedEvents',dest='dumpSelectedEvents',
                  default=False,action='store_true',
                  help='write (run,event,lumi) of selected events to file.')
parser.add_option('--prefix',dest='prefix',
                  default='./',help='The directory we write output to.')


(options,args) = parser.parse_args()

if options.runYear is None:
    print 'need to specify run year: 2011, 2012'
    exit(1)
else:
    if options.runYear != 2011 and options.runYear != 2012:
        raise Exception
if options.runType is None:
    if options.runYear == 2011:
        print 'need to specify run type: data, AB'
    elif options.runYear == 2012:
        print 'need to specify run type: data, ABCD'
    else:
        raise Exception('Run year != 2011 or 2012')
    exit(1)
    
if options.leptonType is None:
    print 'need to specify lepton type: electron,muon'
    exit(1)

if options.datType != 'data' and options.crossSection is None:
    print 'need to specify MC process cross section!'
    exit(1)

#set correction types if none set
if options.leptonCor is None:
    if options.leptonType == 'electron':
        options.leptonCor = 'CorrSmearedNoReg'
    elif options.leptonType == 'muon':
        options.leptonCor = 'RochCor'

for input_file in args:
    print
    print 'Processing: %s\n\tdatType=%s\n\trunType=%s\n\tleptonType=%s'\
          %(input_file,
            options.datType,
            options.runType,
            options.leptonType)
    
    if not options.vanilla:
        print '\tActive corrections: lepton=%s photon=%s'%(options.leptonCor,
                                                            options.gamCor)
    else:
        print '\tNo lepton or photon corrections in use.'
    
    run_analysis(options,input_file)

