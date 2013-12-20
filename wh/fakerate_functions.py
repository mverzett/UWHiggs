import os
from FinalStateAnalysis.StatTools.RooFunctorFromWS import build_roofunctor, make_corrector_from_th2,\
    build_uncorr_2Droofunctor, FunctorFromMVA, MultiFunctorFromMVA, MultiFunctorFromTF1
from FinalStateAnalysis.StatTools.VariableScaler import make_scaler
from FinalStateAnalysis.Utilities.smartdict import SmartDict
from FinalStateAnalysis.MetaData.data_views import read_lumi
import itertools
import optimizer
import re
import glob
from fnmatch import fnmatch
#from optimizer import leading_lepton_iso_tag, subleading_lepton_iso_tag

################################################################################
#### Fitted fake rate functions ################################################
################################################################################

# Get fitted fake rate functions
frfit_dir = os.path.join('results', os.environ['jobid'], 'fakerate_fits')+'/'
is7TeV = '7TeV' in os.environ['jobid']

def get_lumis(pattern):
    lumi_files  = glob.glob( os.path.join('inputs', os.environ['jobid'], '*.lumicalc.sum') )
    ret = 0.
    #print pattern
    for lumi_file in lumi_files:
        sample_name = os.path.basename(lumi_file).split('.')[0]
        #print lumi_file, sample_name
        if fnmatch(sample_name, pattern):
            #print 'Matches!'
            ret += read_lumi(lumi_file)
    return ret

double_e_lumi  = get_lumis('data_DoubleElectron*')
double_mu_lumi = get_lumis('data_DoubleMu*')
mueg_lumi      = get_lumis('data_MuEG*')
wz_lumi        = get_lumis('WZ*ZToTauTau*')
zz_lumi        = get_lumis('ZZ*')


def make_simple_mapper(tomap):
    def f_(string):
        ret = string
        for key, val in tomap.iteritems():
            ret = ret.replace(key, val)
        return ret
    return f_

def make_regex_mapper(tomap):
    def f_(string):
        ret = string
        for key, val in tomap.iteritems():
            ret = re.sub(key, val, ret)
        return ret
    return f_


def build_roofunctor_dict(basename, wsname = 'fit_efficiency', functor_name = 'efficiency', mapper = None):
    ret = {}
    for i in optimizer.lep_id:
        lepid = i
        if mapper:
            lepid = mapper(lepid)
        ret[i] = build_roofunctor(
            basename % lepid,
            wsname, # workspace name
            functor_name
            )
    return ret

def build_2Droofunctor_dict(basename, 
                            xvar, yvar,
                            wsname = 'fit_efficiency', functor_name = 'efficiency', mapper = None):
    ret = {}
    for i in optimizer.lep_id:
        lepid = i
        if mapper:
            lepid = mapper(lepid)
        ret[i] = build_uncorr_2Droofunctor(
            build_roofunctor(
                basename % (lepid, xvar),
                wsname, # workspace name
                functor_name
            ),
            build_roofunctor(
                basename % (lepid, yvar),
                wsname, # workspace name
                functor_name
            ),
            (basename % (lepid, yvar)).replace('.root','.corrected_inputs.root'),
            )
    return ret


def make_corrector_dict(filename, mapname, mapper=None):
    ret = {}
    for i in optimizer.lep_id:
        lepid = i
        if mapper:
            lepid = mapper(lepid)
        ret[i] = make_corrector_from_th2(filename % lepid, mapname)
    return ret
    
def make_scaler_dict(filename, mapname, **kwargs):
    ret = {}
    for i in optimizer.lep_id:
        ret[i] = make_scaler(filename % i, 'mass_scale', **kwargs)
    return ret

def make_mva_functor_dict(template, variables, mapper=None):
    ret = SmartDict() #faster than normal dict
    for i in optimizer.lep_id:
        lepid = i
        if mapper:
            lepid = mapper(lepid)
        ret.book(i, FunctorFromMVA,
            template % lepid,
            template % lepid,
            *variables
        )
    return ret

def null(*args, **kwargs):
    return 0.

def make_null_dict(template, variables, mapper=None):
    ret = {}
    for i in optimizer.lep_id:
        lepid = i
        if mapper:
            lepid = mapper(lepid)
        ret[i] = null
    return ret

SKIP_MMT=eval( os.environ.get('SKIP_MMT','False') )
SKIP_EMT=eval( os.environ.get('SKIP_EMT','False') )
SKIP_EET=eval( os.environ.get('SKIP_EET','False') )

###################
## MMT FUNCTIONS ##
###################

wz_sample = 'WZJetsTo3LNu_ZToTauTau_pythia' if '8TeV' in os.environ['jobid'] else 'WZJetsTo3LNu_ZToTauTau'

if not SKIP_MMT:
    #no changes in muonID in 2013
    mapper    = {'eid1[0-9][A-Z][a-z]+_':'', 'idiso02' : 'pfidiso02'}
    variables = ['muonJetPt', 'muonPt', 'numJets20'] #, 'muonJetCSVBtag']
    lowpt_mu_fr = SmartDict()
    lowpt_mu_fr.book( 'eid12Medium_h2taucuts020', 
                      MultiFunctorFromMVA,
                      'lowpt_mu_fr',
                      (frfit_dir + 'mm_wjets_pt10_h2taucuts020_muonInfo_k50.data.kNN.weights.xml', double_mu_lumi),
                      [ (frfit_dir + 'mm_wjets_pt10_h2taucuts020_muonInfo_k50.%s.kNN.weights.xml' % wz_sample, wz_lumi),
                        (frfit_dir + 'mm_wjets_pt10_h2taucuts020_muonInfo_k50.ZZJetsTo4L_pythia.kNN.weights.xml', zz_lumi)],
                      *variables)
    lowpt_mu_fr.book( 'eid12Medium_h2taucuts', 
                      MultiFunctorFromMVA,
                      'lowpt_mu_fr',
                      (frfit_dir + 'mm_wjets_pt10_h2taucuts_muonInfo_k50.data.kNN.weights.xml', double_mu_lumi),
                      [ (frfit_dir + 'mm_wjets_pt10_h2taucuts_muonInfo_k50.%s.kNN.weights.xml' % wz_sample, wz_lumi),
                        (frfit_dir + 'mm_wjets_pt10_h2taucuts_muonInfo_k50.ZZJetsTo4L_pythia.kNN.weights.xml', zz_lumi)],
                      *variables)
    
    highpt_mu_fr = SmartDict()
    highpt_mu_fr.book( 'eid12Medium_h2taucuts', 
                      MultiFunctorFromMVA,
                      'highpt_mu_fr',
                      (frfit_dir + 'mm_wjets_pt10_h2taucuts_muonInfo_k50.data.kNN.weights.xml', double_mu_lumi),
                      [ (frfit_dir + 'mm_wjets_pt10_h2taucuts_muonInfo_k50.%s.kNN.weights.xml' % wz_sample, wz_lumi),
                        (frfit_dir + 'mm_wjets_pt10_h2taucuts_muonInfo_k50.ZZJetsTo4L_pythia.kNN.weights.xml', zz_lumi)],
                       *variables,
                       phase_space='muonPt > 20'
    )
    
    lowpt_mu_qcd_fr  = SmartDict() 
    lowpt_mu_qcd_fr.book( 'eid12Medium_h2taucuts020', 
                      MultiFunctorFromMVA,
                      'lowpt_mu_fr',
                      (frfit_dir +   'mm_wjets_pt10_h2taucuts020_muonInfo_k20.data.kNN.weights.xml', double_mu_lumi),
                      [ (frfit_dir + 'mm_wjets_pt10_h2taucuts020_muonInfo_k20.%s.kNN.weights.xml' % wz_sample, wz_lumi),
                        (frfit_dir + 'mm_wjets_pt10_h2taucuts020_muonInfo_k20.ZZJetsTo4L_pythia.kNN.weights.xml', zz_lumi)],
                      *variables)

    highpt_mu_qcd_fr = SmartDict() 
    highpt_mu_qcd_fr.book( 'eid12Medium_h2taucuts', 
                      MultiFunctorFromMVA,
                      'highpt_mu_fr',
                      (frfit_dir +   'mm_wjets_pt10_h2taucuts_muonInfo_k20.data.kNN.weights.xml', double_mu_lumi),
                      [ (frfit_dir + 'mm_wjets_pt10_h2taucuts_muonInfo_k20.%s.kNN.weights.xml' % wz_sample, wz_lumi),
                        (frfit_dir + 'mm_wjets_pt10_h2taucuts_muonInfo_k20.ZZJetsTo4L_pythia.kNN.weights.xml', zz_lumi)],
                       *variables,
                       phase_space='muonPt > 20'
    )



    #lowpt_mu_fr  = make_mva_functor_dict(frfit_dir + '/mm_wjets_pt10_%s_muonInfo_k100.kNN.weights.xml', variables, mapper=make_regex_mapper(mapper))
    #highpt_mu_fr = lowpt_mu_fr
    #
    #lowpt_mu_qcd_fr  = make_mva_functor_dict(frfit_dir + '/mm_qcd_pt10_%s_muonInfo_k100.kNN.weights.xml', variables, mapper=make_regex_mapper(mapper))
    #highpt_mu_qcd_fr = lowpt_mu_qcd_fr
    
###################
## EMT FUNCTIONS ##
###################

if not SKIP_EMT:
    variables = ['muonJetPt', 'muonPt', 'numJets20']

    lowpt_mue_fr = SmartDict()
    lowpt_mue_fr.book( 'eid12Medium_h2taucuts',
                      MultiFunctorFromMVA,
                      'lowpt_mue_fr',
                      (frfit_dir + 'em_Mwjets_pt10_h2taucuts_muonInfo_k50.data.kNN.weights.xml', mueg_lumi),
                      [ (frfit_dir + 'em_Mwjets_pt10_h2taucuts_muonInfo_k50.%s.kNN.weights.xml' % wz_sample, wz_lumi),
                        (frfit_dir + 'em_Mwjets_pt10_h2taucuts_muonInfo_k50.ZZJetsTo4L_pythia.kNN.weights.xml', zz_lumi)],
                      *variables)

    highpt_mue_fr = SmartDict()
    highpt_mue_fr.book( 'eid12Medium_h2taucuts',
                        MultiFunctorFromMVA,
                        'lowpt_mue_fr',
                        (frfit_dir + 'em_Mwjets_pt10_h2taucuts_muonInfo_k50.data.kNN.weights.xml', mueg_lumi),
                        [ (frfit_dir + 'em_Mwjets_pt10_h2taucuts_muonInfo_k50.%s.kNN.weights.xml' % wz_sample, wz_lumi),
                          (frfit_dir + 'em_Mwjets_pt10_h2taucuts_muonInfo_k50.ZZJetsTo4L_pythia.kNN.weights.xml', zz_lumi)],
                        *variables,
                       phase_space='muonPt > 20'
    )


    lowpt_mue_qcd_fr = SmartDict()
    lowpt_mue_qcd_fr.book( 'eid12Medium_h2taucuts',
                           MultiFunctorFromMVA,
                           'lowpt_mue_fr',
                           (frfit_dir +   'em_Mwjets_pt10_h2taucuts_muonInfo_k20.data.kNN.weights.xml', mueg_lumi),
                           [ (frfit_dir + 'em_Mwjets_pt10_h2taucuts_muonInfo_k20.%s.kNN.weights.xml' % wz_sample, wz_lumi),
                             (frfit_dir + 'em_Mwjets_pt10_h2taucuts_muonInfo_k20.ZZJetsTo4L_pythia.kNN.weights.xml', zz_lumi)],
                           *variables)

    highpt_mue_qcd_fr = SmartDict()
    highpt_mue_qcd_fr.book( 'eid12Medium_h2taucuts',
                            MultiFunctorFromMVA,
                            'lowpt_mue_fr',
                            (frfit_dir +   'em_Mwjets_pt10_h2taucuts_muonInfo_k20.data.kNN.weights.xml', mueg_lumi),
                            [ (frfit_dir + 'em_Mwjets_pt10_h2taucuts_muonInfo_k20.%s.kNN.weights.xml' % wz_sample, wz_lumi),
                              (frfit_dir + 'em_Mwjets_pt10_h2taucuts_muonInfo_k20.ZZJetsTo4L_pythia.kNN.weights.xml', zz_lumi)],
                            *variables,
                            phase_space='muonPt > 20'
    )


    variables = ['electronJetPt', 'electronPt', 'numJets20']
    lowpt_e_fr  = SmartDict()
    lowpt_e_fr.book( 'eid12Medium_h2taucuts',
                      MultiFunctorFromMVA,
                      'lowpt_e_fr',
                      (frfit_dir + 'em_wjets_pt10_eid12Medium_h2taucuts_electronInfo_k50.data.kNN.weights.xml', mueg_lumi),
                      [ (frfit_dir + 'em_wjets_pt10_eid12Medium_h2taucuts_electronInfo_k50.%s.kNN.weights.xml' % wz_sample, wz_lumi),
                        (frfit_dir + 'em_wjets_pt10_eid12Medium_h2taucuts_electronInfo_k50.ZZJetsTo4L_pythia.kNN.weights.xml', zz_lumi)],
                      *variables)

    highpt_e_fr = SmartDict()
    highpt_e_fr.book( 'eid12Medium_h2taucuts',
                      MultiFunctorFromMVA,
                      'lowpt_e_fr',
                      (frfit_dir + 'em_wjets_pt10_eid12Medium_h2taucuts_electronInfo_k50.data.kNN.weights.xml', mueg_lumi),
                      [ (frfit_dir + 'em_wjets_pt10_eid12Medium_h2taucuts_electronInfo_k50.%s.kNN.weights.xml' % wz_sample, wz_lumi),
                        (frfit_dir + 'em_wjets_pt10_eid12Medium_h2taucuts_electronInfo_k50.ZZJetsTo4L_pythia.kNN.weights.xml', zz_lumi)],
                      *variables,
                       phase_space='electronPt > 20'
    )
    
    lowpt_e_qcd_fr  = SmartDict()
    lowpt_e_qcd_fr.book( 'eid12Medium_h2taucuts',
                         MultiFunctorFromMVA,
                         'lowpt_e_fr',
                         (frfit_dir +   'em_wjets_pt10_eid12Medium_h2taucuts_electronInfo_k20.data.kNN.weights.xml', mueg_lumi),
                         [ (frfit_dir + 'em_wjets_pt10_eid12Medium_h2taucuts_electronInfo_k20.%s.kNN.weights.xml' % wz_sample, wz_lumi),
                           (frfit_dir + 'em_wjets_pt10_eid12Medium_h2taucuts_electronInfo_k20.ZZJetsTo4L_pythia.kNN.weights.xml', zz_lumi)],
                         *variables)

    highpt_e_qcd_fr = SmartDict()
    highpt_e_qcd_fr.book( 'eid12Medium_h2taucuts',
                          MultiFunctorFromMVA,
                          'lowpt_e_fr',
                          (frfit_dir +   'em_wjets_pt10_eid12Medium_h2taucuts_electronInfo_k20.data.kNN.weights.xml', mueg_lumi),
                          [ (frfit_dir + 'em_wjets_pt10_eid12Medium_h2taucuts_electronInfo_k20.%s.kNN.weights.xml' % wz_sample, wz_lumi),
                            (frfit_dir + 'em_wjets_pt10_eid12Medium_h2taucuts_electronInfo_k20.ZZJetsTo4L_pythia.kNN.weights.xml', zz_lumi)],
                          *variables,
                          phase_space='electronPt > 20'
    )

                       
    #lowpt_mue_fr  = make_mva_functor_dict(frfit_dir + '/em_Mwjets_pt10_%s_muonInfo_k100.kNN.weights.xml', variables, mapper=make_regex_mapper(mapper))
    #highpt_mue_fr = lowpt_mue_fr
    #
    #lowpt_mue_qcd_fr  = make_mva_functor_dict(frfit_dir + '/em_Mqcd_pt10_%s_muonInfo_k100.kNN.weights.xml', variables, mapper=make_regex_mapper(mapper))
    #highpt_mue_qcd_fr = lowpt_mue_qcd_fr
    #
    #variables = ['electronJetPt', 'electronPt', 'numJets20']
    #lowpt_e_fr  = make_mva_functor_dict(frfit_dir + '/em_wjets_pt10_%s_electronInfo_k100.kNN.weights.xml', variables)
    #highpt_e_fr = lowpt_e_fr
    #
    #lowpt_e_qcd_fr  = make_mva_functor_dict(frfit_dir + '/em_qcd_pt10_%s_electronInfo_k100.kNN.weights.xml', variables)
    #highpt_e_qcd_fr = lowpt_e_qcd_fr


###################
## EET FUNCTIONS ##
###################

if not SKIP_EET:
    variables = ['electronJetPt', 'electronPt', 'numJets20']
    lowpt_ee_fr  = SmartDict()
    lowpt_ee_fr.book( 'eid12Medium_h2taucuts020',
                      MultiFunctorFromMVA,
                      'lowpt_ee_fr',
                      (frfit_dir + 'ee_wjetsNoZmass_pt10_eid12Medium_h2taucuts020_electronInfo_k50.data.kNN.weights.xml', double_e_lumi),
                      [ (frfit_dir + 'ee_wjetsNoZmass_pt10_eid12Medium_h2taucuts020_electronInfo_k50.%s.kNN.weights.xml' % wz_sample, wz_lumi),
                        (frfit_dir + 'ee_wjetsNoZmass_pt10_eid12Medium_h2taucuts020_electronInfo_k50.ZZJetsTo4L_pythia.kNN.weights.xml', zz_lumi)],
                      *variables)
    lowpt_ee_fr.book( 'eid12Tight_h2taucuts',
                       MultiFunctorFromMVA,
                       'highpt_ee_fr',
                       (frfit_dir + 'ee_wjetsNoZmass_pt10_eid12Tight_h2taucuts_electronInfo_k50.data.kNN.weights.xml', double_e_lumi),
                       [ (frfit_dir + 'ee_wjetsNoZmass_pt10_eid12Tight_h2taucuts_electronInfo_k50.%s.kNN.weights.xml' % wz_sample, wz_lumi),
                         (frfit_dir + 'ee_wjetsNoZmass_pt10_eid12Tight_h2taucuts_electronInfo_k50.ZZJetsTo4L_pythia.kNN.weights.xml', zz_lumi)],
                       *variables
    )
        
    highpt_ee_fr = SmartDict()
    highpt_ee_fr.book( 'eid12Tight_h2taucuts',
                       MultiFunctorFromMVA,
                       'highpt_ee_fr',
                       (frfit_dir + 'ee_wjetsNoZmass_pt10_eid12Tight_h2taucuts_electronInfo_k50.data.kNN.weights.xml', double_e_lumi),
                       [ (frfit_dir + 'ee_wjetsNoZmass_pt10_eid12Tight_h2taucuts_electronInfo_k50.%s.kNN.weights.xml' % wz_sample, wz_lumi),
                         (frfit_dir + 'ee_wjetsNoZmass_pt10_eid12Tight_h2taucuts_electronInfo_k50.ZZJetsTo4L_pythia.kNN.weights.xml', zz_lumi)],
                       *variables,
                       phase_space='electronPt > 20'
    )


    lowpt_ee_qcd_fr  = SmartDict()
    lowpt_ee_qcd_fr.book( 'eid12Medium_h2taucuts020',
                      MultiFunctorFromMVA,
                      'lowpt_ee_fr',
                      (frfit_dir + 'ee_wjetsNoZmass_pt10_eid12Medium_h2taucuts020_electronInfo_k20.data.kNN.weights.xml', double_e_lumi),
                      [ (frfit_dir + 'ee_wjetsNoZmass_pt10_eid12Medium_h2taucuts020_electronInfo_k20.%s.kNN.weights.xml' % wz_sample, wz_lumi),
                        (frfit_dir + 'ee_wjetsNoZmass_pt10_eid12Medium_h2taucuts020_electronInfo_k20.ZZJetsTo4L_pythia.kNN.weights.xml', zz_lumi)],
                      *variables)

    highpt_ee_qcd_fr = SmartDict()
    highpt_ee_qcd_fr.book( 'eid12Tight_h2taucuts',
                       MultiFunctorFromMVA,
                       'highpt_ee_fr',
                       (frfit_dir + 'ee_wjetsNoZmass_pt10_eid12Tight_h2taucuts_electronInfo_k20.data.kNN.weights.xml', double_e_lumi),
                       [ (frfit_dir + 'ee_wjetsNoZmass_pt10_eid12Tight_h2taucuts_electronInfo_k20.%s.kNN.weights.xml' % wz_sample, wz_lumi),
                         (frfit_dir + 'ee_wjetsNoZmass_pt10_eid12Tight_h2taucuts_electronInfo_k20.ZZJetsTo4L_pythia.kNN.weights.xml', zz_lumi)],
                       *variables,
                       phase_space='electronPt > 20'
    )




    #lowpt_ee_fr  = make_mva_functor_dict(frfit_dir + '/ee_wjetsNoZmass_pt10_%s_electronInfo_k100.kNN.weights.xml', variables) 
    #highpt_ee_fr = lowpt_ee_fr
    #
    #lowpt_ee_qcd_fr  = make_mva_functor_dict(frfit_dir + '/ee_qcd_pt10_%s_electronInfo_k100.kNN.weights.xml', variables) 
    #highpt_ee_qcd_fr = lowpt_ee_qcd_fr


##################
## 1D Taus Func ##
##################

filename = '2011.root' if is7TeV else '2012.root'
date     = '2011' if is7TeV else '2012'
tau_fr = MultiFunctorFromTF1(
    'external/%s' % filename,
    [('%s_wmuJetsAll_electronLooseMVA3_muonTight_lepTauOS_fakeable_combinedIsolationLoose3Hits_byEta_projX_barrel/landau'     % date, [-1, 0.8]),
     ('%s_wmuJetsAll_electronLooseMVA3_muonTight_lepTauOS_fakeable_combinedIsolationLoose3Hits_byEta_projX_transition/landau' % date, [0.8, 1.6]),
     ('%s_wmuJetsAll_electronLooseMVA3_muonTight_lepTauOS_fakeable_combinedIsolationLoose3Hits_byEta_projX_forward/landau'    % date, [1.6, 100])]
)

tau_qcd_fr = tau_fr 

############################
## CHARGE FLIP FUNCTIONS  ##
############################

if not (SKIP_EMT or SKIP_EET):
    e_charge_flip       = make_corrector_dict(frfit_dir+"/charge_flip_prob_map_%s.root", "efficiency_map")         
    e_charge_flip_up    = make_corrector_dict(frfit_dir+"/charge_flip_prob_map_%s.root", "efficiency_map_statUp")  
    e_charge_flip_down  = make_corrector_dict(frfit_dir+"/charge_flip_prob_map_%s.root", "efficiency_map_statDown")
    mass_scaler         = make_scaler_dict(frfit_dir+"/charge_flip_prob_map_%s.root", 'mass_scale')
    default_scaler      = mass_scaler[mass_scaler.keys()[0]]
    data_scaler         = None
    if os.path.isfile(frfit_dir+"/charge_flip_prob_map_data.root"):
        data_scaler = make_scaler(frfit_dir+"/charge_flip_prob_map_data.root", 'mass_scale')
