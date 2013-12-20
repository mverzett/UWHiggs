from FinalStateAnalysis.MetaData.data_views import extract_sample, get_best_style, read_lumi
from progressbar import ETA, ProgressBar, FormatLabel, Bar
import rootpy.plotting as plotting
import ROOT
import logging
import os
import sys
from pdb import set_trace
from glob import glob

logging.basicConfig(stream=sys.stderr, level=logging.INFO)
ROOT.gStyle.SetOptTitle(0)
ROOT.gROOT.SetBatch(True)

def mkdir(obj, name):
    if name not in [i.GetName() for i in obj.GetListOfKeys()]:
        ret = obj.mkdir(name)
        ret.cd()
        return ret
    else:
        ret = obj.Get(name)
        ret.cd()
        return ret


######################################
## Get fake rates
######################################
from fakerate_functions import is7TeV, double_e_lumi, \
    double_mu_lumi, mueg_lumi, wz_lumi, zz_lumi
from FinalStateAnalysis.StatTools.RooFunctorFromWS import MultiFunctorFromMVA

jobid = os.environ['jobid']
frfit_dir = os.path.join('results', os.environ['jobid'], 'NeighborsStudy/fakerates')+'/'
wz_sample = 'WZJetsTo3LNu_ZToTauTau_pythia' if '8TeV' in os.environ['jobid'] else 'WZJetsTo3LNu_ZToTauTau'

mm_m_fr = {}
em_m_fr = {}
variables = ['muonJetPt', 'muonPt', 'numJets20']
for i in range(10, 110, 10):
    key1 = ('h2taucuts020', 'k%i' % i)
    key2 = ('h2taucuts'   , 'k%i' % i)
    mm_m_fr[key1] = {}
    mm_m_fr[key2] = {}
    em_m_fr[key2] = {}
    for parity in ['even', 'odd']:
        mm_m_fr[key1][parity] = MultiFunctorFromMVA(
            '-'.join(key1+(parity,)),
            (  frfit_dir + 'mm_wjets_pt10_h2taucuts020_muonInfo_k%i.%s.data.kNN.weights.xml' % (i, parity)             , double_mu_lumi),
            [ (frfit_dir + 'mm_wjets_pt10_h2taucuts020_muonInfo_k%i.%s.%s.kNN.weights.xml' % (i, parity, wz_sample)    , wz_lumi),
              (frfit_dir + 'mm_wjets_pt10_h2taucuts020_muonInfo_k%i.%s.ZZJetsTo4L_pythia.kNN.weights.xml' % (i, parity), zz_lumi)],
            *variables)

        mm_m_fr[key2][parity] = MultiFunctorFromMVA(
            '-'.join(key2+(parity,)),
            (  frfit_dir + 'mm_wjets_pt10_h2taucuts_muonInfo_k%i.%s.data.kNN.weights.xml' % (i, parity)             , double_mu_lumi),
            [ (frfit_dir + 'mm_wjets_pt10_h2taucuts_muonInfo_k%i.%s.%s.kNN.weights.xml' % (i, parity, wz_sample)    , wz_lumi),
              (frfit_dir + 'mm_wjets_pt10_h2taucuts_muonInfo_k%i.%s.ZZJetsTo4L_pythia.kNN.weights.xml' % (i, parity), zz_lumi)],
            *variables)

        em_m_fr[key2][parity] = MultiFunctorFromMVA(
            '-'.join(key2+(parity,)),
            (  frfit_dir + 'em_Mwjets_pt10_h2taucuts_muonInfo_k%i.%s.data.kNN.weights.xml' % (i, parity)             , double_mu_lumi),
            [ (frfit_dir + 'em_Mwjets_pt10_h2taucuts_muonInfo_k%i.%s.%s.kNN.weights.xml' % (i, parity, wz_sample)    , wz_lumi),
              (frfit_dir + 'em_Mwjets_pt10_h2taucuts_muonInfo_k%i.%s.ZZJetsTo4L_pythia.kNN.weights.xml' % (i, parity), zz_lumi)],
            *variables)

em_e_fr = {}
ee_e_fr = {}
variables = ['electronJetPt', 'electronPt', 'numJets20']
for i in range(10, 110, 10):
    key1 = ('eid12Medium_h2taucuts020', 'k%i' % i)
    key2 = ('eid12Medium_h2taucuts'   , 'k%i' % i)
    key3 = ('eid12Tight_h2taucuts'    , 'k%i' % i)
    em_e_fr = {key2 : {}}
    ee_e_fr = {key3 : {}, key1 : {}}
    for parity in ['even', 'odd']:
        ee_e_fr[key1][parity] = MultiFunctorFromMVA(
            '-'.join(key1+(parity,)),
            (  frfit_dir + 'ee_wjetsNoZmass_pt10_eid12Medium_h2taucuts020_electronInfo_k%i.%s.data.kNN.weights.xml' % (i, parity)             , double_mu_lumi),
            [ (frfit_dir + 'ee_wjetsNoZmass_pt10_eid12Medium_h2taucuts020_electronInfo_k%i.%s.%s.kNN.weights.xml' % (i, parity, wz_sample)    , wz_lumi),
              (frfit_dir + 'ee_wjetsNoZmass_pt10_eid12Medium_h2taucuts020_electronInfo_k%i.%s.ZZJetsTo4L_pythia.kNN.weights.xml' % (i, parity), zz_lumi)],
            *variables)

        ee_e_fr[key3][parity] = MultiFunctorFromMVA(
            '-'.join(key3+(parity,)),
            (  frfit_dir + 'ee_wjetsNoZmass_pt10_eid12Tight_h2taucuts_electronInfo_k%i.%s.data.kNN.weights.xml' % (i, parity)             , double_mu_lumi),
            [ (frfit_dir + 'ee_wjetsNoZmass_pt10_eid12Tight_h2taucuts_electronInfo_k%i.%s.%s.kNN.weights.xml' % (i, parity, wz_sample)    , wz_lumi),
              (frfit_dir + 'ee_wjetsNoZmass_pt10_eid12Tight_h2taucuts_electronInfo_k%i.%s.ZZJetsTo4L_pythia.kNN.weights.xml' % (i, parity), zz_lumi)],
            *variables)

        em_e_fr[key2][parity] = MultiFunctorFromMVA(
            '-'.join(key2+(parity,)),
            (  frfit_dir + 'em_wjets_pt10_eid12Medium_h2taucuts_electronInfo_k%i.%s.data.kNN.weights.xml' % (i, parity)             , double_mu_lumi),
            [ (frfit_dir + 'em_wjets_pt10_eid12Medium_h2taucuts_electronInfo_k%i.%s.%s.kNN.weights.xml' % (i, parity, wz_sample)    , wz_lumi),
              (frfit_dir + 'em_wjets_pt10_eid12Medium_h2taucuts_electronInfo_k%i.%s.ZZJetsTo4L_pythia.kNN.weights.xml' % (i, parity), zz_lumi)],
            *variables)

######################################
## Analyzer class
######################################

class NeighborsAnalyzer(object):
    def __init__(self, fakerates, knn_vars):
        self.fakerates = fakerates
        self.knn_vars  = knn_vars
        self.histos    = {}
        self.plot_vars = []
        
    def book(self, variable, *args):
        self.plot_vars.append((variable, args))

    def load_histos(self, tfile):
        for passing in ['pass', 'fail']:
            pass_dir = mkdir(tfile, passing)
            for cut, neighbors in self.fakerates.iterkeys():
                cut_dir = mkdir(pass_dir, cut)
                neighbors_dir = mkdir(cut_dir, neighbors)
                for var, args in self.plot_vars:
                    self.histos[ '/'.join([passing, cut, neighbors, var]) ] = ROOT.TH1F(var, var, *args)

    def run(self, output_tdir, trees):
        logging.info('Running on %s', output_tdir.GetName())

        knn_vars   = self.knn_vars
        plot_vars  = [i[0] for i in self.plot_vars]
        fakerates  = self.fakerates
        cuts       = [i[0] for i in self.fakerates.keys()]
        cuts       = list( set( cuts ) )

        self.load_histos(output_tdir)
        num_entries = trees['even'].GetEntries() + trees['odd'].GetEntries()
        progress= ProgressBar(
            widgets = [
                ETA(),
                Bar('>')],
            maxval = num_entries+1).start()

        counter = 0
        for parity in ['even', 'odd']:
            for row in trees[parity]:
                counter += 1
                if counter % 10 == 0:
                    progress.update(counter+1)

                #compute kNN
                knn_vars_values = dict( [(i, getattr(row, i)) for i in knn_vars] )
                fake_weights    = [ (key, functor[parity](**knn_vars_values)) for key, functor in fakerates.iteritems() ]
                fake_weights    = [ (key, i / (1 - i))  for key, i in fake_weights ]
                
                #compute cut
                selections      = dict([ (i, (getattr(row, i) > 0.5)) for i in cuts ])
                
                #MC weight
                mc_weight = row.weight
                
                variables = [(i, getattr(row, i)) for i in plot_vars]
                #fill histos
                for key, weight in fake_weights:
                    cut_name, neighbors = key
                    pass_key = 'pass' if selections[cut_name] else 'fail'
                    histo_weight = mc_weight
                    if not selections[cut_name]:
                        histo_weight *= weight
                    for var_name, var_value in variables:
                        self.histos[ '/'.join([pass_key, cut_name, neighbors, var_name]) ].Fill(var_value, histo_weight)

        output_tdir.Write()


######################################
## Analyze!
######################################

output_dir = os.path.join('results', os.environ['jobid'], 'NeighborsStudy/AnalyzedFakeRates')+'/'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

def run_analysis(input_file, analyzer, samples_and_files, treeName='muonInfo'):
    for sample, files in samples_and_files.iteritems():
        logging.info("Running on sample %s...", sample)
        sample_dir = mkdir(input_file, sample)
        even_fname = files['even']
        odd_fname  = files['odd']
        logging.info("Using %s as even \nAnd %s as odd files", even_fname, odd_fname)
        
        f_even = ROOT.TFile(even_fname)
        f_odd  = ROOT.TFile(odd_fname)
        trees = {
            'odd'  : f_even.Get(treeName),
            'even' : f_odd.Get( treeName),
        }
        
        analyzer.run(sample_dir, trees)
        f_even.Close()
        f_odd.Close() 
    return

#
# MM FAKES
#
mm_file   = ROOT.TFile(os.path.join(output_dir,'MM.root'), 'recreate')
logging.info("Saving MM channel into %s", mm_file.GetName())

variables = ['muonJetPt', 'muonPt', 'numJets20']
analyzer  = NeighborsAnalyzer(mm_m_fr, variables)

analyzer.book('tagMuonJetMass'  , 200, 0, 200)
analyzer.book('muonPt'     , 200, 0, 200)
analyzer.book('muonJetPt'  , 200, 0, 200)
analyzer.book('probeMuonJetMass', 200, 0, 200)
analyzer.book('LT'              , 100, 0, 500)

samples = {
    'data' : {
        'even' : os.path.join(frfit_dir, 'mm_wjets_pt10_muonInfo_data.even_tree.root'),
        'odd'  : os.path.join(frfit_dir, 'mm_wjets_pt10_muonInfo_data.odd_tree.root') ,
        },
    'wz' : {
        'even' : glob(os.path.join(frfit_dir, 'mm_wjets_pt10*%s.even_tree.root' % wz_sample))[0],
        'odd'  : glob(os.path.join(frfit_dir, 'mm_wjets_pt10*%s.odd_tree.root'  % wz_sample))[0],
        },
    'zz' : {
        'even' : glob(os.path.join(frfit_dir, 'mm_wjets_pt10*ZZ*.even_tree.root'))[0],
        'odd'  : glob(os.path.join(frfit_dir, 'mm_wjets_pt10*ZZ*.odd_tree.root' ))[0],
        },
}

run_analysis(mm_file, analyzer, samples)
mm_file.Close()

#
# EM FAKES
#
em_file   = ROOT.TFile(os.path.join(output_dir,'EM_Muon.root'), 'recreate')
logging.info("Saving EM, muon channel into %s", em_file.GetName())

variables = ['muonJetPt', 'muonPt', 'numJets20']
analyzer  = NeighborsAnalyzer(em_m_fr, variables)

analyzer.book('muonJetMass'    , 200, 0, 200)
analyzer.book('muonPt'         , 200, 0, 200)
analyzer.book('muonJetPt'      , 200, 0, 200)
analyzer.book('electronJetMass', 200, 0, 200)
analyzer.book('LT'             , 100, 0, 500)

samples = {
    'data' : {
        'even' : os.path.join(frfit_dir, 'em_Mwjets_pt10_muonInfo_data.even_tree.root'),
        'odd'  : os.path.join(frfit_dir, 'em_Mwjets_pt10_muonInfo_data.odd_tree.root' ) ,
        },
    'wz' : {
        'even' : glob(os.path.join(frfit_dir, 'em_Mwjets_pt10_muonInfo_%s.even_tree.root' % wz_sample))[0],
        'odd'  : glob(os.path.join(frfit_dir, 'em_Mwjets_pt10_muonInfo_%s.odd_tree.root'  % wz_sample))[0],
        },
    'zz' : {
        'even' : glob(os.path.join(frfit_dir, 'em_Mwjets_pt10_muonInfo_ZZ*.even_tree.root'))[0],
        'odd'  : glob(os.path.join(frfit_dir, 'em_Mwjets_pt10_muonInfo_ZZ*.odd_tree.root' ))[0],
        },
}

run_analysis(em_file, analyzer, samples)
em_file.Close()



em_file   = ROOT.TFile(os.path.join(output_dir,'EM_Electron.root'), 'recreate')
logging.info("Saving EM, electron channel into %s", em_file.GetName())

variables = ['electronJetPt', 'electronPt', 'numJets20']
analyzer  = NeighborsAnalyzer(em_e_fr, variables)

analyzer.book('muonJetMass'    , 200, 0, 200)
analyzer.book('electronPt'     , 200, 0, 200)
analyzer.book('electronJetPt'  , 200, 0, 200)
analyzer.book('electronJetMass', 200, 0, 200)
analyzer.book('LT'             , 100, 0, 500)

samples = {
    'data' : {
        'even' : os.path.join(frfit_dir, 'em_wjets_pt10_electronInfo_data.even_tree.root'),
        'odd'  : os.path.join(frfit_dir, 'em_wjets_pt10_electronInfo_data.odd_tree.root' ) ,
        },
    'wz' : {
        'even' : glob(os.path.join(frfit_dir, 'em_wjets_pt10_electronInfo_%s.even_tree.root' % wz_sample))[0],
        'odd'  : glob(os.path.join(frfit_dir, 'em_wjets_pt10_electronInfo_%s.odd_tree.root'  % wz_sample))[0],
        },
    'zz' : {
        'even' : glob(os.path.join(frfit_dir, 'em_wjets_pt10_electronInfo_ZZ*.even_tree.root'))[0],
        'odd'  : glob(os.path.join(frfit_dir, 'em_wjets_pt10_electronInfo_ZZ*.odd_tree.root' ))[0],
        },
}

run_analysis(em_file, analyzer, samples, 'electronInfo')
em_file.Close()

#
# EE FAKES
#
ee_file   = ROOT.TFile(os.path.join(output_dir,'EE.root'), 'recreate')
logging.info("Saving EE, electron channel into %s", ee_file.GetName())

variables = ['electronJetPt', 'electronPt', 'numJets20']
analyzer  = NeighborsAnalyzer(ee_e_fr, variables)

analyzer.book('tagElectronJetMass'  , 200, 0, 200)
analyzer.book('electronPt'     , 200, 0, 200)
analyzer.book('electronJetPt'  , 200, 0, 200)
analyzer.book('probeElectronJetMass', 200, 0, 200)
analyzer.book('LT'              , 100, 0, 500)

samples = {
    'data' : {
        'even' : os.path.join(frfit_dir, 'ee_wjetsNoZmass_pt10_electronInfo_data.even_tree.root'),
        'odd'  : os.path.join(frfit_dir, 'ee_wjetsNoZmass_pt10_electronInfo_data.odd_tree.root' ) ,
        },
    'wz' : {
        'even' : glob(os.path.join(frfit_dir, 'ee_wjetsNoZmass_pt10_electronInfo_%s.even_tree.root' % wz_sample))[0],
        'odd'  : glob(os.path.join(frfit_dir, 'ee_wjetsNoZmass_pt10_electronInfo_%s.odd_tree.root'  % wz_sample))[0],
        },
    'zz' : {
        'even' : glob(os.path.join(frfit_dir, 'ee_wjetsNoZmass_pt10_electronInfo_ZZ*.even_tree.root'))[0],
        'odd'  : glob(os.path.join(frfit_dir, 'ee_wjetsNoZmass_pt10_electronInfo_ZZ*.odd_tree.root' ))[0],
        },
}

run_analysis(ee_file, analyzer, samples, 'electronInfo')
ee_file.Close()
