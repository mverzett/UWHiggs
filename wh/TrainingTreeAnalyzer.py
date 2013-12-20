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

class TrainingTreeAnalyzer(object):
    def __init__(self, output_dir, *fakerates):
        self.fakerates = fakerates
        self.knn_vars  = fakerates[0].functors_and_weights[0][0].variables

        xml_name  = fakerates[0].functors_and_weights[0][0].xml_filename
        xml_name  = os.path.basename( xml_name )
        xml_split = xml_name.split('_')
        region    = xml_split[1]
        denom     = xml_split[2]
        var       = xml_split[5] \
                    if len(xml_split) == 7 else \
                    xml_split[4]

        self.treepath  = '/'.join([region, denom, var])
        self.cuts = []
        self.neighbors = []
        for fr in fakerates:
            xml_name  = fr.functors_and_weights[0][0].xml_filename
            xml_name  = os.path.basename( xml_name )
            xml_split = xml_name.split('_')
            neigh     = xml_split[-1].split('.')[0]
            cut       = '_'.join(xml_split[3:5]) \
                        if len(xml_split) == 7 else \
                        xml_split[3]
            self.cuts.append(cut)
            self.neighbors.append(neigh)

        self.plot_vars = []
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)
        self.output_dir = output_dir
        self.histos     = {}
        self.plot_vars  = []
        
    def book(self, variable, *args):
        self.plot_vars.append((variable, args))

    @staticmethod
    def mkdir(obj, name):
        if name not in [i.GetName() for i in obj.GetListOfKeys()]:
            return obj.mkdir(name)
        else:
            return obj.Get(name)

    def load_histos(self, tfile):
        for n, cut in zip(self.neighbors, self.cuts):
            neighbors_dir = self.mkdir(tfile, n)
            for passing in ['pass', 'fail']:
                pass_dir = self.mkdir(neighbors_dir, passing)
                pass_dir.cd()

                cut_dir = self.mkdir(pass_dir, cut)
                cut_dir.cd()
                for var, args in self.plot_vars:
                    #print 'loading %s' % os.path.join(neighbors_dir.GetName(), pass_dir.GetName(), cut_dir.GetName(), var)
                    self.histos[ '/'.join([n, passing, cut, var]) ] = ROOT.TH1F(var, var, *args)
        #set_trace()

    def run(self, file_name):
        logging.info('Running on %s', file_name)

        knn_vars   = self.knn_vars
        plot_vars  = [i[0] for i in self.plot_vars]
        cuts       = self.cuts
        fakerates  = self.fakerates
        neightbors = self.neighbors

        tfile = ROOT.TFile.Open(file_name)
        ttree = tfile.Get(self.treepath)
        
        outfile_path = os.path.join(
            self.output_dir, 
            os.path.basename(file_name)
        )

        outfile = ROOT.TFile( outfile_path, 'recreate')
        self.load_histos(outfile)

        progress= ProgressBar(
            widgets = [
                ETA(),
                Bar('>')],
            maxval = ttree.GetEntries() ).start()
        
        for entry, row in enumerate(ttree):
            if entry % 10 == 0:
                progress.update(entry+1)

            #compute kNN
            knn_vars_values = dict( [(i, getattr(row, i)) for i in knn_vars] )
            fake_weights    = [ i(**knn_vars_values) for i in fakerates ]
            fake_weights    = [ i / (1 - i)  for i in fake_weights ]

            #compute cut
            selections      = [ (getattr(row, i) > 0.5) for i in cuts ]
            
            #MC weight
            mc_weight = row.weight

            variables = [(i, getattr(row, i)) for i in plot_vars]
            #fill histos
            for neigh, cut, value, knn in zip(neightbors, cuts, selections, fake_weights):
                pass_key = 'pass' if value else 'fail'
                histo_weight = mc_weight
                if not value:
                    histo_weight *= knn
                for var_name, var_value in variables:
                    self.histos[ '/'.join([neigh, pass_key, cut, var_name]) ].Fill(var_value, histo_weight)

        ## for i in self.histos.itervalues():
        ##     i.Write()
        outfile.Write()
        outfile.Close()
        tfile.Close()

import fakerate_functions as frfits
resultdir = 'results/%s/' % os.environ['jobid']

#
# MM FAKES
#
mm_fakes = TrainingTreeAnalyzer(
    os.path.join(resultdir, 'FakeRatesAnalyzerMM'),
    *(frfits.lowpt_mu_fr.values() + frfits.lowpt_mu_qcd_fr.values())
)

mm_fakes.book('tagMuonJetMass'  , 200, 0, 200)
mm_fakes.book('muonPt'     , 200, 0, 200)
mm_fakes.book('muonJetPt'  , 200, 0, 200)
mm_fakes.book('probeMuonJetMass', 200, 0, 200)
mm_fakes.book('LT'              , 100, 0, 500)

for infile in glob( os.path.join(resultdir, 'FakeRatesMM', '*.root')):
    mm_fakes.run(infile)


#
# EM FAKES
#
em_fakes = TrainingTreeAnalyzer(
    os.path.join(resultdir, 'FakeRatesAnalyzerEM_electron'),
    *(frfits.lowpt_e_fr.values() + frfits.lowpt_e_qcd_fr.values())
)

em_fakes.book('electronJetMass', 200, 0, 200)
em_fakes.book('muonJetMass'    , 200, 0, 200)
em_fakes.book('LT'             , 100, 0, 500)
em_fakes.book('electronPt'     , 200, 0, 200)
em_fakes.book('electronJetPt'  , 200, 0, 200)

for infile in glob( os.path.join(resultdir, 'FakeRatesEM', '*.root')):
    em_fakes.run(infile)


em_fakes = TrainingTreeAnalyzer(
    os.path.join(resultdir, 'FakeRatesAnalyzerEM_muon'),
    *(frfits.lowpt_mue_fr.values() + frfits.lowpt_mue_qcd_fr.values())
)

em_fakes.book('electronJetMass', 200, 0, 200)
em_fakes.book('muonJetMass'    , 200, 0, 200)
em_fakes.book('LT'             , 100, 0, 500)
em_fakes.book('muonPt'     , 200, 0, 200)
em_fakes.book('muonJetPt'  , 200, 0, 200)

for infile in glob( os.path.join(resultdir, 'FakeRatesEM', '*.root')):
    em_fakes.run(infile)

#
# EE FAKES
#
ee_fakes = TrainingTreeAnalyzer(
    os.path.join(resultdir, 'FakeRatesAnalyzerEE'),
    *(frfits.lowpt_ee_fr.values() + frfits.lowpt_ee_qcd_fr.values())
)

ee_fakes.book('tagElectronJetMass'  , 200, 0, 200)
ee_fakes.book('probeElectronJetMass', 200, 0, 200)
ee_fakes.book('LT'                  , 100, 0, 500)
ee_fakes.book('electronPt'     , 200, 0, 200)
ee_fakes.book('electronJetPt'  , 200, 0, 200)

for infile in glob( os.path.join(resultdir, 'FakeRatesEE', '*.root')):
    ee_fakes.run(infile)
