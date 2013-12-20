'''

Measure fake rates in dimuon events.

We measure in QCD (anti-iso mu) and W+jet (iso mu) control regions.

The layout of output is:

    region/denom_tag/var1
    region/denom_tag/var2
    region/denom_tag/num_tag/var1
    region/denom_tag/num_tag/var2

Author: Evan K. Friis, UW

'''

from array import array
import MuMuTree
from FinalStateAnalysis.PlotTools.MegaBase import MegaBase
import mcCorrectors
import baseSelections as selections
import os
import ROOT
from cutflowtracker import cut_flow_tracker
import FinalStateAnalysis.PlotTools.pytree as pytree
from FinalStateAnalysis.Utilities.struct import struct
import math

def inv_mass(pt1,eta1,phi1,pt2,eta2,phi2):
    return math.sqrt(
        2*pt1*pt2*(math.cosh(eta1 - eta2) - math.cos(phi1 - phi2))
    )

#Makes the cut flow histogram
cut_flow_step = ['bare', 'trigger',
                 'm1_m2_SS', 'm1_m2_DR', 'm1_m2_Mass', 
                 'muon Pt', 'muon AbsEta', 'muon PixHits', 
                 'muon JetCSVBtag', 'muon DZ', #'lead mu preselection', 
                 'lead mu pt', 'm1PFIDTight', 'sublead mu preselection',
                 'muon veto', 'bjet veto', 'electron veto',
                 'Jet presence', 'muon isolation', 'MET', 'region assignment'
]


def control_region(row):
    # Figure out what control region we are in.
    if row.m1RelPFIsoDB < 0.15 and row.m1MtToMET > 35 and row.m2MtToMET < 35:
        return 'wjetsLtLow'
    elif row.m1RelPFIsoDB > 0.3 and row.type1_pfMetEt < 25:
        return 'qcd'
    else:
        return None

region_for_event_list = os.environ.get('EVTLIST_REGION','')
SYNC = ('SYNC' in os.environ) and eval(os.environ['SYNC'])

class FakeRatesMM(MegaBase):
    tree = 'mm/final/Ntuple' if not SYNC else 'Ntuple'
    def __init__(self, tree, outfile, **kwargs):
        super(FakeRatesMM, self).__init__(tree, outfile, **kwargs)
        # Use the cython wrapper
        self.tree = MuMuTree.MuMuTree(tree)
        self.out = outfile
        # Histograms for each category
        self.histograms = {}
        self.is7TeV = '7TeV' in os.environ['jobid']
        self.lepIds = ['pfidiso02', 'h2taucuts', 'h2taucuts020'] #, 'h2taucuts025']
        self.pucorrector = mcCorrectors.make_puCorrector('doublemu')

    def begin(self):
        self.book('', "CUT_FLOW", "Cut Flow", len(cut_flow_step), 0, len(cut_flow_step))
        xaxis = self.histograms['CUT_FLOW'].GetXaxis()
        self.cut_flow_histo = self.histograms['CUT_FLOW']
        for i, name in enumerate(cut_flow_step):
            xaxis.SetBinLabel(i+1, name)

        for region in ['wjets', 'wjetsLtLow', 'qcd']:#, 'all']:
            for denom in ['pt10', 'pt20']:
                denom_key = (region, denom)
                denom_histos = {}
                self.histograms[denom_key] = denom_histos
                #SPECIAL NTUPLE!
                denom_histos['muonInfo'] = self.book(
                    os.path.join(region, denom),
                    'muonInfo', 'muonInfo',
                    'run/l:lumi/l:evt/l' +\
                    ':muonPt/D:muonJetPt/D:numJets20/D:weight/D' +\
                    ':muonAbsEta/D:muonJetCSVBtag/D:muonPVDXY/D:numJets40/D' +\
                    ':tagMuonJetMass/D:probeMuonJetMass/D:LT/D:tagMuonPt/D' +\
                    ':'+':'.join(i+'/D' for i in self.lepIds), 
                    type=pytree.PyTree)
                
                for numerator in self.lepIds:
                    num_key = (region, denom, numerator)
                    num_histos = {}
                    self.histograms[num_key] = num_histos

                    def book_histo(name, *args, **kwargs):
                        # Helper to book a histogram
                        if name not in denom_histos:
                            denom_histos[name] = self.book(os.path.join(
                                region, denom), name, *args, **kwargs)
                        num_histos[name] = self.book(os.path.join(
                            region, denom, numerator), name, *args, **kwargs)

                    mu_binning = array(
                        'd', [10, 12, 15, 20, 30, 50, 100, 200])
                    jet_binning = array(
                        'd', [10, 12, 15, 20, 25, 30, 50, 100, 200])
                    eta_binning = array(
                        'd', [0, 0.25, 0.8, 1.3, 1.8, 2.4])


                    book_histo('muonPt', 'Muon Pt', 200, 0, 200)
                    book_histo('muonJetPt', 'Muon Jet Pt', 200, 0, 200)
                    book_histo('muonAbsEta', 'Muon Abs Eta', 100, -2.5, 2.5)
                    book_histo('muonPtRatio', 'Muon Pt', 100, 0., 1.)
                    book_histo('muonPtDiff', 'Muon Pt', 200, 0., 200.)

                    book_histo('m1MtToMET', 'Muon 1 MT', 100, 0, 200)
                    book_histo('m2MtToMET', 'Muon 2 MT', 100, 0, 200)
                    book_histo('m1m2Mass', 'DiMuon Mass', 100, 0, 200)
                    book_histo("m1JetBtag", "Muon 2 Pt", 100, -10, 3.3)
                    book_histo("m2JetBtag", "Muon 2 Pt", 100, -10, 3.3)

                    book_histo('m2JetptD', "", 200, 0, 1)
                    book_histo('m2Jetmult', "", 50, 0, 50)


    def process(self):

        cut_flow_histo = self.cut_flow_histo
        cut_flow_trk   = cut_flow_tracker(cut_flow_histo)
        cut_flow_trk.Fill('bare')

        def preselection(row, cut_flow_trk):
            double_mu_pass =  row.doubleMuPass and \
                row.m1MatchesDoubleMuPaths > 0 and \
                row.m2MatchesDoubleMuPaths > 0
            double_muTrk_pass = row.doubleMuTrkPass and \
                 row.m1MatchesMu17TrkMu8Path > 0 and \
                 row.m2MatchesMu17TrkMu8Path > 0
            if not ( double_mu_pass or double_muTrk_pass ): return False
            cut_flow_trk.Fill('trigger')

            if not row.m1_m2_SS: return False
            cut_flow_trk.Fill('m1_m2_SS')

            if row.m1_m2_DR < 0.5: return False
            cut_flow_trk.Fill('m1_m2_DR')

            if row.m2Pt > row.m1Pt: return False
            if row.m1_m2_Mass < 20: return False
            cut_flow_trk.Fill('m1_m2_Mass')

            if not selections.muSelection(row, 'm1', cut_flow_trk):  return False
            if not row.m1Pt > 20: return False
            cut_flow_trk.Fill('lead mu pt')
            if not row.m1PFIDTight: return False
            cut_flow_trk.Fill('m1PFIDTight')
            #cut_flow_trk.Fill('lead mu preselection')

            if not selections.muSelection(row, 'm2'):  return False
            cut_flow_trk.Fill('sublead mu preselection')

            if not selections.vetos(row, cut_flow_trk):  return False #applies mu bjet e additional tau vetoes
            #if not (row.jetVeto40_DR05 >= 1): return False
            if row.jetVeto20 == 0: return False
            cut_flow_trk.Fill('Jet presence')

            return True

        def fill(the_histos, row, fillNtuple=False):
            # Get PU weight - fix me
            weight = 1
            if row.run == 1:
                weight = self.pucorrector(row.nTruePU) * \
                         mcCorrectors.double_muon_trigger(row,'m1','m2')
                        
            the_histos['muonPt'].Fill(row.m2Pt, weight)
            the_histos['muonJetPt'].Fill(max(row.m2JetPt, row.m2Pt), weight)
            the_histos['muonAbsEta'].Fill(row.m2AbsEta, weight)
            the_histos['muonPtRatio'].Fill(row.m2Pt/max(row.m2JetPt, row.m2Pt), weight)
            the_histos['muonPtDiff'].Fill(max(row.m2JetPt, row.m2Pt) - row.m2Pt, weight)

            the_histos['m1MtToMET'].Fill(row.m1MtToMET, weight)
            the_histos['m1m2Mass' ].Fill(row.m1_m2_Mass, weight)
            the_histos['m2MtToMET'].Fill(row.m2MtToMET, weight)
            the_histos["m1JetBtag"].Fill(row.m1JetBtag, weight)
            the_histos["m2JetBtag"].Fill(row.m2JetBtag, weight)

            the_histos['m2JetptD'].Fill(row.m2JetptD, weight)
            the_histos['m2Jetmult'].Fill(row.m2Jetmult, weight)

            if fillNtuple:
                pfidiso02    = float( row.m2PFIDTight and row.m2RelPFIsoDB < 0.2)
                h2taucuts    = float( row.m2PFIDTight and ((row.m2RelPFIsoDB < 0.15 and row.m2AbsEta < 1.479) or row.m2RelPFIsoDB < 0.1 ))
                h2taucuts020 = float( row.m2PFIDTight and ((row.m2RelPFIsoDB < 0.20 and row.m2AbsEta < 1.479) or row.m2RelPFIsoDB < 0.15))
                muon1_jet_mass = inv_mass(row.m1Pt, row.m1Eta, row.m1Phi, row.leadingJetPt, row.leadingJetEta, row.leadingJetPhi)
                muon2_jet_mass = inv_mass(row.m2Pt, row.m2Eta, row.m2Phi, row.leadingJetPt, row.leadingJetEta, row.leadingJetPhi)
                LT             = row.m2Pt + row.m1Pt + row.leadingJetPt

                entry = struct(**{
                    'run'              : row.run,
                    'lumi'             : row.lumi,
                    'evt'              : row.evt,
                    'muonPt'           : row.m2Pt,
                    'muonJetPt'        : max(row.m2JetPt, row.m2Pt),
                    'numJets20'        : row.jetVeto20,
                    'weight'           : weight,
                    'muonAbsEta'       : row.m2AbsEta,
                    'muonJetCSVBtag'   : max(0, row.m2JetCSVBtag) ,
                    'muonPVDXY'        : abs(row.m2PVDXY),
                    'numJets40'        : row.jetVeto40_DR05,
                    'LT'               : LT,
                    'tagMuonPt'        : row.m1Pt,
                    'pfidiso02'        : pfidiso02,
                    'h2taucuts'        : h2taucuts,
                    'h2taucuts020'     : h2taucuts020,
                    'tagMuonJetMass'   : muon1_jet_mass,
                    'probeMuonJetMass' : muon2_jet_mass,
                })

                the_histos['muonInfo'].Fill( entry ) 
                #array("f", [row.m2Pt, max(row.m2JetPt, row.m2Pt), row.m2AbsEta, max(0, row.m2JetCSVBtag), 
                #            abs(row.m2PVDXY), row.jetVeto20, row.jetVeto40_DR05, weight, muon1_jet_mass,
                #            muon2_jet_mass, LT, pfidiso02, h2taucuts, h2taucuts020] ) )
        
        def fill_region(region, tag):
            # This is a QCD or Wjets
            fill(histos[(region, tag)], row, True)

            if row.m2PFIDTight:
                if row.m2RelPFIsoDB < 0.2:
                    fill(histos[(region, tag, 'pfidiso02')], row)

                if (row.m2RelPFIsoDB < 0.15 and row.m2AbsEta < 1.479) or row.m2RelPFIsoDB < 0.1:
                    fill(histos[(region, tag, 'h2taucuts')], row)

                if (row.m2RelPFIsoDB < 0.20 and row.m2AbsEta < 1.479) or row.m2RelPFIsoDB < 0.15:
                    fill(histos[(region, tag, 'h2taucuts020')], row)

        histos = self.histograms
        for row in self.tree:
            cut_flow_trk.new_row(row.run,row.lumi,row.evt)
            cut_flow_trk.Fill('bare')
            if not preselection(row, cut_flow_trk):
                continue

            region = control_region(row)
            if row.m1RelPFIsoDB > 0.3:
                cut_flow_trk.Fill('muon isolation')
                if row.type1_pfMetEt < 25:
                    cut_flow_trk.Fill('MET')

            if region is None:
                continue
            cut_flow_trk.Fill('region assignment')

            if region_for_event_list and region == region_for_event_list:
                print '%i:%i:%i' % (row.run, row.lumi, row.evt)
                continue

            fill_region(region, 'pt10')
            if region == 'wjetsLtLow' and row.m1MtToMET > 55:
                fill_region('wjets', 'pt10')

            if row.m2Pt > 20:
                fill_region(region, 'pt20')
                if region == 'wjetsLtLow' and row.m1MtToMET > 55:
                    fill_region('wjets', 'pt20')

    def finish(self):
        self.write_histos()
