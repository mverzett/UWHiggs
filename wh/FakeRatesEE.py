'''

Measure fake rates in dielectron events.

NB that the triggers we are using always have the same iso on the leading
and subleading legs.

We measure W+jet (iso electron) control regions.

The layout of output is:

    region/denom_tag/var1
    region/denom_tag/var2
    region/denom_tag/num_tag/var1
    region/denom_tag/num_tag/var2

Author: Evan K. Friis, UW

'''

import EETree
from FinalStateAnalysis.PlotTools.MegaBase import MegaBase
import baseSelections as selections
from FinalStateAnalysis.PlotTools.decorators import decorator
#import ROOT.TMath as math
import mcCorrectors
import os
from array import array
from pprint import pprint
import ROOT
import FinalStateAnalysis.PlotTools.pytree as pytree
from cutflowtracker import cut_flow_tracker
import math

def inv_mass(pt1,eta1,phi1,pt2,eta2,phi2):
    return math.sqrt(
        2*pt1*pt2*(math.cosh(eta1 - eta2) - math.cos(phi1 - phi2))
    )

cut_flow_step = ['bare', 'trigger', 'e1_e2_DR',
                 'e1 presel', 'e2 presel', 'jet requirement',
                 'muon veto', 'bjet veto', 'electron veto', 'tau veto',
                 'region assignment',
]


region_for_event_list = os.environ.get('EVTLIST_REGION','')
zMassCut              = 'NoZmass' in region_for_event_list
if zMassCut:
    cut_flow_step.append('zMassCut')
region_for_event_list = region_for_event_list.replace('NoZmass','')
SYNC = ('SYNC' in os.environ) and eval(os.environ['SYNC'])
print region_for_event_list

class FakeRatesEE(MegaBase):
    tree = 'ee/final/Ntuple' if not SYNC else 'Ntuple'
    def __init__(self, tree, outfile, **kwargs):
        super(FakeRatesEE, self).__init__(tree, outfile, **kwargs)
        # Use the cython wrapper
        self.tree = EETree.EETree(tree)
        self.out = outfile
        # Histograms for each category
        self.histograms = {}
        self.is7TeV = '7TeV' in os.environ['jobid']
        self.pucorrector = mcCorrectors.make_puCorrector('doublee')
        self.defined_eids = selections.electronIds.keys()
        self.iso_points   = ['idiso02', 'h2taucuts', 'h2taucuts020']
        self.lepIds  =  [ '_'.join([i,j])
                          for i in self.defined_eids
                          for j in self.iso_points]

    def begin(self):
        self.book('', "CUT_FLOW", "Cut Flow", len(cut_flow_step), 0, len(cut_flow_step))
        xaxis = self.histograms['CUT_FLOW'].GetXaxis()
        self.cut_flow_histo = self.histograms['CUT_FLOW']
        for i, name in enumerate(cut_flow_step):
            xaxis.SetBinLabel(i+1, name)

        for region in ['wjetsLtLow', 'qcd', 'wjetsNoZmass', 'wjetsLtLowNoZmass', 'qcdNoZmass']:
            for denom in ['pt10', 'pt20']:
                denom_key = (region, denom)
                denom_histos = {}
                self.histograms[denom_key] = denom_histos

                #SPECIAL NTUPLE!
                denom_histos['electronInfo'] = self.book(
                    os.path.join(region, denom),
                    'electronInfo', "electronInfo", 
                    'electronPt:electronJetPt:electronJetCSVBtag:numJets20:numJets40:weight:tagElectronJetMass:probeElectronJetMass:LT:'+':'.join(self.lepIds), 
                    type=ROOT.TNtuple)

                denom_histos['evtInfo'] = self.book(
                    os.path.join(region, denom),
                    'evtInfo', 'evtInfo',
                    'run/l:lumi/l:evt/l:e1Pt/D:e1Eta/D:e1Phi/D:e2Pt/D:e2Eta/D:e2Phi/D:e1ChargeIdTight/D' +\
                    ':e2ChargeIdTight/D:e2RelPFIsoDB/D:e2MtToMET/D:e1MtToMET/D:bjetCSVVetoZHLike/D'+\
                    ':e1RelPFIsoDB/D:e2MVAIDH2TauWP/D:e1MVAIDH2TauWP/D:e1_e2_SS/D:e1MVANonTrig/D:e2MVANonTrig/D',
                    type=pytree.PyTree)

                for numerator in self.lepIds:
                    num_key = (region, denom, numerator)
                    num_histos = {}
                    self.histograms[num_key] = num_histos

                    def book_histo(name, *args, **kwargs):
                        # Helper to book a histogram
                        if name not in denom_histos:
                            denom_histos[name] = self.book(os.path.join(
                                region, denom), name, *args,**kwargs)
                        num_histos[name] = self.book(os.path.join(
                            region, denom, numerator), name, *args,**kwargs)

                    book_histo('electronPt', 'e Pt', 100, 0, 100)
                    book_histo('electronJetPt', 'e Jet Pt', 100, 0, 100)
                    book_histo('electronAbsEta', 'e Abs Eta', 100, -2.5, 2.5)
                    #book_histo('metSignificance', 'MET sig.', 100, 0, 10)
                    book_histo('e1MtToMET', 'e1 MT', 100, 0, 200)
                    book_histo('e2MtToMET', 'e2 MT', 100, 0, 200)
                    book_histo('MET', 'MET', 100, 0, 400)
                    book_histo('OSS', 'SS to OS', 2, 0, 2)
                    book_histo('e1e2Mass', 'DiElectron Mass', 100, 0, 200)
                    book_histo('doubleEPrescale', 'prescale', 10, -0.5, 9.5)

                    book_histo('e2JetArea'        , "", 100, 0, 1)
                    book_histo('e2JetptD', "", 200, 0, 1)
                    book_histo('e2Jetmult', "", 50, 0, 50)

                    book_histo('e2Jetmultvse2JetPt', '', 200, 0, 200, 50, 0, 50,type=ROOT.TH2F)
                    book_histo('e2JetptDvse2JetPt', '' ,  200, 0, 200, 200, 0, 1,type=ROOT.TH2F)



    def process(self):
        cut_flow_histo = self.cut_flow_histo
        cut_flow_trk   = cut_flow_tracker(cut_flow_histo)

        def preselection(row, cut_flow_trk):
            if not row.doubleEPass: return False
            if not (row.e1MatchesDoubleEPath > 0 and \
                row.e2MatchesDoubleEPath > 0): return False 
            cut_flow_trk.Fill('trigger')

            if row.e1_e2_DR < 0.5: return False
            cut_flow_trk.Fill('e1_e2_DR')

            if not row.e1Pt > 20: return False
            if not selections.eSelection(row, 'e1'): return False
            if not row.e1MVAIDH2TauWP: return False
            cut_flow_trk.Fill('e1 presel')

            if not selections.eSelection(row, 'e2'): return False
            cut_flow_trk.Fill('e2 presel')

            #if not (row.jetVeto40_DR05 >= 1):             return False
            if row.jetVeto20 == 0: return False
            cut_flow_trk.Fill('jet requirement')

            if not selections.vetos(row, cut_flow_trk): return False

            return True

        def fill(the_histos, row, fillNtuple=False):
            weight = 1
            if row.run == 1:
                weight = self.pucorrector(row.nTruePU)

            the_histos['electronPt'].Fill(row.e2Pt)
            the_histos['electronJetPt'].Fill(max(row.e2JetPt, row.e2Pt))
            the_histos['electronAbsEta'].Fill(row.e2AbsEta)
            #the_histos['metSignificance'].Fill(row.metSignificance)
            the_histos['e1MtToMET'].Fill(row.e1MtToMET)
            the_histos['e2MtToMET'].Fill(row.e2MtToMET)
            the_histos['MET'].Fill(row.type1_pfMetEt)
            the_histos['OSS'].Fill(int(row.e1_e2_SS))

            the_histos['e2JetptD'].Fill(row.e2JetptD)
            the_histos['e2Jetmult'].Fill(row.e2Jetmult)

            the_histos['e2Jetmultvse2JetPt'].Fill(max(row.e2JetPt, row.e2Pt),row.e2Jetmult)
            the_histos['e2JetptDvse2JetPt'].Fill(max(row.e2JetPt, row.e2Pt),row.e2JetptD) 
                
            the_histos['e1e2Mass'].Fill(row.e1_e2_Mass)
            the_histos['doubleEPrescale'].Fill(row.doubleEPrescale)

            if fillNtuple:
                id_iso_vals = [  selections.lepton_id_iso(row, 'e2', label)  for label in self.lepIds]
                #print id_iso_vals
                #print self.lepIds
                id_iso_vals = [float( i ) for i in id_iso_vals ]
                electron1_jet_mass = inv_mass(row.e1Pt, row.e1Eta, row.e1Phi, row.leadingJetPt, row.leadingJetEta, row.leadingJetPhi)
                electron2_jet_mass = inv_mass(row.e2Pt, row.e2Eta, row.e2Phi, row.leadingJetPt, row.leadingJetEta, row.leadingJetPhi)
                LT = row.e1Pt + row.e2Pt + row.leadingJetPt
                
                the_histos['electronInfo'].Fill( array("f", [row.e2Pt, max(row.e2Pt, row.e2JetPt), max(0, row.e2JetCSVBtag),
                                                             row.jetVeto20, row.jetVeto40_DR05, weight, electron1_jet_mass,
                                                             electron2_jet_mass, LT]+id_iso_vals) )
                the_histos['evtInfo'].Fill( row )


        histos = self.histograms
        #pprint( histos)
        for row in self.tree:
            cut_flow_trk.new_row(row.run,row.lumi,row.evt)
            cut_flow_trk.Fill('bare')

            if not preselection(row, cut_flow_trk):
                continue

            region = selections.control_region_ee(row)
            if region is None:
                continue

            if region == 'zee':
                continue

            if region_for_event_list and region == region_for_event_list:
                cut_flow_trk.Fill('region assignment')
                if zMassCut and not (60 < row.e1_e2_Mass < 120):
                    cut_flow_trk.Fill('zMassCut')

            # This is a QCD or Wjets
            if region_for_event_list and region == region_for_event_list\
               and (not zMassCut or not (60 < row.e1_e2_Mass < 120)) \
               and not SYNC:
                print '%i:%i:%i' % (row.run, row.lumi, row.evt)
                continue

            def make_region_plots(full_region):
                fill(histos[full_region], row, True)

                for idlabel, idfcn in selections.electronIds.iteritems():                
                    if not idfcn(row, 'e2'): #if it does not pass id skip!
                        continue
                    if row.e2RelPFIsoDB < 0.2:
                        fill(histos[full_region + ( idlabel+'_idiso02',)], row)
                    if (row.e2RelPFIsoDB < 0.15 and row.e2AbsEta < 1.479) or row.e2RelPFIsoDB < 0.1:
                        fill(histos[full_region + ( idlabel+'_h2taucuts',)], row)
                    if (row.e2RelPFIsoDB < 0.20 and row.e2AbsEta < 1.479) or row.e2RelPFIsoDB < 0.15:
                        fill(histos[full_region + ( idlabel+'_h2taucuts020',)], row)


            make_region_plots((region, 'pt10'))
            if not (row.e1_e2_Mass > 60 and row.e1_e2_Mass < 120):
                make_region_plots((region+'NoZmass', 'pt10'))
                if region == 'wjetsLtLow' and row.e1MtToMET > 55:
                    make_region_plots(('wjetsNoZmass', 'pt10'))
                    
            if row.e2Pt > 20:
                make_region_plots((region, 'pt20'))
                if not (row.e1_e2_Mass > 60 and row.e1_e2_Mass < 120):
                    make_region_plots((region+'NoZmass', 'pt20'))
                    if region == 'wjetsLtLow' and row.e1MtToMET > 55:
                        make_region_plots(('wjetsNoZmass', 'pt20'))
            

    def finish(self):
        self.write_histos()

