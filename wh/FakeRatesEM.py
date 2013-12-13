'''

Measure fake rates in the E+Mu channel

We measure in QCD (anti-iso mu) and W+jet (iso mu) control regions.

The layout of output is:
    region/denom_tag/var1
    region/denom_tag/var2
    region/denom_tag/num_tag/var1
    region/denom_tag/num_tag/var2

Author: Evan K. Friis, UW

'''

import ROOT
import EMuTree
from FinalStateAnalysis.PlotTools.MegaBase import MegaBase
import baseSelections as selections
import mcCorrectors
import os
from array import array
import FinalStateAnalysis.PlotTools.pytree as pytree
from cutflowtracker import cut_flow_tracker
import math

def inv_mass(pt1,eta1,phi1,pt2,eta2,phi2):
    return math.sqrt(
        2*pt1*pt2*(math.cosh(eta1 - eta2) - math.cos(phi1 - phi2))
    )

cut_flow_step = ['bare', #'trigger',
                 'e_m_SS', 'e_m_DR', #'e_m_Mass', 
                 'mu presel', 'e presel', 'jet presence',
                 'muon veto', 'bjet veto', 'electron veto', 'tau veto',
                 'trigger', 'tag ID ISO', 'tag MT', 'probe MT'
]

REGION = os.environ['REGION'] if ('REGION' in os.environ) else ''
env_lepton = REGION[0]  if REGION else ''
env_region = REGION[1:] if REGION else ''

def control_region(row):
    # Figure out what control region we are in.
    if row.mPFIDTight and row.mRelPFIsoDB < 0.15 and row.mMtToMET > 35 and row.eMtToMET < 35 and row.bjetCSVVetoZHLike == 0:
        return 'wjetsLtLow'
    elif row.mPFIDTight and row.mRelPFIsoDB < 0.15 and row.mMtToMET > 35 and row.eMtToMET < 35 and row.bjetCSVVetoZHLike > 0:
        return 'ttbar'
    elif row.mPFIDTight and row.mRelPFIsoDB > 0.3 and row.type1_pfMetEt < 25:
        return 'qcd'
    else:
        return None

def control_region_muon(row):
    # Figure out what control region we are in.
    if selections.lepton_id_iso(row, 'e', 'eid12Medium_h2taucuts') and row.eMtToMET > 35 and row.mMtToMET < 35 and row.bjetCSVVetoZHLike == 0:
        return 'MwjetsLtLow'
    elif selections.lepton_id_iso(row, 'e', 'eid12Medium_h2taucuts') and row.eMtToMET > 35 and row.mMtToMET < 35 and row.bjetCSVVetoZHLike > 0:
        return 'Mttbar'
    elif row.eRelPFIsoDB > 0.3 and row.type1_pfMetEt < 25: # selections.electronIds['eid12Medium'](row, 'e') and
        return 'Mqcd'
    else:
        return None


class FakeRatesEM(MegaBase):
    tree = 'em/final/Ntuple' if 'SYNC' not in os.environ else 'Ntuple'
    def __init__(self, tree, outfile, **kwargs):
        super(FakeRatesEM, self).__init__(tree, outfile, **kwargs)
        # Use the cython wrapper
        self.tree = tree #EMuTree.EMuTree(tree)
        self.out = outfile
        # Histograms for each category
        self.histograms = {}
        self.is7TeV = '7TeV' in os.environ['jobid']
        self.pucorrector = mcCorrectors.make_puCorrector('mueg')
        self.defined_eids = selections.electronIds.keys()
        self.iso_points   = ['idiso02', 'h2taucuts', 'h2taucuts020']
        self.lepIds  =  [ '_'.join([i,j])
                          for i in self.defined_eids
                          for j in self.iso_points]
        self.muon_lepIds = ['pfidiso02', 'h2taucuts', 'h2taucuts020'] #, 'h2taucuts025']

    def begin(self):
        #Electron Fake Rates
        self.book('', "CUT_FLOW", "Cut Flow", len(cut_flow_step), 0, len(cut_flow_step))
        xaxis = self.histograms['CUT_FLOW'].GetXaxis()
        self.cut_flow_histo = self.histograms['CUT_FLOW']
        for i, name in enumerate(cut_flow_step):
            xaxis.SetBinLabel(i+1, name)

        for region in ['wjets', 'wjetsLtLow', 'qcd', 'ttbar']:
            for denom in ['pt10','pt20']:
                denom_key = (region, denom)
                denom_histos = {}
                self.histograms[denom_key] = denom_histos

                #SPECIAL NTUPLE!
                denom_histos['electronInfo'] = self.book(
                    os.path.join(region, denom),
                    'electronInfo', "electronInfo", 
                    'electronPt:electronJetPt:electronJetCSVBtag:numJets20:numJets40:weight:electronJetMass:muonJetMass:LT:'+':'.join(self.lepIds), 
                    type=ROOT.TNtuple)

                denom_histos['evtInfo'] = self.book(
                    os.path.join(region, denom),
                    'evtInfo', 'evtInfo',
                    'run/l:lumi/l:evt/l' +\
                    ':ePt/D:eEta/D:ePhi/D:mPt/D:mEta/D:mPhi/D' +\
                    ':eChargeIdTight/D:mPFIDTight/D:mRelPFIsoDB/D' +\
                    ':mMtToMET/D:eMtToMET/D:bjetCSVVetoZHLike/D' +\
                    ':eRelPFIsoDB/D:eMVAIDH2TauWP/D:e_m_SS/D',
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

                    book_histo('ePt', 'e Pt', 100, 0, 100)
                    book_histo('eJetPt', 'e Jet Pt', 100, 0, 100)
                    book_histo('eAbsEta', 'e Abs Eta', 100, -2.5, 2.5)
                    #book_histo('metSignificance', 'MET sig.', 100, 0, 10)
                    book_histo('mMtToMET', 'm MT', 100, 0, 200)
                    book_histo('eJetptD', "", 200, 0, 1)
                    book_histo('eJetmult', "", 50, 0, 50)

                    book_histo('eJetmultvseJetPt', '', 200, 0, 200, 50, 0, 50,type=ROOT.TH2F)
                    book_histo('eJetptDvseJetPt', '' ,  200, 0, 200, 200, 0, 1,type=ROOT.TH2F)

        #Muon Fake Rates
        for region in ['Mwjets', 'MwjetsLtLow', 'Mqcd', 'Mttbar']:
            for denom in ['pt10', 'pt20']:
                denom_key = (region, denom)
                denom_histos = {}
                self.histograms[denom_key] = denom_histos
                
                #SPECIAL NTUPLE!
                denom_histos['muonInfo'] = self.book(
                    os.path.join(region, denom),
                    'muonInfo', "muonInfo", 
                    'muonPt:muonJetPt:muonJetCSVBtag:numJets20:numJets40:weight:electronJetMass:muonJetMass:LT:'+':'.join(self.muon_lepIds), 
                    type=ROOT.TNtuple)
                
                denom_histos['evtInfo'] = self.book(
                    os.path.join(region, denom),
                    'evtInfo', 'evtInfo',
                    'run/l:lumi/l:evt/l' +\
                    ':ePt/D:eEta/D:ePhi/D' +\
                    ':mPt/D:mEta/D:mPhi/D' +\
                    ':eChargeIdTight/D:mPFIDTight/D:mRelPFIsoDB/D' +\
                    ':mMtToMET/D:eMtToMET/D:bjetCSVVetoZHLike/D' +\
                    ':eRelPFIsoDB/D:eMVAIDH2TauWP/D:e_m_SS/D',
                    type=pytree.PyTree)

                for numerator in self.muon_lepIds:
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
                    book_histo('eMtToMET', 'e MT', 100, 0, 200)



    def process(self):
        cut_flow_histo = self.cut_flow_histo
        cut_flow_trk   = cut_flow_tracker(cut_flow_histo)

        def preselection(row):
            if not row.e_m_SS: return False
            cut_flow_trk.Fill('e_m_SS')

            if row.e_m_DR < 0.5:                      return False
            cut_flow_trk.Fill('e_m_DR')

            if not selections.muSelection(row, 'm'):  return False #applies basic selection (eta, pt > 10, DZ, pixHits, jetBTag)
            cut_flow_trk.Fill('mu presel')

            if not selections.eSelection(row, 'e'):   return False #applies basic selection (eta, pt > 10, DZ, missingHits, jetBTag, HasConversion and chargedIdTight)
            if not row.eChargeIdTight:                return False
            cut_flow_trk.Fill('e presel')

            #if not (row.jetVeto40_DR05 >= 1):         return False
            if row.jetVeto20 == 0: return False
            cut_flow_trk.Fill('jet presence')

            if not selections.vetos(row, cut_flow_trk, False):  return False #applies mu bjet e additional tau vetoes
            #DO NOT APPLY BJet veto. It is used for ttbar CR
            return True
        #if self.is7TeV:
            #base_selection = 'mu17ele8Pass && ' + base_selection

        def fill(the_histos, row, fillNtuple=False):
            weight = 1
            if row.run == 1:
                weight = self.pucorrector(row.nTruePU) *\
                         mcCorrectors.correct_mueg_mu(row.mPt, row.mAbsEta) * \
                         mcCorrectors.correct_mueg_e(row.ePt, row.eAbsEta)
            
            the_histos['ePt'].Fill(row.ePt)
            the_histos['eJetPt'].Fill(max(row.eJetPt, row.ePt))
            the_histos['eAbsEta'].Fill(row.eAbsEta)
            #the_histos['metSignificance'].Fill(row.metSignificance)
            the_histos['mMtToMET'].Fill(row.mMtToMET)
            the_histos['eJetptD'].Fill(row.mJetptD)
            the_histos['eJetmult'].Fill(row.mJetmult)
            the_histos['eJetmultvseJetPt'].Fill(max(row.eJetPt, row.ePt),row.eJetmult)
            the_histos['eJetptDvseJetPt'].Fill(max(row.eJetPt, row.ePt),row.eJetptD) 
            if fillNtuple:
                id_iso_vals = [ float( selections.lepton_id_iso(row, 'e', label) ) for label in self.lepIds]
                electron_jet_mass = inv_mass(row.ePt, row.eEta, row.ePhi, row.leadingJetPt, row.leadingJetEta, row.leadingJetPhi)
                muon_jet_mass     = inv_mass(row.mPt, row.mEta, row.mPhi, row.leadingJetPt, row.leadingJetEta, row.leadingJetPhi)
                LT                = row.ePt + row.mPt + row.leadingJetPt
                
                the_histos['electronInfo'].Fill( array("f", [row.ePt, max(row.eJetPt, row.ePt), max(0, row.eJetCSVBtag), 
                                                             row.jetVeto20, row.jetVeto40_DR05, weight, electron_jet_mass,
                                                             muon_jet_mass, LT]+id_iso_vals) )

                the_histos['evtInfo'].Fill(row)

        def fill_muon(the_histos, row, fillNtuple=False):
            weight = 1.
            if row.run == 1:
                weight = self.pucorrector(row.nTruePU) *\
                         mcCorrectors.correct_mueg_mu(row.mPt, row.mAbsEta) * \
                         mcCorrectors.correct_mueg_e(row.ePt, row.eAbsEta)
                        
            the_histos['muonPt'].Fill(row.mPt, weight)
            the_histos['eMtToMET'].Fill(row.eMtToMET)
            the_histos['muonJetPt'].Fill(max(row.mJetPt, row.mPt), weight)
            the_histos['muonAbsEta'].Fill(row.mAbsEta, weight)
            the_histos['muonPtRatio'].Fill(row.mPt/max(row.mJetPt, row.mPt), weight)
            the_histos['muonPtDiff'].Fill(max(row.mJetPt, row.mPt) - row.mPt, weight)
            if fillNtuple:
                pfidiso02    = float( row.mPFIDTight and row.mRelPFIsoDB < 0.2)
                h2taucuts    = float( row.mPFIDTight and ((row.mRelPFIsoDB < 0.15 and row.mAbsEta < 1.479) or row.mRelPFIsoDB < 0.1 ))
                h2taucuts020 = float( row.mPFIDTight and ((row.mRelPFIsoDB < 0.20 and row.mAbsEta < 1.479) or row.mRelPFIsoDB < 0.15))
                electron_jet_mass = inv_mass(row.ePt, row.eEta, row.ePhi, row.leadingJetPt, row.leadingJetEta, row.leadingJetPhi)
                muon_jet_mass     = inv_mass(row.mPt, row.mEta, row.mPhi, row.leadingJetPt, row.leadingJetEta, row.leadingJetPhi)
                LT                = row.ePt + row.mPt + row.leadingJetPt
                
                the_histos['muonInfo'].Fill( array("f", [row.mPt, max(row.mJetPt, row.mPt), max(0, row.mJetCSVBtag), 
                                                         row.jetVeto20, row.jetVeto40_DR05, weight, electron_jet_mass,
                                                         muon_jet_mass, LT, pfidiso02, h2taucuts, h2taucuts020] ) )

                the_histos['evtInfo'].Fill(row)
                
        def fill_region(region,pt_cut):
            if region is None:
                return None
            fill(histos[(region, pt_cut)], row, True)

            for idlabel, idfcn in selections.electronIds.iteritems():                
                if not idfcn(row, 'e'): #if it does not pass id skip!
                    continue

                if row.eRelPFIsoDB < 0.2:
                    fill(histos[(region, pt_cut, idlabel+'_idiso02')], row)

                if (row.eRelPFIsoDB < 0.15 and row.eAbsEta < 1.479) or row.eRelPFIsoDB < 0.1:
                    fill(histos[(region, pt_cut, idlabel+'_h2taucuts')], row)

                if (row.eRelPFIsoDB < 0.20 and row.eAbsEta < 1.479) or row.eRelPFIsoDB < 0.15:
                    fill(histos[(region, pt_cut, idlabel+'_h2taucuts020')], row)

        def fill_muon_region(region, tag):
            if region is None:
                return None
            # This is a QCD or Wjets
            fill_muon(histos[(region, tag)], row, True)
            
            if row.mPFIDTight:
                if row.mRelPFIsoDB < 0.2:
                    fill_muon(histos[(region, tag, 'pfidiso02')], row)
                    
                if (row.mRelPFIsoDB < 0.15 and row.mAbsEta < 1.479) or row.mRelPFIsoDB < 0.1:
                    fill_muon(histos[(region, tag, 'h2taucuts')], row)

                if (row.mRelPFIsoDB < 0.20 and row.mAbsEta < 1.479) or row.mRelPFIsoDB < 0.15:
                    fill_muon(histos[(region, tag, 'h2taucuts020')], row)

                
        histos = self.histograms
        for row in self.tree:
            cut_flow_trk.new_row(row.run,row.lumi,row.evt)
            cut_flow_trk.Fill('bare')
            if not preselection(row):
                continue
            region = control_region(row)
            muon_region = control_region_muon(row)

            is7TeV = bool('7TeV' in os.environ['jobid'])
            use_iso_trigger = not is7TeV
            mu17e8 = (row.mu17ele8isoPass and row.mPt >= 20) if use_iso_trigger else (row.mu17ele8Pass and row.mPt >= 20)
            mu8e17 = (row.mu8ele17isoPass and row.ePt >= 20) if use_iso_trigger else (row.mu8ele17Pass and row.ePt >= 20)

            if env_lepton and env_lepton == 'e':
                if mu17e8:
                    cut_flow_trk.Fill('trigger')
                    if row.mPFIDTight and row.mRelPFIsoDB < 0.15:
                        cut_flow_trk.Fill('tag ID ISO')
                        if row.mMtToMET > 35:
                            cut_flow_trk.Fill('tag MT')
                            if row.eMtToMET < 35:
                                cut_flow_trk.Fill('probe MT')

            elif env_lepton and env_lepton == 'm':
                if mu8e17:
                    cut_flow_trk.Fill('trigger')
                    if selections.lepton_id_iso(row, 'e', 'eid12Medium_h2taucuts'):
                        cut_flow_trk.Fill('tag ID ISO')
                        if row.eMtToMET > 35:
                            cut_flow_trk.Fill('tag MT')
                            if row.mMtToMET < 35:
                                cut_flow_trk.Fill('probe MT')

            cut_flow_trk.flush()

            if mu17e8:
                fill_region(region,'pt10')
                if region == 'wjetsLtLow' and row.mMtToMET > 55:
                    fill_region('wjets', 'pt10')
            if mu8e17:
                fill_region(region,'pt20')
                if region == 'wjetsLtLow' and row.mMtToMET > 55:
                    fill_region('wjets', 'pt20')
                if muon_region:
                    fill_muon_region(muon_region, 'pt10')
                    if muon_region == 'MwjetsLtLow' and row.eMtToMET > 55:
                        fill_muon_region('Mwjets', 'pt10')
                    if row.mPt > 20:
                        fill_muon_region(muon_region, 'pt20')
                        if muon_region == 'MwjetsLtLow' and row.eMtToMET > 55:
                            fill_muon_region('Mwjets', 'pt20')

    def finish(self):
        self.write_histos()
