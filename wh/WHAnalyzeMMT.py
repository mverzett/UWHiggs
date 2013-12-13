'''

Analyze MMT events for the WH analysis

'''

from MuMuTauTree import MuMuTauTree
import ROOT
import os
from WHAnalyzerBase import WHAnalyzerBase, quad, inv_mass
import mcCorrectors
import baseSelections as selections
import fakerate_functions as frfits
import math
from array import array
from FinalStateAnalysis.PlotTools.decorators import memo_last
import optimizer
import FinalStateAnalysis.PlotTools.pytree as pytree

#initialize FRFits
optimizer_keys   = [ i for i in optimizer.grid_search.keys() if i.startswith('MMT') ]
print optimizer_keys
grid_search = {}
if len(optimizer_keys) > 1:
        grid_search[key] = optimizer.grid_search[key]

NO_HISTO_FILL = ('NO_HISTO_FILL' in os.environ) and eval(os.environ['NO_HISTO_FILL']) 

if not NO_HISTO_FILL:
    for key in optimizer_keys:
        frfits.highpt_mu_fr.__activate__(optimizer.grid_search[key]['leading_iso'])
        frfits.highpt_mu_qcd_fr.__activate__(optimizer.grid_search[key]['leading_iso'])
        frfits.lowpt_mu_fr.__activate__(optimizer.grid_search[key]['subleading_iso'])
        frfits.lowpt_mu_qcd_fr.__activate__(optimizer.grid_search[key]['subleading_iso'])

SYNC = ('SYNC' in os.environ) and eval(os.environ['SYNC'])

################################################################################
#### Analysis logic ############################################################
################################################################################

class WHAnalyzeMMT(WHAnalyzerBase):
    tree = 'mmt/final/Ntuple' if not SYNC else 'Ntuple'
    def __init__(self, tree, outfile, **kwargs):
        self.channel = 'MMT'
        super(WHAnalyzeMMT, self).__init__(tree, outfile, MuMuTauTree, **kwargs)

	def attr_getter(attribute):
            def f(row, weight):
                return (getattr(row,attribute), weight)
            return f

        def merge_functions(fcn_1, fcn_2):
            def f(row, weight):
                r1, w1 = fcn_1(row, weight)
                r2, w2 = fcn_2(row, weight)
                w = w1 if w1 and w2 else None
                return ((r1, r2), w)
            return f

        lead_iso = self.grid_search['']['leading_iso']
        sublead_iso = self.grid_search['']['subleading_iso']
        
        #@memo_last
        #def f_par_prob(m1Pt, m1JetPt,
        #               m2Pt, m2JetPt,
        #               tPt):
        #    p_m1 = (( frfits.highpt_mu_fr[lead_iso](muonJetPt=max(m1JetPt, m1Pt), muonPt=m1Pt) +\
        #              frfits.highpt_mu_qcd_fr[lead_iso](muonJetPt=max(m1JetPt, m1Pt), muonPt=m1Pt) )/2)
        #    p_m2 = (( frfits.lowpt_mu_fr[sublead_iso](muonJetPt=max(m2JetPt, m2Pt), muonPt=m2Pt) + \
        #              frfits.lowpt_mu_qcd_fr[sublead_iso](muonJetPt=max(m2JetPt, m2Pt), muonPt=m2Pt))/2)
        #    p_t  = frfits.tau_fr(tPt)
        #    return (p_m1 + p_m2*(1 - p_m1) + p_t*(1 - p_m1)*(1 - p_m2))
        #    
        #def f_prob(row, weight):
        #    val = f_par_prob(row.m1Pt, row.m1JetPt,
        #                     row.m2Pt, row.m2JetPt,
        #                     row.tPt)
        #    return val, weight
        #
        #def log_prob(row, weight):
        #    prob, weight = f_prob(row, weight)
        #    return ROOT.TMath.Log10(prob), weight
        
        #self.hfunc['faking_prob'] = f_prob
        #self.hfunc['log_prob']    = log_prob
        #self.hfunc["m2_t_Mass#faking_prob"] = merge_functions( attr_getter('m2_t_Mass'), f_prob  )
        #self.hfunc["m2_t_Mass#log_prob"   ] = merge_functions( attr_getter('m2_t_Mass'), log_prob)

        self.hfunc['subMTMass'] = lambda row, weight: (row.m2_t_Mass, weight) if row.m1MtToMET > row.m2MtToMET else (row.m1_t_Mass, weight) #maps the name of non-trivial histograms to a function to get the proper value, the function MUST have two args (evt and weight). Used in WHAnalyzerBase.fill_histos later
        self.hfunc['pt_ratio' ] = lambda row, weight: (row.m2Pt/row.m1Pt, weight)
        self.hfunc['SYNC'] = lambda row, weight: (row, None)#((row.run, row.lumi, row.evt, row.m1Pt, row.m1Eta, row.m1Phi, row.m2Pt, row.m2Eta, row.m2Phi, row.tPt, row.tEta, row.tPhi, weight), None)
        
        self.pucorrector = mcCorrectors.make_puCorrector('doublemu')

    def book_histos(self, folder):
        LTBinning = array('d',[0] + range(50,210,10) + [600])
        nLTBins   = len(LTBinning) -1
        for key in self.grid_search:
            prefix = key+'$' if key else ''
            #self.book(folder, prefix+"m2_t_Mass", "subleadingMass", 200, 0, 200)
            self.book(folder, prefix+"LT" ,  "LT" , 100, 0., 500)
            self.book(folder, prefix+"m2_t_Pt", "subleadingPt", 400, 0, 400)

        if len(self.grid_search.keys()) == 1:

            self.book(folder, 'SYNC', 'SYNC', 
                      'run/l:lumi/l:evt/l' +\
                      ':m1Pt/D:m1Eta/D:m1Phi/D:m1PFIDTight/D:m1RelPFIsoDB/D' +\
                      ':m2Pt/D:m2Eta/D:m2Phi/D:m2PFIDTight/D:m2RelPFIsoDB/D' +\
                      ':tPt/D:tEta/D:tPhi/D:tAntiMuonTight/D:tDecayMode/D:tLooseIso3Hits/D:tAntiElectronMVA3Loose/D' +\
                      ':m1_m2_Mass/D:m2_t_Mass/D:LT/D:m1_m2_SS/D:m1_t_SS/D',
                      type=pytree.PyTree)

            self.book(folder, "m2_t_Mass#LT" , "subleadingMass", 300, 0, 300, nLTBins, LTBinning, type=ROOT.TH2F)
            #self.book(folder, "Event_ID", "Event ID", 'run:lumi:evt1:evt2', type=ROOT.TNtuple)
            #self.book(folder, "DEBUG", "DEBUG", 'run:lumi:m2_t_Mass:LT:m1Pt:m2Pt:tPt', type=ROOT.TNtuple)
            #Pt
            self.book(folder, "m1Pt#LT" , "subleadingMass", 150, 0, 150, nLTBins, LTBinning, type=ROOT.TH2F)
            self.book(folder, "tPt#LT"  , "subleadingMass", 150, 0, 150, nLTBins, LTBinning, type=ROOT.TH2F)
            self.book(folder, "m2Pt#LT" , "subleadingMass", 150, 0, 150, nLTBins, LTBinning, type=ROOT.TH2F)
            self.book(folder, "m1JetPt#LT" , "subleadingMass", 150, 0, 150, nLTBins, LTBinning, type=ROOT.TH2F)
            self.book(folder, "m2JetPt#LT" , "subleadingMass", 150, 0, 150, nLTBins, LTBinning, type=ROOT.TH2F)

            #eta
            self.book(folder, "m1AbsEta#LT" , "subleadingMass", 100, 0, 2.5, nLTBins, LTBinning, type=ROOT.TH2F)
            self.book(folder, "m2AbsEta#LT" , "subleadingMass", 100, 0, 2.5, nLTBins, LTBinning, type=ROOT.TH2F)
            self.book(folder, "tAbsEta#LT"  , "subleadingMass", 100, 0, 2.5, nLTBins, LTBinning, type=ROOT.TH2F)

            #DR
            self.book(folder, "m1_t_DR#LT" , "subleadingMass", 100, 0, 10, nLTBins, LTBinning, type=ROOT.TH2F)
            self.book(folder, "m2_t_DR#LT" , "subleadingMass", 100, 0, 10, nLTBins, LTBinning, type=ROOT.TH2F)
	    self.book(folder, "m1_m2_DR#LT", "subleadingMass", 100, 0, 10, nLTBins, LTBinning, type=ROOT.TH2F)
	   
            #Jet BTag
            self.book(folder, "m2JetBtag#LT", "Muon 2 Pt", 100, -100, 100, nLTBins, LTBinning, type=ROOT.TH2F)
            self.book(folder, "m1JetBtag#LT", "Muon 2 Pt", 100, -100, 100, nLTBins, LTBinning, type=ROOT.TH2F)
            self.book(folder, "m1JetCSVBtag#LT", "Muon 2 Pt", 120, -5, 1 , nLTBins, LTBinning, type=ROOT.TH2F)
            self.book(folder, "m2JetCSVBtag#LT", "Muon 2 Pt", 120, -5, 1 , nLTBins, LTBinning, type=ROOT.TH2F)

            #Mt To MET
            self.book(folder, "m1MtToMET#LT", "Muon 2 Pt", 150, 0, 150, nLTBins, LTBinning, type=ROOT.TH2F)
            self.book(folder, "m2MtToMET#LT", "Muon 2 Pt", 150, 0, 150, nLTBins, LTBinning, type=ROOT.TH2F)

            #tau ISO 
            self.book(folder, "tRawIso3Hits#LT", "Muon 2 Pt", 200, 0, 200, nLTBins, LTBinning, type=ROOT.TH2F)
            self.book(folder, "tRawIsoMVA2#LT", "Muon 2 Pt" , 200, -1, -1, nLTBins, LTBinning, type=ROOT.TH2F)

            self.book(folder, "m2RelPFIsoDB", "m2Iso", 100, 0, 0.3)
            self.book(folder, "m1_t_Mass", "leadingMass", 200, 0, 200)
            self.book(folder, "m1_m2_Mass", "Muon 1-2 Mass", 120, 0, 120)

            # Rank muons by less MT to MET, for WZ control region
            self.book(folder, "weight", "Event weight", 100, 0, 5)
            self.book(folder, "rho", "Fastjet #rho", 100, 0, 25)
            self.book(folder, "nvtx", "Number of vertices", 31, -0.5, 30.5)
            self.book(folder, "subMTMass", "subMTMass", 200, 0, 200)
            self.book(folder, "nTruePU", "NPU", 62, -1.5, 60.5)
            self.book(folder, "m1DZ",  "m1DZ", 100, 0., 1)
            self.book(folder, "m2DZ",  "m2DZ", 100, 0., 1)
            self.book(folder, "tDZ" ,  "tDZ" , 100, 0., 1)
            self.book(folder, "doubleMuPrescale", "HLT prescale", 26, -5.5, 20.5)

    def preselection(self, row, cut_flow_trk = None, LT_threshold = 80., taupt_thr = 0.):
        ''' Preselection applied to events.

        Excludes FR object IDs and sign cut.
        '''

        #
        # NEW EXPLICIT CUT, BEFORE IT WAS IMPLICIT AND MADE WHILE PLOTTING!
        #
        if row.m1_t_DR < 0.5 or row.m2_t_DR < 0.5 or row.m1_m2_DR < 0.5: return False
        cut_flow_trk.Fill('DR separation')

        if row.m2_t_Mass < 20: return False
        cut_flow_trk.Fill('sub mass cut')

        double_mu_pass =  row.doubleMuPass and \
            row.m1MatchesDoubleMuPaths > 0 and \
            row.m2MatchesDoubleMuPaths > 0
        double_muTrk_pass = row.doubleMuTrkPass and \
             row.m1MatchesMu17TrkMu8Path  > 0 and \
             row.m2MatchesMu17TrkMu8Path  > 0
        if not ( double_mu_pass or double_muTrk_pass ): return False
        cut_flow_trk.Fill('trigger')

        if row.m1Pt < row.m2Pt:                   return False
        if row.m1Pt < 20:                         return False
        if not selections.muSelection(row, 'm1'): return False #applies basic selection (eta, pt > 10, DZ, pixHits, jetBTag)
        cut_flow_trk.Fill('obj1 Presel')

        if not selections.muSelection(row, 'm2'): return False #applies basic selection (eta, pt > 10, DZ, pixHits, jetBTag)
        cut_flow_trk.Fill('obj2 Presel')

        if not selections.tauSelection(row, 't'): return False #applies basic selection (eta, pt > 20, DZ)
        if not row.tAntiElectronMVA3Loose:        return False
        if row.tPt < taupt_thr: return False
        cut_flow_trk.Fill('obj3 Presel')

        if row.LT < LT_threshold: return False
        cut_flow_trk.Fill('LT')

        if row.m1_m2_Mass < 20:          return False
        cut_flow_trk.Fill('dilepton mass cut')

        if not selections.vetos(row, cut_flow_trk):    return False #applies mu bjet e additional tau vetoes
        
        return True

    @staticmethod
    def tau_sign_cut(row):
        ''' Returns true if muons are SS '''
        return not bool(row.m1_t_SS)

    @staticmethod
    def sign_cut(row):
        ''' Returns true if muons are SS '''
        return bool(row.m1_m2_SS)

    @staticmethod
    def obj1_id(row, leadleptonId='h2taucuts', subleadleptonId=None):
        return selections.lepton_id_iso(row, 'm1', leadleptonId)

    @staticmethod
    def obj2_id(row, leadleptonId=None, subleadleptonId='h2taucuts'):
        return selections.lepton_id_iso(row, 'm2', subleadleptonId)

    @staticmethod
    def obj3_id(row, tauId=None, LT_threshold = 80., taupt_thr = 0.):
        retval = False
        if row.LT >= LT_threshold and row.tPt >= taupt_thr:
            retval = bool(row.tLooseIso3Hits)
        else:
            retval = bool( getattr(row, tauId) )
        return retval

    @staticmethod
    def obj1_matches_gen(row):
        return row.m1GenPdgId == -1*row.m1Charge*13

    @staticmethod
    def obj2_matches_gen(row):
        return row.m2GenPdgId == -1*row.m2Charge*13

    @staticmethod
    def obj3_matches_gen(row):
        return t.genDecayMode != -2 
    
    @staticmethod
    def anti_wz(row):
        return bool(row.tAntiMuonTight and (( row.tDecayMode == 0 and row.tLeadChargeCandEMFraction > 0.2 ) or row.tDecayMode <> 0)) # and not row.tMuOverlap

    def enhance_wz(self, row):
        # Require the "tau" to be a muon, and require the third muon
        # to have M_Z +- 20
        # Cut on m2 PT > 20
        #if row.m2Pt < 20:
            #return False
        # Make sure any Z is from m1
        return row.tMuOverlap

    def event_weight(self, row):
        if row.run > 2:
            return 1.
        return self.pucorrector(row.nTruePU) * \
            mcCorrectors.get_muon_corrections(row,'m1','m2') * \
            mcCorrectors.double_muon_trigger(row,'m1','m2')

    def obj1_weight(self, row, ledleptonId='h2taucuts', subledleptonId=None):
        return frfits.highpt_mu_fr[ledleptonId](muonJetPt=max(row.m1JetPt, row.m1Pt), muonPt=row.m1Pt, numJets20=row.jetVeto20+1) #

    def obj2_weight(self, row, ledleptonId=None, subledleptonId='h2taucuts'):
        return frfits.lowpt_mu_fr[subledleptonId](muonJetPt=max(row.m2JetPt, row.m2Pt), muonPt=row.m2Pt, numJets20=row.jetVeto20+1) #

    def obj3_weight(self, row, notUsed1=None, notUsed2=None):
        return frfits.tau_fr(row.tPt, row.tAbsEta)

    def obj1_qcd_weight(self, row, ledleptonId='h2taucuts', subledleptonId=None):
        return frfits.highpt_mu_qcd_fr[ledleptonId](muonJetPt=max(row.m1JetPt, row.m1Pt), muonPt=row.m1Pt, numJets20=row.jetVeto20+1) #

    def obj2_qcd_weight(self, row, ledleptonId=None, subledleptonId='h2taucuts'):
        return frfits.lowpt_mu_qcd_fr[subledleptonId](muonJetPt=max(row.m2JetPt, row.m2Pt), muonPt=row.m2Pt, numJets20=row.jetVeto20+1) #

    def obj3_qcd_weight(self, row, notUsed1=None, notUsed2=None):
        return frfits.tau_qcd_fr(row.tPt, row.tAbsEta)

    # For measuring charge flip probability
    # Not really used in this channel
    def obj1_obj3_SS(self, row):
        return row.m1_t_SS

    ## def obj1_charge_flip(self, row):
    ##     return 0

    ## def obj2_charge_flip(self, row):
    ##     return 0

