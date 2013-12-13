'''

Analyze EET events for the WH analysis

'''

from EETauTree import EETauTree
import os
from WHAnalyzerBase import WHAnalyzerBase, quad, inv_mass
import ROOT
import mcCorrectors
import baseSelections as selections
import fakerate_functions as frfits
import optimizer
import math
from array import array
from chargeflipcuts import charge_flip_funcs
from FinalStateAnalysis.PlotTools.decorators import memo_last
import FinalStateAnalysis.PlotTools.pytree as pytree

#initialize FRFits
optimizer_keys   = [ i for i in optimizer.grid_search.keys() if i.startswith('EET') ]
grid_search = {}
if len(optimizer_keys) > 1:
        grid_search[key] = optimizer.grid_search[key]

NO_HISTO_FILL = ('NO_HISTO_FILL' in os.environ) and eval(os.environ['NO_HISTO_FILL']) 

if not NO_HISTO_FILL:
    for key in optimizer_keys:
        frfits.highpt_ee_fr.__activate__(optimizer.grid_search[key]['leading_iso'])
        frfits.highpt_ee_qcd_fr.__activate__(optimizer.grid_search[key]['leading_iso'])
        frfits.lowpt_ee_fr.__activate__(optimizer.grid_search[key]['subleading_iso'])
        frfits.lowpt_ee_qcd_fr.__activate__(optimizer.grid_search[key]['subleading_iso'])

SYNC = ('SYNC' in os.environ) and eval(os.environ['SYNC'])

#mtr = frfits.mt_likelihood_ratio
################################################################################
#### Analysis logic ############################################################
################################################################################

class WHAnalyzeEET(WHAnalyzerBase):
    tree = 'eet/final/Ntuple' if not SYNC else 'Ntuple'
    def __init__(self, tree, outfile, **kwargs):
        self.channel = 'EET'	
        super(WHAnalyzeEET, self).__init__(tree, outfile, EETauTree, **kwargs)
        self.hfunc['subMTMass'] = lambda row, weight: (row.e2_t_Mass, weight) if row.e1MtToMET > row.e2MtToMET else (row.e1_t_Mass, weight) #maps the name of non-trivial histograms to a function to get the proper value, the function MUST have two args (evt and weight). Used in WHAnalyzerBase.fill_histos later
	def make_both_barrel(fcn, negval=(0.,0.)):
		def f_(row, weight):
			return fcn(row, weight) if row.e1AbsEta < 1.48 and row.e2AbsEta < 1.48 else negval
		return f_

	def make_both_endcap(fcn, negval=(0.,0.)):
		def f_(row, weight):
			return fcn(row, weight) if row.e1AbsEta >= 1.48 and row.e2AbsEta >= 1.48 else negval
		return f_

	def make_mixed(fcn, negval=(0.,0.)):
		def f_(row, weight):
			return fcn(row, weight) \
			    if (row.e1AbsEta >= 1.48 and row.e2AbsEta < 1.48) or \
			    (row.e1AbsEta < 1.48 and row.e2AbsEta >= 1.48) \
		            else negval
		return f_

	def attr_getter(attribute):
            def f(row, weight):
                return (getattr(row,attribute), weight)
            return f

	def double_attr_getter(attr1, attr2):
            def f(row, weight):
                return ( (getattr(row,attr1), getattr(row,attr1)), weight)
            return f

        def charge_selector(fcn, chargeIds, negval=(0.,0.)):
            def f(row, weight):
                ok = all( getattr(row, chid) for chid in chargeIds )
                return fcn(row, weight) if ok else negval
            return f

        def mass_scaler(fcn, scaler=frfits.data_scaler):
            def f(row, weight):
                val, w = fcn(row, weight)
                return scaler(val), w
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
        
        @memo_last
        def f_par_prob(e1Pt, e1JetPt,
                       e2Pt, e2JetPt,
                       tPt):
            p_e1 = (( frfits.highpt_ee_fr[lead_iso](electronJetPt=max(e1JetPt, e1Pt), electronPt=e1Pt) +\
                      frfits.highpt_ee_qcd_fr[lead_iso](electronJetPt=max(e1JetPt, e1Pt), electronPt=e1Pt) )/2)
            p_e2 = (( frfits.lowpt_ee_fr[sublead_iso](electronJetPt=max(e2JetPt, e2Pt), electronPt=e2Pt) + \
                      frfits.lowpt_ee_qcd_fr[sublead_iso](electronJetPt=max(e2JetPt, e2Pt), electronPt=e2Pt))/2)
            p_t  = frfits.tau_fr(tPt)
            return (p_e1 + p_e2*(1 - p_e1) + p_t*(1 - p_e1)*(1 - p_e2))
            

        def f_prob(row, weight):
            val = f_par_prob(row.e1Pt, row.e1JetPt,
                             row.e2Pt, row.e2JetPt,
                             row.tPt)
            return val, weight

        def log_prob(row, weight):
            prob, weight = f_prob(row, weight)
            return ROOT.TMath.Log10(prob), weight
        
        self.hfunc['faking_prob'] = f_prob
        self.hfunc['log_prob']    = log_prob
        self.hfunc["e*2_t_Mass#LT" ] = merge_functions( mass_scaler( attr_getter('e2_t_Mass')), attr_getter('LT'))
        self.hfunc["e*2_t_Mass#tPt"] = merge_functions( mass_scaler( attr_getter('e2_t_Mass')), attr_getter('tPt'))
        self.hfunc["e*1_e2_Mass"] = mass_scaler( attr_getter('e1_e2_Mass'))
        self.hfunc["e*1_t_Mass" ] = mass_scaler( attr_getter('e1_t_Mass')) 
        self.hfunc["e1_e*2_Mass"] = mass_scaler( attr_getter('e1_e2_Mass'))
        self.hfunc["e*2_t_Mass" ] = mass_scaler( attr_getter('e2_t_Mass')) 

        self.hfunc["e*+2_t_Mass#LT" ] = merge_functions( mass_scaler( attr_getter('e2_t_Mass'), scaler=frfits.default_scaler), attr_getter('LT'))
        self.hfunc["e*+2_t_Mass#tPt"] = merge_functions( mass_scaler( attr_getter('e2_t_Mass'), scaler=frfits.default_scaler), attr_getter('tPt'))
        self.hfunc["e*+1_e2_Mass"] = mass_scaler( attr_getter('e1_e2_Mass'), scaler=frfits.default_scaler)
        self.hfunc["e*+1_t_Mass" ] = mass_scaler( attr_getter('e1_t_Mass') , scaler=frfits.default_scaler) 
        self.hfunc["e1_e*+2_Mass"] = mass_scaler( attr_getter('e1_e2_Mass'), scaler=frfits.default_scaler)
        self.hfunc["e*+2_t_Mass" ] = mass_scaler( attr_getter('e2_t_Mass') , scaler=frfits.default_scaler) 

        self.hfunc['pt_ratio' ] = lambda row, weight: (row.e2Pt/row.e1Pt, weight)
        self.hfunc['SYNC'] = lambda row, weight: (row, None) #((row.run, row.lumi, row.evt, row.e1Pt, row.e1Eta, row.e1Phi, row.e2Pt, row.e2Eta, row.e2Phi, row.tPt, row.tEta, row.tPhi, weight), None)

        #self.hfunc['evt_info'] = lambda row, weight: (array.array("f", [row.e1Pt, row.e2Pt, row.tPt, row.LT, weight] ), None)
        self.pucorrector = mcCorrectors.make_puCorrector('doublee')

    @staticmethod
    def anti_charge_flip(row, rejection_power=80):
        #print 'running anti_charge_flip %s' % rejection_power
        return charge_flip_funcs[rejection_power](row)

    def book_histos(self, folder):
	#PLOTS TO FILL IN ANY CASE
        LTBinning = array('d',[0] + range(50,210,10) + [600])
        nLTBins   = len(LTBinning) -1
	for key in self.grid_search:
	    prefix = key+'$' if key else ''
	    self.book(folder, prefix+"LT", "L_T", 100, 0, 500)
            self.book(folder, prefix+"e2_t_Mass", "subleadingMass", 300, 0, 300)
            self.book(folder, prefix+"e2_t_Pt",   "subleadingPt", 400, 0, 400)
            #Charge mis-id special histograms
            if 'c2' in folder:
                self.book(folder, prefix+"e*2_t_Mass", "subleadingMass with misid sclaing correction", 300, 0, 300)
                self.book(folder, prefix+"e*+2_t_Mass", "subleadingMass with misid sclaing correction", 300, 0, 300)

        if len(self.grid_search.keys()) == 1:
            if 'c1' in folder:
                self.book(folder, "e*1_e2_Mass", "E 1-2 Mass with misid sclaing correction", 120, 0, 120)
                self.book(folder, "e*1_t_Mass", "leadingMass with misid sclaing correction", 200, 0, 200)
                self.book(folder, "e*+1_e2_Mass", "E 1-2 Mass with misid sclaing correction", 120, 0, 120)
                self.book(folder, "e*+1_t_Mass", "leadingMass with misid sclaing correction", 200, 0, 200)
            elif 'c2' in folder:
                self.book(folder, "e1_e*2_Mass", "E 1-2 Mass with misid sclaing correction", 120, 0, 120)
                self.book(folder, "e*2_t_Mass#LT"         , '', 300, 0, 300, nLTBins, LTBinning, type=ROOT.TH2F)
                self.book(folder, "e1_e*+2_Mass", "E 1-2 Mass with misid sclaing correction", 120, 0, 120)
                self.book(folder, "e*+2_t_Mass#LT"         , '', 300, 0, 300, nLTBins, LTBinning, type=ROOT.TH2F)


            self.book(folder, prefix+"e2_t_Mass#LT" , "subleadingMass", 300, 0, 300, nLTBins, LTBinning, type=ROOT.TH2F)

            self.book(folder, 'SYNC', 'SYNC', 
                      'run/l:lumi/l:evt/l' +\
                      ':e1Pt/D:e1Eta/D:e1Phi/D:e1RelPFIsoDB/D:e1MVANonTrig/D' +\
                      ':e2Pt/D:e2Eta/D:e2Phi/D:e2RelPFIsoDB/D:e2MVANonTrig/D' +\
                      ':tPt/D:tEta/D:tPhi/D:tAntiMuonTight/D:tDecayMode/D:tLooseIso3Hits/D:tAntiElectronMVA3Loose/D' +\
                      ':tAntiMuonLoose/D:tAntiElectronMVA3Medium/D:tAntiElectronMVA3Tight/D' +\
                      ':e1_e2_Mass/D:e2_t_Mass/D:LT/D:e1_e2_SS/D:e1_t_SS/D',
                      type=pytree.PyTree)

            #Pt
            self.book(folder, prefix+"e1Pt#LT" , "subleadingMass", 150, 0, 150, nLTBins, LTBinning, type=ROOT.TH2F)
            self.book(folder, prefix+"e2Pt#LT" , "subleadingMass", 150, 0, 150, nLTBins, LTBinning, type=ROOT.TH2F)
            self.book(folder, prefix+"tPt#LT" , "subleadingMass", 150, 0, 150, nLTBins, LTBinning, type=ROOT.TH2F)
            self.book(folder, prefix+"e1JetPt#LT" , "subleadingMass", 150, 0, 150, nLTBins, LTBinning, type=ROOT.TH2F)
            self.book(folder, prefix+"e2JetPt#LT" , "subleadingMass", 150, 0, 150, nLTBins, LTBinning, type=ROOT.TH2F)

            #eta
            self.book(folder, prefix+"e1AbsEta#LT" , "subleadingMass", 100, 0, 2.5, nLTBins, LTBinning, type=ROOT.TH2F)
            self.book(folder, prefix+"e2AbsEta#LT" , "subleadingMass", 100, 0, 2.5, nLTBins, LTBinning, type=ROOT.TH2F)
            self.book(folder, prefix+"tAbsEta#LT" , "subleadingMass", 100, 0, 2.5, nLTBins, LTBinning, type=ROOT.TH2F)

            #DR
            self.book(folder, prefix+"e1_t_DR#LT" , "subleadingMass", 100, 0, 10, nLTBins, LTBinning, type=ROOT.TH2F)
            self.book(folder, prefix+"e2_t_DR#LT" , "subleadingMass", 100, 0, 10, nLTBins, LTBinning, type=ROOT.TH2F)

            #Jet BTag
            self.book(folder, "e1JetBtag#LT", "Muon 2 Pt", 30, -10, 5, 60, 0, 600, type=ROOT.TH2F)
            self.book(folder, "e2JetBtag#LT", "Muon 2 Pt", 30, -10, 5, 60, 0, 600, type=ROOT.TH2F)

            #tau ISO 
            self.book(folder, "tRawIso3Hits#LT", "Muon 2 Pt", 200, 0, 200, nLTBins, LTBinning, type=ROOT.TH2F)
            self.book(folder, "tRawIsoMVA2#LT", "Muon 2 Pt" , 200, -1, -1, nLTBins, LTBinning, type=ROOT.TH2F)

            self.book(folder, "e2RelPFIsoDB", "e2RelPFIsoDB", 30, 0, 0.3)
            self.book(folder, "e1RelPFIsoDB", "e1RelPFIsoDB", 30, 0, 0.3)
            self.book(folder, "e1_e2_Mass", "E 1-2 Mass", 120, 0, 120)
            self.book(folder, "e1_t_Mass", "leadingMass", 200, 0, 200)
            self.book(folder, "weight", "Event weight", 100, 0, 5)
            self.book(folder, "rho", "Fastjet #rho", 100, 0, 25)
            self.book(folder, "nvtx", "Number of vertices", 31, -0.5, 30.5)
            self.book(folder, "e1Pt", "E 1 Pt", 100, 0, 100)
            self.book(folder, "e2Pt", "E 2 Pt", 100, 0, 100)
            self.book(folder, "e1JetPt", "E 1 Pt", 100, 0, 100)
            self.book(folder, "e2JetPt", "E 2 Pt", 100, 0, 100)
            self.book(folder, "e1AbsEta", "Muon 1 AbsEta", 100, 0, 2.5)
            self.book(folder, "e2AbsEta", "Muon 2 AbsEta", 100, 0, 2.5)
	    self.book(folder, "tPt", "tPt", 100, 0,100)
	    self.book(folder, "tAbsEta", "tAbsEta", 100, 0, 2.3)
            self.book(folder, "e1_e2_DR" , "e1_e2_DR" , 500, 0, 10)
            self.book(folder, "e1_t_DR"  , "e1_e2_DR" , 500, 0, 10)
            self.book(folder, "e2_t_DR"  , "e1_e2_DR" , 500, 0, 10)

            #let's look for osme other possible selections
            self.book(folder, "Mass"          , "Mass"      , 100, 0, 1)
            self.book(folder, "pt_ratio"      , "pt_ratio"      , 100, 0, 1)
            self.book(folder, "tToMETDPhi"    , "tToMETDPhi"    , 100, 0, 4)
            self.book(folder, "e1_e2_Pt"       , "lepRecoil"     , 600, 0, 8000)
            self.book(folder, "e1_e2_CosThetaStar", "e1_e2_CosThetaStar", 110, 0., 1.1)
            self.book(folder, "e1_t_CosThetaStar" , "e1_t_CosThetaStar" , 110, 0., 1.1)
            self.book(folder, "e2_t_CosThetaStar" , "e2_t_CosThetaStar" , 110, 0., 1.1)
            self.book(folder, "type1_pfMetEt"         , "metEt"         , 300, 0, 300)
            self.book(folder, "mva_metEt"         , "metEt"         , 300, 0, 300)
            #SPECIAL NTUPLE!
            #self.book(folder, 'evt_info', 'evt_info', 'e1Pt:e2Pt:tPt:LT:weight', type=ROOT.TNtuple)

            
    #There is no call to self, so just promote it to statucmethod, to allow usage by other dedicated analyzers
    def preselection(self, row, cut_flow_trk = None, LT_threshold = 80., taupt_thr = 0.):
        ''' Preselection applied to events.

        Excludes FR object IDs and sign cut.
        '''
        if row.e1_t_DR < 0.5 or row.e2_t_DR < 0.5 or row.e1_e2_DR < 0.5: return False
        cut_flow_trk.Fill('DR separation')

        if row.e2_t_Mass < 20: return False
        cut_flow_trk.Fill('sub mass cut')

        if not row.doubleEPass:                   return False
	if not (row.e1MatchesDoubleEPath > 0 and \
		row.e2MatchesDoubleEPath > 0): return False 
        cut_flow_trk.Fill('trigger')

	if row.e1Pt < 20:           return False
        if not selections.eSelection(row, 'e1'):  return False
	cut_flow_trk.Fill('obj1 Presel')

        if not selections.eSelection(row, 'e2'):  return False
        cut_flow_trk.Fill('obj2 Presel')

        if not selections.tauSelection(row, 't'): return False
        if row.tPt < taupt_thr: return False
        if not row.tAntiMuonLoose:   return False
        cut_flow_trk.Fill('obj3 Presel')

        if row.LT < LT_threshold: return False
        cut_flow_trk.Fill('LT')

        #if row.e1_e2_SS and row.e1_t_SS: return False #remove three SS leptons
        if row.e1_e2_Mass < 20:          return False
        cut_flow_trk.Fill('dilepton mass cut')

        if not selections.vetos(row, cut_flow_trk):    return False #applies mu bjet e additional tau vetoes

        return True

    #There is no call to self, so just promote it to statucmethod, to allow usage by other dedicated analyzers
    @staticmethod
    def sign_cut(row):
        ''' Returns true if muons are SS '''
        return bool(row.e1_e2_SS)

    @staticmethod
    def tau_sign_cut(row):
        ''' Returns true if muons are SS '''
        return not bool(row.e1_t_SS)

    #There is no call to self, so just promote it to statucmethod, to allow usage by other dedicated analyzers
    @staticmethod
    def obj1_id(row, leadleptonId='eid13Looseh2taucuts', subleadleptonId=None):
        return selections.lepton_id_iso(row, 'e1', leadleptonId)

    #There is no call to self, so just promote it to statucmethod, to allow usage by other dedicated analyzers
    @staticmethod
    def obj2_id(row, leadleptonId=None, subleadleptonId='eid13Looseh2taucuts'):
        return selections.lepton_id_iso(row, 'e2', subleadleptonId)

    #There is no call to self, so just promote it to statucmethod, to allow usage by other dedicated analyzers
    @staticmethod
    def obj3_id(row, tauId=None, LT_threshold = 80., taupt_thr = 0.):
        if row.LT >= LT_threshold and row.tPt >= taupt_thr:
            return bool(row.tLooseIso3Hits)
        else:
            return bool( getattr(row, tauId) )

    @staticmethod
    def obj1_matches_gen(row):
        return row.e1GenPdgId == -1*row.e1Charge*11

    @staticmethod
    def obj2_matches_gen(row):
        return row.e2GenPdgId == -1*row.e2Charge*11

    @staticmethod
    def obj3_matches_gen(row):
        return t.genDecayMode != -2 

    #There is no call to self, so just promote it to statucmethod, to allow usage by other dedicated analyzers
    @staticmethod
    def anti_wz(row):
        if row.e2_t_Zcompat < 10 or row.e1_t_Zcompat < 10 :
            if not row.tAntiElectronMVA3Tight:
                return False
        if row.e2_t_Zcompat < 20 or row.e1_t_Zcompat < 20 :
            if not row.tAntiElectronMVA3Medium:
                return False
        return row.tAntiElectronMVA3Loose

    def enhance_wz(self, row):
        # Require the "tau" to be a electron, and require the third electron
        # to have M_Z +- 20 
        return row.tCiCTightElecOverlap

    def event_weight(self, row):
        if row.run > 2: #FIXME! add tight ID correction
            return 1.
        return self.pucorrector(row.nTruePU) * \
                mcCorrectors.get_electron_corrections(row,'e2') *\
                mcCorrectors.electron_tight_corrections(row.e1Pt, row.e1AbsEta) *\
                mcCorrectors.double_electron_trigger(row)
   
    def obj1_weight(self, row, leadleptonId='eid13Looseh2taucuts', subleadleptonId=None):
        return frfits.highpt_ee_fr[leadleptonId](electronJetPt=max(row.e1JetPt, row.e1Pt), electronPt=row.e1Pt, numJets20=row.jetVeto20+1)

    def obj2_weight(self, row, leadleptonId=None, subleadleptonId='eid13Looseh2taucuts'):
	return frfits.lowpt_ee_fr[subleadleptonId](electronJetPt=max(row.e2JetPt, row.e2Pt), electronPt=row.e2Pt, numJets20=row.jetVeto20+1)

    def obj3_weight(self, row, notUsed1=None, notUsed2=None):
        return frfits.tau_fr(row.tPt, row.tAbsEta)

    def obj1_qcd_weight(self, row, leadleptonId='eid13Looseh2taucuts', subleadleptonId=None):
        return frfits.highpt_ee_qcd_fr[leadleptonId](electronJetPt=max(row.e1JetPt, row.e1Pt), electronPt=row.e1Pt, numJets20=row.jetVeto20+1)

    def obj2_qcd_weight(self, row, leadleptonId=None, subleadleptonId='eid13Looseh2taucuts'):
        return frfits.lowpt_ee_qcd_fr[subleadleptonId](electronJetPt=max(row.e2JetPt, row.e2Pt), electronPt=row.e2Pt, numJets20=row.jetVeto20+1)

    def obj3_qcd_weight(self, row, notUsed1=None, notUsed2=None):
        return frfits.tau_qcd_fr(row.tPt, row.tAbsEta)

    # For measuring charge flip probability
    # Not really used in this channel
    def obj1_obj3_SS(self, row):
        return not row.e1_t_SS

    def obj1_charge_flip(self, row, leadleptonId='eid13Looseh2taucuts', subleadleptonId=None):
        return frfits.e_charge_flip[leadleptonId](row.e1AbsEta,row.e1Pt) #highpt_e_charge_flip

    def obj2_charge_flip(self, row, leadleptonId=None, subleadleptonId='eid13Looseh2taucuts'):
        return frfits.e_charge_flip[subleadleptonId](row.e2AbsEta,row.e2Pt) #lowpt_e_charge_flip

    def obj1_charge_flip_sysup(self, row, leadleptonId='eid13Looseh2taucuts', subleadleptonId=None):
        return frfits.e_charge_flip_up[leadleptonId](row.e1AbsEta,row.e1Pt) #highpt_e_charge_flip

    def obj2_charge_flip_sysup(self, row, leadleptonId=None, subleadleptonId='eid13Looseh2taucuts'):
        return frfits.e_charge_flip_up[subleadleptonId](row.e2AbsEta,row.e2Pt) #lowpt_e_charge_flip







