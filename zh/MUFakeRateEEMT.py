'''
Base analyzer for muon id fake-rate estimation: Z->e e
'''

import EMUFakeRatesBase
import baseSelections as selections
from EEMuTauTree import EEMuTauTree

class MUFakeRateEEMT(EMUFakeRatesBase.EMUFakeRatesBase):
    tree = 'eemt/final/Ntuple'
    def __init__(self, tree, outfile, **kwargs):
        super(MUFakeRateEEMT, self).__init__(tree, outfile, EEMuTauTree, **kwargs)
        self.lepton   = 'muon'
        self.branchId = 'm'
        
    def zSelection(self, row):
        if not selections.ZEESelectionNoVetos(row): return False
        if bool(row.muGlbIsoVetoPt10):  return False
        if selections.overlap(row, 'e1','e2','m','t') : return False
        return selections.signalMuonSelection(row,'m')

    def lepton_passes_iso(self, row):
        return bool(row.mRelPFIsoDB < 0.15) ##THIS SEEMS too low        
    
