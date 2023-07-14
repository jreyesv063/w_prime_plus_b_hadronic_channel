from functions import *
from variables import *
from topsXFinder import *


class declareVariables(variables):
    def __init__(self, name):
        super(declareVariables, self).__init__(name)
        
        
    
    def beginJob(self):
        print("Here is beginJob")
        
        
    def endJob(self):
        print("Here is endJob")
        
        
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        print("Here is beginFile")
        
        self.sumNumEvt = 0
        self.sumgenWeight = 0
        self.out = declareVariables(inputFile)
        
        
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        print("Here is endFile")
        
        self.out.sumNumEvt[0] = self.sumNumEvt
        self.out.sumgenWeight[0] = self.sumgenWeight
        self.out.evtree.Fill()
        self.out.outputfile.Write()
        self.out.outputfile.Close()
       
       
###########3

    def analyze(self, events, selected_idx):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        
        # For all events
        if (self.sumNumEvt > self.maxNumEvt and self.maxNumEvt != -1): 
        
            return False
                
        self.sumNumEvt = self.sumNumEvt + 1
        if not self.isData: self.sumgenWeight = self.sumgenWeight + (events.genWeight / abs(events.genWeight))
        if not self.sumNumEvt % self.prescaleEvt == 0: return False

        selected_idx = {'muon': [], 'electron': [],
                        'jetsTop_L': [], 'jetsTop_M': [], 'jetsTop_T': [], 
                        'jetsW_L': [], 'jetsW_M': [], 'jetsW_T': [], 
                        'jets_cleaned': [], 'jets_noCleaned_against_boostedJets': [],
                        'jetsB_L': [], 'jetsB_M': [], 'jetsB_T': []}
                        
        muons = events.Muon
        select_muons(events)

        eles = events.Electron
        select_electrons(events)

        fatjets = events.FatJet
        select_top_jets(events, self.year, mode)
        select_w_jets(events, self.year, mode)

        jets = events.Jet
        select_jets(eventS, self.year)

        if not ((len(selected_idx['muon']) + len(selected_idx['electron'])) == 1):
            return False

        # top reconstruction
        # these are for the total results
        results = {'top': [], 
                    'w': [], 
                    'top_topology_decay': [],
                    'chi2': []}
        # Top to keep all the Top(TLorentzVector), W to keep all the W(TLorentzVector) (boosted top part,just using (0,0,0,0) instead)
        # topTopologyDecay to define the kind of top (0: boosted, 1: partially-boosted, 2: resolved hadronic, 3: resolved leptonic)
        
        top_type = {'boosted_top': True, 'partially_boosted_top': True, 'hadronic_top': True, 'leptonic_top': True}
        topsXFinder(event, selected_idx, top_type, results)

        # Event
        self.out.run[0] = events.run
        self.out.luminosityBlock[0] = events.luminosityBlock
        self.out.event[0] = events.events  # & 0xffffffffffffffff
        self.out.num_top_total[0] = len(results['top'])

        # Save tree
        self.out.Events.Fill()
        
        return True
