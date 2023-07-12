# -------------------------
#   Top tagger 
# -------------------------
"""
  It is requiered: 
        * Muons:
                  pt >= 20:                      events.Muon.pt >= 20
                 |eta| <= 2.4                    np.abs(events.Muon.eta) <= 2.4
                  Tight ID                       events.Muon.tightId                     (True or False)
                  Iso <= 0.15                    events.Muon.pfRelIso04_all <= 0.15

        * Electrons:
                  pt >= 20                      events.Electron.pt >= 20
                |eta| <= 2.5                    np.abs(events.Electron.eta) <= 2.5
                Electron_cutBased >= 1          events.Electron.cutBased >= 1
                overlap(muons)                  events.Electron.delta_r(events.Muon)

"""


def select_muons(events: ak.Array):

    good_muons = (
        (events.Muon.pt >= 20)
        & (np.abs(events.Muon) <= 2.4)
        & (events.Muon.tightID)
        & (events.Muon.pfRelIso04_all <= 0.15)
    )
    
    n_good_muons = ak.sum(good_muons, axis=1)
    muons = ak.firsts(events.Muon[good_muons])
       
    return muons
    
    
    
        

def select_electrons(events: ak.Array):

    good_electrons = (
        (events.Electron.pt >= 20)
        & (np.abs(events.Electron) <= 2.5)
        & (events.Electron.cutBased >= 1)
        & (is_overlap(events,0.4))
    )
    
    n_good_electrons = ak.sum(good_electrons, axis=1)
    electrons = ak.firsts(events.Muon[good_electrons])
    
    
    return electrons


    electron_muon_dr = select_muons(events).delta_r(select_electrons(events))    # Overlap electron with

		

def is_overlap(events, value):
    if select_muons(events).delta_r(select_electrons(events) > value: return True
    else: return False




def select_top_jets(events: ak.Array, year: str, mod: str = ""):


    types_topjets = {"jetsTop_L": [], "jetsTop_M": [], "jetsTop_T": []}
    
    
    jet_id = 10000 #Recommendation in https://twiki.cern.ch/twiki/bin/view/CMS/JetID#Recommendations_for_13_TeV_data and https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVUL ("Please note: For AK8 jets, the corresponding (CHS or PUPPI) AK4 jet ID should be used.")
    pNet_id = 10000 #Recommendation in https://twiki.cern.ch/twiki/bin/view/CMS/ParticleNetTopWSFs
    
    if year + mod = '2018': 
        jet_id = 2
        pNet_id_L = 0.58
        pNet_id_M = 0.80
        pNet_id_T = 0.97
    elif year + mod= '2017': 
        jet_id = 2
        pNet_id_L = 0.58
        pNet_id_M = 0.80
        pNet_id_T = 0.97
    elif year + mod = '2016': 
        jet_id = 3
        pNet_id_L = 0.50
        pNet_id_M = 0.73 
        pNet_id_T = 0.96
    elif year + mod = '2016APV':
        jet_id = 3
        pNet_id_L = 0.49
        pNet_id_M = 0.74
        pNet_id_T = 0.96
    else: print('Check the year parameter entered.')

    
    good_fatjets_L(
        (events.FatJet.pt >= 300)
        & (np.abs(events.FatJet.eta) < 2.4)
        & (events.FatJet.jetId >= jet_id)   
 #       & (events.FatJet.delta_r(select_muons(events)) > 0.8)
 #       & (events.FatJet.delta_r(select_electrons(events)) > 0.8)
        & (events.FatJet.particleNet_TvsQCD >= pNet_id_L)
        
    )
        
    n_good_jetsTop_L = ak.sum(good_fatjets_L, axis=1)
    jetsTop_L = ak.firsts(events.FatJet[good_fatjets_L])
    types_topjets['jetsTop_L'].append(jetsTop_L)



####

    good_fatjets_M(
        (events.FatJet.pt >= 300)
        & (np.abs(events.FatJet.eta) < 2.4)
        & (events.FatJet.jetId >= jet_id)   
 #       & (events.FatJet.delta_r(select_muons(events)) > 0.8)
 #       & (events.FatJet.delta_r(select_electrons(events)) > 0.8)
        & (events.FatJet.particleNet_TvsQCD >= pNet_id_M)
        
    )
        
    n_good_jetsTop_M = ak.sum(good_fatjets_M, axis=1)
    jetsTop_M = ak.firsts(events.FatJet[good_fatjets_M])
    types_topjets['jetsTop_M'].append(jetsTop_M)
    
    
####

    good_fatjets_T(
        (events.FatJet.pt >= 300)
        & (np.abs(events.FatJet.eta) < 2.4)
        & (events.FatJet.jetId >= jet_id)   
 #       & (events.FatJet.delta_r(select_muons(events)) > 0.8)
 #       & (events.FatJet.delta_r(select_electrons(events)) > 0.8)
        & (events.FatJet.particleNet_TvsQCD >= pNet_id_T)
        
    )
        
    n_good_jetsTop_T = ak.sum(good_fatjets_T, axis=1)
    jetsTop_T = ak.firsts(events.FatJet[good_fatjets_T])
    types_topjets['jetsTop_T'].append(jetsTop_T)
    
    return types_topjets
    
    
def select_w_jets(event, year, selected_idx):

    types_wjets = {"jetsW_L": [], "jetsW_M": [], "jetsW_T": []}


    jet_id = 10000 #Recommendation in https://twiki.cern.ch/twiki/bin/view/CMS/JetID#Recommendations_for_13_TeV_data and https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVUL ("Please note: For AK8 jets, the corresponding (CHS or PUPPI) AK4 jet ID should be used.")
    pNet_id = 10000 #Recommendation in https://twiki.cern.ch/twiki/bin/view/CMS/ParticleNetTopWSFs

    if year + mod = '2018': 
        jet_id = 2
        pNet_id_L = 0.70
        pNet_id_M = 0.94
        pNet_id_T = 0.98
    elif year + mod= '2017': 
        jet_id = 2
        pNet_id_L = 0.71
        pNet_id_M = 0.94
        pNet_id_T = 0.98
    elif year + mod = '2016': 
        jet_id = 3
        pNet_id_L = 0.67
        pNet_id_M = 0.93 
        pNet_id_T = 0.97
    elif year + mod = '2016APV':
        jet_id = 3
        pNet_id_L = 0.68
        pNet_id_M = 0.94
        pNet_id_T = 0.97
    else: print('Check the year parameter entered.')
    
 
     good_wjets_L(
        (events.FatJet.pt >= 200)
        & (np.abs(events.FatJet.eta) < 2.4)
        & (events.FatJet.jetId >= jet_id)   
 #       & (events.FatJet.delta_r(select_muons(events)) > 0.8)
 #       & (events.FatJet.delta_r(select_electrons(events)) > 0.8)
 #       & (events.FatJet.delta_r(select_top_jets(events, year, mod)) > 0.8)   
       
    )
        
    n_good_wjetsTop_L = ak.sum(good_wjets_L, axis=1)
    wjetsTop_L = ak.firsts(events.FatJet[good_wjets_L])
    types_wjets['jetsW_L'].append(wjetsTop_L)



####

    good_wjets_M(
        (events.FatJet.pt >= 300)
        & (np.abs(events.FatJet.eta) < 2.4)
        & (events.FatJet.jetId >= jet_id)   
 #       & (events.FatJet.delta_r(select_muons(events), 0.8))
 #       & (events.FatJet.delta_r(select_electrons(events), 0.8))  
 #       & (events.FatJet.delta_r(select_top_jets(events, year, mod)) > 0.8)  
    )
        
    n_good_wjetsTop_M = ak.sum(good_wjets_M, axis=1)
    wjetsTop_M = ak.firsts(events.FatJet[good_fatjets_M])
    types_wjets['jetsW_M'].append(wjetsTop_M)
    
    
####

    good_wjets_T(
        (events.FatJet.pt >= 300)
        & (np.abs(events.FatJet.eta) < 2.4)
        & (events.FatJet.jetId >= jet_id)   
 #       & (events.FatJet.delta_r(select_muons(events), 0.8))
 #       & (events.FatJet.delta_r(select_electrons(events), 0.8))     
 #       & (events.FatJet.delta_r(select_top_jets(events, year, mod)) > 0.8)  
    )
        
    n_good_wjetsTop_T = ak.sum(good_wjets_T, axis=1)
    wjetsTop_T = ak.firsts(events.FatJet[good_wjets_T])
    types_wjets['jetsW_T'].append(wjetsTop_T)
    
    return types_wjets
    
   
 
    
def select_jets(event, year, selected_idx):
    types_jets = {'jets_noCleaned_against_boostedJets': [],  'jets_cleaned': [], 'jetsB_L': [], 'jetsB_M': [], 'jetsB_T': [] }
    
    jet_id = 10000 #Recommendation in https://twiki.cern.ch/twiki/bin/view/CMS/JetID#Recommendations_for_13_TeV_data and https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVUL ("Please note: For AK8 jets, the corresponding (CHS or PUPPI) AK4 jet ID should be used.")
    bjet_id = 10000 #Recommendation in https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation
    bjet_eta = -1
    jet_eta = -1

    
    if year + mod = '2018': 
        jet_id = 2
        jet_eta = 2.5
        
        bjet_eta = 2.5
        bjet_id = 0.0490
        bjet_id = 0.2783
        bjet_id = 0.7100
        
    elif year + mod= '2017': 
        jet_id = 2
        jet_eta = 2.5
        
        bjet_eta = 2.5
        bjet_id = 0.0532
        bjet_id = 0.3040
        bjet_id = 0.7476

    elif year + mod = '2016': 
        jet_id = 3
        jet_eta = 2.4
        
        bjet_eta = 2.4
        bjet_id = 0.0480
        bjet_id = 0.2489
        bjet_id = 0.6377

    elif year + mod = '2016APV':
        jet_id = 3
        jet_eta = 2.4
        
        bjet_eta = 2.4
        bjet_id = 0.0508
        bjet_id = 0.2598
        bjet_id = 0.6502

    else: print('Check the year parameter entered.')
    
   
     good_jets_no_cleaned(
        (events.Jet.pt >= 20)
        & (np.abs(events.Jet.eta) < jet_eta)
        & (events.Jet.jetId >= jet_id)   
 #       & (events.Jet.delta_r(select_muons(events)) > 0.4)
 #       & (events.Jet.delta_r(select_electrons(events)) > 0.4)      
    )
        
    n_good_no_jets_cleaned = ak.sum(good_jetsv_cleaned, axis=1)
    jets_no_cleaned = ak.firsts(events.FatJet[good_jets_no_cleaned])
    types_jets_no_cleaned['jets_noCleaned_against_boostedJets'].append(jets_no_cleaned)



    
     good_jets_cleaned(
        (events.Jet.pt >= 20)
        & (np.abs(events.Jet.eta) < jet_eta)
        & (events.Jet.jetId >= jet_id)   
 #       & (events.Jet.delta_r(select_muons(events)) > 0.4)
 #       & (events.Jet.delta_r(select_electrons(events)) > 0.4)     
 #       & (events.FatJet.delta_r(select_top_jets(events, year, mod)) > 0.8)  
 #       & (events.FatJet.delta_r(select_w_jets(events, year, mod)) > 0.8)   
    )
        
    n_good_jets_cleaned = ak.sum(good_jets_cleaned, axis=1)
    jets_cleaned = ak.firsts(events.FatJet[good_jets_cleaned])
    types_jets_cleaned['jets_cleaned'].append(jets_cleaned)


     good_bjets_L(
        (events.Jet.pt >= 20)
        & (np.abs(events.Jet.eta) < bjet_eta)
        & (events.Jet.btagDeepFlavB >= bjet_id)   
 #       & (events.Jet.delta_r(select_muons(events)) > 0.4)
 #       & (events.Jet.delta_r(select_electrons(events)) > 0.4)     
 #       & (events.FatJet.delta_r(select_top_jets(events, year, mod)) > 0.8)  
 #       & (events.FatJet.delta_r(select_w_jets(events, year, mod)) > 0.8)   
    )
    n_good_bjets_L = ak.sum(good_bjets_cleaned, axis=1)
    bjets_L = ak.firsts(events.FatJet[good_bjets_L])
    types_bjets_cleaned['jetsB_L'].append(bjets_L)



     good_bjets_M(
        (events.Jet.pt >= 20)
        & (np.abs(events.Jet.eta) < bjet_eta)
        & (events.Jet.btagDeepFlavB >= bjet_id)   
 #       & (events.Jet.delta_r(select_muons(events)) > 0.4)
 #       & (events.Jet.delta_r(select_electrons(events)) > 0.4)     
 #       & (events.FatJet.delta_r(select_top_jets(events, year, mod)) > 0.8)  
 #       & (events.FatJet.delta_r(select_w_jets(events, year, mod)) > 0.8)   
    )
    n_good_bjets_M = ak.sum(good_bjets_cleaned, axis=1)
    bjets_M = ak.firsts(events.FatJet[good_bjets_M])
    types_bjets_cleaned['jetsB_M'].append(bjets_M)
    
  
     good_bjets_T(
        (events.Jet.pt >= 20)
        & (np.abs(events.Jet.eta) < bjet_eta)
        & (events.Jet.btagDeepFlavB >= bjet_id)   
 #       & (events.Jet.delta_r(select_muons(events)) > 0.4)
 #       & (events.Jet.delta_r(select_electrons(events)) > 0.4)     
 #       & (events.FatJet.delta_r(select_top_jets(events, year, mod)) > 0.8)  
 #       & (events.FatJet.delta_r(select_w_jets(events, year, mod)) > 0.8)   
    )
    n_good_bjets_T = ak.sum(good_bjets_cleaned, axis=1)
    bjets_T = ak.firsts(events.FatJet[good_bjets_T])
    types_bjets_cleaned['jetsB_T'].append(bjets_T)  
    
    

    return types_jets
