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


def select_muons(events: ak.Array, selected_idx):
    
    selections_muons = PackedSelection()
    
    Muon = events.Muon
    
    good_muons = (
        (Muon.pt < 20)
        & (np.abs(Muon.eta) > 2.4)
        & (Muon.tightId)
        & (Muon.pfRelIso04_all > 0.15)
    )
    
    n_good_muons = ak.sum(good_muons, axis=1)
    muons = ak.firsts(Muon[good_muons])
    
    
    selections_muons.add("one_muon", n_good_muons == 1)  
        
    # define the event selection criteria
    regions = {
        "mu": ["one_muon"],
    }
    
    Filtered = selections_muons.all(*regions["mu"])                                 
    selected_idx['muon'].append(Filtered)
                                
    return selected_idx['muon']
    
    
    
        

def select_electrons(events: ak.Array, selected_idx):

    selections_electrons = PackedSelection()
    
    Muon = events.Muon
    Electron = events.Electron


    good_electrons = (
        (Electron.pt < 20)
        & (np.abs(Electron.eta) > 2.5)
        & (Electron.cutBased < 1)
    )

        
    n_good_electrons = ak.sum(good_electrons, axis=1) 
    overlap_e_mu = ak.firsts(Muon).delta_r(ak.firsts(Electron))
    
    
    selections_electrons.add("one_electron", n_good_electrons == 1)  
    selections_electrons.add("electron_muon_dr", overlap_e_mu > 0.4)
    
    
    # define the event selection criteria
    regions = {
        "ele": ["one_electron", "electron_muon_dr"],
    }
        
    Filtered = selections_electrons.all(*regions["ele"])                                 
    selected_idx['electron'].append(Filtered)
    
    return selected_idx['electron']

##########################################


def select_top_jets(events: ak.Array, year: str, mod: str, selected_idx):
   
    selections_top_jets = PackedSelection()
 

    FatJet = events.FatJet
    Muon = events.Muon
    Electron = events.Electron
    
    
    jet_id = 10000 #Recommendation in https://twiki.cern.ch/twiki/bin/view/CMS/JetID#Recommendations_for_13_TeV_data and https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVUL ("Please note: For AK8 jets, the corresponding (CHS or PUPPI) AK4 jet ID should be used.")
    pNet_id = 10000 #Recommendation in https://twiki.cern.ch/twiki/bin/view/CMS/ParticleNetTopWSFs
    
    # Loose: L subindex.
    # Medium: M subindex.
    # Tight: T subindex.
    # Nomenclature: "year": [jetd_id, pNet_id_L, pNet_id_M, pNet_id_T]
    
    constants_as_fuction_wp = {
                                "2016" : [3, 0.50, 0.73, 0.96], 
                                "2016APV" : [3, 0.49, 0.74, 0.96], 
                                "2017": [2, 0.58, 0.80, 0.97], 
                                "2018": [2, 0.58, 0.80, 0.97]
                            }
    

    jet_id = constants_as_fuction_wp[year + mod][0]
    pNet_id_L = constants_as_fuction_wp[year + mod][1]
    pNet_id_M = constants_as_fuction_wp[year + mod][2]
    pNet_id_T = constants_as_fuction_wp[year + mod][3]
    
    
    
    good_fatjets = (
        (FatJet.pt < 300)
        & (np.abs(FatJet.eta) >= 2.4) #This is the recommendation for all the fat jets (there are not reconstructed forward fat jets)
        & (FatJet.jetId < jet_id)   
    )
    
    n_good_fatjets = ak.sum(good_fatjets, axis=1)
    
    
    leading_fatjets = ak.firsts(FatJet[good_fatjets])
    overlap_fatjets_electron = ak.firsts(FatJet).delta_r(ak.firsts(Electron))
    overlap_fatjets_muon = ak.firsts(FatJet).delta_r(ak.firsts(Muon))
    
    
    selections_top_jets.add("fatjets_electron_dr", overlap_fatjets_electron > 0.8)
    selections_top_jets.add("fatjets_muon_dr", overlap_fatjets_muon > 0.8)              
    selections_top_jets.add("leading_fatjets", n_good_fatjets == 1)
    selections_top_jets.add("id_L", FatJet.particleNet_TvsQCD < pNet_id_L)
    selections_top_jets.add("id_M", FatJet.particleNet_TvsQCD < pNet_id_M)
    selections_top_jets.add("id_T", FatJet.particleNet_TvsQCD < pNet_id_T)
    
        
    # define selection regions for each channel
    regions = {  
        "Loose": ["leading_fatjets", "fatjets_electron_dr", "fatjets_muon_dr", "id_L"],
        "Medium": ["leading_fatjets", "fatjets_electron_dr", "fatjets_muon_dr", "id_M"],
        "Tight": ["leading_fatjets", "fatjets_electron_dr", "fatjets_muon_dr", "id_T"],
    }
    
    Filtered_L = selections_top_jets.all(*regions["Loose"])                                 
    selected_idx['jetsTop_L'].append(Filtered_L)
    
    Filtered_M = selections_top_jets.all(*regions["Medium"])                                 
    selected_idx['jetsTop_M'].append(Filtered_M)
    
    Filtered_T = selections_top_jets.all(*regions["Tight"])                                 
    selected_idx['jetsTop_T'].append(Filtered_T)
    
    return selected_idx['jetsTop_L'], selected_idx['jetsTop_M'], selected_idx['jetsTop_T']




def select_w_jets(events: ak.Array, year: str, mod: str, selected_idx):
    
    selections_w_jets = PackedSelection()
 
    FatJet = events.FatJet
    Muon = events.Muon
    Electron = events.Electron
    

    jet_id = 10000 #Recommendation in https://twiki.cern.ch/twiki/bin/view/CMS/JetID#Recommendations_for_13_TeV_data and https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVUL ("Please note: For AK8 jets, the corresponding (CHS or PUPPI) AK4 jet ID should be used.")
    pNet_id = 10000 #Recommendation in https://twiki.cern.ch/twiki/bin/view/CMS/ParticleNetTopWSFs
    
    # Loose: L subindex.
    # Medium: M subindex.
    # Tight: T subindex.
    # Nomenclature: "year": [jetd_id, pNet_id_L, pNet_id_M, pNet_id_T]
    
    constants_as_fuction_wp = {
                                "2016" : [3, 0.68, 0.94, 0.97], 
                                "2016APV" : [3, 0.67, 0.93, 0.97], 
                                "2017": [2, 0.71, 0.94, 0.98], 
                                "2018": [2, 0.70, 0.94, 0.98]
                            }


    jet_id = constants_as_fuction_wp[year + mod][0]
    pNet_id_L = constants_as_fuction_wp[year + mod][1]
    pNet_id_M = constants_as_fuction_wp[year + mod][2]
    pNet_id_T = constants_as_fuction_wp[year + mod][3]
    
   
    good_wjets = (
        (FatJet.pt < 200)
        & (np.abs(FatJet.eta) >= 2.4) #This is the recommendation for all the fat jets (there are not reconstructed forward fat jets)
        & (FatJet.jetId < jet_id)   
    )
    
    n_good_wjets = ak.sum(good_wjets, axis=1) 
    leading_wjets = ak.firsts(FatJet[good_wjets])
    
    
    overlap_wjets_electron = ak.firsts(FatJet).delta_r(ak.firsts(Electron)) # Overlap with electrons
    overlap_wjets_muon = ak.firsts(FatJet).delta_r(ak.firsts(Muon))         # Overlap with muons
    overlap_wjets_fatjet = ak.firsts(FatJet).delta_r(ak.firsts(Muon))       # Overlap with fatjet  (pendiente)
    
   
   
    selections_w_jets.add("leading_wjets", n_good_wjets == 1)
    selections_w_jets.add("wjets_electron_dr", overlap_wjets_electron > 0.8)
    selections_w_jets.add("wjets_muon_dr", overlap_wjets_muon > 0.8)
    selections_w_jets.add("wjets_fatjet_dr", overlap_wjets_fatjet > 0.8)
    selections_w_jets.add("id_L", FatJet.particleNet_TvsQCD < pNet_id_L)
    selections_w_jets.add("id_M", FatJet.particleNet_TvsQCD < pNet_id_M)
    selections_w_jets.add("id_T", FatJet.particleNet_TvsQCD < pNet_id_T)
    

    
    # define selection regions for each channel
    regions = {
        "Loose": ["leading_wjets", "wjets_electron_dr", "wjets_muon_dr", "wjets_fatjet_dr", "id_L"],
        "Medium": ["leading_wjets", "wjets_electron_dr", "wjets_muon_dr", "wjets_fatjet_dr", "id_M"], 
        "Tight": ["leading_wjets", "wjets_electron_dr", "wjets_muon_dr",  "wjets_fatjet_dr", "id_T"],
    }
    
    
    
    Filtered_L = selections_w_jets.all(*regions["Loose"])                                 
    selected_idx['jetsW_L'].append(Filtered_L)
    
    Filtered_M = selections_w_jets.all(*regions["Medium"])                                 
    selected_idx['jetsW_M'].append(Filtered_M)
    
    Filtered_T = selections_w_jets.all(*regions["Tight"])                                 
    selected_idx['jetsW_T'].append(Filtered_T)
    
    return selected_idx['jetsW_L'], selected_idx['jetsW_M'], selected_idx['jetsW_T']



def select_jets(events: ak.Array, year: str, mod: str, selected_idx):
    
    selections_jets = PackedSelection()
 
    FatJet = events.FatJet
    Muon = events.Muon
    Electron = events.Electron
    Jet = events.Jet
    

    pNet_id = 10000 #Recommendation in https://twiki.cern.ch/twiki/bin/view/CMS/ParticleNetTopWSFs
   
    jet_eta = -1
    jet_id = 10000 #Recommendation in https://twiki.cern.ch/twiki/bin/view/CMS/JetID#Recommendations_for_13_TeV_data and https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVUL ("Please note: For AK8 jets, the corresponding (CHS or PUPPI) AK4 jet ID should be used.")
    b_jet_eta = -1
    b_jet_id = 10000 #Recommendation in https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation
        
    # Loose: L subindex.
    # Medium: M subindex.
    # Tight: T subindex.
    # Nomenclature: "year": [jetd_eta, jet_id_T, b_jet_eta, b_jet_id_L, b_jet_id_M, b_jet_id_T]
    
    constants_as_fuction_wp = {
                                "2016" : [2.4, 3, 2.4, 0.0480, 0.2489, 0.6377], 
                                "2016APV" : [2.4, 3, 2.4, 0.0508, 0.2598, 0.6502], 
                                "2017": [2.5, 2, 2.5, 0.0532, 0.3040, 0.7476], 
                                "2018": [2.5, 2, 2.5, 0.0490, 0.2783, 0.7100]
                            } 

    jet_eta = constants_as_fuction_wp[year + mod][0]
    jet_id_T = constants_as_fuction_wp[year + mod][1]
    b_jet_eta = constants_as_fuction_wp[year + mod][2]
    b_jet_id_L = constants_as_fuction_wp[year + mod][3]
    b_jet_id_M = constants_as_fuction_wp[year + mod][4]
    b_jet_id_T = constants_as_fuction_wp[year + mod][5]
    
    good_jets = (
        (FatJet.pt < 20)
        & (np.abs(Jet.eta) > jet_eta) 
        & (Jet.jetId < jet_id)   
    )
      
    n_good_jets = ak.sum(good_jets, axis=1) 
    leading_jets = ak.firsts(Jet[good_jets])
    

    overlap_jets_electron = ak.firsts(Jet).delta_r(ak.firsts(Electron)) # Overlap with electrons
    overlap_jets_muon = ak.firsts(Jet).delta_r(ak.firsts(Muon))         # Overlap with muons
    overlap_jets_fatjet = ak.firsts(Jet).delta_r(ak.firsts(Muon))       # Overlap with fatjet (pendiente)
    overlap_jets_w = ak.firsts(Jet).delta_r(ak.firsts(Electron)) # Overlap with wjets (pendiente)
    overlap_jets_top = ak.firsts(Jet).delta_r(ak.firsts(Muon))         # Overlap with topjet (pendiente)
    
    
    
    selections_jets.add("leading_jets", n_good_jets == 1)
    selections_jets.add("jets_electron_dr", overlap_jets_electron > 0.4)
    selections_jets.add("jets_muon_dr", overlap_jets_muon > 0.4)
    selections_jets.add("jets_fatjet_dr", overlap_jets_fatjet > 0.8)
    selections_jets.add("jets_wjetM_dr", overlap_jets_w > 0.8)
    selections_jets.add("jets_TopM_dr", overlap_jets_top > 0.8)
    
    selections_jets.add("b_jet_eta", np.abs(Jet.eta) > b_jet_eta)

    selections_jets.add("id_L", Jet.btagDeepFlavB < b_jet_id_L)
    selections_jets.add("id_M", Jet.btagDeepFlavB < b_jet_id_M)
    selections_jets.add("id_T", Jet.btagDeepFlavB < b_jet_id_T)
    
    # define selection regions for each channel
    regions = {
        # all the AK4 jets contain boosted jet        
        "noCleaned": ["leading_jets", "jets_electron_dr", "jets_muon_dr"],
        "Cleaned": ["leading_jets", "jets_electron_dr", "jets_muon_dr", "jets_fatjet_dr", "jets_wjetM_dr", "jets_TopM_dr" ],
        "jetsB_L": ["leading_jets", "jets_electron_dr", "jets_muon_dr", "jets_fatjet_dr", "jets_wjetM_dr", "jets_TopM_dr", "b_jet_eta", "id_L"],
        "jetsB_M": ["leading_jets", "jets_electron_dr", "jets_muon_dr", "jets_fatjet_dr", "jets_wjetM_dr", "jets_TopM_dr", "b_jet_eta", "id_M"],
        "jetsB_T": ["leading_jets", "jets_electron_dr", "jets_muon_dr", "jets_fatjet_dr", "jets_wjetM_dr", "jets_TopM_dr", "b_jet_eta", "id_T"],   
    }
    
    #No cleaned
    Filtered_NC = selections_jets.all(*regions["noCleaned"])                                 
    selected_idx['jets_noCleaned_against_boostedJets'].append(Filtered_NC)
    
    Filtered_C = selections_jets.all(*regions["Cleaned"])                                 
    selected_idx['jets_cleaned'].append(Filtered_C)
    
    Filtered_L = selections_jets.all(*regions["jetsB_L"])                                 
    selected_idx['jetsB_L'].append(Filtered_L)
    
    Filtered_M = selections_jets.all(*regions["jetsB_M"])                                 
    selected_idx['jetsB_M'].append(Filtered_M)
    
    Filtered_T = selections_jets.all(*regions["jetsB_T"])                                 
    selected_idx['jetsB_T'].append(Filtered_T)



