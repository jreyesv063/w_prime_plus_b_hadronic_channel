import hist
import pandas as pd
import numpy as np
import awkward as ak
from coffea.analysis_tools import PackedSelection
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema



 selected_object = {'muon': [], 'electron': [],                       # Leptons
                    'jetsTop': [],                                    # FatJet
                    'jetsW': [],                                      # WJet
                    'jets_noCleaned_against_boostedJets': [],         # No overlap with Electrons and Muons
                    'jets_cleaned': [],                               # No overlap with FatJet, WJet, Electrons and Muons
                    'jetsB': []                                       # BJet
                   }
    
selections = PackedSelection()



# ------------
# Electrons
# ------------

good_electrons = (
    (events.Electron.pt >= 30)
    & (np.abs(events.Electron.eta) < 2.4)
    & ((np.abs(events.Electron.eta) < 1.44) | (np.abs(events.Electron.eta) > 1.57))
    & (events.Electron.cutBased >= 1)
)

leading_electron = ak.firsts(events.Electron[good_electrons])

# Events
n_good_electrons = ak.sum(good_electrons, axis=1)
electrons = events.Electron[(good_electrons & n_good_electrons == 1)]


# Selections
selections.add("one_electron", n_good_electrons == 1)
    
    
    

# ------------
# Muons
# ------------

good_muons = (
    (events.Muon.pt >= 20)
    & (np.abs(events.Muon.eta) < 2.4)
    & (events.Muon.tightId)
    & (events.Muon.pfRelIso04_all > 0.15)
)

# Check that muons don't overlap with electrons (DeltaR > 0.4)
leading_muon = ak.firsts(events.Muon[good_muons])
overlap_muon_electron = leading_muon.delta_r(leading_electron) > 0.4


# Events
n_good_muons = ak.sum(good_muons, axis=1)
muons = events.Muon[(good_muons & overlap_muon_electron & n_good_muons == 1)]


# Selections
selections.add("one_muon", n_good_muons == 1)
selections.add("muon_electron_dr", overlap_muon_electron)



# -----
# Fat Jets
# -----

# Loose: L subindex.
# Medium: M subindex.
# Tight: T subindex.
# Nomenclature: "year": [jetd_id, pNet_id_L, pNet_id_M, pNet_id_T]
    
constants_as_fuction_wp_fatjet = {
                            "2016" : [3, 0.50, 0.73, 0.96], 
                            "2016APV" : [3, 0.49, 0.74, 0.96], 
                            "2017": [2, 0.58, 0.80, 0.97], 
                            "2018": [2, 0.58, 0.80, 0.97]
                        }


jet_id = constants_as_fuction_wp_fatjet["2017"][0]
pNet_id = constants_as_fuction_wp_fatjet["2017"][2]


# Select good fat jets
good_fatjets = (
    (events.FatJet.pt >= 300)
    & (np.abs(events.FatJet.eta) <= 2.4)
    & (events.FatJet.particleNet_TvsQCD >= pNet_id)   # Medium   (Top vs QCD)
    & (events.FatJet.jetId >= jet_id)   # Medium
)

leading_top = ak.firsts(events.FatJet[good_fatjets])
# check that top don't overlap with electrons and muons (DeltaR > 0.8)
overlap_top_electron = leading_top.delta_r(leading_electron) > 0.8
overlap_top_muon = leading_top.delta_r(leading_muon) > 0.8 


# Events
n_good_fatjets = ak.sum(good_fatjets, axis=1)
tops = events.FatJet[(good_fatjets & overlap_top_electron & overlap_top_muon)]



# Selections
selections.add("one_top", n_good_fatjets == 1)
selections.add("top_electron_dr", overlap_top_electron)
selections.add("top_muon_dr", overlap_top_muon)


# -----
# Top Jets
# -----
    
# Loose: L subindex.
# Medium: M subindex.
# Tight: T subindex.
# Nomenclature: "year": [jetd_id, pNet_id_L, pNet_id_M, pNet_id_T]
    
constants_as_fuction_wp_wjet = {
                            "2016" : [3, 0.68, 0.94, 0.97], 
                            "2016APV" : [3, 0.67, 0.93, 0.97], 
                            "2017": [2, 0.71, 0.94, 0.98], 
                            "2018": [2, 0.70, 0.94, 0.98]
                        }

jet_id = constants_as_fuction_wp_wjet["2017"][0]
pNet_id = constants_as_fuction_wp_wjet["2017"][2]

good_wjets = (
    (events.FatJet.pt >= 200)
    & (np.abs(events.FatJet.eta) <= 2.4) #This is the recommendation for all the fat jets (there are not reconstructed forward fat jets)
    & (events.FatJet.jetId >= jet_id)   
    & (events.FatJet.particleNet_WvsQCD >= pNet_id)   # Medium  (W vs QCD)
)


leading_wjet = ak.firsts(events.FatJet[good_wjets])
# Leading wjet don't overlap with electrons, muons (DeltaR > 0.8) and topJets (Delta R > 0.8) 
overlap_w_electron = leading_wjet.delta_r(leading_electron) > 0.8
overlap_w_muon = leading_wjet.delta_r(leading_muon) > 0.8 
overlap_w_top = leading_wjet.delta_r(leading_top) > 0.8 

# Events
n_good_wjets = ak.sum(good_wjets, axis=1) 
wjets = events.FatJet[(good_wjets & overlap_w_electron & overlap_w_muon & overlap_w_top & n_good_wjets == 1)]


# Selections
selections.add("one_w", n_good_wjets == 1) 
selections.add("two_w", n_good_wjets == 2)   # we have two tops, so we will have 2 W (see decay squema)
selections.add("w_electron_dr", overlap_w_electron)
selections.add("w_muon_dr", overlap_w_muon)
selections.add("w_top_dr", overlap_w_top)



# -----
# Jets and bjets
# -----
    
# Nomenclature: "year": [jetd_eta, jet_id, b_jet_id_L, b_jet_id_M, b_jet_id_T]

constants_as_fuction_wp_jet = {
                            "2016" : [2.4, 3,  0.0480, 0.2489, 0.6377], 
                            "2016APV" : [2.4, 3, 0.0508, 0.2598, 0.6502], 
                            "2017": [2.5, 2, 0.0532, 0.3040, 0.7476], 
                            "2018": [2.5, 2, 0.0490, 0.2783, 0.7100]
                        } 



jet_eta = constants_as_fuction_wp_jet["2017"][0]
jet_id = constants_as_fuction_wp_jet["2017"][1]


good_jets = (
        (events.Jet.pt >= 20)
        & (np.abs(events.Jet.eta) <= jet_eta) 
        & (events.Jet.jetId >= jet_id)   
)     



bjet_wp = constants_as_fuction_wp_jet["2017"][3]   # wp deppFlavour (medium)

good_bjets = (
                events.Jet.btagDeepFlavB >= bjet_wp
            )


leading_jet = ak.firsts(events.Jet[good_jets])
# Leading jet don't overlap with electrons, muons (DeltaR > 0.4), topJets and jets (Delta R > 0.8) 
overlap_jet_electron = leading_jet.delta_r(leading_electron) > 0.4
overlap_jet_muon = leading_jet.delta_r(leading_muon) > 0.4
overlap_jet_top = leading_jet.delta_r(leading_top) > 0.8
overlap_jet_w = leading_jet.delta_r(leading_wjet) > 0.8



# Events
n_good_jets = ak.sum(good_jets, axis=1) 
n_good_bjets = ak.sum(good_bjets, axis=1) 
jets_NC = events.Jet[(good_jets & overlap_jet_electron & overlap_jet_muon)]
jets = events.Jet[(good_jets & overlap_jet_electron & overlap_jet_muon & overlap_jet_top & overlap_jet_w)]
bjets = events.Jet[(good_jets & good_bjets & overlap_jet_electron & overlap_jet_muon & overlap_jet_top & overlap_jet_w & n_good_bjets == 1)]



# Selections
selections.add("two_jets", n_good_jets == 2)   # When one top decay we have two jets (see decay squema)
selections.add("one_bjet", n_good_bjets == 1)   
selections.add("jet_electron_dr", overlap_jet_electron)
selections.add("jet_muon_dr", overlap_jet_muon)
selections.add("jet_top_dr", overlap_jet_top)
selections.add("jet_wjet_dr", overlap_jet_w)



def number_entries(object):  
   
    num_none = ak.count_nonzero(ak.num(object) == 0)  # Number of empty entries.
    total = len(object)
    n = total - num_none
    
    return "passed the criteria = {}".format(n)



objects = {
        "electron": ["one_electron"],
        "muon": ["one_muon", "muon_electron_dr"],
        "jetsTop": ["one_top", "top_electron_dr", "top_muon_dr"],
        "jetsW": ["one_w", "w_electron_dr", "w_muon_dr", "w_top_dr"],
        "jets_NC": ["two_jets", "jet_electron_dr", "jet_muon_dr"],
        "jets_cleaned": ["two_jets", "jet_electron_dr", "jet_muon_dr", "jet_top_dr", "jet_wjet_dr"],
        "jetsB": ["one_bjet", "jet_electron_dr", "jet_muon_dr", "jet_top_dr", "jet_wjet_dr"]
    }

print("Total events {}".format(len(events)))

selected_object['electron'].append(electrons)
selected_object['electron'].append(selections.all(*objects["electron"]))
#print("Number of electrons that " + number_entries(selected_object['electron'][0]))


selected_object['muon'].append(muons)
selected_object['muon'].append(selections.all(*objects["muon"]))
#print("Number of muons that " + number_entries(selected_object['muon'][0]))


selected_object['jetsTop'].append(tops)
selected_object['jetsTop'].append(selections.all(*objects["jetsTop"]))
#print("Number of tops that " + number_entries(selected_object['jetsTop'][0]))


selected_object['jetsW'].append(wjets)
selected_object['jetsW'].append(selections.all(*objects["jetsW"]))
#print("Number of ws that " + number_entries(selected_object['jetsW'][0]))


selected_object['jets_noCleaned_against_boostedJets'].append(jets_NC)
selected_object['jets_noCleaned_against_boostedJets'].append(selections.all(*objects["jets_NC"]))
#print("Number of jets_NC that " + number_entries(selected_object['jets_noCleaned_against_boostedJets'][0]))


selected_object['jets_cleaned'].append(jets)
selected_object['jets_cleaned'].append(selections.all(*objects["jets_cleaned"]))
#print("Number of jets that " + number_entries(selected_object['jets_cleaned'][0]))


selected_object['jetsB'].append(bjets)
selected_object['jetsB'].append(selections.all(*objects["jetsB"]))
#print("Number of bjets that " + number_entries(selected_object['jetsB'][0]))


# --------------------
# TLorentz 
# --------------------
def TLorentzVector(E, px, py, pz):
    
    TLorentz = np.array([E, px, py, pz])
   
    return TLorentz


# -------------------
# calculate chi2
# -------------------

def chi2(top_mass, w_mass, masses_sigmas_chi2, pdg_masses):
    
    # PDG values
    top_mass_pdg = pdg_masses[0]
    w_mass_pdg = pdg_masses[1]
    
    
    # Sigmas
    top_sigma = masses_sigmas_chi2["Partially-boosted top"][0]
    w_sigma = masses_sigmas_chi2["Partially-boosted top"][3]
    
    
    t = (top_mass - top_mass_pdg) / top_sigma
    w = (w_mass - w_mass_pdg) / w_sigma
    
    
    chi2 = t**2 + w**2 

    
    return chi2


# ---------------
# add_results
# ---------------
def add_results(chi2, top_AK4jets_list, w, b, all_combinations_results, topology):
    all_combinations_results['chi2'].append(chi2)
    all_combinations_results['top_AK4jets'].append(top_AK4jets_list)
    all_combinations_results['w'].append(w)
    all_combinations_results['top'].append(w + b)
    all_combinations_results['top_topology_decay'].append(topology)



# PDG values

# Nomenclature: [top_mass_pdg, w_mass_pdg, muon_mass_pdg, electron_mass_pdg]

pdg_masses = [173.1, 80.4, 0.10566, 0.511e-3]


# Nomenclature: [top_sigma, top_low, top_down, w_sigma, w_low, w_up, chi2]
masses_sigmas_chi2 = {
                        "boosted top" : [37.14, 100, 400],
                        "Partially-boosted top" : [37.14, 100, 300, 20.09, 37.14, 20.09, 5],
                        "Resolved hadronic top" : [35.02, 100, 300, 24.98, 40, 200, 5],
                        "Resolved leptonic top" : [120.8, 100, 600, 127.2, 40, 600, 5]
                    }


# Results will be stored in results{}
results ={
          'top':[], 
          'w':[], 
          'top_topology_decay':[], 
          'chi2':[]
         }


# reconstruct the boosted top
def boosted_top(jetsTop, 
                masses_sigmas_chi2,
                results):
                    
        
    top_mass_low = masses_sigmas_chi2["boosted top"][1]
    top_mass_up =  masses_sigmas_chi2["boosted top"][2]
    

    good_top_mass = (
        (jetsTop.mass > top_mass_low)
        & (jetsTop.mass < top_mass_up)
    )
    
    
    boosted_top = jetsTop[good_top_mass]   
    
    results['top'].append(boosted_top)
    results['w'].append(TLorentzVector(0,0,0,0))    # append(ROOT.TLorentzVector(0, 0, 0, 0))  # https://coffeateam.github.io/coffea/api/coffea.nanoevents.methods.vector.LorentzVector.html
    results['top_topology_decay'].append(0)
    results['chi2'].append(0.0000001)
