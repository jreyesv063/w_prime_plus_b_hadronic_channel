import json
import numpy as np
import awkward as ak
import importlib.resources
from coffea.nanoevents.methods.nanoaod import JetArray
from coffea.nanoevents.methods.base import NanoEventsArray


"""
R determines the size of the clustering cone, which in CMS is either 0.4 (AK4 jets) or 0.8 (AK8 jets)
https://ml4physicalsciences.github.io/2019/files/NeurIPS_ML4PS_2019_53.pdf
https://indico.cern.ch/event/608530/contributions/2464546/attachments/1416703/2169289/LHCC_vlq2ht_tholen.pdf

"""

def select_good_Fatjets(
    events: NanoEventsArray, 
    year="2017", 
    working_point="T") -> ak.highlevel.Array:
    
    """
    ak8 PFJets with Puppi and JECs applied, after pT > 175 GeV are stored
    
    """
    # Wps top tagger, jet_id and jet_eta
    with open("wprime_plus_b/jsons/topWps.json", "r") as f: 
        Wps = json.load(f)
 
    pNet_id = Wps[year]["TvsQCD"][working_point]                      
    jet_id = Wps[year]['jet_id']  
    
    return (
                (events.FatJet.pt >= 300)
                & (np.abs(events.FatJet.eta) <= 2.4)
                & (events.FatJet.particleNet_TvsQCD >= pNet_id)   # Top vs QCD (tight)
                & (events.FatJet.jetId >= jet_id)                 # tight ID
            )
  

def select_good_Wjets(
    events: NanoEventsArray, 
    year="2017", 
    working_point="T") -> ak.highlevel.Array:
    """
    ak8 PFJets with Puppi and JECs applied, after pT > 175 GeV are stored
    
    """
    # Wps top tagger, jet_id and jet_eta
    with open("wprime_plus_b/jsons/topWps.json", "r") as f: 
        Wps = json.load(f)
        
    pNet_id = Wps[year]["WvsQCD"][working_point]                      
    jet_id = Wps[year]['jet_id']  
    
        
    return (
                (events.FatJet.pt >= 200)
                & (np.abs(events.FatJet.eta) <= 2.4)
                & (events.FatJet.particleNet_WvsQCD >= pNet_id)   # W vs QCD (tight) 
                & (events.FatJet.jetId >= jet_id)   
            )
   

    
def select_good_jets(
    jets, 
    year="2017") -> ak.highlevel.Array:
    """
    More information about the Jet flags:
    https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookNanoAOD
    
    ak4 jets
    
    """            
    # Wps top tagger, jet_id and jet_eta
    with open("wprime_plus_b/jsons/topWps.json", "r") as f: 
        Wps = json.load(f)

        
    jet_eta = Wps[year]['jet_eta']
    jet_id = Wps[year]['jet_id']

    return (
                (jets.pt >= 30)           
                & (np.abs(jets.eta) <= jet_eta) 
                & (jets.jetId >= jet_id)                    #  pass tight and tightLepVeto ID.         
                & (jets.puId == 7)                          #  pass loose, medium, tight ID
            )   
    

def select_good_bjets_prev(
    jets,
    year: str = "2017",
    btag_working_point: str = "M",
) -> ak.highlevel.Array:
    """
    Selects and filters 'good' b-jets from a collection of jets based on specified criteria

    Parameters:
    -----------
    events:
        A collection of events represented using the NanoEventsArray class.

    year: {'2016', '2017', '2018'}
        Year for which the data is being analyzed. Default is '2017'.

    btag_working_point: {'L', 'M', 'T'}
        Working point for b-tagging. Default is 'M'.

    jet_id: https://twiki.cern.ch/twiki/bin/view/CMS/JetID#Run_II
        Jet ID flags {1, 2, 3, 6, 7}
        For 2016 samples:
            1 means: pass loose ID, fail tight, fail tightLepVeto
            3 means: pass loose and tight ID, fail tightLepVeto
            7 means: pass loose, tight, tightLepVeto ID.
        For 2017 and 2018 samples:
            2 means: pass tight ID, fail tightLepVeto
            6 means: pass tight and tightLepVeto ID.

    jet_pileup_id: https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJetID
        Pileup ID flags for pre-UL trainings {0, 4, 6, 7}. Should be applied only to AK4 CHS jets with pT < 50 GeV
        0 means 000: fail all PU ID;
        4 means 100: pass loose ID, fail medium, fail tight;
        6 means 110: pass loose and medium ID, fail tight;
        7 means 111: pass loose, medium, tight ID.

    Returns:
    --------
        An Awkward Array mask containing the selected "good" b-jets that satisfy the specified criteria.
    """
    # open and load btagDeepFlavB working point
    with importlib.resources.open_text("wprime_plus_b.data", "btagWPs.json") as file:
        btag_threshold = json.load(file)["deepJet"][year][btag_working_point]


    return (
        (jets.pt > 20)
        & (np.abs(jets.eta) < 2.4)
        & (jets.jetId == 6)
        & (jets.puId == 7)
        & (jets.btagDeepFlavB > btag_threshold)      
    )
    
    