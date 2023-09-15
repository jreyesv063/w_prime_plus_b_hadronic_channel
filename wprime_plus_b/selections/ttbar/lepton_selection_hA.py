import numpy as np
import awkward as ak
from coffea.nanoevents.methods.base import NanoEventsArray


def select_good_electrons(events: NanoEventsArray) -> ak.highlevel.Array:
    
    return (
                (events.Electron.pt >= 30)
                & (np.abs(events.Electron.eta) <= 2.4)
                & (
                    (np.abs(events.Electron.eta) <= 1.44)
                    | (np.abs(events.Electron.eta) >= 1.57)
                )
                & (events.Electron.cutBased >= 1)
                & (events.Electron.mvaFall17V2Iso_WP90)
            )
    


def select_good_muons(events: NanoEventsArray) -> ak.highlevel.Array:
    
    return (
                (events.Muon.pt >= 35)
                & (np.abs(events.Muon.eta) <= 2.4)
                & (events.Muon.tightId)
                & (events.Muon.pfRelIso04_all <= 0.15)
            )
    
    

def select_good_taus(events: NanoEventsArray) -> ak.highlevel.Array:
    """
    https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePFTauID
    https://github.com/cms-tau-pog/TauIDSFs
    
    The number of prongs remains to be applied.
    
    Example: https://github.com/schaefes/hh2bbtautau/blob/da6d47a7ddb2b1e7ffda06b8a96c6ddead2824b8/hbt/production/tau.py#L87   
    
    Recommendations: pT > 20 GeV, |eta| < 2.3, |dz| < 0.2
    https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendationForRun2
    
    """
    return (
                (events.Tau.pt >= 20)
                & (np.abs(events.Tau.eta) <= 2.3)
                # github.com/bonanomi/hh2bbww/blob/388efda4e9a6a207d4e983c7f7528acb3a4c374f/hbw/selection/default.py#L295
                & (events.Tau.idDeepTau2017v2p1VSjet >= 8)  # VVVLoose,VVLoose,VLoose,Loose,Medium,Tight,VTight,VVTight
                & (events.Tau.idDeepTau2017v2p1VSe >= 8)    # VVVLoose,VVLoose,VLoose,Loose,Medium,Tight,VTight,VVTight
                & (events.Tau.idDeepTau2017v2p1VSmu >= 1)   # VLoose,Loose,Medium,Tight
                & (np.abs(events.Tau.dz) < 0.2)
                & (events.Tau.idDecayModeNewDMs)    
                & (
                    (events.Tau.decayMode == 0)             # 0 (tau->pi)
                    | (events.Tau.decayMode == 1)           # 1 (tau->rho->pi+pi0)
                    | (events.Tau.decayMode == 1)           # 2 (tau->a1->pi+2pi0)
                    | (events.Tau.decayMode == 10)          # 10 (tau->a1->3pi)
                    | (events.Tau.decayMode == 11)          # 11 (tau->3pi+pi0)
                )
            )
