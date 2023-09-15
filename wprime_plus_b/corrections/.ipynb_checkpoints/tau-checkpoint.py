import correctionlib
import numpy as np
import awkward as ak
import importlib.resources
from typing import Type
from coffea.analysis_tools import Weights
from wprime_plus_b.corrections.utils import pog_years, get_pog_json

# ----------------------------------
# Tau scale factors
"""
TauID corrections

https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendationForRun2
https://github.com/cms-tau-pog/TauIDSFs
https://github.com/uhh-cms/hh2bbtautau/blob/7666ed0426c87baa8d143ec26a216c3cadde513b/hbt/calibration/tau.py#L60

Good example:
https://github.com/schaefes/hh2bbtautau/blob/da6d47a7ddb2b1e7ffda06b8a96c6ddead2824b8/hbt/production/tau.py#L108


* DeepTau2017v2p1VSe = eta [0.0, 2.3); genmatch (0, 1); wp (Loose, Medium, Tight, VLoose, VTight, VVLoose, VVTight; syst (down, nom, up)

* DeepTau2017v2p1VSmu (eta [0.0, 2.3); genmatch (0, 2); wp (Loose, Medium, Tight, VLoose); syst (down, nom, up))

* DeepTau2017v2p1VSjet (pt [-inf, inf); dm (0, 1, 2, 10, 11) ; genmatch (0, 1, 2, 3, 4, 5, 6 ); wp (Loose, Medium, Tight, VTight); wp_VSe (Tight, VVLoose );  syst ; flag (dm, pt))
        
* tau_energy_scale (pt [-inf, inf), eta [0.0, 2.5), dm (0, 1, 2, 10, 11), genmatch (1, 2, 5, 6), id (DeepTau2017v2p1), syst (down, nom, up))



"""
# -----------------------------------
class TauCorrector:
    def __init__(
        self,
        taus: ak.Array,
        weights: Type[Weights],
        year: str = "2017",
        year_mod: str = "",
        tag: str = "tau",
        variation: str = "nom"
    ) -> None:

        self.variation = variation
        
        # tau array
        self.taus = taus

        # tau transverse momentum and pseudorapidity
        self.tau_pt = np.array(ak.fill_none(self.taus.pt, 0.0))
        self.tau_eta = np.array(ak.fill_none(self.taus.eta, 0.0))
        self.tau_dm = ak.to_numpy(self.taus.decayMode, allow_missing=True)
        self.tau_genPart = ak.to_numpy(self.taus.genPartFlav, allow_missing=True)

        # weights container
        self.weights = weights

        # define correction set_id
        self.cset = correctionlib.CorrectionSet.from_file(
            get_pog_json(json_name="tau", year=year + year_mod)
        )
                        
        
        self.year = year
        self.year_mod = year_mod
        self.pog_year = pog_years[year + year_mod]

        self.tag = tag


    def add_id_weight_DeepTau2017v2p1VSe(self, working_point: str = "Tight", systematic: str = "nom"):
        # tau gen particle. We only need to consider values: 1, 3. 0 is unmached 
        e_mask = (
            (self.tau_genPart == 1) | 
            (self.tau_genPart == 3)
        )
        
        # tau pseudorapidity range: [0.0, 2.3)
        tau_eta = np.clip(self.tau_eta.copy(), 0.0, 2.3)
        
        
        # genmatch
        tau_gen = ak.fill_none(ak.mask(self.tau_genPart, e_mask),0)
        
        
        # get scale factors
        values = {}    
        
        
        """
        Sf is called with:
        
        evaluate(eta (real),  genmatch (int) , wp (string), syst (string))
        
        """       
        values["nominal"] = self.cset["DeepTau2017v2p1VSe"].evaluate(tau_eta, 
                                                                     tau_gen, 
                                                                     "Medium", 
                                                                     "nom")
        
        # -------------------
        # Systematics        
        # -------------------    
        if self.variation == "nominal":
            values["up"] = self.cset["DeepTau2017v2p1VSe"].evaluate(tau_eta, 
                                                                    tau_gen, 
                                                                    "Medium", 
                                                                    "up")
            
            values["down"] = self.cset["DeepTau2017v2p1VSe"].evaluate(tau_eta, 
                                                                      tau_gen, 
                                                                      "Medium", 
                                                                      "down")
            
            # add scale factors to weights container
            self.weights.add(
                name=f"{self.tag}_id",
                weight=values["nominal"],
                weightUp=values["up"],
                weightDown=values["down"]
            )
        else:
            self.weights.add(
                name=f"{self.tag}_id",
                weight=values["nominal"]
            )
        
        
    def add_id_weight_DeepTau2017v2p1VSmu(self, working_point: str = "Tight", systematic: str = "nom"):
        # tau gen particle. We only need to consider values: 2, 4. 0 is unmached 
        mu_mask = (
            (self.tau_genPart == 2) | 
            (self.tau_genPart == 4)
        )

        # tau pseudorapidity range: [0.0, 2.3)
        tau_eta = np.clip(self.tau_eta.copy(), 0.0, 2.3)
        
        
        # genmatch
        tau_gen = ak.fill_none(ak.mask(self.tau_genPart, mu_mask),0)
            
        
        # get scale factors
        values = {}    
        
        
        """
        Sf is called with:
        
        evaluate(eta (real),  genmatch (int) , wp (string), syst (string))
        
        """       
        values["nominal"] = self.cset["DeepTau2017v2p1VSmu"].evaluate(tau_eta, 
                                                                      tau_gen, 
                                                                      working_point, 
                                                                      "nom")

        # -------------------
        # Systematics        
        # -------------------    
        if self.variation == "nominal":
            values["up"] = self.cset["DeepTau2017v2p1VSmu"].evaluate(tau_eta, 
                                                                     tau_gen, 
                                                                     working_point, 
                                                                     "up")
            
            values["down"] = self.cset["DeepTau2017v2p1VSmu"].evaluate(tau_eta, 
                                                                       tau_gen, 
                                                                       working_point, 
                                                                       "down")
            
            # add scale factors to weights container
            self.weights.add(
                name=f"{self.tag}_id",
                weight=values["nominal"],
                weightUp=values["up"],
                weightDown=values["down"]
        )
        else:
            self.weights.add(
                name=f"{self.tag}_id",
                weight=values["nominal"]
        )
            
        
    def add_id_weight_DeepTau2017v2p1VSjet(self, 
                                           working_point: str = "Tight", 
                                           working_point_VSe: str = "Tight", 
                                           systematic: str = "nom", 
                                           flag: str = "pt"
                                          ):
        
        # tau gen particle. We only need to consider values: 5 in genmatch. 0 is unmached
        # For dm, the possible values will be 0, 1, 2, 10, 11   
        tau_mask_gm = (self.tau_genPart == 5)
        tau_mask_dm = (
            (self.tau_dm == 0) |
            (self.tau_dm == 1) |
            (self.tau_dm == 2) |
            (self.tau_dm == 10) |
            (self.tau_dm == 11)
        )

        # tau pt
        tau_pt = self.tau_pt
        

        # tau decay mode
        tau_dm = ak.fill_none(ak.mask(self.tau_dm, tau_mask_dm),0)
        
        
        # genmatch
        tau_gen = ak.fill_none(ak.mask(self.tau_genPart, tau_mask_gm),0)

       
        # get scale factors
        values = {}
        
        
        """
        https://github.com/LEAF-HQ/LEAF/blob/d22cc55594a4b16d061c25dbf7ecdec04eedbc34/Analyzer/src/TauScaleFactorApplicatorJson.cc#L28
        
        Sf is called with:
        
        evaluate(pt (real),  dm (int), genmatch (int), wp (string), wp_VSe (string), syst (string), flag (string))
        
         - dm (decay mode): 0 (tau->pi); 1 (tau->rho->pi+pi0); 2 (tau->a1->pi+2pi0); 10 (tau->a1->3pi); 11 (tau->3pi+pi0)
         - getmatch: 0 or 6 = unmatched or jet, 1 or 3 = electron, 2 or 4 = muon, 5 = real tau
         - flag: We have worked in 'pt' = pT-dependent
         
        """
  
        values["nominal"] = self.cset["DeepTau2017v2p1VSjet"].evaluate(tau_pt, 
                                                                       tau_dm, 
                                                                       tau_gen, 
                                                                       working_point, 
                                                                       working_point_VSe, 
                                                                       "nom", 
                                                                       flag)
        
        # -------------------
        # Systematics        
        # -------------------    
        if self.variation == "nominal":
            values["up"] = self.cset["DeepTau2017v2p1VSjet"].evaluate(tau_pt, 
                                                                      tau_dm, 
                                                                      tau_gen, 
                                                                      working_point, 
                                                                      working_point_VSe, 
                                                                      "up", 
                                                                      flag)
            
            values["down"] = self.cset["DeepTau2017v2p1VSjet"].evaluate(tau_pt, 
                                                                        tau_dm, 
                                                                        tau_gen, 
                                                                        working_point, 
                                                                        working_point_VSe, 
                                                                        "down", 
                                                                        flag)
        
            
            # add scale factors to weights container
            self.weights.add(
                name=f"{self.tag}_id",
                weight=values["nominal"],
                weightUp=values["up"],
                weightDown=values["down"]
        )
        else:
            self.weights.add(
                name=f"{self.tag}_id",
                weight=values["nominal"]
        )


    def add_id_weight_tau_energy_scale(self, id: str = "DeepTau2017v2p1" , systematic: str = "nom"):
        mask_gm = (
            (self.tau_genPart == 1) | 
            (self.tau_genPart == 2) | 
            (self.tau_genPart == 5) | 
            (self.tau_genPart == 6)
        )
        mask_dm = (
            (self.tau_dm == 0) |
            (self.tau_dm == 1) |
            (self.tau_dm == 2) |
            (self.tau_dm == 10) |
            (self.tau_dm == 11)
        )
        
        # tau pt
        tau_pt = self.tau_pt
        
        # tau pseudorapidity range: [0.0, 2.3)
        tau_eta = np.clip(self.tau_eta.copy(), 0.0, 2.5)
        
        # tau decay mode
        tau_dm = ak.fill_none(ak.mask(self.tau_dm, mask_dm),0)
        
        
        # genmatch
        tau_gen = ak.fill_none(ak.mask(self.tau_genPart, mask_gm),0)
           
        
        # get scale factors
        values = {}
        
                
        """
        Sf is called with:
        
        evaluate(pt (real); eta (real);  dm (int); genmatch (int); id (string); syst (string))
        
        """
        
        values["nominal"] = self.cset["tau_energy_scale"].evaluate(tau_pt, 
                                                                   tau_eta, 
                                                                   tau_dm, 
                                                                   tau_gen, 
                                                                   id, 
                                                                   "nom")

        # -------------------
        # Systematics        
        # -------------------    
        if self.variation == "nominal":
            values["up"] = self.cset["tau_energy_scale"].evaluate(tau_pt, 
                                                                  tau_eta, 
                                                                  tau_dm, 
                                                                  tau_gen, 
                                                                  id, 
                                                                  "up")
            
            values["down"] = self.cset["tau_energy_scale"].evaluate(tau_pt, 
                                                                    tau_eta, 
                                                                    tau_dm, 
                                                                    tau_gen, 
                                                                    id, 
                                                                    "down")
        
            
            # add scale factors to weights container
            self.weights.add(
                name=f"{self.tag}_id",
                weight=values["nominal"],
                weightUp=values["up"],
                weightDown=values["down"]
        )
        else:
            self.weights.add(
                name=f"{self.tag}_id",
                weight=values["nominal"]
        )        
        
        
 