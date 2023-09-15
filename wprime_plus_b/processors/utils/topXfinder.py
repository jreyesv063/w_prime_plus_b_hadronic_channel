import numpy as np
import awkward as ak
import json
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
from wprime_plus_b.processors.utils.analysis_utils import chi2_test, delta_r_mask, entries, true_sum

class topXfinder:
    def __init__(
        self,
        year: str = "2017",
        year_mod: str = "",
        b_wp: str = "M",
    ) -> None:
        

        deepJet = {
                "2016APV": {
                            "L": 0.0508,
                            "M": 0.2598,
                            "T": 0.6502
                        },
        "2016": {
            "L": 0.048,
            "M": 0.2489,
            "T": 0.6377
        },
        "2017": {
            "L": 0.0532,
            "M": 0.304,
            "T": 0.7476
        },
        "2018": {
            "L": 0.049,
            "M": 0.2783,
            "T": 0.71
        }
        }
       
        # year
        self.year = year
        self.year_mod = year_mod
        self._btag_wp = deepJet["2017"][b_wp]
        

        
    ######################################
    #########  1 Jet #####################
    ######################################
    # ----------------------------------------
    #  Scenario I (FatJet = 1)
    # ----------------------------------------
    def Scenario_I(self, tops, top_low_mass, top_up_mass):

        jet_top = tops.mask[ak.num(tops) == 1]

        good_top_mass = (
            (jet_top.mass > top_low_mass) & 
            (jet_top.mass < top_up_mass)
        )

        
        top_scenario_I = jet_top.mask[good_top_mass]  
        top_I_mass = ak.fill_none(ak.firsts(top_scenario_I.mass), 0)
        mask_scenario_I = ak.fill_none(ak.firsts(good_top_mass), False)


        #print("Top tagger_mask_I (True): ", ak.sum(mask_scenario_I), mask_scenario_I.ndim, mask_scenario_I[mask_scenario_I == True])

        
        return top_scenario_I, mask_scenario_I, top_I_mass, ak.sum(mask_scenario_I)
    
    
    ######################################
    #########  2 Jet #####################
    ######################################
        
    # ----------------------------------------
    #  Scenario II (FatJet = 1  + b = 1)
    # ----------------------------------------
    def Scenario_II(self, tops, bjets, 
                    top_low_mass, top_up_mass):
        
        jet_top = tops.mask[ak.num(tops) == 1]
        jet_b = bjets.mask[ak.num(bjets) == 1]
          

        good_top_mass = (
            (jet_top.mass > top_low_mass) & 
            (jet_top.mass < top_up_mass) & 
            (delta_r_mask(jet_top, jet_b, threshold=0.8))
        )    

    
        top_scenario_II = jet_top.mask[good_top_mass]
        top_II_mass = ak.fill_none(ak.firsts(top_scenario_II.mass), 0)
        mask_scenario_II = ak.fill_none(ak.firsts(good_top_mass), False)

        #print("Top tagger_mask_II (True): ", ak.sum(mask_scenario_II), mask_scenario_II.ndim, mask_scenario_II[mask_scenario_II == True])
        
        
        return top_scenario_II, mask_scenario_II, top_II_mass, ak.sum(mask_scenario_II)
    
    # ----------------------------------------
    #  Scenario III (W = 1  + b = 1)
    # ----------------------------------------
    def Scenario_III(self, wjets, bjets, 
                     top_sigma, top_low_mass, top_up_mass, top_pdg, 
                     w_sigma, w_low_mass, w_up_mass, w_pdg, 
                     chi2):
        
        jet_b = bjets.mask[ak.num(bjets) == 1]
        jet_w = wjets.mask[ak.num(wjets) == 1]
       

        good_w_mass = (
            (jet_w.mass > top_low_mass) &
            (jet_w.mass < top_up_mass)  &
            (delta_r_mask(jet_b, jet_w, threshold=0.8))
        )

        jet_top = jet_w.mask[good_w_mass] + jet_b
    
    
        good_top_mass = (
            (jet_top.mass >  w_low_mass) &
            (jet_top.mass < w_up_mass)
        )    
    
        tops = jet_top.mask[good_top_mass]
        w = jet_top.mask[good_top_mass]
        
        # ----------------------
        # chi2 criteria
        # ----------------------
        chi2_cal = chi2_test(tops,  w, 
                            top_sigma, 
                            w_sigma, 
                            top_pdg, 
                            w_pdg)

        good_chi2 = (chi2_cal < chi2)
            
            
        top_scenario_III = tops.mask[good_chi2]
        top_III_mass = ak.fill_none(ak.firsts(top_scenario_III.mass), 0)
        mask_scenario_III = ak.fill_none(ak.firsts(good_chi2), False) 
    

        #print("Top tagger_mask_III (True): ", ak.sum(mask_scenario_III), mask_scenario_III.ndim, mask_scenario_III[mask_scenario_III == True])
        
        return top_scenario_III, mask_scenario_III, top_III_mass, ak.sum(mask_scenario_III) 
    
    
    # ----------------------------------------
    #  Scenario IV (W = 1  + b = 1)
    # ----------------------------------------              
    def Scenario_IV(self, wjets, bjets,
                    top_low_mass, top_up_mass,
                    w_low_mass, w_up_mass
                   ):
        jet_b = bjets.mask[ak.num(bjets) == 1]
        jet_w = wjets.mask[ak.num(wjets) == 1]
       

        good_w_mass = (
            (jet_w.mass > w_low_mass) &
            (jet_w.mass < w_up_mass)  &
            (delta_r_mask(jet_b, jet_w, threshold=0.8))
        )

        jet_top = jet_w.mask[good_w_mass] + jet_b
    
    
        good_top_mass = (
            (jet_top.mass > top_low_mass) &
            (jet_top.mass < top_up_mass)
        )    
    
            
        """
        We want to invert the filter of the mass of the top so as not to overwrite events of scenario 3.
        We do not apply chi2 since we cannot reconstruct the top.
            
        """

        top_scenario_IV = jet_top.mask[~good_top_mass]
        top_IV_mass = ak.fill_none(ak.firsts(top_scenario_IV.mass), 0)
        mask_scenario_IV = ak.fill_none(ak.firsts(~good_top_mass), False) # The mask has been inverted
    

        #print("Top tagger_mask_IV (True): ", ak.sum(mask_scenario_IV), mask_scenario_IV.ndim, mask_scenario_IV[mask_scenario_IV == True])
        
        return top_scenario_IV, mask_scenario_IV, top_IV_mass, ak.sum(mask_scenario_IV)

    
    
    ######################################
    #########  3 Jet #####################
    ######################################
        
    # ----------------------------------------
    #  Scenario V (W = 1  +  b = 2)
    # ----------------------------------------
    def Scenario_V(self, bjets, wjets,
                    top_sigma, top_low_mass, top_up_mass, top_pdg, 
                     w_sigma, w_low_mass, w_up_mass, w_pdg, 
                     chi2):
        
        
        jet_b = bjets.mask[ak.num(bjets) == 2]
        jet_w = wjets.mask[ak.num(wjets) == 1]
            
 
        jet_b1 = ak.pad_none(jet_b, 2)[:, 0]
        jet_b2 = ak.pad_none(jet_b, 2)[:, 1]
            
            
        # Cross cleaning bjets
        cross_cleaning = (jet_b2.delta_r(jet_b1) > 0.8)

        b1_cc = jet_b1.mask[cross_cleaning]            
        b2_cc = jet_b2.mask[cross_cleaning]

        
        
        good_w_mass = (
            (jet_w.mass > w_low_mass) &
            (jet_w.mass < w_up_mass) &
            (delta_r_mask(jet_b1, jet_w, threshold=0.8)) &
            (delta_r_mask(jet_b2, jet_w, threshold=0.8))
        )
            
        w = jet_w.mask[good_w_mass]
   

        top = ak.where(((w + b1_cc).mass > top_low_mass) & ((w + b1_cc).mass < top_up_mass),
                       w + b1_cc,
                       w + b2_cc)

            
        # ----------------------
        # chi2 criteria
        # ----------------------
        chi2_cal = chi2_test(top, w,
                                top_sigma, 
                                w_sigma, 
                                top_pdg, 
                                w_pdg)
            
        good_chi2 = (chi2_cal < chi2)
            
            
        top_scenario_V = top.mask[good_chi2]
        top_V_mass = ak.fill_none(ak.firsts(top_scenario_V.mass), 0)
        mask_scenario_V = ak.fill_none(ak.firsts(good_chi2), False)
  

        #print("Top tagger_mask_V (True): ", ak.sum(mask_scenario_V), mask_scenario_V.ndim, mask_scenario_V[mask_scenario_V == True])
    
        return top_scenario_V, mask_scenario_V, top_V_mass, ak.sum(mask_scenario_V)
    
    
        
    def Scenario_VI(self, jets,
                    top_sigma, top_low_mass, top_up_mass, top_pdg, 
                     w_sigma, w_low_mass, w_up_mass, w_pdg, 
                     chi2, btag_wp):
        
        multi_jets = jets.mask[ak.num(jets) == 3]
        
            
        trijet = ak.combinations(multi_jets, 3, fields=["j1", "j2", "j3"])
            
        # Cross cleaning
        cross_cleaning = (
            (trijet.j1.delta_r(trijet.j2) > 0.8) &
            (trijet.j1.delta_r(trijet.j3) > 0.8) &
            (trijet.j2.delta_r(trijet.j3) > 0.8)
        )

        trijet = trijet.mask[cross_cleaning]
            
        dijet = trijet.j1 + trijet.j2
        

        
        trijet["p4"] = dijet + trijet.j3
  
               
        """
        We identify the maximum b_wp of the jets, so 
        we can guarantee that at least one of the 3 jets can be considered a b-jet.
            
        """
        trijet["max_btag"] = np.maximum(trijet.j1.btagDeepFlavB,
                                        np.maximum(trijet.j2.btagDeepFlavB, 
                                                    trijet.j3.btagDeepFlavB),
        )
        
        # at least one-btag in bjj candidates
        trijet = trijet[trijet.max_btag > self._btag_wp]
            
            
            
        # w (jj) candidates mass
        dijet = ak.where(
                trijet.j1.btagDeepFlavB > self._btag_wp,             # Condition 1
                trijet.j2 + trijet.j3,                         # True
                (ak.where(trijet.j2.btagDeepFlavB > self._btag_wp,   # False: Condition 2
                    trijet.j1 + trijet.j3,                     # True
                    trijet.j1 + trijet.j2,                     # False
                        )
                ),
        )


        good_w_mass = (
            (dijet.mass > w_low_mass) &
            (dijet.mass < w_up_mass)
        )  
        
        jet_w = dijet.mask[good_w_mass]
        jet_top = trijet["p4"].mask[good_w_mass]
        

        good_top_mass = (
            (jet_top.mass > top_low_mass) &
            (jet_top.mass < top_up_mass)
        )  
        
        tops = jet_top.mask[good_top_mass]
        w = jet_w.mask[good_top_mass]
            
            
        # ----------------------
        # chi2 criteria
        # ----------------------
        chi2_cal = chi2_test(tops,  w, 
                            top_sigma, 
                            w_sigma, 
                            top_pdg, 
                            w_pdg)          
            
            
        good_chi2 = (chi2_cal < chi2)

        
        top_scenario_VI = ak.mask(tops,good_chi2)
        top_VI_mass = ak.fill_none(ak.firsts(top_scenario_VI.mass), 0)
        mask_scenario_VI = ak.fill_none(ak.any(good_chi2, axis = 1), False) 
        

        #print("Top tagger_mask_VI (True): ", ak.sum(mask_scenario_VI), mask_scenario_VI.ndim, mask_scenario_VI[mask_scenario_VI == True])
    
        return top_scenario_VI, mask_scenario_VI, top_VI_mass, ak.sum(mask_scenario_VI)   
    
    # ----------------------------------------
    #  Scenario VII (b = 1  + light_jets = 2)
    # ----------------------------------------
    def Scenario_VII(self, jets,
                    top_low_mass, top_up_mass,
                    w_low_mass, w_up_mass,
                    btag_wp):
        
        multi_jets = jets.mask[ak.num(jets) == 3]
            
        trijet = ak.combinations(multi_jets, 3, fields=["j1", "j2", "j3"])
            
        # Cross cleaning
        cross_cleaning = (
            (trijet.j1.delta_r(trijet.j2) > 0.8) &
            (trijet.j1.delta_r(trijet.j3) > 0.8) &
            (trijet.j2.delta_r(trijet.j3) > 0.8)
        )

        trijet = trijet.mask[cross_cleaning]
            
        dijet = trijet.j1 + trijet.j2
        trijet["p4"] = dijet + trijet.j3
            
               
        """
        We identify the maximum b_wp of the jets, so 
        we can guarantee that at least one of the 3 jets can be considered a b-jet.
            
        """
        trijet["max_btag"] = np.maximum(trijet.j1.btagDeepFlavB,
                                        np.maximum(trijet.j2.btagDeepFlavB, 
                                                    trijet.j3.btagDeepFlavB),
        )
        
        # at least one-btag in bjj candidates
        trijet = trijet[trijet.max_btag > self._btag_wp]
            
        # We verify that there are no more b-jets
            
            
        # w (jj) candidates mass
        dijet = ak.where(
                trijet.j1.btagDeepFlavB > btag_wp,             # Condition 1
                trijet.j2 + trijet.j3,                         # True
                (ak.where(trijet.j2.btagDeepFlavB > btag_wp,   # False: Condition 2
                    trijet.j1 + trijet.j3,                     # True
                    trijet.j1 + trijet.j2,                     # False
                        )
                ),
        )
       
    
        good_w_mass = (
            (dijet.mass > w_low_mass) &
            (dijet.mass < w_up_mass)
        )  
                    
            
        jet_w = dijet.mask[good_w_mass]
        jet_top = trijet["p4"].mask[good_w_mass]
        

        good_top_mass = (
            (jet_top.mass > top_low_mass) &
            (jet_top.mass < top_up_mass)
        )  
            
            
        top_scenario_VII = jet_top.mask[~good_top_mass]
        top_VII_mass = ak.fill_none(ak.firsts(top_scenario_VII.mass), 0)
        mask_scenario_VII = ak.fill_none(ak.any(~good_chi2, axis = 1), False) # The mask has been inverted
    
    
        #print("Top tagger_mask_VII (True): ", ak.sum(mask_scenario_VII), mask_scenario_VII.ndim, mask_scenario_VII[mask_scenario_VII == True])
        
        return top_scenario_VII, mask_scenario_VII, top_VII_mass, ak.sum(mask_scenario_VII)
    
    
    
    ######################################
    #########  4 Jet #####################
    ######################################
        
    # ----------------------------------------
    #  Scenario IX (b=2  + light_jets = 2)
    # ----------------------------------------        
    def Scenario_IX(self, jets,
                      top_sigma, top_low_mass, top_up_mass, top_pdg,
                      w_sigma, w_low_mass, w_up_mass, w_pdg,
                      chi2, btag_wp):


        multi_jets = jets.mask[ak.num(jets) == 4]


        trijet = ak.combinations(multi_jets, 3, fields=["j1", "j2", "j3"])

        # -------------------------------------------
        # Cross cleaning
        cross_cleaning = (
                (trijet.j1.delta_r(trijet.j2) > 0.8)
                & (trijet.j1.delta_r(trijet.j3) > 0.8)
                & (trijet.j2.delta_r(trijet.j3) > 0.8)
                & (trijet.j2.delta_r(trijet.j3) > 0.8)
        )

        trijet = trijet.mask[cross_cleaning]
        # -------------------------------------------

        dijet = trijet.j1 + trijet.j2
        trijet["p4"] = dijet + trijet.j3

        trijet["max_btag"] = np.maximum(
            trijet.j1.btagDeepFlavB,
            np.maximum(trijet.j2.btagDeepFlavB, trijet.j3.btagDeepFlavB),
        )

        # at least one-btag in bjj candidates
        trijet = trijet[trijet.max_btag > self._btag_wp]


        # w (jj) candidates mass
        dijet = ak.where(
            trijet.j1.btagDeepFlavB > self._btag_wp,
            trijet.j2 + trijet.j3,
            (
                ak.where(
                    trijet.j2.btagDeepFlavB > self._btag_wp,
                    trijet.j1 + trijet.j3,
                    trijet.j1 + trijet.j2,
                )
            ),
        )


        good_w = (
                (dijet.mass > w_low_mass)
                & (dijet.mass < w_up_mass)
        )  

        dijet = dijet.mask[good_w]
        trijet = trijet.mask[good_w]


        good_tops = (
                (trijet["p4"].mass > top_low_mass)
                & (trijet["p4"].mass < top_up_mass)
        )  

        good_trijet = trijet.mask[good_tops]
        good_dijet  = dijet.mask[good_tops]


        # ----------------------
        # chi2 criteria
        # ----------------------
        chi2_cal = chi2_test(good_trijet["p4"], good_dijet, top_sigma, w_sigma, top_pdg, w_pdg)


        good_chi2 = chi2_cal < chi2


        top_chi2 = ak.mask(good_trijet["p4"],good_chi2)
        w_chi2 = ak.mask(good_dijet,good_chi2)





        top_scenario_IX = trijet.mask[good_chi2]["p4"]
        top_IX_mass = ak.fill_none(ak.firsts(top_scenario_IX.mass), 0)
        w_scenario_IX = dijet.mask[good_chi2]
        mask_scenario_IX = ak.fill_none(ak.any(good_chi2, axis=1), False)

        #print("Top tagger_mask_IX (True): ", ak.sum(mask_scenario_IX), mask_scenario_IX.ndim, mask_scenario_IX[mask_scenario_IX == True])

        return top_scenario_IX, mask_scenario_IX, top_IX_mass, ak.sum(mask_scenario_IX)
    
    ######################################
    #########  >= 5 Jet ##################
    ######################################
        
    # ----------------------------------------
    #  Scenario X (b=2  + light_jets >= 3)
    # ---------------------------------------- 
    def Scenario_X(self, jets,
                    top_sigma, top_low_mass, top_up_mass, top_pdg, 
                     w_sigma, w_low_mass, w_up_mass, w_pdg, 
                     chi2, btag_wp):
        
        multi_jets = jets.mask[ak.num(jets) >= 5]
    
        trijet = ak.combinations(multi_jets, 3, fields=["j1", "j2", "j3"])
    
        # Cross cleaning
        cross_cleaning = (
            (trijet.j1.delta_r(trijet.j2) > 0.8) &
            (trijet.j1.delta_r(trijet.j3) > 0.8) &
            (trijet.j2.delta_r(trijet.j3) > 0.8) &
            (trijet.j2.delta_r(trijet.j3) > 0.8)
        )

        trijet = trijet.mask[cross_cleaning]

            
        dijet = trijet.j1 + trijet.j2
        trijet["p4"] = dijet + trijet.j3
    
        """
        We identify the maximum b_wp of the jets, so 
        we can guarantee that at least one of the 3 jets can be considered a b-jet.
            
        """
        trijet["max_btag"] = np.maximum(trijet.j1.btagDeepFlavB,
                                        np.maximum(trijet.j2.btagDeepFlavB, 
                                                    trijet.j3.btagDeepFlavB),
        )
    
    
        # at least one-btag in bjj candidates
        trijet = trijet[trijet.max_btag > self._btag_wp]
    
    
        # w (jj) candidates mass
        dijet = ak.where(
                trijet.j1.btagDeepFlavB > self._btag_wp,             # Condition 1
                trijet.j2 + trijet.j3,                         # True
                (ak.where(trijet.j2.btagDeepFlavB > self._btag_wp,   # False: Condition 2
                    trijet.j1 + trijet.j3,                     # True
                    trijet.j1 + trijet.j2,                     # False
                        )
                ),
        )
            
    
        good_w_mass = (
            (dijet.mass > w_low_mass) &
            (dijet.mass < w_up_mass)
        )  
        
        jet_w = dijet.mask[good_w_mass]
        jet_top = trijet["p4"].mask[good_w_mass]
      

    
        good_top_mass = (
            (trijet["p4"].mass > top_low_mass) &
            (trijet["p4"].mass < top_up_mass)
        )  
        
        
        tops = jet_top.mask[good_top_mass]
        w = jet_w.mask[good_top_mass]
            
            
        # ----------------------
        # chi2 criteria
        # ----------------------
        chi2_cal = chi2_test(tops,  w, 
                            top_sigma, 
                            w_sigma, 
                            top_pdg, 
                            w_pdg)          
            
            
        good_chi2 = (chi2_cal < chi2)
            
            
        top_scenario_X = ak.mask(tops,good_chi2)
        mask_scenario_X = ak.fill_none(ak.any(good_chi2, axis=1), False) 
            
        #print("Top tagger_mask_X (True): ", ak.sum(mask_scenario_X))
        
        
        return top_scenario_X, mask_scenario_X
