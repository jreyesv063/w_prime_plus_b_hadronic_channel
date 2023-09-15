import json
import copy
import pickle
import numpy as np
import awkward as ak
import importlib.resources
from coffea import processor
from coffea.analysis_tools import Weights, PackedSelection

# Utils
from wprime_plus_b.processors.utils import histograms, topXfinder
from wprime_plus_b.processors.utils.analysis_utils import delta_r_mask, normalize, chi2_test, topTagger_mask, chi2_test, topTagger_mask, pdg_masses, tagger_constants, entries, true_sum


# Corrections
from wprime_plus_b.corrections.btag import BTagCorrector
from wprime_plus_b.corrections.jec import jet_corrections
from wprime_plus_b.corrections.met import met_phi_corrections
from wprime_plus_b.corrections.pileup import add_pileup_weight
from wprime_plus_b.corrections.lepton import ElectronCorrector, MuonCorrector
from wprime_plus_b.corrections.tau import TauCorrector

# Objects
from wprime_plus_b.selections.ttbar.jet_selection_hA import (
    select_good_Fatjets,
    select_good_Wjets,
    select_good_jets,
    select_good_bjets_prev 
)
from wprime_plus_b.selections.ttbar.lepton_selection_hA import (
    select_good_electrons,
    select_good_muons,
    select_good_taus  
)

# Top tagger
from wprime_plus_b.processors.utils.topXfinder import topXfinder


class TtbarAnalysis_h(processor.ProcessorABC):
    def __init__(
        self,
        year: str = "2017",
        yearmod: str = " ",
        btag_wp: str = "L",
        TvsQCD: str = "T",
        WvsQCD: str = "T",
        syst: str = "nominal",
        output_type="hist",
    ):
        
        self._year = year
        self._yearmod = yearmod
        self._btag_wp = btag_wp
        self._TvsQCD = TvsQCD
        self._WvsQCD = WvsQCD
        self._syst = syst
        self._output_type = output_type
        
        
        # initialize dictionary of hists for control regions
        self.hist_dict = {
            "tau_kin": histograms.ttbar_h_tau_hist,
            "met_kin": histograms.ttbar_h_met_hist,
            "b_kin": histograms.ttbar_h_b_hist,
            "tau_met_kin": histograms.ttbar_h_tau_met_hist,
            "tau_b_kin": histograms.ttbar_h_tau_b_hist,
            "top_mass": histograms.ttbar_h_top_hist
        }


        # define dictionary to store analysis variables
        self.features = {}   

        # initialize dictionary of arrays
        self.array_dict = {}
        
        
    def add_feature(self, name: str, var: ak.Array) -> None:
        """add a variable array to the out dictionary"""
        self.features = {**self.features, name: var}
        
        
    def process(self, events):
        # get dataset name
        dataset = events.metadata["dataset"]

        # get number of events before selection
        nevents = len(events)

        # check if sample is MC
        self.is_mc = hasattr(events, "genWeight")

        # create copies of histogram objects
        hist_dict = copy.deepcopy(self.hist_dict)

        # create copy of array dictionary
        array_dict = copy.deepcopy(self.array_dict)
        

        # define systematic variations
        syst_variations = ["nominal"]
        if self.is_mc:
            jet_jec_syst_variations = ["JESUp", "JESDown"]
            jet_jer_syst_variations = ["JERUp", "JERDown"]
            met_obj_syst_variations = ["UEUp", "UEDown"]

            if self._syst == "jec":
                syst_variations.extend(jet_jec_syst_variations)
            elif self._syst == "jer":
                syst_variations.extend(jet_jer_syst_variations)
            elif self._syst == "jec":
                syst_variations.extend(jet_jec_syst_variations)
                syst_variations.extend(jet_jer_syst_variations)
            elif self._syst == "met":
                syst_variations.extend(met_obj_syst_variations)
            elif self._syst == "full":
                syst_variations.extend(jet_jec_syst_variations)
                syst_variations.extend(jet_jer_syst_variations)
                syst_variations.extend(met_obj_syst_variations)
                
        for syst_var in syst_variations:
            # ------------------
            # event preselection
            # leptons
            # ------------------
            good_electrons = (select_good_electrons(events))
            electrons = events.Electron[good_electrons]
            
            good_muons = (select_good_muons(events) &
                         delta_r_mask(events.Muon, electrons, threshold=0.4))
            muons = events.Muon[good_muons]
            
            # select good taus
            good_taus = (select_good_taus(events)
                        & (delta_r_mask(events.Tau, electrons, threshold=0.4))
                        & (delta_r_mask(events.Tau, muons, threshold=0.4))
            )
            taus = events.Tau[good_taus]
            
            
            # ------------------
            # Fatjets
            # -------------------            
            # select good Fatjets
            good_fatjets = (select_good_Fatjets(events, 
                                                self._year + self._yearmod,           
                                                working_point=self._TvsQCD)
                            & (delta_r_mask(events.FatJet, electrons, threshold=0.4))
                            & (delta_r_mask(events.FatJet, muons, threshold=0.4))
            )
            fatjets = events.FatJet[good_fatjets]

            # select good W jets
            good_wjets = (select_good_Wjets(events, 
                                            self._year + self._yearmod, 
                                            working_point=self._WvsQCD)
                            & (delta_r_mask(events.FatJet, electrons, threshold=0.4))
                            & (delta_r_mask(events.FatJet, muons, threshold=0.4))
                            & (delta_r_mask(events.FatJet, fatjets, threshold=0.8))
                         )
            wjets = events.FatJet[good_wjets]
            
 
            # ------------------
            # jets
            # -------------------
            # apply JEC/JER corrections to jets (in data, the corrections are already applied)
            if self.is_mc:
                corrected_jets, met = jet_corrections(
                    events, self._year + self._yearmod
                )
                # jet JEC/JER shift
                if syst_var == "JESUp":
                    corrected_jets = corrected_jets.JES_Total.up
                elif syst_var == "JESDown":
                    corrected_jets = corrected_jets.JES_Total.down
                elif syst_var == "JERUp":
                    corrected_jets = corrected_jets.JER.up
                elif syst_var == "JERDown":
                    corrected_jets = corrected_jets.JER.down
                # MET UnclusteredEnergy shift
                elif syst_var == "UEUp":
                    met = met.MET_UnclusteredEnergy.up
                elif syst_var == "UEDown":
                    met = met.MET_UnclusteredEnergy.down
            else:
                corrected_jets, met = events.Jet, events.MET
                

            # select good jets
            good_jets = (select_good_jets(corrected_jets, 
                                          self._year + self._yearmod)
                            & (delta_r_mask(events.Jet, electrons, threshold=0.4))
                            & (delta_r_mask(events.Jet, muons, threshold=0.4))
                            & (delta_r_mask(events.Jet, fatjets, threshold=0.8))
                            & (delta_r_mask(events.Jet, wjets, threshold=0.8))
            )    
            jets = events.Jet[good_jets]
            

            # select good bjets
            good_bjets = (select_good_bjets_prev(corrected_jets, 
                                          self._year + self._yearmod,
                                           "M")
                            & (delta_r_mask(events.Jet, electrons, threshold=0.4))
                            & (delta_r_mask(events.Jet, muons, threshold=0.4))
                            & (delta_r_mask(events.Jet, fatjets, threshold=0.8))
                            & (delta_r_mask(events.Jet, wjets, threshold=0.8))
                            & (delta_r_mask(events.Jet, jets, threshold=0.4))
            )    
            bjets = events.Jet[good_bjets]
            
            
            # apply MET phi corrections
            met_pt, met_phi = met_phi_corrections(
                met_pt=met.pt,
                met_phi=met.phi,
                npvs=events.PV.npvs,
                is_mc=self.is_mc,
                year=self._year,
                year_mod=self._yearmod,
            )
            met["pt"], met["phi"] = met_pt, met_phi  
 
            # ---------------
            # event selection
            # ---------------
            # make a PackedSelection object to store selection masks
            self.selections = PackedSelection()   

            # add luminosity calibration mask (only to data)
            with importlib.resources.path(
                "wprime_plus_b.data", "lumi_masks.pkl"
            ) as path:
                with open(path, "rb") as handle:
                    self._lumi_mask = pickle.load(handle)
            if not self.is_mc:
                lumi_mask = self._lumi_mask[self._year](
                    events.run, events.luminosityBlock
                )
            else:
                lumi_mask = np.ones(len(events), dtype="bool")
            
            self.selections.add("lumi", lumi_mask)
            
            
            # add lepton triggers masks
            with importlib.resources.path(
                "wprime_plus_b.data", "triggers.json"
            ) as path:
                with open(path, "r") as handle:
                    self._triggers = json.load(handle)[self._year]
            
            trigger = {}       
            for ch in ["tau"]: 
                trigger[ch] = np.zeros(nevents, dtype="bool")
                for t in self._triggers[ch]:
                    #print(t)
                    if t in events.HLT.fields:
                        trigger[ch] = trigger[ch] | events.HLT[t]
                        
            self.selections.add("trigger_tau", trigger["tau"])  
            
            
            # add MET filters mask
            # open and load met filters
            with importlib.resources.path(
                "wprime_plus_b.data", "metfilters.json"
            ) as path:
                with open(path, "r") as handle:
                    self._metfilters = json.load(handle)[self._year]
            metfilters = np.ones(nevents, dtype="bool")
            metfilterkey = "mc" if self.is_mc else "data"
            for mf in self._metfilters[metfilterkey]:
                if mf in events.Flag.fields:
                    metfilters = metfilters & events.Flag[mf]
                    
            self.selections.add("metfilters", metfilters)
            
            # check that there be a minimum MET greater than 250 GeV
            self.selections.add("met_pt", met.pt > 250)

            # add number of leptons and jets
            self.selections.add("electron_veto", ak.num(electrons) == 0)
            self.selections.add("muon_veto", ak.num(muons) == 0)
            self.selections.add("one_tau", ak.num(taus) == 1)
            one_or_two_bjets = topTagger_mask(ak.num(bjets) == 1, ak.num(bjets) == 2)
            self.selections.add("2or1bjets", one_or_two_bjets)

            # define selection regions for each channel
            region_selection = {
                    "tt_CR": [
                        "lumi",
                        "metfilters",
                        "trigger_tau",
                        "muon_veto",
                        "electron_veto",
                        "one_tau",
                        "met_pt",
                        "2or1bjets",
                    ]
            }
            
            # ---------------
            # event variables
            # ---------------
            region_selection =  self.selections.all(*region_selection["tt_CR"])
            
            # if there are no events left after selection cuts continue to the next .root file
            nevents_after = ak.sum(region_selection)
            

            if nevents_after > 0:
                # select region objects
                region_bjets = bjets[region_selection]
                region_taus = taus[region_selection]
                region_met = met[region_selection]
                
                region_fatjets = fatjets[region_selection]
                region_jets = jets[region_selection]
                region_wjets = wjets[region_selection]
                
                # -------------------------------------------------------
                # -- Top Tagger: resolve, unresolve, partially_resolve --
                # -------------------------------------------------------
                
                self.topXfinder = PackedSelection()  # We will store all the mask used.
                
                top_pdg, w_pdg = pdg_masses()
                
                top_sigma, top_low, top_up, w_sigma, w_low, w_up, chi2 = tagger_constants("hadronic")  
                
                # Object and cases
                topX = topXfinder(self._year, self._yearmod, self._btag_wp) 

                top_I, mask_I, mass_I, nevents_caseI = topX.Scenario_I(region_fatjets, 
                                                                       top_low, top_up)

                
                top_II, mask_II, mass_II, nevents_caseII = topX.Scenario_II(region_fatjets, 
                                                                            region_bjets, 
                                                                            top_low, top_up)


                top_III, mask_III, mass_III, nevents_caseIII = topX.Scenario_III(region_wjets,        
                                                                                 region_bjets, 
                                                                        top_sigma, top_low, top_up, top_pdg,
                                                                        w_sigma, w_low, w_up, w_pdg,
                                                                        chi2)
                           
                nevents_caseIV = 0
                
                top_V, mask_V, mass_V, nevents_caseV = topX.Scenario_V(region_wjets, 
                                                                       region_bjets, 
                                                                      top_sigma, top_low, top_up, top_pdg,
                                                                      w_sigma, w_low, w_up, w_pdg,
                                                                      chi2)                                

                
                top_VI, mask_VI, mass_VI, nevents_caseVI = topX.Scenario_VI(region_jets, 
                                                                      top_sigma, top_low, top_up, top_pdg,
                                                                      w_sigma, w_low, w_up, w_pdg,
                                                                      chi2, self._btag_wp)
                
                nevents_caseVII = 0             
                                
                top_IX, mask_IX, mass_IX, nevents_caseIX = topX.Scenario_IX(region_jets, 
                                                                  top_sigma, top_low, top_up, top_pdg,
                                                                  w_sigma, w_low, w_up, w_pdg,
                                                                  chi2, self._btag_wp)
               
                top_mass = mass_I + mass_II + mass_III + mass_V + mass_VI + mass_IX 
                
                mask_topXfinder = mask_I
                mask_topXfinder = topTagger_mask(mask_topXfinder, mask_II)
                mask_topXfinder = topTagger_mask(mask_topXfinder, mask_III)                
                mask_topXfinder = topTagger_mask(mask_topXfinder, mask_V)
                mask_topXfinder = topTagger_mask(mask_topXfinder, mask_VI)
                mask_topXfinder = topTagger_mask(mask_topXfinder, mask_IX)
                
                
                self.topTagger = PackedSelection()
                self.topTagger.add("top_tagger", mask_topXfinder) 
                
                final_selection = {
                    "top_tagger",
                }
                
                
                final_selection =  self.topTagger.all(*final_selection)
            
                region_bjets = region_bjets[final_selection]
                region_taus = region_taus[final_selection]
                region_met = region_met[final_selection]
                top_mass = top_mass[final_selection]

                nevents_after = ak.sum(final_selection)
                
                if nevents_after > 0 :

                    # ----------------------------------------------------
                    # -------------------- Plots -------------------------
                    # ----------------------------------------------------

                    tau_MET_dphi = region_taus.delta_phi(region_met) 
                    tau_met_transMass = np.sqrt(
                        2.0
                        * region_taus.pt
                        * region_met.pt
                        * (
                            ak.ones_like(region_met.pt)
                            - np.cos(region_taus.delta_phi(region_met))
                        )
                    )


                     # leading bjets
                    leading_bjets = ak.firsts(region_bjets)
                    tau_bjet_dr = leading_bjets.delta_r(region_taus)
                    tau_bjet_dphi = leading_bjets.delta_phi(region_taus)


                    # -- Tau histograms --
                    self.add_feature("tau_pt", region_taus.pt)
                    self.add_feature("tau_eta", region_taus.eta) 
                    self.add_feature("tau_phi", region_taus.phi) 

                    # -- MET histograms --
                    self.add_feature("met_pt", region_met.pt)
                    self.add_feature("met_phi", region_met.phi) 

                    # -- bjet histograms --
                    self.add_feature("b_pt", region_bjets.pt)
                    self.add_feature("b_eta", region_bjets.eta) 
                    self.add_feature("b_phi", region_bjets.phi) 

                    # -- Tau + MET histograms --
                    self.add_feature("tau_MET_deltaPhi", region_bjets.pt)
                    self.add_feature("tau_met_transverseMass", tau_met_transMass)                     

                    # -- Tau + bjet histograms --
                    self.add_feature("tau_b_deltaPhi", tau_bjet_dphi)
                    self.add_feature("tau_b_deltaR", tau_bjet_dr) 
                    
                    # -- Top histogram --
                    self.add_feature("top_mrec", top_mass)




                    # ----------------------------------------------------
                    # ------------------ event weights -------------------
                    # ----------------------------------------------------
                    weights_container = Weights(
                        len(events[final_selection]), storeIndividual=True
                    )
                    if self.is_mc:
                        # add gen weigths
                        gen_weight = events.genWeight[final_selection]
                        weights_container.add("genweight", gen_weight)

                        # add L1prefiring weights
                        if self._year in ("2016", "2017"):
                            if syst_var == "nominal":
                                weights_container.add(
                                    "L1Prefiring",
                                    weight=events.L1PreFiringWeight.Nom[final_selection],
                                    weightUp=events.L1PreFiringWeight.Up[final_selection],
                                    weightDown=events.L1PreFiringWeight.Dn[final_selection],
                                )
                            else:
                                weights_container.add(
                                    "L1Prefiring",
                                    weight=events.L1PreFiringWeight.Nom[region_selection],
                                )
                        # add pileup reweighting
                        add_pileup_weight(
                            n_true_interactions=ak.to_numpy(
                                events.Pileup.nPU[final_selection]
                            ),
                            weights=weights_container,
                            year=self._year,
                            year_mod=self._yearmod,
                            variation=syst_var,
                        )

                        # b-tagging corrector
                        njets = 2
                        btag_corrector = BTagCorrector(
                            jets=region_bjets,
                            njets=njets,
                            weights=weights_container,
                            sf_type="comb",
                            worging_point="M",
                            tagger="deepJet",
                            year=self._year,
                            year_mod=self._yearmod,
                            full_run=False,
                            variation=syst_var,
                        )
                        # add b-tagging weights
                        btag_corrector.add_btag_weights(flavor="bc")
                        btag_corrector.add_btag_weights(flavor="light")


                        # Tau corrector
                        tau_corrector = TauCorrector(
                                            taus = ak.firsts(region_taus), 
                                            weights=weights_container, 
                                            year=self._year, 
                                            year_mod=self._yearmod,
                                            tag="tau", 
                                            variation=syst_var
                        )
                        tau_corrector.add_id_weight_DeepTau2017v2p1VSe("Tight", "nom")
                        tau_corrector.add_id_weight_DeepTau2017v2p1VSmu("Tight", "nom")
                        tau_corrector.add_id_weight_DeepTau2017v2p1VSjet("Tight", "Tight", "nom", "pt")
                        tau_corrector.add_id_weight_tau_energy_scale("DeepTau2017v2p1", "nom")


                    # ------------------
                    # histogram filling
                    # ------------------
                    if self._output_type == "hist":       

                        # break up the histogram filling for event-wise variations and object-wise variations
                        # apply event-wise variations only for nominal
                        if self.is_mc and syst_var == "nominal":
                            # get event weight systematic variations for MC samples
                            event_weights = [
                                weight
                                for weight in weights_container.weightStatistics
                                if "genweight" not in weight
                            ]
                            event_weight_syst_variations_up = [
                                f"{event_weight}Up" for event_weight in event_weights
                            ]
                            event_weight_syst_variations_down = [
                                f"{event_weight}Down" for event_weight in event_weights
                            ]
                            event_weight_syst = ["nominal"]
                            event_weight_syst.extend(event_weight_syst_variations_up)
                            event_weight_syst.extend(event_weight_syst_variations_down)

                            for variation in event_weight_syst:
                                # get weight
                                if variation == "nominal":
                                    syst_weight = weights_container.weight()
                                else:
                                    syst_weight = weights_container.weight(variation)

                                for kin in hist_dict:
                                    # get filling arguments
                                    fill_args = {
                                        feature: normalize(self.features[feature])
                                        for feature in hist_dict[
                                            kin
                                        ].axes.name[:-1]
                                        if "dataset" not in feature
                                    }
                                    # fill histograms
                                    hist_dict[kin].fill(
                                        **fill_args,
                                        dataset=dataset,
                                        variation=variation,
                                        weight=syst_weight
                                    )

                        # object-wise variations
                        syst_weight = weights_container.weight()
                        for kin in hist_dict:
                            # get filling arguments
                            fill_args = {
                                feature: normalize(self.features[feature])
                                for feature in hist_dict[kin].axes.name[:-1]
                                if "dataset" not in feature
                            }
                            # fill histograms
                            hist_dict[kin].fill(
                                **fill_args,
                                dataset=dataset,
                                variation=syst_var,
                                weight=syst_weight,
                            )

                    elif self._output_type == "array":
                        self.add_feature("weights", weights_container.weight())
                        self.add_feature("genweights", weights_container.partial_weight("genweight"))
                        # select variables and put them in column accumulators
                        array_dict = {
                            feature_name: processor.column_accumulator(
                                normalize(feature_array)
                            )
                            for feature_name, feature_array in self.features.items()
                        }  


            else:
                nevents_caseI  = 0
                nevents_caseII = 0
                nevents_caseIII= 0
                nevents_caseIV = 0
                nevents_caseV  = 0
                nevents_caseVI = 0
                nevents_caseVII= 0
                nevents_caseIX = 0
                
        # define output dictionary accumulator
        output = {}
        if self._output_type == "hist":
            output["histograms"] = hist_dict
        elif self._output_type == "array":
            output["arrays"] = array_dict
        # save metadata
        output["metadata"] = {
            "events_before": nevents,
            "events_after": nevents_after,
            "events_case_I": nevents_caseI,
            "events_case_II": nevents_caseII,
            "events_case_III": nevents_caseIII,
            "events_case_IV": nevents_caseIV,
            "events_case_V": nevents_caseV,
            "events_case_VI": nevents_caseVI,
            "events_case_VII": nevents_caseVII,
            "events_case_IX": nevents_caseIX
        }
        # save sumw for MC samples
        if self.is_mc:
            output["metadata"].update({"sumw": ak.sum(events.genWeight)})
        return output
        
    def postprocess(self, accumulator):
        return accumulator