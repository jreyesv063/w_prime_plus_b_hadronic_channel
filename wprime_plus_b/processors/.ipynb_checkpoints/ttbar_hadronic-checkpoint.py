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
    select_good_bjets,
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
        yearmod: str = "",
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
#            "top_kin": histograms.ttbar_h_top_hist
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
            # leptons
            # -------------------
            # select good electrons
            good_electrons = select_good_electrons(events)
            electrons = events.Electron.mask[good_electrons]

            # select good muons
            good_muons = (select_good_muons(events) 
                        & delta_r_mask(events.Muon, electrons, threshold=0.4)
            )
            muons = events.Muon.mask[good_muons]

            # select good taus
            good_taus = (select_good_taus(events)
                        & (delta_r_mask(events.Tau, electrons, threshold=0.4))
                        & (delta_r_mask(events.Tau, muons, threshold=0.4))
            )
            taus = events.Tau.mask[good_taus]
            
       
            # ------------------
            # Fatjets
            # -------------------            
            # select good Fatjets
            good_fatjets = (select_good_Fatjets(events.FatJet, self._year + self._yearmod, working_point=self._TvsQCD)
                            & (delta_r_mask(events.FatJet, electrons, threshold=0.4))
                            & (delta_r_mask(events.FatJet, muons, threshold=0.4))
            )
            fatjets = events.FatJet.mask[good_fatjets]

            # select good W jets
            good_wjets = (select_good_Wjets(events.FatJet, self._year + self._yearmod, working_point=self._WvsQCD)
                            & (delta_r_mask(events.FatJet, electrons, threshold=0.4))
                            & (delta_r_mask(events.FatJet, muons, threshold=0.4))
                            & (delta_r_mask(events.FatJet, fatjets, threshold=0.8))
                         )
            wjets = events.FatJet.mask[good_wjets]
                      
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
            good_jets = (select_good_jets(corrected_jets, self._year)
                            & (delta_r_mask(corrected_jets, electrons, threshold=0.4))
                            & (delta_r_mask(corrected_jets, muons, threshold=0.4))
                            & (delta_r_mask(corrected_jets, fatjets, threshold=0.8))
                            & (delta_r_mask(corrected_jets, wjets, threshold=0.8))
            )    
            jets = corrected_jets.mask[good_jets]
                
                
            good_bjets = select_good_bjets_prev(
                events=events,
                year=self._year,
                btag_working_point= "M",
                jet_pt_threshold= 20,
                jet_id= 6,
                jet_pileup_id= 7
            )
            good_bjets = (
                good_bjets
                & (delta_r_mask(corrected_jets, electrons, threshold=0.4))
                & (delta_r_mask(corrected_jets, muons, threshold=0.4))
            )
            bjets = corrected_jets[good_bjets]
                
            """
            # select good bjets
            good_bjets = (select_good_bjets(corrected_jets, year=self._year, working_point= self._btag_wp)
                            & (delta_r_mask(corrected_jets, electrons, threshold=0.4))
                            & (delta_r_mask(corrected_jets, muons, threshold=0.4))
                            & (delta_r_mask(corrected_jets, fatjets, threshold=0.8))
                            & (delta_r_mask(corrected_jets, wjets, threshold=0.8))
            )
            bjets = corrected_jets.mask[good_bjets]
            """
            
            
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
            
            # add number of bjets
            self.selections.add("one_bjet", ak.num(bjets) == 1)

            
            region_selection = {
                "tt_tTauNuB": [
                                "lumi",
                                "metfilters",
                                "trigger_tau",
                                "muon_veto",
                                "electron_veto",
                                "one_tau",
                                "met_pt"
                ]
            }            
            
            # ---------------
            # event variables
            # ---------------
            region_selection =  self.selections.all(*region_selection["tt_tTauNuB"])
            
            # if there are no events left after selection cuts continue to the next .root file
            nevents_after = ak.sum(region_selection)
    
        
            if nevents_after > 0:
                # select region objects
                preselection_electrons = electrons[region_selection]
                preselection_muons = muons[region_selection]
                preselection_taus = taus[region_selection]
                preselection_met = met[region_selection]
                preselection_fatjets = fatjets[region_selection]
                preselection_jets = jets[region_selection]
                preselection_bjets = bjets[region_selection]
                preselection_wjets = wjets[region_selection]
      
                # -----------------
                # Top Tagger
                # resolve, unresolve, partially_resolve
                # -----------------
                self.b_top_selection = PackedSelection()  # We will store all the mask used.
                               
                top_pdg, w_pdg = pdg_masses()
                
                top_sigma, top_low, top_up, w_sigma, w_low, w_up, chi2 = tagger_constants("hadronic")  
                
                # Object and cases
                topX = topXfinder(self._year, self._yearmod, self._btag_wp) 
                
                top_I, mask_I, mass_I = topX.Scenario_I(preselection_fatjets, top_low, top_up)

                
                top_II, mask_II, mass_II = topX.Scenario_II(preselection_fatjets, preselection_bjets, 
                                                   top_low, top_up)


                top_III, mask_III, mass_III = topX.Scenario_III(preselection_wjets, preselection_bjets, 
                                                      top_sigma, top_low, top_up, top_pdg,
                                                      w_sigma, w_low, w_up, w_pdg,
                                                      chi2)
                            
                
                
                top_V, mask_V, mass_V = topX.Scenario_V(preselection_wjets, preselection_bjets, 
                                                      top_sigma, top_low, top_up, top_pdg,
                                                      w_sigma, w_low, w_up, w_pdg,
                                                      chi2)                                

                
                top_VI, mask_VI, mass_VI = topX.Scenario_VI(preselection_jets, 
                                                      top_sigma, top_low, top_up, top_pdg,
                                                      w_sigma, w_low, w_up, w_pdg,
                                                      chi2, self._btag_wp)
                
                               
                                
                top_IX, mask_IX, mass_IX = topX.Scenario_IX(preselection_jets, 
                                  top_sigma, top_low, top_up, top_pdg,
                                  w_sigma, w_low, w_up, w_pdg,
                                  chi2, self._btag_wp)
               
                
                mask_topXfinder = mask_I
                mask_topXfinder = topTagger_mask(mask_topXfinder, mask_II)
                mask_topXfinder = topTagger_mask(mask_topXfinder, mask_III)                
                mask_topXfinder = topTagger_mask(mask_topXfinder, mask_V)
                mask_topXfinder = topTagger_mask(mask_topXfinder, mask_VI)
                mask_topXfinder = topTagger_mask(mask_topXfinder, mask_IX)
                
                
                # Sumamos todas las contribucciones para ver como se recontruyen los tops
                top_jet_mass = mass_I + mass_II + mass_III + mass_V + mass_VI + mass_IX
                

                self.b_top_selection.add("topTagger", mask_topXfinder) 
                self.b_top_selection.add("one_bjet", ak.num(preselection_bjets) == 1) 
                
                
                final_selection = {
                    "topTagger",
                    "one_bjet"
                }
                
                
                final_selection =  self.b_top_selection.all(*final_selection)                                
     
            
                region_bjets = preselection_bjets[final_selection]
                region_electrons = preselection_electrons[final_selection]
                region_muons = preselection_muons[final_selection]
                region_taus = preselection_taus[final_selection]
                region_met = preselection_met[final_selection]
                region_fatjets = preselection_fatjets[final_selection]
                region_jets = preselection_jets[final_selection]
                
  
#                self.add_feature("top_mrec", top_jet_mass)
    
                self.add_feature("tau_pt", region_taus.pt)
                self.add_feature("tau_eta", region_taus.eta) 
                self.add_feature("tau_phi", region_taus.phi) 

                
                self.add_feature("met_pt", region_met.pt)
                self.add_feature("met_phi", region_met.phi) 
                
                self.add_feature("b_pt", region_bjets.pt)
                self.add_feature("b_eta", region_bjets.eta) 
                self.add_feature("b_phi", region_bjets.phi) 
                
                
                
                
                nevents_after = ak.sum(final_selection)
 
                """
                print("mask_topXfinder: ",ak.sum(mask_topXfinder),  mask_topXfinder.ndim,  mask_topXfinder[mask_topXfinder == True]) 
                print("Number of events after the cuts: ", nevents_after)
                print("\n")  
                """
                

                if nevents_after > 0:
                    # -------------
                    # event weights
                    # -------------
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
                                    weightDown=events.L1PreFiringWeight.Dn[final_selection]
                                )
                            else:
                                weights_container.add(
                                    "L1Prefiring",
                                    weight=events.L1PreFiringWeight.Nom[final_selection]
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
                        njets = 1
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
                                a = 3
                                #syst_weight = weights_container.weight(variation)

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
                                    weight=syst_weight,
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
        }
        # save sumw for MC samples
        if self.is_mc:
            output["metadata"].update({"sumw": ak.sum(events.genWeight)})
        return output
        
    def postprocess(self, accumulator):
        return accumulator