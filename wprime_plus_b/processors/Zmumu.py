import os
import json
import pickle
import correctionlib
import numpy as np
import pandas as pd
import awkward as ak
import hist as hist2
from datetime import datetime
from typing import List, Union
from typing import Type
from coffea import util
from coffea import processor
from coffea.nanoevents.methods import candidate, vector
from coffea.analysis_tools import Weights, PackedSelection
from .utils import normalize, build_p4
from .corrections import (
    add_pileup_weight,
    add_electronID_weight,
    add_electronReco_weight,
    add_electronTrigger_weight,
    add_muon_weight,
    add_muonTriggerIso_weight,
    get_met_corrections,
)

class CandleProcessor(processor.ProcessorABC):
    def __init__(
        self,
        year: str = "2017",
        yearmod: str = "",
        channel: str = "ele",
    ):
        self._year = year
        self._yearmod = yearmod
        self._channel = channel


        # open triggers
        with open("wprime_plus_b/data/triggers.json", "r") as f:
            self._triggers = json.load(f)[self._year]

            
        # open met filters
        # https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2
        with open("wprime_plus_b/data/metfilters.json", "rb") as handle:
            self._metfilters = json.load(handle)[self._year]

            
        # open lumi masks
        with open("wprime_plus_b/data/lumi_masks.pkl", "rb") as handle:
            self._lumi_mask = pickle.load(handle)


        # output histograms
        self.make_output = lambda: {
            "sumw": 0,
            "mass": hist2.Hist(
                hist2.axis.Regular(
                    40, 50, 200, name="invariant_mass", label="$m_{ll}$ [GeV]"
                ),
                hist2.storage.Weight(),
            ),
        }



    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        dataset = events.metadata["dataset"]
        nevents = len(events)
        self.isMC = hasattr(events, "genWeight")
        self.output = self.make_output()

        # luminosity
        if not self.isMC:
            lumi_mask = self._lumi_mask[self._year](events.run, events.luminosityBlock)
        else:
            lumi_mask = np.ones(len(events), dtype="bool")
            
            
        # MET filters
        metfilters = np.ones(nevents, dtype="bool")
        metfilterkey = "mc" if self.isMC else "data"
        
        
        for mf in self._metfilters[metfilterkey]:
            if mf in events.Flag.fields:
                metfilters = metfilters & events.Flag[mf]
                
                
        # triggers
        trigger = {}
        for ch in ["ele", "mu", "tau"]:
            trigger[ch] = np.zeros(nevents, dtype="bool")
            for t in self._triggers[ch]:
                if t in events.HLT.fields:
                    trigger[ch] = trigger[ch] | events.HLT[t]
                    
                    
        # electrons
        good_electrons = (
            (events.Electron.pt >= 10)
            & (np.abs(events.Electron.eta) < 2.1)
            & events.Electron.mvaFall17V2Iso_WP90
        )
        
        n_good_electrons = ak.sum(good_electrons, axis=1)
        electrons = events.Electron[good_electrons]


        # taus: https://github.com/cms-sw/cmssw/blob/43cb25ff8489ff15adccb8016fb048bdc32d2ce5/PhysicsTools/NanoAOD/python/nanoDQM_cfi.py#L802
        good_taus = (
            (events.Electron.pt >= 20)
            & (np.abs(events.Electron.eta) < 2.1)
            & (events.Tau.idDeepTau2017v2p1VSjet > 6)          # int 1 = VVVLoose, 2 = VVLoose, 3 = VLoose, 4 = Loose, 5 = Medium, 6 = Tight, 7 = VTight, 8 = VVTight'
            & (events.Tau.idDeepTau2017v2p1VSe > 5)            # int 1 = VVVLoose, 2 = VVLoose, 3 = VLoose, 4 = Loose, 5 = Medium, 6 = Tight, 7 = VTight, 8 = VVTight
            & (events.Tau.idDeepTau2017v2p1VSmu > 4)           # int 1 = VLoose, 2 = Loose, 3 = Medium, 4 = Tight
        )
        
        n_good_taus = ak.sum(good_electrons, axis=1)
        taus = events.Electron[good_electrons]
        
        
        # muons
        good_muons = (
            (events.Muon.pt >= 35)
            & (np.abs(events.Muon.eta) < 2.1)
            & (events.Muon.tightId)
            & (
                ak.pad_none(events.Muon, 2).charge[:, 0]
                * ak.pad_none(events.Muon, 2).charge[:, 1]    # OS: charge[:,0] = First muon charge; charge[:,1] = Second muon charge
            )
        )
        n_good_muons = ak.sum(good_muons, axis=1)
        muons = events.Muon[good_muons]

        leading_muon = build_p4(ak.pad_none(muons, 2)[:, 0])
        subleading_muon = build_p4(ak.pad_none(muons, 2)[:, 1])


        # invariant mass
        invariant_mass = {
            "mu": (leading_muon + subleading_muon).mass,
        }
        
        
        # weights
        self.weights = Weights(nevents, storeIndividual=True)
        
        if self.isMC:
        
            # genweight
            self.output["sumw"] = ak.sum(events.genWeight)
            self.weights.add("genweight", events.genWeight)
            
            # L1prefiring
            if self._year in ("2016", "2017"):
                self.weights.add(
                    "L1Prefiring",
                    weight=events.L1PreFiringWeight.Nom,
                    weightUp=events.L1PreFiringWeight.Up,
                    weightDown=events.L1PreFiringWeight.Dn,
                )
            # pileup
            add_pileup_weight(
                weights=self.weights,
                year=self._year,
                mod=self._yearmod,
                nPU=ak.to_numpy(events.Pileup.nPU),
            )
            
            # muon weights
            if self._channel == "mu":
                add_muon_weight(
                    weights=self.weights,
                    muons=muons,
                    sf_type="id",
                    year=self._year,
                    mod=self._yearmod,
                    wp="tight",
                )
                add_muon_weight(
                    weights=self.weights,
                    muons=muons,
                    sf_type="iso",
                    year=self._year,
                    mod=self._yearmod,
                    wp="tight",
                )
                add_muonTriggerIso_weight(
                    weights=self.weights,
                    muons=muons,
                    year=self._year,
                    mod=self._yearmod,
                )
                
        # selections
        n_good_muons = n_good_muons == 2
        veto_electron = n_good_electrons == 0
        veto_tau = n_good_taus == 0
        self.selections = PackedSelection()
        self.selections.add("trigger_mu", trigger["mu"])
        self.selections.add("lumi", lumi_mask)
        self.selections.add("metfilters", metfilters)
        self.selections.add("two_muons", n_good_muons)
        self.selections.add("veto_electron", veto_electron)
        self.selections.add("veto_muon", veto_tau)
        self.selections.add(
            "mass_range",
            (75 < invariant_mass[self._channel])
            & (invariant_mass[self._channel] < 105),
        )

        # regions
        regions = {
            "mu": [
                "lumi",
                "metfilters",
                "trigger_mu",
                "mass_range",
                "two_muons",
                "veto_electron",
                "veto_tau",
            ],
        }

        selections = regions[self._channel]
        cut = self.selections.all(*selections)
        region_weight = self.weights.weight()[cut]

        self.output["mass"].fill(
            invariant_mass=normalize(invariant_mass[self._channel], cut),
            weight=region_weight,
        )

        return {dataset: self.output}

    def postprocess(self, accumulator):
        return accumulator
