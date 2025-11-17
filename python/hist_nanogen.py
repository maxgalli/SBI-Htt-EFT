import argparse
from typing import Callable
import os
import json

import awkward as ak
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
from coffea import processor
import hist
import uproot
import numpy as np

NanoAODSchema.warn_missing_crossrefs = False

def check_valid_floats(arr : ak.Array):
    arr = ak.to_numpy(arr)
    if np.any(np.isnan(arr)):
        raise ValueError("NaN detected")
    if np.any(np.isinf(arr)):
        raise ValueError("inf detected")

class SingleObjectProcessor(processor.ProcessorABC):
    def __init__(self, extract_candidates : Callable, candidate_name : str = "", extract_reweights : dict[str, Callable] | None = None):
        self.extract = extract_candidates
        self.name = candidate_name
        self.extract_reweights = extract_reweights if extract_reweights is not None else {}

    def process(self, events):
        assert self.extract_reweights is not None
        candidates = self.extract(events)
        candidates = candidates[:,:2] # At most 2 of each candidate

        weight_sets = {rw_name : self.extract_reweights[rw_name](events) for rw_name in self.extract_reweights}
        jagged_weight_sets = {rw_name : ak.broadcast_arrays(weight_sets[rw_name], candidates)[0] for rw_name in weight_sets}

        h_pt = hist.Hist.new.Reg(
            50, 0, 200, name="pt", label="{} $p_T$ [GeV]".format(self.name)
        ).Weight()
        h_eta = hist.Hist.new.Reg(
            50, -5.0, 5.0, name="eta", label="{} $\\eta$".format(self.name)
        ).Weight()
        h_phi = hist.Hist.new.Reg(
            50, -np.pi, np.pi, name="phi", label="{} $\\phi$".format(self.name)
        ).Weight()
        
        flat_pt = ak.flatten(candidates.pt)
        flat_eta = ak.flatten(candidates.eta)
        flat_phi = ak.flatten(candidates.phi)
        h_pt.fill(pt=flat_pt)
        h_eta.fill(eta=flat_eta)
        h_phi.fill(phi=flat_phi)

        histograms = {
            "norw_{}_pt".format(self.name) : h_pt,
            "norw_{}_eta".format(self.name) : h_eta,
            "norw_{}_phi".format(self.name) : h_phi
        }

        for rw_name in jagged_weight_sets:
            flat_reweights = ak.flatten(jagged_weight_sets[rw_name])
            assert len(flat_reweights) == len(flat_pt)

            rw_h_pt = hist.Hist.new.Reg(
            50, 0, 200, name="pt", label="{} {} $p_T$ [GeV]".format(rw_name, self.name)
            ).Weight()
            rw_h_eta = hist.Hist.new.Reg(
                50, -5.0, 5.0, name="eta", label="{} {} $\\eta$".format(rw_name, self.name)
            ).Weight()
            rw_h_phi = hist.Hist.new.Reg(
                50, -np.pi, np.pi, name="phi", label="{} {} $\\phi$".format(rw_name, self.name)
            ).Weight()

            rw_h_pt.fill(pt=flat_pt, weight=flat_reweights)
            rw_h_eta.fill(eta=flat_eta, weight=flat_reweights)
            rw_h_phi.fill(phi=flat_phi, weight=flat_reweights)

            histograms["{}_{}_pt".format(rw_name, self.name)] = rw_h_pt
            histograms["{}_{}_eta".format(rw_name, self.name)] = rw_h_eta
            histograms["{}_{}_phi".format(rw_name, self.name)] = rw_h_phi

        return histograms

    def postprocess(self, accumulator):
        pass

def make_4vector(cand, replace_mass=False):
    if replace_mass and hasattr(cand, "pdgId"):
        masses = ak.where(abs(cand.pdgId) == 15, 1.77686, cand.mass)
    else:
        masses = cand.mass

    return ak.zip({"pt": cand.pt, "eta": cand.eta, "phi": cand.phi, "mass": masses}, with_name="Momentum4D")

class DiobjectProcessor(processor.ProcessorABC):
    def __init__(self, 
        extract_candidates1 : Callable, 
        extract_candidates2 : Callable = None, 
        candidate_name : str = "",
        extract_reweights : dict[str, Callable] | None = None
    ):
        self.extract1 = extract_candidates1
        self.extract2 = extract_candidates2
        self.name = candidate_name
        self.extract_reweights = extract_reweights if extract_reweights is not None else {}

    def process(self, events):
        if self.extract2 is None:
            diobject = self.extract1(events)
            both_present = ak.num(diobject) >= 2
            diobject = diobject[both_present]
            candidates1 = ak.singletons(diobject[:,0])
            candidates2 = ak.singletons(diobject[:,1])
        else:
            candidates1 = self.extract1(events)
            candidates2 = self.extract2(events)
            both_present = (ak.num(candidates1) == 1) & (ak.num(candidates2) == 1)
            candidates1 = candidates1[both_present]
            candidates2 = candidates2[both_present]

        events = events[both_present]
        weight_sets = {rw_name : self.extract_reweights[rw_name](events) for rw_name in self.extract_reweights}
        vec1 = make_4vector(candidates1)
        vec2 = make_4vector(candidates2)
        
        vec_sum = vec1 + vec2
        mass = vec_sum.mass
        dR = vec1.deltaR(vec2)
        pt = vec_sum.pt
        absdeltaphi = abs(vec1.deltaphi(vec2))

        h_mass = hist.Hist.new.Reg(
            50, 0, 250, name="mass", label="{} $m$ [GeV]".format(self.name)
        ).Weight()
        h_dR = hist.Hist.new.Reg(
            50, 0, 5, name="dR", label="{} $\\Delta R$".format(self.name)
        ).Weight()
        h_pt = hist.Hist.new.Reg(
            50, 0, 200, name="pt", label="{} $p_T$".format(self.name)
        ).Weight()
        h_absdeltaphi = hist.Hist.new.Reg(
            50, 0, np.pi, name="absdeltaphi", label="{} $|\\Delta \\phi|$".format(self.name)
        ).Weight()

        h_mass.fill(mass=ak.flatten(mass))
        h_dR.fill(dR=ak.flatten(dR))
        h_pt.fill(pt=ak.flatten(pt))
        h_absdeltaphi.fill(absdeltaphi=ak.flatten(absdeltaphi))

        histograms = {}
        histograms["norw_{}_mass".format(self.name)] = h_mass
        histograms["norw_{}_dR".format(self.name)] = h_dR
        histograms["norw_{}_pt".format(self.name)] = h_pt
        histograms["norw_{}_pt".format(self.name)]

        for rw_name in weight_sets:
            h_mass_rw = hist.Hist.new.Reg(
                50, 0, 250, name="mass", label="{} $m$ [GeV]".format(self.name)
            ).Weight()
            h_dR_rw = hist.Hist.new.Reg(
                50, 0, 5, name="dR", label="{} $\\Delta R$".format(self.name)
            ).Weight()
            h_pt_rw = hist.Hist.new.Reg(
                50, 0, 200, name="pt", label="{} $p_T$".format(self.name)
            ).Weight()
            h_absdeltaphi_rw = hist.Hist.new.Reg(
                50, 0, np.pi, name="absdeltaphi", label="{} $|\\Delta \\phi|$".format(self.name)
            ).Weight()

            h_mass_rw.fill(mass=ak.flatten(mass), weight=weight_sets[rw_name])
            h_dR_rw.fill(dR=ak.flatten(dR), weight=weight_sets[rw_name])
            h_pt_rw.fill(pt=ak.flatten(pt), weight=weight_sets[rw_name])
            h_absdeltaphi_rw.fill(absdeltaphi=ak.flatten(absdeltaphi), weight=weight_sets[rw_name])

            histograms["{}_{}_mass".format(rw_name, self.name)] = h_mass_rw
            histograms["{}_{}_dR".format(rw_name, self.name)] = h_dR_rw
            histograms["{}_{}_pt".format(rw_name, self.name)] = h_pt_rw
            histograms["{}_{}_absdeltaphi".format(rw_name, self.name)] = h_absdeltaphi_rw

        return histograms

    def postprocess(self, accumulator):
        pass

class PhiCPProcessor(processor.ProcessorABC):
    def __init__(self, 
        extract_candidates1 : Callable, 
        extract_candidates2 : Callable = None, 
        candidate1_name : str = "",
        candidate2_name : str = "",
        extract_reweights : dict[str,Callable] | None = None,
        verbosity : int = 1
    ):
        self.extract1 = extract_candidates1
        self.extract2 = extract_candidates2
        self.name1 = candidate1_name
        self.name2 = candidate2_name
        self.name = self.name1 + "_" + self.name2
        self.extract_reweights = extract_reweights if extract_reweights is not None else {}
        self.verbosity = verbosity

    def match_gentau(self, events, cand):
        if hasattr(cand, "genPartIdxMother"):
            # Handles GenVisTau
            taus = events.GenPart[cand.genPartIdxMother]
            return taus

        elif hasattr(cand, "pdgId"):
            # Handles GenDressedLepton. Have to iteratively follow ancestry.
            pdg_matches = events.GenPart[events.GenPart.pdgId == ak.flatten(cand.pdgId)]
            p_pdg_matches = make_4vector(pdg_matches)
            p_cand = make_4vector(cand)

            dR = ak.flatten(p_cand).deltaR(p_pdg_matches)

            pointer = pdg_matches[ak.singletons(ak.argmin(dR, axis=1))]
            not_tau = (pointer.genPartIdxMother != -1) & (abs(events.GenPart[pointer.genPartIdxMother].pdgId) == 15)
            niter = 0
            while ak.any(not_tau):
                niter += 1
                if niter > 20:
                    raise Exception("Maximum depth reached searching for candidate ancestry")

                pointer = ak.where(not_tau, events.GenPart[pointer.genPartIdxMother], pointer)
                not_tau = (pointer.genPartIdxMother != -1) & (abs(events.GenPart[pointer.genPartIdxMother].pdgId) == 15)

            return pointer

        raise ValueError("Invalid object to match to taus")

    def process(self, events) -> dict[str, hist.Hist]:
        if self.extract2 is None:
            diobject = self.extract1(events)
            both_present = ak.num(diobject) == 2
            diobject = diobject[both_present]
            candidates1 = ak.singletons(diobject[:,0])
            candidates2 = ak.singletons(diobject[:,1])
        else:
            candidates1 = self.extract1(events)
            candidates2 = self.extract2(events)
            both_present = (ak.num(candidates1) == 1) & (ak.num(candidates2) == 1)
            candidates1 = candidates1[both_present]
            candidates2 = candidates2[both_present]

        events = events[both_present]
        assert len(events) == len(candidates1) == len(candidates2)
        assert ak.all(ak.num(candidates1) == 1) and ak.all(ak.num(candidates2) == 1)

        tau1 = self.match_gentau(events, candidates1)
        tau2 = self.match_gentau(events, candidates2)
        assert len(tau1) == len(tau2) == len(events)
        assert ak.max(ak.num(tau1)) <= 1 and ak.max(ak.num(tau2)) <= 1

        both_found = (ak.num(tau1) == 1) & (ak.num(tau2) == 1)
        both_valid = (abs(tau1.pdgId) == 15) & (abs(tau2.pdgId) == 15)
        opposite_sign = tau1.pdgId + tau2.pdgId == 0
        assert len(both_found) == len(both_valid) == len(opposite_sign)
        both_matched = ak.flatten(both_found & both_valid & opposite_sign)

        match_rate = ak.sum(both_matched) / len(both_matched)
        if self.verbosity >= 1:
            print("Match rate:", match_rate)
        assert match_rate > 0.9

        tau1 = tau1[both_matched]
        tau2 = tau2[both_matched]
        candidates1 = candidates1[both_matched]
        candidates2 = candidates2[both_matched]

        events = events[both_matched]
        weight_sets = {rw_name : self.extract_reweights[rw_name](events) for rw_name in self.extract_reweights}

        cand2_is_pos = tau2.pdgId > 0

        ptau1 = make_4vector(tau1)
        ptau2 = make_4vector(tau2)
        pvis1 = make_4vector(candidates1)
        pvis2 = make_4vector(candidates2)

        pH = ptau1 + ptau2

        ptau1_rf = ptau1.boostCM_of_p4(pH).to_3D()
        ptau2_rf = ptau2.boostCM_of_p4(pH).to_3D()
        pvis1_rf = pvis1.boostCM_of_p4(pH).to_3D()
        pvis2_rf = pvis2.boostCM_of_p4(pH).to_3D()
        
        ptau_pos_rf = ak.where(cand2_is_pos, ptau2_rf, ptau1_rf)
        ptau_neg_rf = ak.where(cand2_is_pos, ptau1_rf, ptau2_rf)
        pvis_pos_rf = ak.where(cand2_is_pos, pvis2_rf, pvis1_rf)
        pvis_neg_rf = ak.where(cand2_is_pos, pvis1_rf, pvis2_rf)

        npos = ptau_pos_rf.cross(pvis_pos_rf).unit()
        nneg = ptau_neg_rf.cross(pvis_neg_rf).unit()

        ptau_pos_norm = ptau_pos_rf.unit()
        num = nneg.cross(ptau_pos_norm).dot(npos)
        den = npos.dot(nneg)

        phicp = np.arctan2(num, den)

        defangle1 = np.arccos(ptau1_rf.unit().dot(pvis1_rf.unit()))
        defangle2 = np.arccos(ptau2_rf.unit().dot(pvis2_rf.unit()))

        h_phicp = hist.Hist.new.Reg(
            50, -np.pi, np.pi, name="phicp", label="{} phi CP".format(self.name)
        ).Weight()
        h_defangle1 = hist.Hist.new.Reg(
            50, 0, np.pi, name="theta", label="{} deflection angle".format(self.name1)
        ).Weight()
        h_defangle2 = hist.Hist.new.Reg(
            50, 0, np.pi, name="theta", label="{} deflection angle".format(self.name2)
        ).Weight()

        h_phicp.fill(phicp=ak.flatten(phicp))
        h_defangle1.fill(theta=ak.flatten(defangle1))
        h_defangle2.fill(theta=ak.flatten(defangle2))

        histograms = {
            "norw_{}_phicp".format(self.name) : h_phicp,
            "norw_{}_defangle1".format(self.name1) : h_defangle1,
            "norw_{}_defangle2".format(self.name2) : h_defangle2,
        }
        
        for rw_name in weight_sets:
            h_phicp_rw = hist.Hist.new.Reg(
                50, -np.pi, np.pi, name="phicp", label="{} {} phi CP".format(rw_name, self.name)
            ).Weight()
            h_defangle1_rw = hist.Hist.new.Reg(
                50, 0, np.pi, name="theta", label="{} {} deflection angle".format(rw_name, self.name1)
            ).Weight()
            h_defangle2_rw = hist.Hist.new.Reg(
                50, 0, np.pi, name="theta", label="{} {} deflection angle".format(rw_name, self.name2)
            ).Weight()

            h_phicp_rw.fill(phicp=ak.flatten(phicp), weight=weight_sets[rw_name])
            h_defangle1_rw.fill(theta=ak.flatten(defangle1), weight=weight_sets[rw_name])
            h_defangle2_rw.fill(theta=ak.flatten(defangle2), weight=weight_sets[rw_name])

            histograms["{}_{}_phicp".format(rw_name, self.name)] = h_phicp_rw
            histograms["{}_{}_defangle1".format(rw_name, self.name1)] = h_defangle1_rw
            histograms["{}_{}_defangle2".format(rw_name, self.name2)] = h_defangle2_rw

        return histograms

    def postprocess(self, accumulator):
        pass

class LHETauProcessor(processor.ProcessorABC):
    def __init__(self, 
        extract_candidates : Callable,
        candidate_name : str = "LHEtau",
        extract_reweights : dict[str,Callable] | None = None,
        verbosity : int = 1
    ):
        self.extract = extract_candidates
        self.name = candidate_name
        self.extract_reweights = extract_reweights if extract_reweights is not None else {}
        self.verbosity = verbosity

    def process(self, events) -> dict[str, hist.Hist]:
        candidates = self.extract(events)
        both_present = ak.num(candidates) == 2
        candidates = candidates[both_present]
        events = events[both_present]

        lhetau1 = candidates[:,0]
        lhetau2 = candidates[:,1]
        weight_sets = {rw_name : self.extract_reweights[rw_name](events) for rw_name in self.extract_reweights}

        up_up = (lhetau1.spin > 0) & (lhetau2.spin > 0)
        up_down = (lhetau1.spin > 0) & (lhetau2.spin < 0)
        down_up = (lhetau1.spin < 0) & (lhetau2.spin > 0)
        down_down = (lhetau1.spin < 0) & (lhetau2.spin < 0)

        spin_idx = ak.Array(np.zeros(len(lhetau1)))
        spin_idx = ak.where(up_up, 0, spin_idx)
        spin_idx = ak.where(up_down, 1, spin_idx)
        spin_idx = ak.where(down_up, 2, spin_idx)
        spin_idx = ak.where(down_down, 3, spin_idx)

        check_valid_floats(spin_idx)

        h_spin = hist.Hist.new.Reg(
            4, -0.5, 3.5, name="spin", label="{} spin".format(self.name)
        ).Weight()
        h_spin.fill(spin=spin_idx)
        histograms = {"norw_{}_spin".format(self.name) : h_spin}

        for rw_name in weight_sets:
            h_spin_rw = hist.Hist.new.Reg(
                4, -0.5, 3.5, name="spin", label="{} {} spin".format(rw_name, self.name)
            ).Weight()
            h_spin_rw.fill(spin=spin_idx, weight=weight_sets[rw_name])
            histograms["{}_{}_spin".format(rw_name, self.name)] = h_spin_rw

        return histograms

    def postprocess(self, accumulator):
        pass

class WeightProcessor(processor.ProcessorABC):
    def __init__(self, extract_reweights : Callable, weight_ranges_path : str):
        with open(weight_ranges_path, "r") as f:
            self.weight_ranges = json.load(f)
        self.extract_reweights = extract_reweights

    def process(self, events) -> dict[str, hist.Hist]:
        weight_sets = {rw_name : self.extract_reweights[rw_name](events) for rw_name in self.extract_reweights}
        histograms = {}
        
        for rw_name in weight_sets:
            h_weights = hist.Hist.new.Reg(
                50, 0, self.weight_ranges[rw_name][1] + 0.1, name="weight_values", label="{} weight values".format(rw_name)
            ).Weight()
            h_weights.fill(weight_values=weight_sets[rw_name])
            histograms[rw_name + "_weight_values"] = h_weights

        return histograms

    def postprocess(self, accumulator):
        pass

def serial(infile, reweight_names, histtypes):
    # python hist_nanogen.py data/uncorrelated/cpodd3.root histograms/central_uncorr_hists.root rw0000,rw0021,rw0023
    assert infile.endswith(".root")
    print("Opening:", infile)
    events = NanoEventsFactory.from_root(
        {infile : "Events"},
        schemaclass=NanoAODSchema
    ).events()

    extract_gen_taus = lambda ev: ev.GenPart[(abs(ev.GenPart.pdgId) == 15) & (ev.GenPart.hasFlags(['isPrompt', 'isLastCopy']))]
    extract_dressed_elec = lambda ev: ev.GenDressedLepton[ev.GenDressedLepton.hasTauAnc & (abs(ev.GenDressedLepton.pdgId) == 11)]
    extract_dressed_mu = lambda ev: ev.GenDressedLepton[ev.GenDressedLepton.hasTauAnc & (abs(ev.GenDressedLepton.pdgId) == 13)]
    extract_tauh = lambda ev: ev.GenVisTau
    # extract_tauh = lambda ev: ev.GenVisTau[ev.GenVisTau.status <= 3] # 1-prong
    extract_gen_jets = lambda ev: ev.GenJet
    extract_lhetaus = lambda ev: ev.LHEPart[abs(ev.LHEPart.pdgId) == 15]
    extract_reweights = {rw_name : (lambda ev, name=rw_name: getattr(ev.LHEWeight, name)) for rw_name in reweight_names}

    histograms = {}

    # Single object kinematics
    if "single" in histtypes:
        gen_tau_processor = SingleObjectProcessor(extract_gen_taus, candidate_name="GenTau", extract_reweights=extract_reweights)
        dressed_elec_processor = SingleObjectProcessor(extract_dressed_elec, candidate_name="GenDressedElectron", extract_reweights=extract_reweights)
        dressed_mu_processor = SingleObjectProcessor(extract_dressed_mu, candidate_name="GenDressedMu", extract_reweights=extract_reweights)
        tauh_processor = SingleObjectProcessor(extract_tauh, candidate_name="GenVisTau", extract_reweights=extract_reweights)
        gen_jet_processor = SingleObjectProcessor(extract_gen_jets, candidate_name="2_leading_GenJet", extract_reweights=extract_reweights)

        histograms |= gen_tau_processor.process(events)
        histograms |= dressed_elec_processor.process(events)
        histograms |= dressed_mu_processor.process(events)
        histograms |= tauh_processor.process(events)
        histograms |= gen_jet_processor.process(events)

    # Diobject kinematics
    if "diobject" in histtypes:
        gengen_processor = DiobjectProcessor(extract_gen_taus, candidate_name="GenTau_GenTau", extract_reweights=extract_reweights)
        mutau_processor = DiobjectProcessor(extract_dressed_mu, extract_tauh, candidate_name="GenDressedMu_GenVisTau", extract_reweights=extract_reweights)
        etau_processor = DiobjectProcessor(extract_dressed_elec, extract_tauh, candidate_name="GenDressedElectron_GenVisTau", extract_reweights=extract_reweights)
        tautau_processor = DiobjectProcessor(extract_tauh, candidate_name="GenVisTau_GenVisTau", extract_reweights=extract_reweights)
        dijet_processor = DiobjectProcessor(extract_gen_jets, candidate_name="GenJet_GenJet", extract_reweights=extract_reweights)

        histograms |= gengen_processor.process(events)
        histograms |= mutau_processor.process(events)
        histograms |= etau_processor.process(events)
        histograms |= tautau_processor.process(events)
        histograms |= dijet_processor.process(events)

    # Phi_CP
    if "phicp" in histtypes:
        mutau_CP_processor = PhiCPProcessor(extract_dressed_mu, extract_tauh, candidate1_name="GenDressedMu", candidate2_name="GenVisTau", extract_reweights=extract_reweights)
        etau_CP_processor = PhiCPProcessor(extract_dressed_elec, extract_tauh, candidate1_name="GenDressedElectron", candidate2_name="GenVisTau", extract_reweights=extract_reweights)
        tautau_CP_processor = PhiCPProcessor(extract_tauh, candidate1_name="GenVisTau", candidate2_name="GenVisTau", extract_reweights=extract_reweights)

        histograms |= mutau_CP_processor.process(events)
        histograms |= etau_CP_processor.process(events)
        histograms |= tautau_CP_processor.process(events)

    # LHE Tau
    if "lhetau" in histtypes:
        lhetau_spin_processor = LHETauProcessor(extract_lhetaus, extract_reweights=extract_reweights)

        histograms |= lhetau_spin_processor.process(events)

    # Weight distribution
    if "weight" in histtypes:
        weight_processor = WeightProcessor(extract_reweights, "weight_ranges.json")

        histograms |= weight_processor.process(events)

    return histograms

def main(args):
    infile = args.infile
    outfile = args.outfile
    reweight_names = args.reweights.split(",") if args.reweights is not None else []
    histtypes = set(args.histtypes.split(","))

    histograms = serial(infile, reweight_names, histtypes)

    print("Creating:", outfile)
    with uproot.recreate(outfile) as fout:
        for hist_name in histograms:
            # values = histograms[hist_name].values()
            # edges = histograms[hist_name].axes[0].edges
            
            # fout[hist_name] = (values.astype(np.float64), edges.astype(np.float64))
            fout[hist_name] = histograms[hist_name]

    return 0

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("infile")
    parser.add_argument("outfile")
    parser.add_argument("-w", "--reweights")
    parser.add_argument("-h", "--histtypes", default="single,diobject,phicp,lhetau,weight")
    parser.add_argument("-d", "--decaymodes", default="all")
    print("\nFinished with exit code:", main(parser.parse_args()))
