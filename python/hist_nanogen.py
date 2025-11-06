import sys
from typing import Callable
import os

import awkward as ak
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
from coffea import processor
import hist
import uproot
import numpy as np

NanoAODSchema.warn_missing_crossrefs = False

class SingleObjectProcessor(processor.ProcessorABC):
    def __init__(self, extract_candidates : Callable, candidate_name : str = "", extract_reweights : Callable = None):
        self.extract = extract_candidates
        self.name = candidate_name
        self.extract_reweights = extract_reweights

    def process(self, events):
        candidates = self.extract(events)

        h_pt = hist.Hist.new.Reg(
            50, 0, 200, name="pt", label="{} $p_T$ [GeV]".format(self.name)
        ).Double()
        h_eta = hist.Hist.new.Reg(
            50, -5.0, 5.0, name="eta", label="{} $\\eta$".format(self.name)
        ).Double()
        h_phi = hist.Hist.new.Reg(
            50, -np.pi, np.pi, name="phi", label="{} $\\phi$".format(self.name)
        ).Double()
        
        h_pt.fill(pt=ak.flatten(candidates.pt))
        h_eta.fill(eta=ak.flatten(candidates.eta))
        h_phi.fill(phi=ak.flatten(candidates.phi))

        return {
            "norw_{}_pt".format(self.name) : h_pt,
            "norw_{}_eta".format(self.name) : h_eta,
            "norw_{}_phi".format(self.name) : h_phi
        }

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
        extract_reweights : Callable = None
    ):
        self.extract1 = extract_candidates1
        self.extract2 = extract_candidates2
        self.name = candidate_name
        self.extract_reweights = extract_reweights

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

        mass = (vec1 + vec2).mass
        dR = vec1.deltaR(vec2)

        h_mass = hist.Hist.new.Reg(
            50, 0, 250, name="mass", label="{} $m$ [GeV]".format(self.name)
        ).Double()
        h_dR = hist.Hist.new.Reg(
            50, 0, 5, name="dR", label="{} $\\Delta R$".format(self.name)
        ).Double()

        h_mass.fill(mass=ak.flatten(mass))
        h_dR.fill(dR=ak.flatten(dR))

        histograms = {}
        histograms["norw_{}_mass".format(self.name)] = h_mass
        histograms["norw_{}_dR".format(self.name)] = h_dR

        for rw_name in weight_sets:
            h_mass_rw = hist.Hist.new.Reg(
                50, 0, 250, name="mass", label="{} $m$ [GeV]".format(self.name)
            ).Double()
            h_dR_rw = hist.Hist.new.Reg(
                50, 0, 5, name="dR", label="{} $\\Delta R$".format(self.name)
            ).Double()

            h_mass_rw.fill(mass=ak.flatten(mass), weight=weight_sets[rw_name])
            h_dR_rw.fill(dR=ak.flatten(dR), weight=weight_sets[rw_name])

            histograms["{}_{}_mass".format(rw_name, self.name)] = h_mass_rw
            histograms["{}_{}_dR".format(rw_name, self.name)] = h_dR_rw

        return histograms

    def postprocess(self, accumulator):
        pass

class PhiCPProcessor(processor.ProcessorABC):
    def __init__(self, 
        extract_candidates1 : Callable, 
        extract_candidates2 : Callable = None, 
        candidate_name : str = "",
        extract_reweights : dict[str,Callable] = None
    ):
        self.extract1 = extract_candidates1
        self.extract2 = extract_candidates2
        self.name = candidate_name
        self.extract_reweights = extract_reweights if extract_reweights is not None else {}

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
        assert len(both_found) == len(both_valid)
        both_matched = ak.flatten(both_found & both_valid)

        assert ak.sum(both_matched) / len(both_matched) > 0.8

        tau1 = tau1[both_matched]
        tau2 = tau2[both_matched]
        candidates1 = candidates1[both_matched]
        candidates2 = candidates2[both_matched]
        weight_sets = {rw_name : self.extract_reweights[rw_name](events) for rw_name in self.extract_reweights}

        if not hasattr(candidates2, "charge"):
            raise ValueError("Candidates 2 should be tauh")
        cand2_is_pos = candidates2.charge > 0

        tau_pos = ak.where(cand2_is_pos, tau2, tau1)
        tau_neg = ak.where(cand2_is_pos, tau1, tau2)
        candidates_pos = ak.where(cand2_is_pos, candidates2, candidates1)
        candidates_neg = ak.where(cand2_is_pos, candidates1, candidates2)

        pvis_pos = make_4vector(candidates_pos)
        pvis_neg = make_4vector(candidates_neg)
        ptau_pos = make_4vector(tau_pos)
        ptau_neg = make_4vector(tau_neg)

        pH = ptau_neg + ptau_pos

        pvis_pos_rf = pvis_pos.boostCM_of_p4(pH).to_3D()
        pvis_neg_rf = pvis_neg.boostCM_of_p4(pH).to_3D()
        ptau_pos_rf = ptau_pos.boostCM_of_p4(pH).to_3D()
        ptau_neg_rf = ptau_neg.boostCM_of_p4(pH).to_3D()

        npos = ptau_pos_rf.cross(pvis_pos_rf).unit()
        nneg = ptau_neg_rf.cross(pvis_neg_rf).unit()

        ptau_pos_norm = ptau_pos_rf.unit()
        num = nneg.cross(ptau_pos_norm).dot(npos)
        den = npos.dot(nneg)

        phicp = np.arctan2(num, den)

        histograms = {}

        h_phicp = hist.Hist.new.Reg(
            50, -np.pi, np.pi, name="phicp", label="{} phi CP".format(self.name)
        ).Double()
        h_phicp.fill(phicp=ak.flatten(phicp))
        histograms["norw_{}_phicp".format(self.name)] = h_phicp
        
        for rw_name in weight_sets:
            h_phicp_rw = hist.Hist.new.Reg(
                50, -np.pi, np.pi, name="phicp", label="{} {} phi CP".format(rw_name, self.name)
            ).Double()
            h_phicp_rw.fill(phicp=ak.flatten(phicp), weight=weight_sets[rw_name])
            histograms["{}_{}_phicp".format(rw_name, self.name)] = h_phicp_rw

        return histograms

    def postprocess(self, accumulator):
        pass

def main(argv):
    # python hist_nanogen.py data/uncorrelated/cpodd3.root histograms/central_uncorr_hists.root rw0000,rw0021,rw0023
    inpath = argv[1]
    outfile = argv[2]
    reweight_names = argv[3].split(",") if len(argv) >= 4 else []

    infiles = [inpath] if inpath.endswith(".root") else [inpath + fname for fname in os.listdir(inpath) if fname.endswith(".root")]
    assert outfile.endswith(".root")
    print("Opening:", infiles)

    events = ak.concatenate([
        NanoEventsFactory.from_root(
        {fname : "Events"},
        schemaclass=NanoAODSchema
        ).events() for fname in infiles
    ])

    histograms = {}

    # Single object kinematics
    extract_gen_taus = lambda ev: ev.GenPart[(abs(ev.GenPart.pdgId) == 15) & (ev.GenPart.hasFlags(['isPrompt', 'isLastCopy']))]
    extract_dressed_elec = lambda ev: ev.GenDressedLepton[ev.GenDressedLepton.hasTauAnc & (abs(ev.GenDressedLepton.pdgId) == 11)]
    extract_dressed_mu = lambda ev: ev.GenDressedLepton[ev.GenDressedLepton.hasTauAnc & (abs(ev.GenDressedLepton.pdgId) == 13)]
    extract_tauh = lambda ev: ev.GenVisTau
    extract_gen_jets = lambda ev: ev.GenJet[ak.num(ev.GenJet) >= 2]
    extract_reweights = {rw_name : (lambda ev, name=rw_name: getattr(ev.LHEWeight, name)) for rw_name in reweight_names}

    gen_tau_processor = SingleObjectProcessor(extract_gen_taus, candidate_name="GenTau")
    dressed_elec_processor = SingleObjectProcessor(extract_dressed_elec, candidate_name="GenDressedElectron")
    dressed_mu_processor = SingleObjectProcessor(extract_dressed_mu, candidate_name="GenDressedMu")
    tauh_processor = SingleObjectProcessor(extract_tauh, candidate_name="GenVisTau")
    gen_jet_processor = SingleObjectProcessor(extract_gen_jets, candidate_name="2_leading_GenJet")

    histograms |= gen_tau_processor.process(events)
    histograms |= dressed_elec_processor.process(events)
    histograms |= dressed_mu_processor.process(events)
    histograms |= tauh_processor.process(events)
    histograms |= gen_jet_processor.process(events)

    # Diobject kinematics
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
    mutau_CP_processor = PhiCPProcessor(extract_dressed_mu, extract_tauh, "GenDressedMu_GenVisTau", extract_reweights=extract_reweights)
    etau_CP_processor = PhiCPProcessor(extract_dressed_elec, extract_tauh, "GenDressedElectron_GenVisTau", extract_reweights=extract_reweights)
    tautau_CP_processor = PhiCPProcessor(extract_tauh, candidate_name="GenVisTau_GenVisTau", extract_reweights=extract_reweights)

    histograms |= mutau_CP_processor.process(events)
    histograms |= etau_CP_processor.process(events)
    histograms |= tautau_CP_processor.process(events)

    with uproot.recreate(outfile) as fout:
        for hist_name in histograms:
            values = histograms[hist_name].values()
            edges = histograms[hist_name].axes[0].edges
            
            fout[hist_name] = (values.astype(np.float64), edges.astype(np.float64))

    return 0

if __name__ == "__main__":
    print("\nFinished with exit code:", main(sys.argv))
