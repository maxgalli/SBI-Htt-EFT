import sys
from typing import Callable

import awkward as ak
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
from coffea import processor
import hist
import uproot
import numpy as np

NanoAODSchema.warn_missing_crossrefs = False

class SingleObjectProcessor(processor.ProcessorABC):
    def __init__(self, extract_candidates : Callable, candidate_name : str = "", extract_weights : Callable = None):
        self.extract = extract_candidates
        self.name = candidate_name
        self.extract_weights = extract_weights

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
            "{}_pt".format(self.name) : h_pt,
            "{}_eta".format(self.name) : h_eta,
            "{}_phi".format(self.name) : h_phi
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
        extract_weights : Callable = None
    ):
        self.extract1 = extract_candidates1
        self.extract2 = extract_candidates2
        self.name = candidate_name
        self.extract_weights = extract_weights

    def process(self, events):
        if self.extract2 is None:
            diobject = self.extract1(events)
            diobject = diobject[ak.num(diobject) == 2]
            candidates1 = diobject[:,0]
            candidates2 = diobject[:,1]
        else:
            candidates1 = self.extract1(events)
            candidates2 = self.extract2(events)
            both_present = (ak.num(candidates1) == 1) & (ak.num(candidates2) == 1)
            candidates1 = candidates1[both_present]
            candidates2 = candidates2[both_present]

        vec1 = make_4vector(candidates1)
        vec2 = make_4vector(candidates2)

        h_mass = hist.Hist.new.Reg(
            50, 0, 250, name="mass", label="{} $m$ [GeV]".format(self.name)
        ).Double()
        h_dR = hist.Hist.new.Reg(
            50, 0, 5, name="dR", label="{} $\\Delta R$".format(self.name)
        ).Double()

        mass = (vec1 + vec2).mass
        dR = vec1.deltaR(vec2)
        if self.extract2 is not None:
            mass = ak.flatten(mass)
            dR = ak.flatten(dR)

        h_mass.fill(mass=mass)
        h_dR.fill(dR=dR)

        return {
            "{}_mass".format(self.name) : h_mass,
            "{}_dR".format(self.name) : h_dR
        }

    def postprocess(self, accumulator):
        pass

class PhiCPProcessor(processor.ProcessorABC):
    def __init__(self, 
        extract_gentaus : Callable, 
        extract_candidates1 : Callable, 
        extract_candidates2 : Callable = None, 
        candidate_name : str = "",
        extract_weights : Callable = None
    ):
        self.extract1 = extract_candidates1
        self.extract2 = extract_candidates2
        self.extract_taus = extract_gentaus
        self.name = candidate_name
        self.extract_weights = extract_weights

    def match_gentau(self, events, cand):
        if hasattr(cand, "genPartIdxMother"):
            taus = events.GenPart[ak.singletons(cand.genPartIdxMother)]
            taus = ak.mask(taus, cand.genPartIdxMother >= 0)
            if not ak.all(abs(taus.pdgId) == 15):
                # print(cand.genPartIdxMother)
                # print(ak.singletons(cand.genPartIdxMother))
                # print()
                # for elem in events.GenPart[0].pdgId:
                #     print(elem)
                raise Exception("Hadronig tau index mismatch")
            return taus

        gentaus = self.extract_taus(events)

        pvis = make_4vector(cand)
        ptau = make_4vector(gentaus)

        dR = pvis.deltaR(ptau)
        return gentaus[ak.singletons(ak.argmin(dR, axis=1))]

    def process(self, events):
        if self.extract2 is None:
            diobject = self.extract1(events)
            diobject = diobject[ak.num(diobject) == 2]
            candidates1 = diobject[:,0]
            candidates2 = diobject[:,1]
            events = events[ak.num(diobject) == 2]
        else:
            candidates1 = self.extract1(events)
            candidates2 = self.extract2(events)
            both_present = (ak.num(candidates1) == 1) & (ak.num(candidates2) == 1)
            candidates1 = ak.flatten(candidates1[both_present])
            candidates2 = ak.flatten(candidates2[both_present])
            events = events[both_present]

        print(candidates1)
        print(candidates2)

        tau1 = self.match_gentau(events, candidates1)

        print

        tau2 = self.match_gentau(events, candidates2)
        # print(ak.num(tau1))
        # print(ak.num(tau2))
        tau_mismatch_filter = ak.flatten(tau1.deltaR(tau2) > 0.01) # Both candidates matched to the same GenPart tau

        vis_p1 = make_4vector(candidates1)[tau_mismatch_filter]
        vis_p2 = make_4vector(candidates2)[tau_mismatch_filter]
        tau_p1 = make_4vector(tau1)[tau_mismatch_filter]
        tau_p2 = make_4vector(tau2)[tau_mismatch_filter]

        assert ak.all(ak.num(tau_p1) == 1) and ak.all(ak.num(tau_p2) == 1)
        tau_p1 = ak.flatten(tau_p1)
        tau_p2 = ak.flatten(tau_p2)
    
        H_p = tau_p1 + tau_p2

        vis1_rf = vis_p1.boostCM_of_p4(H_p).to_3D()
        vis2_rf = vis_p2.boostCM_of_p4(H_p).to_3D()
        tau1_rf = tau_p1.boostCM_of_p4(H_p).to_3D()
        tau2_rf = tau_p2.boostCM_of_p4(H_p).to_3D()

        n1 = tau1_rf.cross(vis1_rf).unit()
        n2 = tau2_rf.cross(vis2_rf).unit()

        if hasattr(candidates2, "charge"):
            tau_pos_norm = ak.where(candidates2.charge > 0, tau2_rf, tau1_rf).unit()
        else:
            raise ValueError("Candidates 2 should be tauh")

        num = n2.cross(tau_pos_norm).dot(n1)
        den = n1.dot(n2)

        phicp = np.arctan2(num, den)

        h_phicp = hist.Hist.new.Reg(
            50, -np.pi, np.pi, name="phicp", label="{} phi CP".format(self.name)
        ).Double()

        if self.extract_weights is None:
            h_phicp.fill(phicp=phicp)
        else:
            h_phicp.fill(phicp=phicp, weights=self.extract_weights(events))

        return {"{}_phicp".format(self.name) : h_phicp}

    def postprocess(self, accumulator):
        pass

def main(argv):
    infile = argv[1]
    outfile = argv[2]

    events = NanoEventsFactory.from_root(
        {infile : "Events"},
        schemaclass=NanoAODSchema
    ).events()

    events = events[:100]

    histograms = {}

    # Single object kinematics
    extract_gen_taus = lambda ev: ev.GenPart[(abs(ev.GenPart.pdgId) == 15) & (ev.GenPart.hasFlags(['isPrompt', 'isLastCopy']))]
    extract_dressed_elec = lambda ev: ev.GenDressedLepton[ev.GenDressedLepton.hasTauAnc & (abs(ev.GenDressedLepton.pdgId) == 11)]
    extract_dressed_mu = lambda ev: ev.GenDressedLepton[ev.GenDressedLepton.hasTauAnc & (abs(ev.GenDressedLepton.pdgId) == 13)]
    extract_tauh = lambda ev: ev.GenVisTau
    extract_gen_jets = lambda ev: ev.GenJet[ak.num(ev.GenJet) >= 2][:,:2]

    tautau_CP_processor = PhiCPProcessor(extract_gen_taus, extract_tauh, candidate_name="GenDressedMu_GenVisTau")
    histograms |= tautau_CP_processor.process(events)
    return 1

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
    gengen_processor = DiobjectProcessor(extract_gen_taus, candidate_name="GenTau_GenTau")
    mutau_processor = DiobjectProcessor(extract_dressed_mu, extract_tauh, candidate_name="GenDressedMu_GenVisTau")
    etau_processor = DiobjectProcessor(extract_dressed_elec, extract_tauh, candidate_name="GenDressedElectron_GenVisTau")
    tautau_processor = DiobjectProcessor(extract_tauh, candidate_name="GenVisTau_GenVisTau")
    dijet_processor = DiobjectProcessor(extract_gen_jets, candidate_name="GenJet_GenJet")

    histograms |= gengen_processor.process(events)
    histograms |= mutau_processor.process(events)
    histograms |= etau_processor.process(events)
    histograms |= tautau_processor.process(events)
    histograms |= dijet_processor.process(events)

    # Phi_CP
    mutau_CP_processor = PhiCPProcessor(extract_gen_taus, extract_dressed_mu, extract_tauh, "GenDressedMu_GenVisTau")
    etau_CP_processor = PhiCPProcessor(extract_gen_taus, extract_dressed_elec, extract_tauh, "GenDressedMu_GenVisTau")
    tautau_CP_processor = PhiCPProcessor(extract_gen_taus, extract_tauh, candidate_name="GenDressedMu_GenVisTau")

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
