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
    def __init__(self, extract_candidates : Callable, candidate_name : str = ""):
        self.extract = extract_candidates
        self.name = candidate_name

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

class DiobjectProcessor(processor.ProcessorABC):
    def __init__(self, extract_candidates1 : Callable, extract_candidates2 : Callable = None, candidate_name : str = ""):
        self.extract1 = extract_candidates1
        self.extract2 = extract_candidates2
        self.name = candidate_name

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

        make_4vector = lambda x: ak.zip({"pt": x.pt, "eta": x.eta, "phi": x.phi, "mass": x.mass}, with_name="Momentum4D")
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

def main(argv):
    infile = argv[1]
    outfile = argv[2]

    events = NanoEventsFactory.from_root(
        {infile : "Events"},
        schemaclass=NanoAODSchema
    ).events()

    histograms = {}

    # Single object kinematics
    extract_gen_taus = lambda ev: ev.GenPart[(abs(ev.GenPart.pdgId) == 15) & (ev.GenPart.hasFlags(['isPrompt', 'isLastCopy']))]
    extract_dressed_elec = lambda ev: ev.GenDressedLepton[ev.GenDressedLepton.hasTauAnc & (abs(ev.GenDressedLepton.pdgId) == 11)]
    extract_dressed_mu = lambda ev: ev.GenDressedLepton[ev.GenDressedLepton.hasTauAnc & (abs(ev.GenDressedLepton.pdgId) == 13)]
    extract_tauh = lambda ev: ev.GenVisTau
    extract_gen_jets = lambda ev: ev.GenJet[ak.num(ev.GenJet) >= 2][:,:2]

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

    with uproot.recreate(outfile) as fout:
        for hist_name in histograms:
            values = histograms[hist_name].values()
            edges = histograms[hist_name].axes[0].edges
            
            fout[hist_name] = (values.astype(np.float64), edges.astype(np.float64))

    return 0

if __name__ == "__main__":
    print("\nFinished with exit code:", main(sys.argv))
