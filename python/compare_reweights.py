import sys
from typing import Callable

import awkward as ak
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
from coffea import processor

NanoAODSchema.warn_missing_crossrefs = False

class DiffProcessor(processor.ProcessorABC):
    def __init__(self, extract1 : Callable, extract2 : Callable):
        self.extract1 = extract1
        self.extract2 = extract2
    
    def process(self, events):
        total = ak.sum(abs(self.extract1(events) - self.extract2(events)))
        return total

    def postprocess(self):
        pass

def main(argv):
    nreweights = 25
    infile = argv[1]
    assert infile.endswith(".root")

    events = NanoEventsFactory.from_root(
        {infile : "Events"},
        schemaclass=NanoAODSchema
    ).events()

    for ireweight in range(nreweights):
        for jreweight in range(ireweight + 1, nreweights):
            rw_name1 = "rw{:04}".format(ireweight)
            rw_name2 = "rw{:04}".format(jreweight)

            extract1 = lambda ev: getattr(ev.LHEWeight, rw_name1)
            extract2 = lambda ev: getattr(ev.LHEWeight, rw_name2)
            rw_processor = DiffProcessor(extract1, extract2)

            print("{} {} Difference: {:.6f}".format(rw_name1, rw_name2, rw_processor.process(events)))

    return 0

if __name__ == "__main__":
    print("\nFinished with exit code:", main(sys.argv))