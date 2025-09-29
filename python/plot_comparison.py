import sys

import uproot
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np

def main(argv):
    infile1 = argv[1]
    infile2 = argv[2]
    proc_name1 = argv[3]
    proc_name2 = argv[4]

    f1 = uproot.open(infile1)
    f2 = uproot.open(infile2)

    histograms = {name : (f1[name].to_numpy(), f2[name].to_numpy()) for name in f1.keys()}

    f1.close()
    f2.close()

    for hist_name in histograms:
        name = hist_name[:-2]
        values1, edges1 = histograms[hist_name][0]
        values2, edges2 = histograms[hist_name][1]

        # Apply CMS style
        plt.style.use(hep.style.CMS)

        fig, ax = plt.subplots()

        # Plot histograms
        hep.histplot(
            values1,
            edges1,
            ax=ax,
            histtype="step",
            label=proc_name1,
            color="C0",
            density=True
        )
        hep.histplot(
            values2,
            edges2,
            ax=ax,
            histtype="step",
            label=proc_name2,
            color="C1",
            density=True
        )

        # Labels and legend
        label = name.replace("_", " ")

        label = label.replace("GenTau GenTau", r"$\tau\tau$")
        label = label.replace("GenDressedElectron GenVisTau", r"$e\tau_h$")
        label = label.replace("GenDressedMu GenVisTau", r"$\mu\tau_h$")
        label = label.replace("GenVisTau GenVisTau", r"$\tau_h \tau_h$")
        label = label.replace("GenJet GenJet", r"$jj$")

        label = label.replace("GenTau", r"$\tau$")
        label = label.replace("GenDressedElectron", r"$e$")
        label = label.replace("GenDressedMu", r"$e$")
        label = label.replace("GenVisTau", r"$\tau_h$")
        label = label.replace("GenJet", "jet")
        
        label = label.replace("pt", r"$p_T$ [GeV]")
        label = label.replace("eta", r"$\eta$")
        label = label.replace("phi", r"$\phi$")
        label = label.replace("mass", r"$m$")
        label = label.replace("dR", r"$\Delta R$")

        ax.set_xlabel(label)
        ax.set_ylabel("Events Normalized")
        ax.legend()

        # Add CMS text
        hep.cms.label(ax=ax, label="Internal", data=True, year=2023, com=13.6)

        # Save and show
        plt.savefig("plots/{}_{}_{}.png".format(proc_name1, proc_name2, name))

    return 0


if __name__ == "__main__":
    print("\nFinished with exit code:", main(sys.argv))