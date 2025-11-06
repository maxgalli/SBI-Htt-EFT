import sys

import uproot
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np

def modify_label(label):
    label = label.replace("_", " ")
    label = label.replace("norw", "")
    label = label.replace("ms", "MS")
    label = label.replace("eft", "EFT")
    label = label.replace("stride5", "")
    label = label.replace("rw0000", "RW:0,0,0,0")
    label = label.replace("rw0023", "RW:0,0,0,10")
    label = label.replace("rw0022", "RW:100,100,100,100")

    return label

def main(argv):
    # python plot_comparison.py histograms/central_VBF_hists_stride5.root,norw histograms/central_uncorr_hists_stride5.root,norw
    # python plot_comparison.py histograms/central_VBF_hists_stride5.root,norw histograms/eft_noms_hists_stride5.root,rw0000
    # python plot_comparison.py histograms/eft_noms_hists_stride5.root,rw0000 histograms/eft_noms_hists_stride5.root,rw0023
    histograms = {}

    for arg in argv[1:]:
        # norw for no reweight applied
        # rw#### otherwise
        file_name, rw_name = arg.split(",")
        assert file_name.endswith(".root")

        sample_name = file_name.split("/")[-1] + "_" + rw_name
        sample_name = sample_name.replace(".root", "")
        sample_name = sample_name.replace("_hists", "")

        infile = uproot.open(file_name)

        for hist_name in infile.keys():
            if hist_name.startswith(rw_name):
                hist_cat = hist_name.replace(rw_name + "_", "")
                if hist_cat not in histograms:
                    histograms[hist_cat] = {}
                
                assert sample_name not in histograms[hist_cat]
                histograms[hist_cat][sample_name] = infile[hist_name].to_numpy()

        infile.close()

    for hist_cat in histograms:
        name = hist_cat[:-2]

        if not len(histograms[hist_cat]) == 2:
            continue

        print("Plotting:", name)
        samples = list(histograms[hist_cat].keys())
        
        values1, edges1 = histograms[hist_cat][samples[0]]
        values2, edges2 = histograms[hist_cat][samples[1]]

        plt.style.use(hep.style.CMS)

        fig, ax = plt.subplots()

        hep.histplot(
            values1,
            edges1,
            ax=ax,
            histtype="step",
            label=modify_label(samples[0]),
            color="C0",
            density=True
        )
        hep.histplot(
            values2,
            edges2,
            ax=ax,
            histtype="step",
            label=modify_label(samples[1]),
            color="C1",
            density=True
        )

        if "phicp" in name:
            ax.set_ylim(0.1592 - 0.03, 0.1592 + 0.03)

        label = name.replace("_", " ")

        label = label.replace("GenTau GenTau", r"$\tau\tau$")
        label = label.replace("GenDressedElectron GenVisTau", r"$e\tau_h$")
        label = label.replace("GenDressedMu GenVisTau", r"$\mu\tau_h$")
        label = label.replace("GenVisTau GenVisTau", r"$\tau_h \tau_h$")
        label = label.replace("GenJet GenJet", r"$jj$")

        label = label.replace("GenTau", r"$\tau$")
        label = label.replace("GenDressedElectron", r"$e$")
        label = label.replace("GenDressedMu", r"$\mu$")
        label = label.replace("GenVisTau", r"$\tau_h$")
        label = label.replace("GenJet", "jet")
        
        label = label.replace("pt", r"$p_T$ [GeV]")
        label = label.replace("eta", r"$\eta$")
        label = label.replace("phi", r"$\phi$")
        label = label.replace("mass", r"$m$")
        label = label.replace("dR", r"$\Delta R$")

        label = label.replace(r"$cp", r"_{CP}$")

        ax.set_xlabel(label)
        ax.set_ylabel("Events Normalized")
        ax.legend()

        hep.cms.label(ax=ax, label="Internal", data=True, year=2023, com=13.6)

        plt.savefig("plots/{}_{}_{}.png".format(samples[0], samples[1], name))
        plt.close(fig)

    return 0


if __name__ == "__main__":
    print("\nFinished with exit code:", main(sys.argv))