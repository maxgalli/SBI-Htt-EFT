import sys

import uproot
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np

def modify_label(label):
    label = label.replace("-", " ")
    label = label.replace("norw", "")
    label = label.replace("noms", "")
    label = label.replace("eft", "")
    label = label.replace("stride5", "")
    
    label = label.replace(" 0p0 ", "=0p0,")
    label = label.replace(" 10p0", "=10p0")
    label = label.replace("p0", "")

    return label.strip()

def main(argv):
    # python plot_comparison.py histograms/central_VBF_hists.root,norw histograms/eft_noms_hists.root,cHbox_0p0_cHDD_0p0_ceHRe_0p0_ceHIm_0p0 histograms/eft_noms_hists.root,cHbox_0p0_cHDD_0p0_ceHRe_0p0_ceHIm_10p0
    # python plot_comparison.py histograms/central_VBF_hists_stride5.root,norw histograms/eft_noms_hists_stride5.root,cHbox_0p0_cHDD_0p0_ceHRe_0p0_ceHIm_0p0 histograms/eft_noms_hists_stride5.root,cHbox_0p0_cHDD_0p0_ceHRe_0p0_ceHIm_10p0
    # python plot_comparison.py histograms/eft_noms_hists.root,cHbox_0p0_cHDD_0p0_ceHRe_0p0_ceHIm_0p0
    # python plot_comparison.py histograms/eft_noms_hists.root,cHbox_0p0_cHDD_0p0_ceHRe_0p0_ceHIm_10p0
    histograms = {}

    for arg in argv[1:]:
        # norw for no reweight applied
        # rw#### otherwise
        file_name, rw_name = arg.split(",")
        assert file_name.endswith(".root")

        sample_name = file_name.split("/")[-1] + "-" + rw_name
        sample_name = sample_name.replace(".root", "")
        sample_name = sample_name.replace("_hists", "")
        sample_name = sample_name.replace("_", "-")

        infile = uproot.open(file_name)

        for hist_name in infile.keys():
            if hist_name.startswith(rw_name):
                hist_cat = hist_name.replace(rw_name + "_", "")
                if hist_cat not in histograms:
                    histograms[hist_cat] = {}
                
                assert sample_name not in histograms[hist_cat]
                histograms[hist_cat][sample_name] = infile[hist_name]

        infile.close()

    plt.style.use(hep.style.CMS)

    for hist_cat in histograms:
        name = hist_cat[:-2]
        print("Plotting:", name)

        fig, ax = plt.subplots()
        for isamp, sample_name in enumerate(histograms[hist_cat]):
            if histograms[hist_cat][sample_name].to_numpy()[0].sum() == 0:
                print("Empty histogram: {} - {}".format(sample_name, hist_cat))
                continue
            hep.histplot(
                histograms[hist_cat][sample_name],
                histtype="step",
                ax=ax,
                label=modify_label(sample_name),
                color="C{}".format(isamp),
                density=True,
            )

        if "phicp" in name:
            ax.set_ylim(0.1592 - 0.03, 0.1592 + 0.03)
        if "LHEtau_spin" in name:
            ax.set_xticks([0, 1, 2, 3], ["++", "+-", "-+", "--"])

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

        label = label.replace(r"absdelta$\phi$", r"$|\Delta\phi|$")
        label = label.replace(r"$cp", r"_{CP}$")
        label = label.replace("defangle1", "deflection angle")
        label = label.replace("defangle2", "deflection angle")

        ax.set_xlabel(label)
        ax.set_ylabel("Events Normalized")
        ax.legend()

        hep.cms.label(ax=ax, label="Internal", data=True, year=2023, com=13.6)

        outpath = "plots/"
        for sample_name in histograms[hist_cat]:
            outpath += sample_name + "_"
        outpath += name + ".png"

        plt.savefig(outpath)
        plt.close(fig)

        if False:
            samples = list(histograms[hist_cat].keys())
        
            # values1, edges1 = histograms[hist_cat][samples[0]]
            # values2, edges2 = histograms[hist_cat][samples[1]]
            hist1 = histograms[hist_cat][samples[0]]

            plt.style.use(hep.style.CMS)

            fig, ax = plt.subplots()

            hep.histplot(
                histograms[hist_cat][samples[0]],
                histtype="step",
                ax=ax,
                label=modify_label(samples[0]),
                color="C0",
                density=True,
            )
            hep.histplot(
                histograms[hist_cat][samples[1]],
                histtype="step",
                ax=ax,
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