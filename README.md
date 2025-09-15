# SBI $H \rightarrow \tau \tau$ EFT

## Setup

To setup the environment and the necessary software stack (on lxplus) run:
```
. setup.sh
```

## Generation

This part describes the instructions to generate gridpacks, LHE, nanogen and
nanoaod files for a certain process. We will consider as an example ```ppH_Htt_SMEFTsim_topU3l_quadratic```, but it will be the same for all of them.

### Creating gridpack

The creation of gridpacks follows the standard procedure explained
[here](https://cms-generators.docs.cern.ch/how-to-produce-gridpacks/mg5-amcnlo/). The MG cards are thus stored inside ```madgraph-cards```, divided per process, with the naming scheme required by the genproduction scripts.
First, access the correct path:
```
cd genproductions_scripts/bin/MadGraph5_aMCatNLO
```

As exaplined in the documentation, gridpacks for Run 3 should be produced using
```cmssw-el8```. Following the instructions
[here](https://cms-sw.github.io/singularity.html), on lxplus we just activate it
by typing ```cmssw-el8```.
To generate the gridpack we then run:
```
./gridpack_generation.sh ppH_Htt_SMEFTsim_topU3l_quadratic ../../../madgraph-cards/ppH_Htt_SMEFTsim_topU3l_quadratic
```
This will produce a gridpack called
```ppH_Htt_SMEFTsim_topU3l_quadratic_el8_amd64_gcc10_CMSSW_12_4_8_tarball.tar.xz```.
Note that the procedure runs locally; in order to run this procedure on
HTCondor, check out the documentation for
```submit_condor_gridpack_generation.sh```.
