# SBI $H \rightarrow \tau \tau$ EFT samples production

## Setup

To setup the environment and the necessary software stack (on lxplus) run:
```
cmssw-el8
. setup.sh
```

## Generation

This part describes the instructions to generate gridpacks, LHE, nanogen and
nanoaod files for a certain process. We will consider as an example ```ppH_Htt_SMEFTsim_topU3l_quadratic```, but it will be the same for all of them.

### MG cards and EFT (TODO)

### Creating gridpack

The creation of gridpacks follows the standard procedure explained [here](https://cms-generators.docs.cern.ch/how-to-produce-gridpacks/mg5-amcnlo/). The MG cards are thus stored inside ```madgraph-cards```, divided per process, with the naming scheme required by the genproduction scripts.
First, access the correct path:
```
cd genproductions_scripts/bin/MadGraph5_aMCatNLO
```

As exaplined in the documentation, gridpacks for Run 3 should be produced using ```cmssw-el8```. Following the instructions [here](https://cms-sw.github.io/singularity.html), on lxplus we just activate it by typing ```cmssw-el8```.
To generate the gridpack we then run:
```
./gridpack_generation.sh ppH_Htt_SMEFTsim_topU3l_quadratic ../../../madgraph-cards/ppH_Htt_SMEFTsim_topU3l_quadratic
```
This will produce a gridpack called ```ppH_Htt_SMEFTsim_topU3l_quadratic_el8_amd64_gcc10_CMSSW_12_4_8_tarball.tar.xz```.
Note that the procedure runs locally; in order to run this procedure on HTCondor, check out the documentation for ```submit_condor_gridpack_generation.sh```.

### Creating EDM GEN files

This part will require to exit from ```cmssw-el8```. First of all, change the path inside the pythia fragment at ```CMSSW_13_2_9/src/Configuration/GenProduction/python/pythia_fragment.py``` to the one that was just created. 
We can then create the config file by running (taken from the [cmseft repo](https://github.com/FNALLPC/cmseft/tree/main)):
```
cmsDriver.py Configuration/GenProduction/python/pythia_fragment.py \
    --mc \
    --python_filename ppH_Htt_SMEFTsim_topU3l_quadratic_gen_cfg.py \
    --eventcontent RAWSIM,LHE \
    --datatier GEN,LHE \
    --conditions 130X_mcRun3_2023_realistic_postBPix_v5 \
    --beamspot Realistic25ns13p6TeVEarly2023Collision \
    --step LHE,GEN \
    --nThreads 1 \
    --geometry DB:Extended \
    --era Run3_2023 \
    --customise Configuration/DataProcessing/Utils.addMonitoring \
    --customise_commands process.RandomNumberGeneratorService.externalLHEProducer.initialSeed=123 \
    --fileout file:gen.root \
    --no_exec -n 100
```
and run it with:
```
cmsRun ppH_Htt_SMEFTsim_topU3l_quadratic_gen_cfg.py
```

### Creating NanoGEN files

The NanoGEN format can be used for exploratory studies (not for the full analysis setup, as in this case we need reco-level quantities). A NanoGEN config file can be produced with:
```
cmsDriver.py Configuration/GenProduction/python/pythia_fragment.py \
    --python_filename ppH_Htt_SMEFTsim_topU3l_quadratic_nanogen_cfg.py --eventcontent NANOAODGEN \
    --customise Configuration/DataProcessing/Utils.addMonitoring --datatier NANOAOD \
    --customise_commands process.RandomNumberGeneratorService.externalLHEProducer.initialSeed=123 \
    --fileout file:nanogen.root --conditions 130X_mcRun3_2023_realistic_postBPix_v5 --beamspot Realistic25ns13p6TeVEarly2023Collision \
    --step LHE,GEN,NANOGEN --geometry DB:Extended --era Run3_2023 --no_exec --mc -n 100
```

The different weights, taken from the ```reweight_card.dat```, can be saved by adding lines like
```
named_weights = [
    "rw0000",
    "rw0001",
]
process.genWeightsTable.namedWeightIDs = named_weights
process.genWeightsTable.namedWeightLabels = named_weights
```
(TODO script to aumomate will be written).

#### Scale out (TODO)

The two previous steps can be unified into one as follows:
```
cmsDriver.py Configuration/GenProduction/python/pythia_fragment.py \
    --python_filename ppH_Htt_SMEFTsim_topU3l_quadratic_gennanogen_cfg.py \
    --eventcontent NANOAODGEN \
    --datatier NANOAOD \
    --conditions 130X_mcRun3_2023_realistic_postBPix_v5 \
    --beamspot Realistic25ns13p6TeVEarly2023Collision \
    --step LHE,GEN,NANOGEN \
    --geometry DB:Extended \
    --era Run3_2023 \
    --customise Configuration/DataProcessing/Utils.addMonitoring \
    --fileout file:gen.root \
    --no_exec --mc -n 100
```
Also in this case, one has to add the columns that contain the weights. Moreover, the two following modifications are necessary:
- move the gridpacks to ```/eos/user/g/gallim/www/EFTstudies/gridpacks``` and change the path to the first step to ```root://eosuser.cern.ch//eos/user/g/gallim/www/EFTstudies/gridpacks/<gpackname>``` when submitting on crab
- for this to work, also change the line ```scriptName = cms.FileInPath('GeneratorInterface/LHEInterface/data/run_generic_tarball_cvmfs.sh')``` to ```scriptName = cms.FileInPath('GeneratorInterface/LHEInterface/data/run_generic_tarball_xrootd.sh')```

A ```crab_config.py``` file then has to be written (it is alread present in the repo for the dummy example) and production at scale can be submitted with:
```
crab submit -c crab_config.py
```

### Creating NanoAOD files

In this part we summarize the procedure to produce nanoAOD files for a sample (i.e. with reconstruction included). The procedure can be find in different twikis, for convenience we base ours on [this](https://indico.cern.ch/event/1500035/contributions/6575125/attachments/3091743/5476084/IITH_GEN_Tutorial.pdf).
Note that, in order for the commands to work and have tags and CMSSW versions that are compatible with each other, the best idea is to take a NanoAOD sample similar to what we want to generate and copy the config commands.
This example is made using [this workflow](https://cms-pdmv-prod.web.cern.ch/mcm/chained_requests?prepid=HIG-chain_Run3Summer23BPixwmLHEGS_flowRun3Summer23BPixDRPremix_flowRun3Summer23BPixMiniAODv4_flowRun3Summer23BPixNanoAODv13-00001&page=0&shown=15), but could be different for other processes.

GENSIM step:
```
cmsDriver.py Configuration/GenProduction/python/pythia_fragment.py \
    --python_filename ppH_Htt_SMEFTsim_topU3l_quadratic_gensim_cfg.py \
    --eventcontent RAWSIM,LHE \
    --customise Configuration/DataProcessing/Utils.addMonitoring --datatier GEN-SIM,LHE \
    --fileout file:gensim.root --conditions 130X_mcRun3_2023_realistic_postBPix_v5 --beamspot Realistic25ns13p6TeVEarly2023Collision \
    --step LHE,GEN,SIM --geometry DB:Extended --era Run3_2023 --no_exec --mc -n -1
```
For the scaling to work, it is also convenient to make the number of events processed customizable via command line (this is taken from Davide's samples production for EFT ttHbb). We thus add these lines:
```
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('analysis')
options.register('jobNum', 0, VarParsing.multiplicity.singleton,VarParsing.varType.int,"jobNum")
options.register('nEvents', 0, VarParsing.multiplicity.singleton,VarParsing.varType.int,"nEvents")
options.parseArguments()
```
and change the following lines:
```
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.nEvents))
```
```
# Input source
process.source = cms.Source(
        "EmptySource",
        firstLuminosityBlock  = cms.untracked.uint32(options.jobNum+1),
        numberEventsInLuminosityBlock = cms.untracked.uint32(options.nEvents)
)
```
```nEvents = cms.untracked.uint32(options.nEvents)``` in ```process.externalLHEProducer```.

One can then start the production with:
```
cmsRun ppH_Htt_SMEFTsim_topU3l_quadratic_gensim_cfg.py
```

DIGI-RAW step (TODO resources for careful and meaningful inclusion of pileup and HLT have to be found):
```
cmsDriver.py \
    --python_filename ppH_Htt_SMEFTsim_topU3l_quadratic_digi_cfg.py \
    --eventcontent PREMIXRAW \
    --customise Configuration/DataProcessing/Utils.addMonitoring --datatier GEN-SIM-RAW \
    --fileout file:digi.root \
    --pileup_input "dbs:/Neutrino_E-10_gun/Run3Summer21PrePremix-Summer23BPix_130X_mcRun3_2023_realistic_postBPix_v1-v1/PREMIX" \
    --conditions 130X_mcRun3_2023_realistic_postBPix_v5 \
    --step DIGI,DATAMIX,L1,DIGI2RAW,HLT:2023v12 \
    --procModifiers premix_stage2,siPixelQualityRawToDigi --geometry DB:Extended --filein file:gensim.root \
    --datamix PreMix --era Run3_2023 --no_exec --mc -n -1
```
```
cmsRun ppH_Htt_SMEFTsim_topU3l_quadratic_digi_cfg.py
```

RECO/AOD step:
```
cmsDriver.py \
    --python_filename ppH_Htt_SMEFTsim_topU3l_quadratic_reco_cfg.py \
    --eventcontent AODSIM --customise Configuration/DataProcessing/Utils.addMonitoring --datatier AODSIM \
    --fileout file:aod.root --conditions 130X_mcRun3_2023_realistic_postBPix_v5 --step RAW2DIGI,L1Reco,RECO,RECOSIM \
    --geometry DB:Extended --filein file:digi.root --era Run3_2023 --no_exec --mc -n -1
```
```
cmsRun ppH_Htt_SMEFTsim_topU3l_quadratic_reco_cfg.py
```

MINIAOD step:
```
cmsDriver.py \
    --python_filename ppH_Htt_SMEFTsim_topU3l_quadratic_miniaod_cfg.py --eventcontent MINIAODSIM \
    --customise Configuration/DataProcessing/Utils.addMonitoring --datatier MINIAODSIM \
    --fileout file:miniaod.root --conditions 130X_mcRun3_2023_realistic_postBPix_v5 \
    --step PAT --geometry DB:Extended --filein file:aod.root --era Run3_2023 \
    --no_exec --mc -n -1
```
```
cmsRun ppH_Htt_SMEFTsim_topU3l_quadratic_miniaod_cfg.py
```

NANOAOD step:
```
cmsDriver.py  \
    --eventcontent NANOAODSIM --customise Configuration/DataProcessing/Utils.addMonitoring --datatier NANOAODSIM \
    --conditions 133X_mcRun3_2023_realistic_postBPix_ForNanov13_v2 --step NANO \
    --scenario pp --era Run3_2023,run3_nanoAOD_124 --python_filename ppH_Htt_SMEFTsim_topU3l_quadratic_nanoaod_cfg.py \
    --fileout file:nanoaod.root --filein file:miniaod.root --no_exec --mc -n -1
```
Also in this case, the different weights, taken from the ```reweight_card.dat```, can be saved by adding lines like
```
named_weights = [
    "rw0000",
    "rw0001",
]
process.genWeightsTable.namedWeightIDs = named_weights
process.genWeightsTable.namedWeightLabels = named_weights
```

Remember, for the purpose of having the dataset published on DBS, to add the following line to the config file:
```
fakeNameForCrab = cms.untracked.bool(True)
```
right underneath the line ```outputCommands = process.NANOAODSIMEventContent.outputCommands```.

Finally, run:
```
cmsRun ppH_Htt_SMEFTsim_topU3l_quadratic_nanoaod_cfg.py
```

#### Scale out (TODO)

Here we show how to use CRAB to produce NanoAOD with the relevant EFT weights. Relevant documentation can be found [here](https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3AdvancedTopic). As it can be seen in the link, the procedure consists in using a `scriptExe.sh` file where the different stages are chained. In order to control the number of events produced, extra lines to configure via command line are added to the config file of the first stage, while the others simpli take `-1` (i.e., hopefully, all the events coming from the previous steps).

Important notes:

- remember to move the gridpacks to ```/eos/user/g/gallim/www/EFTstudies/gridpacks``` and change the path to the first gensim step to ```root://eosuser.cern.ch//eos/user/g/gallim/www/EFTstudies/gridpacks/<gpackname>``` when submitting on crab
- for this to work, also change the line ```scriptName = cms.FileInPath('GeneratorInterface/LHEInterface/data/run_generic_tarball_cvmfs.sh')``` to ```scriptName = cms.FileInPath('GeneratorInterface/LHEInterface/data/run_generic_tarball_xrootd.sh')```

## Other useful links

- generators short exercise at CMS DAS ([link](https://fnallpc.github.io/generators/))
- presentation for EFT samples production at LPC school ([link](https://indico.fnal.gov/event/68174/timetable/?view=standard#26-mc-generation-tutorial))
- Ilaria Brivio's SMEFTsim tutorial ([link](https://github.com/IlariaBrivio/SMEFTsim_tutorial_June23))# SBI $H \rightarrow \tau \tau$ EFT samples production

## Setup

To setup the environment and the necessary software stack (on lxplus) run:
```
cmssw-el8
. setup.sh
```

## Generation

This part describes the instructions to generate gridpacks, LHE, nanogen and
nanoaod files for a certain process. We will consider as an example ```ppH_Htt_SMEFTsim_topU3l_quadratic```, but it will be the same for all of them.

### MG cards and EFT (TODO)

### Creating gridpack

The creation of gridpacks follows the standard procedure explained [here](https://cms-generators.docs.cern.ch/how-to-produce-gridpacks/mg5-amcnlo/). The MG cards are thus stored inside ```madgraph-cards```, divided per process, with the naming scheme required by the genproduction scripts.
First, access the correct path:
```
cd genproductions_scripts/bin/MadGraph5_aMCatNLO
```

As exaplined in the documentation, gridpacks for Run 3 should be produced using ```cmssw-el8```. Following the instructions [here](https://cms-sw.github.io/singularity.html), on lxplus we just activate it by typing ```cmssw-el8```.
To generate the gridpack we then run:
```
./gridpack_generation.sh ppH_Htt_SMEFTsim_topU3l_quadratic ../../../madgraph-cards/ppH_Htt_SMEFTsim_topU3l_quadratic
```
This will produce a gridpack called ```ppH_Htt_SMEFTsim_topU3l_quadratic_el8_amd64_gcc10_CMSSW_12_4_8_tarball.tar.xz```.
Note that the procedure runs locally; in order to run this procedure on HTCondor, check out the documentation for ```submit_condor_gridpack_generation.sh```.

### Creating EDM GEN files

This part will require to exit from ```cmssw-el8```. First of all, change the path inside the pythia fragment at ```CMSSW_13_2_9/src/Configuration/GenProduction/python/pythia_fragment.py``` to the one that was just created. 
We can then create the config file by running (taken from the [cmseft repo](https://github.com/FNALLPC/cmseft/tree/main)):
```
cmsDriver.py Configuration/GenProduction/python/pythia_fragment.py \
    --mc \
    --python_filename ppH_Htt_SMEFTsim_topU3l_quadratic_gen_cfg.py \
    --eventcontent RAWSIM,LHE \
    --datatier GEN,LHE \
    --conditions 130X_mcRun3_2023_realistic_postBPix_v5 \
    --beamspot Realistic25ns13p6TeVEarly2023Collision \
    --step LHE,GEN \
    --nThreads 1 \
    --geometry DB:Extended \
    --era Run3_2023 \
    --customise Configuration/DataProcessing/Utils.addMonitoring \
    --customise_commands process.RandomNumberGeneratorService.externalLHEProducer.initialSeed=123 \
    --fileout file:gen.root \
    --no_exec -n 100
```
and run it with:
```
cmsRun ppH_Htt_SMEFTsim_topU3l_quadratic_gen_cfg.py
```

### Creating NanoGEN files

The NanoGEN format can be used for exploratory studies (not for the full analysis setup, as in this case we need reco-level quantities). A NanoGEN config file can be produced with:
```
cmsDriver.py Configuration/GenProduction/python/pythia_fragment.py \
    --python_filename ppH_Htt_SMEFTsim_topU3l_quadratic_nanogen_cfg.py --eventcontent NANOAODGEN \
    --customise Configuration/DataProcessing/Utils.addMonitoring --datatier NANOAOD \
    --customise_commands process.RandomNumberGeneratorService.externalLHEProducer.initialSeed=123 \
    --fileout file:nanogen.root --conditions 130X_mcRun3_2023_realistic_postBPix_v5 --beamspot Realistic25ns13p6TeVEarly2023Collision \
    --step LHE,GEN,NANOGEN --geometry DB:Extended --era Run3_2023 --no_exec --mc -n 100
```

The different weights, taken from the ```reweight_card.dat```, can be saved by adding lines like
```
named_weights = [
    "rw0000",
    "rw0001",
]
process.genWeightsTable.namedWeightIDs = named_weights
process.genWeightsTable.namedWeightLabels = named_weights
```
(TODO script to aumomate will be written).

#### Scale out (TODO)

The two previous steps can be unified into one as follows:
```
cmsDriver.py Configuration/GenProduction/python/pythia_fragment.py \
    --python_filename ppH_Htt_SMEFTsim_topU3l_quadratic_gennanogen_cfg.py \
    --eventcontent NANOAODGEN \
    --datatier NANOAOD \
    --conditions 130X_mcRun3_2023_realistic_postBPix_v5 \
    --beamspot Realistic25ns13p6TeVEarly2023Collision \
    --step LHE,GEN,NANOGEN \
    --geometry DB:Extended \
    --era Run3_2023 \
    --customise Configuration/DataProcessing/Utils.addMonitoring \
    --fileout file:gen.root \
    --no_exec --mc -n 100
```
Also in this case, one has to add the columns that contain the weights. Moreover, the two following modifications are necessary:
- move the gridpacks to ```/eos/user/g/gallim/www/EFTstudies/gridpacks``` and change the path to the first step to ```root://eosuser.cern.ch//eos/user/g/gallim/www/EFTstudies/gridpacks/<gpackname>``` when submitting on crab
- for this to work, also change the line ```scriptName = cms.FileInPath('GeneratorInterface/LHEInterface/data/run_generic_tarball_cvmfs.sh')``` to ```scriptName = cms.FileInPath('GeneratorInterface/LHEInterface/data/run_generic_tarball_xrootd.sh')```

A ```crab_config.py``` file then has to be written (it is alread present in the repo for the dummy example) and production at scale can be submitted with:
```
crab submit -c crab_config.py
```

### Creating NanoAOD files

In this part we summarize the procedure to produce nanoAOD files for a sample (i.e. with reconstruction included). The procedure can be find in different twikis, for convenience we base ours on [this](https://indico.cern.ch/event/1500035/contributions/6575125/attachments/3091743/5476084/IITH_GEN_Tutorial.pdf).
Note that, in order for the commands to work and have tags and CMSSW versions that are compatible with each other, the best idea is to take a NanoAOD sample similar to what we want to generate and copy the config commands.
This example is made using [this workflow](https://cms-pdmv-prod.web.cern.ch/mcm/chained_requests?prepid=HIG-chain_Run3Summer23BPixwmLHEGS_flowRun3Summer23BPixDRPremix_flowRun3Summer23BPixMiniAODv4_flowRun3Summer23BPixNanoAODv13-00001&page=0&shown=15), but could be different for other processes.

GENSIM step:
```
cmsDriver.py Configuration/GenProduction/python/pythia_fragment.py \
    --python_filename ppH_Htt_SMEFTsim_topU3l_quadratic_gensim_cfg.py \
    --eventcontent RAWSIM,LHE \
    --customise Configuration/DataProcessing/Utils.addMonitoring --datatier GEN-SIM,LHE \
    --fileout file:gensim.root --conditions 130X_mcRun3_2023_realistic_postBPix_v5 --beamspot Realistic25ns13p6TeVEarly2023Collision \
    --step LHE,GEN,SIM --geometry DB:Extended --era Run3_2023 --no_exec --mc -n -1
```
For the scaling to work, it is also convenient to make the number of events processed customizable via command line (this is taken from Davide's samples production for EFT ttHbb). We thus add these lines:
```
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('analysis')
options.register('jobNum', 0, VarParsing.multiplicity.singleton,VarParsing.varType.int,"jobNum")
options.register('nEvents', 0, VarParsing.multiplicity.singleton,VarParsing.varType.int,"nEvents")
options.parseArguments()
```
and change the following lines:
```
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.nEvents))
```
```
# Input source
process.source = cms.Source(
        "EmptySource",
        firstLuminosityBlock  = cms.untracked.uint32(options.jobNum+1),
        numberEventsInLuminosityBlock = cms.untracked.uint32(options.nEvents)
)
```
```nEvents = cms.untracked.uint32(options.nEvents)``` in ```process.externalLHEProducer```.

One can then start the production with:
```
cmsRun ppH_Htt_SMEFTsim_topU3l_quadratic_gensim_cfg.py
```

DIGI-RAW step (TODO resources for careful and meaningful inclusion of pileup and HLT have to be found):
```
cmsDriver.py \
    --python_filename ppH_Htt_SMEFTsim_topU3l_quadratic_digi_cfg.py \
    --eventcontent PREMIXRAW \
    --customise Configuration/DataProcessing/Utils.addMonitoring --datatier GEN-SIM-RAW \
    --fileout file:digi.root \
    --pileup_input "dbs:/Neutrino_E-10_gun/Run3Summer21PrePremix-Summer23BPix_130X_mcRun3_2023_realistic_postBPix_v1-v1/PREMIX" \
    --conditions 130X_mcRun3_2023_realistic_postBPix_v5 \
    --step DIGI,DATAMIX,L1,DIGI2RAW,HLT:2023v12 \
    --procModifiers premix_stage2,siPixelQualityRawToDigi --geometry DB:Extended --filein file:gensim.root \
    --datamix PreMix --era Run3_2023 --no_exec --mc -n -1
```
```
cmsRun ppH_Htt_SMEFTsim_topU3l_quadratic_digi_cfg.py
```

RECO/AOD step:
```
cmsDriver.py \
    --python_filename ppH_Htt_SMEFTsim_topU3l_quadratic_reco_cfg.py \
    --eventcontent AODSIM --customise Configuration/DataProcessing/Utils.addMonitoring --datatier AODSIM \
    --fileout file:aod.root --conditions 130X_mcRun3_2023_realistic_postBPix_v5 --step RAW2DIGI,L1Reco,RECO,RECOSIM \
    --geometry DB:Extended --filein file:digi.root --era Run3_2023 --no_exec --mc -n -1
```
```
cmsRun ppH_Htt_SMEFTsim_topU3l_quadratic_reco_cfg.py
```

MINIAOD step:
```
cmsDriver.py \
    --python_filename ppH_Htt_SMEFTsim_topU3l_quadratic_miniaod_cfg.py --eventcontent MINIAODSIM \
    --customise Configuration/DataProcessing/Utils.addMonitoring --datatier MINIAODSIM \
    --fileout file:miniaod.root --conditions 130X_mcRun3_2023_realistic_postBPix_v5 \
    --step PAT --geometry DB:Extended --filein file:aod.root --era Run3_2023 \
    --no_exec --mc -n -1
```
```
cmsRun ppH_Htt_SMEFTsim_topU3l_quadratic_miniaod_cfg.py
```

NANOAOD step:
```
cmsDriver.py  \
    --eventcontent NANOAODSIM --customise Configuration/DataProcessing/Utils.addMonitoring --datatier NANOAODSIM \
    --conditions 133X_mcRun3_2023_realistic_postBPix_ForNanov13_v2 --step NANO \
    --scenario pp --era Run3_2023,run3_nanoAOD_124 --python_filename ppH_Htt_SMEFTsim_topU3l_quadratic_nanoaod_cfg.py \
    --fileout file:nanoaod.root --filein file:miniaod.root --no_exec --mc -n -1
```
Also in this case, the different weights, taken from the ```reweight_card.dat```, can be saved by adding lines like
```
named_weights = [
    "rw0000",
    "rw0001",
]
process.genWeightsTable.namedWeightIDs = named_weights
process.genWeightsTable.namedWeightLabels = named_weights
```

Remember, for the purpose of having the dataset published on DBS, to add the following line to the config file:
```
fakeNameForCrab = cms.untracked.bool(True)
```
right underneath the line ```outputCommands = process.NANOAODSIMEventContent.outputCommands```.

Finally, run:
```
cmsRun ppH_Htt_SMEFTsim_topU3l_quadratic_nanoaod_cfg.py
```

#### Scale out (TODO)

Here we show how to use CRAB to produce NanoAOD with the relevant EFT weights. Relevant documentation can be found [here](https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3AdvancedTopic). As it can be seen in the link, the procedure consists in using a `scriptExe.sh` file where the different stages are chained. In order to control the number of events produced, extra lines to configure via command line are added to the config file of the first stage, while the others simpli take `-1` (i.e., hopefully, all the events coming from the previous steps).

Important notes:

- remember to move the gridpacks to ```/eos/user/g/gallim/www/EFTstudies/gridpacks``` and change the path to the first gensim step to ```root://eosuser.cern.ch//eos/user/g/gallim/www/EFTstudies/gridpacks/<gpackname>``` when submitting on crab
- for this to work, also change the line ```scriptName = cms.FileInPath('GeneratorInterface/LHEInterface/data/run_generic_tarball_cvmfs.sh')``` to ```scriptName = cms.FileInPath('GeneratorInterface/LHEInterface/data/run_generic_tarball_xrootd.sh')```

## Other useful links

- generators short exercise at CMS DAS ([link](https://fnallpc.github.io/generators/))
- presentation for EFT samples production at LPC school ([link](https://indico.fnal.gov/event/68174/timetable/?view=standard#26-mc-generation-tutorial))
- Ilaria Brivio's SMEFTsim tutorial ([link](https://github.com/IlariaBrivio/SMEFTsim_tutorial_June23))
