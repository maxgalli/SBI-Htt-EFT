#!/bin/bash

if [ ! -d genproductions_scripts ]; then
  git clone https://gitlab.cern.ch/cms-gen/genproductions_scripts.git
fi

# in case this is not already done, setup cms packaging commands
. /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=el9_amd64_gcc12
if [ ! -d CMSSW_13_2_9 ]; then
  cmsrel CMSSW_13_2_9
  pushd CMSSW_13_2_9/src
  export SCRAM_ARCH=el9_amd64_gcc12
  cmsenv
  git cms-init
  git cms-addpkg PhysicsTools/NanoAOD
  mkdir -p Configuration/GenProduction/python/
  cp ../../fragments/pythia_fragment.py Configuration/GenProduction/python/

  scram b -j 4
  popd
fi

# ensure environment if running again
pushd CMSSW_13_2_9/src && cmsenv && popd
