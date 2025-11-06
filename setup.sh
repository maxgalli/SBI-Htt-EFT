#!/bin/bash

# in case this is not already done, setup cms packaging commands
. /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=el8_amd64_gcc11
if [ ! -d CMSSW_13_0_14 ]; then
  cmsrel CMSSW_13_0_14
  pushd CMSSW_13_0_14/src
  export SCRAM_ARCH=el8_amd64_gcc11
  cmsenv
  git cms-init
  git cms-addpkg PhysicsTools/NanoAOD
  mkdir -p Configuration/GenProduction/python/
  cp ../../fragments/pythia_fragment.py Configuration/GenProduction/python/

  scram b -j 4
  popd
fi

# ensure environment if running again
pushd CMSSW_13_0_14/src && cmsenv && popd
