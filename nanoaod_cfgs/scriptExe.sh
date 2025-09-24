#!/bin/bash
NEVENTS=7000

echo "================= CMSRUN starting jobNum= ====================" | tee -a job.log

export SCRAM_ARCH=el8_amd64_gcc11
source /cvmfs/cms.cern.ch/cmsset_default.sh

echo "================= CMSRUN starting GEN-SIM step ====================" | tee -a job.log
cmsRun -j ppH_Htt_SMEFTsim_topU3l_quadratic_gensim_cfg.log ppH_Htt_SMEFTsim_topU3l_quadratic_gensim_cfg.py jobNum=$1 nEvents=${NEVENTS}

echo "================= CMSRUN starting DIGI-RAW step ====================" | tee -a job.log

cmsRun -j ppH_Htt_SMEFTsim_topU3l_quadratic_digi_cfg.log ppH_Htt_SMEFTsim_topU3l_quadratic_digi_cfg.py

echo "================= CMSRUN starting RECO-AOD step ====================" | tee -a job.log

cmsRun -j ppH_Htt_SMEFTsim_topU3l_quadratic_reco_cfg.log ppH_Htt_SMEFTsim_topU3l_quadratic_reco_cfg.py

echo "================= CMSRUN starting MiniAOD step  ====================" | tee -a job.log

cmsRun -j ppH_Htt_SMEFTsim_topU3l_quadratic_miniaod_cfg.log ppH_Htt_SMEFTsim_topU3l_quadratic_miniaod_cfg.py

echo "================= CMSRUN starting NanoAOD step  ====================" | tee -a job.log

cmsRun -e -j FrameworkJobReport.xml ppH_Htt_SMEFTsim_topU3l_quadratic_nanoaod_cfg.py

echo "================= CMSRUN finished ====================" | tee -a job.log
