from CRABClient.UserUtilities import config
config = config()

config.General.requestName     = 'VBF_Htt_SMEFTsim_SM_MS_50kevents'
config.General.workArea        = 'crab_jobs'
config.General.transferOutputs = True
config.General.transferLogs    = True

config.JobType.pluginName = 'PrivateMC'
config.JobType.psetName = 'VBF_Htt_SMEFTsim_SM_gennanogen_cfg.py'
config.JobType.allowUndistributedCMSSW = True
config.JobType.maxMemoryMB = 3000
config.JobType.numCores = 1
config.JobType.inputFiles = [
    '/eos/user/z/zhaom/htautau/samples/SBI-Htt-EFT/genproductions_mg35x_gh/bin/MadGraph5_aMCatNLO/VBF_Htt_SMEFTsim_SM_MS_el8_amd64_gcc10_CMSSW_12_4_8_tarball.tar.xz'
]

config.Data.splitting = 'EventBased'
config.Data.unitsPerJob = 1000
NJOBS = 50 # This is not a configuration parameter, but an auxiliary variable that we use in the next line.
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
# config.Data.outLFNDirBase = '/store/user/gallim/SBI-Htt-EFT/nanogen'
config.Data.publication = False
config.Data.outputPrimaryDataset = 'qqHtoTauTau_SM'
config.Data.outputDatasetTag     = '130X_mcRun3_2023_realistic_postBPix_v5_SM_MS'

config.Site.storageSite = 'T3_CH_CERNBOX'
