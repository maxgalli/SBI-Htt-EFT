from CRABClient.UserUtilities import config
config = config()

config.General.requestName     = 'VBF_Htt_SMEFTsim_topU3l_quadratic_gennanogen'
config.General.workArea        = 'crab_jobs'
config.General.transferOutputs = True
config.General.transferLogs    = False

config.JobType.pluginName = 'PrivateMC'
config.JobType.psetName = 'VBF_Htt_SMEFTsim_topU3l_quadratic_gennanogen_cfg.py'
config.JobType.allowUndistributedCMSSW = True
config.JobType.maxMemoryMB = 4500

config.JobType.numCores = 1
#process.options.numberOfThreads=2
config.Data.splitting = 'EventBased'
config.Data.unitsPerJob = 2500
NJOBS = 2000  # This is not a configuration parameter, but an auxiliary variable that we use in the next line.
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.outLFNDirBase = '/store/user/gallim/SBI-Htt-EFT/nanogen'
config.Data.publication = False
config.Data.outputPrimaryDataset = 'qqHtoTauTau'
config.Data.outputDatasetTag     = '130X_mcRun3_2023_realistic_postBPix_v5'

config.Site.storageSite = 'T3_CH_PSI'
