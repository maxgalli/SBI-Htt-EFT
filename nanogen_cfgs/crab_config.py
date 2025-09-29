from CRABClient.UserUtilities import config
config = config()

config.General.requestName     = 'ppH_Htt_SMEFTsim_topU3l_quadratic_nanogen_2'
config.General.workArea        = 'crab_jobs'
config.General.transferOutputs = True
config.General.transferLogs    = False

config.JobType.pluginName = 'PrivateMC'
config.JobType.psetName = 'ppH_Htt_SMEFTsim_topU3l_quadratic_gennanogen_cfg.py'
config.JobType.allowUndistributedCMSSW = True

config.JobType.numCores = 1
config.Data.splitting = 'EventBased'
config.Data.unitsPerJob = 5000
NJOBS = 2  # This is not a configuration parameter, but an auxiliary variable that we use in the next line.
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.outLFNDirBase = '/store/user/gallim/SBI-Htt-EFT/nanogen_test2'
config.Data.publication = False
config.Data.outputPrimaryDataset = 'PrimaryDatasetTest2'
config.Data.outputDatasetTag     = 'DatasetTagTest2'

config.Site.storageSite = 'T3_CH_PSI'
