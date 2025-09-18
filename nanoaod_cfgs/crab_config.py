from CRABClient.UserUtilities import config
config = config()

config.General.requestName     = 'ppH_Htt_SMEFTsim_topU3l_quadratic'
config.General.workArea        = 'crab_jobs'
config.General.transferOutputs = True
config.General.transferLogs    = False

config.JobType.pluginName  = 'PrivateMC'
# we use the 'fake' one to avoid the following error
# Found a 'PoolSource' in the parameter set. But that's not compatible with PrivateMC. Please switch to either 'EmptySource' or 'LHESource' for event generation, or set 'JobType.pluginName' to 'Analysis'.
# seemed to have worked for Davide https://github.com/valsdav/ttHbb-EFT/blob/master/nanoaod_cfgs/RunIISummer20UL18/ttHbb_p1j_EFTcenter_5F_tbarlnutqq_madspin/crabConfig.py
config.JobType.psetName    = 'ppH_Htt_SMEFTsim_topU3l_quadratic_nanoaod_cfg_fake.py'
config.JobType.inputFiles  = ['scriptExe.sh', 'ppH_Htt_SMEFTsim_topU3l_quadratic_gensim_cfg.py', 'ppH_Htt_SMEFTsim_topU3l_quadratic_digi_cfg.py', 'ppH_Htt_SMEFTsim_topU3l_quadratic_reco_cfg.py', 'ppH_Htt_SMEFTsim_topU3l_quadratic_miniaod_cfg.py', 'ppH_Htt_SMEFTsim_topU3l_quadratic_nanoaod_cfg.py']
config.JobType.scriptExe   ='scriptExe.sh'
config.JobType.allowUndistributedCMSSW = True

config.JobType.numCores    = 1
# as only 2500 are guaranteed to be available
config.JobType.maxMemoryMB = 2500

config.Data.splitting   = 'EventBased'
config.Data.unitsPerJob = 100
NJOBS = 2
config.Data.totalUnits  = config.Data.unitsPerJob * NJOBS
config.Data.outLFNDirBase = '/store/user/gallim/SBI-Htt-EFT/test'
config.Data.publication = True
config.Data.outputPrimaryDataset = 'PrimaryDatasetTest'
config.Data.outputDatasetTag     = 'DatasetTagTest'

config.Site.storageSite = 'T3_CH_PSI'
