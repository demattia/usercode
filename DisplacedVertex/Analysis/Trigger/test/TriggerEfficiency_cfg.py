import FWCore.ParameterSet.Config as cms

process = cms.Process("TriggerEfficiency")

# # initialize MessageLogger and output report
# process.load("FWCore.MessageLogger.MessageLogger_cfi")
# # process.MessageLogger.cerr.threshold = 'INFO'
# # process.MessageLogger.categories.append('Demo')
# # process.MessageLogger.cerr.INFO = cms.untracked.PSet(
# #     limit = cms.untracked.int32(-1)
# # )

process.load("FWCore.MessageService.test.Services_cff")
# Here is the configuration of the MessgeLogger Service:
process.MessageLogger = cms.Service("MessageLogger",
    destinations = cms.untracked.vstring('Message'),
    Message = cms.untracked.PSet(
        threshold = cms.untracked.string('INFO')
    )
)

process.load("Configuration.StandardSequences.Services_cff")
process.load('Configuration.StandardSequences.Reconstruction_cff')

process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")

# Careful, this needs to be changed for the data
process.GlobalTag.globaltag = 'START52_V2::All'

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )



readFiles = cms.untracked.vstring()
readFiles.extend([
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_10.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_100.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_101.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_102.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_103.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_104.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_105.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_106.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_107.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_108.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_109.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_11.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_110.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_111.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_112.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_113.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_114.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_115.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_116.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_117.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_118.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_119.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_12.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_120.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_121.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_122.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_123.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_124.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_125.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_126.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_127.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_128.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_129.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_13.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_130.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_131.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_132.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_133.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_134.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_135.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_136.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_137.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_138.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_139.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_14.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_140.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_141.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_142.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_143.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_144.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_145.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_146.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_147.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_148.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_149.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_15.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_150.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_151.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_152.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_153.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_154.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_155.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_156.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_157.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_158.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_159.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_16.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_160.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_161.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_162.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_163.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_164.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_165.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_166.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_167.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_168.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_169.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_17.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_170.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_171.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_172.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_173.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_174.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_175.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_176.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_177.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_178.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_179.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_18.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_180.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_181.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_182.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_183.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_184.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_185.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_186.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_187.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_188.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_189.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_19.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_190.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_191.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_192.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_193.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_194.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_195.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_196.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_197.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_198.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_199.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_2.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_20.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_200.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_21.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_22.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_23.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_24.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_25.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_26.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_27.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_28.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_29.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_3.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_30.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_31.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_32.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_33.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_34.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_35.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_36.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_37.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_38.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_39.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_4.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_40.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_41.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_42.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_43.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_44.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_45.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_46.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_47.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_48.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_49.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_5.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_50.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_51.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_52.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_53.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_54.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_55.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_56.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_57.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_58.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_59.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_6.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_60.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_61.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_62.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_63.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_64.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_65.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_66.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_67.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_68.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_69.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_7.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_70.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_71.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_72.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_73.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_74.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_75.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_76.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_77.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_78.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_79.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_8.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_80.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_81.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_82.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_83.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_84.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_85.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_86.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_87.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_88.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_89.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_9.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_90.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_91.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_92.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_93.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_94.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_95.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_96.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_97.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_98.root",
"castor:/castor/cern.ch/user/d/demattia/DisplacedVertex/MC/MuonGun/UniformD0/5_2_0_pre6/HLT/Refit/HLTOutput_new_99.root"]);

process.source = cms.Source('PoolSource', fileNames = readFiles)

# process.source = cms.Source("PoolSource",
#     fileNames = cms.untracked.vstring(
#     "file:/afs/cern.ch/user/d/demattia/scratch0/DisplacedVertex/HLTOutput_new.root"
#     )
# )

process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

# For the propagator
process.MaterialPropagator = cms.ESProducer("PropagatorWithMaterialESProducer",
                                            MaxDPhi = cms.double(1.6),
                                            ComponentName = cms.string('PropagatorWithMaterial'),
                                            Mass = cms.double(0.105),
                                            PropagationDirection = cms.string('alongMomentum'),
                                            useRungeKutta = cms.bool(False),
                                            # If ptMin > 0, uncertainty in reconstructed momentum will be taken into account when estimating rms scattering angle.
                                            # (By default, it is neglected). However, it will also be assumed that the track pt can't be below specified value,
                                            # to prevent this scattering angle becoming too big.
                                            ptMin = cms.double(-1.)
                                            )
process.SteppingHelixPropagatorAny = cms.ESProducer("SteppingHelixPropagatorESProducer",
                                                    ComponentName = cms.string('SteppingHelixPropagatorAny'),
                                                    NoErrorPropagation = cms.bool(False),
                                                    PropagationDirection = cms.string('anyDirection'),
                                                    useTuningForL2Speed = cms.bool(False),
                                                    useIsYokeFlag = cms.bool(True),
                                                    endcapShiftInZNeg = cms.double(0.0),
                                                    SetVBFPointer = cms.bool(False),
                                                    AssumeNoMaterial = cms.bool(False),
                                                    endcapShiftInZPos = cms.double(0.0),
                                                    useInTeslaFromMagField = cms.bool(False),
                                                    VBFName = cms.string('VolumeBasedMagneticField'),
                                                    useEndcapShiftsInZ = cms.bool(False),
                                                    sendLogWarning = cms.bool(False),
                                                    useMatVolumes = cms.bool(True),
                                                    debug = cms.bool(False),
                                                    #This sort of works but assumes a measurement at propagation origin
                                                    ApplyRadX0Correction = cms.bool(True),
                                                    useMagVolumes = cms.bool(True),
                                                    returnTangentPlane = cms.bool(True)
                                                    )

# Smart propagator with IP
process.smartPropagatorWithIPESProducer = cms.ESProducer("SmartPropagatorWithIPESProducer",
                                                         ComponentName = cms.string('SmartPropagatorWithIP'),
                                                         TrackerPropagator = cms.string('PropagatorWithMaterial'),
                                                         # TrackerPropagator = cms.string('RungeKuttaTrackerPropagator'),
                                                         MuonPropagator = cms.string('SteppingHelixPropagatorAny'),
                                                         PropagationDirection = cms.string('alongMomentum'),
                                                         Epsilon = cms.double(10.0) # the standard one has 5., but uses 10 hardcoded internally...
                                                         )



process.demo = cms.EDAnalyzer('TriggerEfficiency',
                              UseMCtruth = cms.bool(True),
                              TrackCollection = cms.InputTag("hltL2Muons"),
                              )

process.p = cms.Path(process.demo)
