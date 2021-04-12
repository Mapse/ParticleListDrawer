import FWCore.ParameterSet.Config as cms
# Include varParsing
from FWCore.ParameterSet.VarParsing import VarParsing

process = cms.Process("GENLevelDump")
import FWCore.Utilities.FileUtils as FileUtils

options = VarParsing('analysis')

options.register ('channel',
				  '',
				  VarParsing.multiplicity.singleton,
				  VarParsing.varType.string,
				  "Channel Chosen")

options.parseArguments()

## 2017 geometry 
#from Configuration.Eras.Era_Run2_2018_cff import Run2_2018

## Phase 2
from Configuration.Eras.Era_Phase2C8_timing_layer_bar_cff import Phase2C8_timing_layer_bar
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')

## Global tag for 10.6 phase2 mc
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '102X_mc2017_realistic_v8', '')

# Messagelogger
""" process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('Demo')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit = cms.untracked.int32(-1)
) """

#process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# Choose the right channel to draw the particle list.

if options.channel == 'Jpsi_D0':
  path = 'Jpsi_OpenCharm/files_path_D0.txt'
  print '-------------- Particle List for Jpsi + D0 fragment --------------'
  
elif options.channel == 'Jpsi_D*':
  path = 'Jpsi_OpenCharm/files_path_D*.txt'
  print '-------------- Particle List for Jpsi + D* fragment --------------'

elif options.channel == 'Jpsi_D+':
  path = 'Jpsi_OpenCharm/files_path_D+.txt'
  print '-------------- Particle List for Jpsi + D+ fragment --------------'

elif options.channel == 'Jpsi_Ds+':
  path = 'Jpsi_OpenCharm/files_path_Ds+.txt'
  print '-------------- Particle List for Jpsi + Ds+ fragment --------------'

elif options.channel == 'Jpsi_Lambdac+':
  path = 'Jpsi_OpenCharm/files_path_Lambdac+.txt'
  print '-------------- Particle List for Jpsi + Lambdac+ fragment --------------'

elif options.channel == 'Psi_D0':
  path = 'Psi_OpenCharm/files_path_D0.txt'
  print '-------------- Particle List for Psi + D0 fragment --------------'
  
elif options.channel == 'Psi_D*':
  path = 'Psi_OpenCharm/files_path_D*.txt'
  print '-------------- Particle List for Psi + D* fragment --------------'

elif options.channel == 'Psi_D+':
  path = 'Psi_OpenCharm/files_path_D+.txt'
  print '-------------- Particle List for Psi + D+ fragment --------------'

elif options.channel == 'Psi_Ds+':
  path = 'Psi_OpenCharm/files_path_Ds+.txt'
  print '-------------- Particle List for Psi + Ds+ fragment --------------'

elif options.channel == 'Psi_Lambdac+':
  path = 'Psi_OpenCharm/files_path_Lambdac+.txt'
  print '-------------- Particle List for Psi + Lambdac+ fragment --------------'


# File source
mylist = FileUtils.loadListFromFile(path) 
process.source = cms.Source(
    "PoolSource",
    #fileNames  = cms.untracked.vstring(options.inputFiles),
    fileNames  = cms.untracked.vstring(*mylist),
    duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
)


## Show GenParticles
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.printTree = cms.EDAnalyzer("ParticleListDrawer",
  maxEventsToPrint = cms.untracked.int32(10),
  printVertex = cms.untracked.bool(True),
  printOnlyHardInteraction = cms.untracked.bool(False), # Print only status=3 particles. This will not work for Pythia8, which does not have any such particles.
  src = cms.InputTag("genParticles"), #genParticles
     
)


process.genlevel = cms.EDAnalyzer('GenLevelStudies',
 src = cms.InputTag("genParticles")
)

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '106X_upgrade2023_realistic_v2', '')

process.TFileService = cms.Service("TFileService",
	fileName = cms.string("mc_DZero.root"),
	closeFileFast = cms.untracked.bool(True))

# Path and endpath to run the producer and output modules
process.p =  cms.Path(process.printTree) #process.genlevel+

