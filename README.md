# ParticleListDrawer
This repository contains all information necessary to run the particle list drawer for a generator level Monte Carlo.


## Setting the environment:

export SCRAM_ARCH=slc7_amd64_gcc700

cmsrel CMSSW_10_6_4

cd CMSSW_10_6_4/src

cmsenv

mkdir PartDrawer
cd PartDrawer

git clone git@github.com:mapse/ParticleListDrawer .

scram b
