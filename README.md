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

## To run

With the current configuration the config file is set to run the following fragments:

J/psi + D0; J/psi + D*; J/psi + Ds+; J/psi + Lambdac+ as well as Psi(2S) + same open charms.

In order to add more options go to test/Config.py and add another **elif** with the path of the files.

**example:**

```
elif options.channel == 'my_new_channel':
   path = 'my_new_channel/my_new_math.txt'
   print '-------------- Particle List for channel fragment --------------'  
```

Finally run it:

cmsRun test/Config.py channel=my_new_channel
