#!/bin/sh

export CMSSW_PROJECT_SRC="/afs/cern.ch/work/n/namin/public/ring/test/new/CMSSW_5_3_11_patch6/src"

export CFG_FILE="hitanalyzerrefit_A1_cfg.py"

cd $CMSSW_PROJECT_SRC
eval `scramv1 run -sh`
cmsRun $CMSSW_PROJECT_SRC/$CFG_FILE
