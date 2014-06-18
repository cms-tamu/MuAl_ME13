# Summary
This package allows for the collection of layer information (hit positions, trajectory positions, residuals, etc.)
TODO for CSC right now, but will add DT

# Instructions
Set up working area:
``` bash
mkdir layer_info
cd layer_info
git clone https://github.com/cms-tamu/MuAl_ME13.git
cmsrel CMSSW_5_3_6_patch1
cd CMSSW_5_3_6_patch1/src
cp -rp ../../MuAl_ME13/alignment_layer_ntuple/* .
cmsenv
scram b -j8
```

Now make some soft links:
``` bash
ln -s Alignment/MuonAlignmentAlgorithms/scripts/createJobs.py
ln -s Alignment/MuonAlignmentAlgorithms/interface/CSCTTree.h
ln -s Alignment/MuonAlignmentAlgorithms/python/gather_cfg.py
ln -s Alignment/MuonAlignmentAlgorithms/plugins/MuonAlignmentFromReference.cc
ln -s Alignment/MuonAlignmentAlgorithms/python/align_cfg.py
ln -s Alignment/MuonAlignmentAlgorithms/src/MuonResidualsFromTrack.cc
ln -s Alignment/MuonAlignmentAlgorithms/scripts/layerPlots.py
ln -s Alignment/MuonAlignmentAlgorithms/python/MuonAlignmentFromReference_cfi.py
```

Soft links for the `*.cc` files are for convenience, as they are the files (along with their corresponding header files) that I edited for layer information collection. `CSCTTree.h` was created to contain the layer information structure.

`createJobs.py` and `gather_cfg.py` were both modified to allow for easy switching between data and MC alignment using a `--isMC` flag in the `createJobs.py` line (example provided later on). Note: at minimum, `MuonAlignmentFromReference_cfi.py`, `createJobs.py`, `gather_cfg.py`, `MuonAlignmentFromReference.cc` need to be modified when adding another createJobs syntax parameter that gets used within a gather job.

A sample filelist (one file) generated using the singleMuonGun is provided in ` Cert_singleMuonGun_MC_MEp_1_3_17_ABRIDGED_V2.py`.
`idealGeometry.db` and `inertGlobalPositionRcd.db` are used for MC.


Now run using this example line:
``` bash
./createJobs.py MuonFilter_2012_MEp_1_3_17_MC_TEST_ 1 idealGeometry.db Cert_singleMuonGun_MC_MEp_1_3_17_ABRIDGED_V2.py \
-s MuonFilter_2012_MEp_1_3_17_MC_TEST.sh --validationLabel MuonFilter_2012_MEp_1_3_17_MC_TEST \
--user_mail 157866test@gmail.com --minTrackPt 30 --maxTrackPt 200 --maxDxy 0.2 --minNCrossedChambers 1 \
--residualsModel pureGaussian --peakNSigma 2. --station123params 000010 --station4params 000010 --cscparams 100001 \
--useResiduals 1100 --noDT --mapplots --curvatureplots --segdiffplots --extraPlots --globalTag MC_53_V14::All \
--gprcd inertGlobalPositionRcd --gprcdconnect sqlite_file:inertGlobalPositionRcd.db  --createAlignNtuple -j 1 \
--isMC --createLayerNtuple --layerPlots
```
And, finally, submit via
``` bash
. MuonFilter_2012_MEp_1_3_17_MC_TEST.sh
```

This will run over 20k events and should take no more than 15 minutes to run.

Note that the `--isMC` flag was added to accomodate MC inputs. `--createLayerNtuple` collects layer information and places it in the *_plotting.root file inside of the alignment job's directory. `--layerPlots` will take this layer information and produce plots in a tree-like directory structure. Each folder will contain a PHP file allowing for easy viewing of the plots in the directories if they are copied over to a website. If a chamber receives less than 200 tracks, it will not be considered for plotting.

Also note that the `--b` tag is used to tell `createJobs.py` that jobs will be *b*ig, thus forcing them to use the daily queue instead of the hourly queue by default.

In this example, `MuonFilter_2012_MEp_1_3_17_MC_TEST_01/MuonFilter_2012_MEp_1_3_17_MC_TEST_01_plotting.root` will contain a TTree `csc_layer_ttree` with layer information. The .tgz file in the `CMSSW_5_3_6_patch1/src` directory that normally contains the plots will be accompanied by another .tgz file with the layer plots, for example (timestamp will be different),  `MuonFilter_2012_MEp_1_3_17_MC_TEST_20140528235539_layer_plots.tgz`. Doing
``` bash
tar xzf MuonFilter_2012_MEp_1_3_17_MC_TEST_20140528235539_layer_plots.tgz
```
produces a folder called `layer_plots/` which can then be copied to a web area and browsed using the directory and index PHP pages.

# Cuts
The optional `--cutTypes` flag allows for one of six cuts to be implemented during the gather jobs. For example, `--cutTypes SIDEBOX` will implement side cuts and box cuts. Cuts are described in detail in https://twiki.cern.ch/twiki/bin/viewauth/CMS/MuAlME13CutComparisons. The tags associated with the various cuts are:
* SIDE
* SIDEBOX
* SIDEBOXFID
* BOX
* BOXFID
* SIDESYMBOXFID

# COMPLETE COPY&PASTE TEST
The below commands can be completely copied and pasted as a whole to obtain the code, copy it, compile, create and submit a job with a cut ("BOX" in this case).
```bash
mkdir layer_info
cd layer_info
git clone https://github.com/cms-tamu/MuAl_ME13.git
cmsrel CMSSW_5_3_6_patch1
cd CMSSW_5_3_6_patch1/src
cp -rp ../../MuAl_ME13/alignment_layer_ntuple/* .
cmsenv
scram b -j8
ln -s Alignment/MuonAlignmentAlgorithms/scripts/createJobs.py
ln -s Alignment/MuonAlignmentAlgorithms/interface/CSCTTree.h
ln -s Alignment/MuonAlignmentAlgorithms/python/gather_cfg.py
ln -s Alignment/MuonAlignmentAlgorithms/plugins/MuonAlignmentFromReference.cc
ln -s Alignment/MuonAlignmentAlgorithms/python/align_cfg.py
ln -s Alignment/MuonAlignmentAlgorithms/src/MuonResidualsFromTrack.cc
ln -s Alignment/MuonAlignmentAlgorithms/scripts/layerPlots.py
ln -s Alignment/MuonAlignmentAlgorithms/python/MuonAlignmentFromReference_cfi.py
./createJobs.py MuonFilter_2012_MEp_1_3_17_MC_TEST_ 1 idealGeometry.db Cert_singleMuonGun_MC_MEp_1_3_17_ABRIDGED_V2.py \
-s MuonFilter_2012_MEp_1_3_17_MC_TEST.sh --validationLabel MuonFilter_2012_MEp_1_3_17_MC_TEST \
--user_mail 157866test@gmail.com --minTrackPt 30 --maxTrackPt 200 --maxDxy 0.2 --minNCrossedChambers 1 \
--residualsModel pureGaussian --peakNSigma 2. --station123params 000010 --station4params 000010 --cscparams 100001 \
--useResiduals 1100 --noDT --mapplots --curvatureplots --segdiffplots --extraPlots --globalTag MC_53_V14::All \
--gprcd inertGlobalPositionRcd --gprcdconnect sqlite_file:inertGlobalPositionRcd.db  --createAlignNtuple -j 2 --cutTypes BOX \
--isMC --createLayerNtuple --layerPlots
. MuonFilter_2012_MEp_1_3_17_MC_TEST.sh
bjobs -w
```
Afterwards, you may extract the output `*_layer_plots.tgz` file and copy it to a web directory. Navigate to `layer_plots/summary.php` and click on a few plots to ensure everything is working correctly. The 2D occupancy plots should confirm that the cut was implemented correctly.
