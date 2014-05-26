# Summary
This package allows for the collection of layer information (hit positions, trajectory positions, residuals, etc.)
TODO for CSC right now, but will add DT

# Instructions
Set up working area:
``` bash
mkdir layer_info
cd layer_info
git clone http://github.com/cms-tamu/MuAl_ME13.git
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
```


Soft links for the `*.cc` files are for convenience, as they are the only files (along with their corresponding header files) that I edited for layer information collection. `CSCTTree.h` was created to contain the layer information structure.

`createJobs.py` and `gather_cfg.py` were both modified to allow for easy switching between data and MC alignment using a `--isMC` flag in the `createJobs.py` line (example provided later on).


A sample filelist (one file) generated using the singleMuonGun is provided in ` Cert_singleMuonGun_MC_MEp_1_3_17_ABRIDGED_V2.py`.
`idealGeometry.db` and `inertGlobalPositionRcd.db` are used for MC.


Now run using an example line below:
``` bash
./createJobs.py MuonFilter_2012_MEp_1_3_17_MC_TEST_ 1 idealGeometry.db Cert_singleMuonGun_MC_MEp_1_3_17_ABRIDGED_V2.py -s MuonFilter_2012_MEp_1_3_17_MC_TEST.sh --validationLabel MuonFilter_2012_MEp_1_3_17_MC_TEST --user_mail 157866test@gmail.com --minTrackPt 30 --maxTrackPt 200 --maxDxy 0.2 --minNCrossedChambers 1 --residualsModel pureGaussian --peakNSigma 2. --station123params 000010 --station4params 000010 --cscparams 100001 --useResiduals 1100 --noDT --mapplots --curvatureplots --segdiffplots --extraPlots --globalTag MC_53_V14::All --gprcd inertGlobalPositionRcd --gprcdconnect sqlite_file:inertGlobalPositionRcd.db  --createAlignNtuple -j 1 --isMC
```

This has 10k events and should take no more than 12 minutes to run.
