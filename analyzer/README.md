Muon Tracks Analyzer
----

##### Pre-requisites
CMSSW_5_3_11_patch6, SLC5

##### Included example configuration files, scripts

* formattedfilespyA.py
* hitanalyzerrefit_A1_cfg.py
* hitanalyzerrefit_MC1_cfg.py
* runA1.sh

`formattedfilespyA.py` contains an example filelist of length ~300 .root files stored in the `/eos/store/caf/user/` area. 

You can use CRAB to submit jobs, but as I found this to be an unresolved issue, I opted to make scripts to run this via local batch submission to LSF.

Submitting one job is done with:
```bsub -R "type=SLC5_64" -q 1nd -J jobrunA1 -u 12098asdftest@gmail.com < runA1.sh```
Note that you will need to edit the `CMSSW_PROJECT_SRC` variable within `runA1.sh`, as well as the python cfg it references.

The python cfg `hitanalyzerrefit_A1_cfg.py` is configured for data, and `hitanalyzerrefit_MC1_cfg.py` is configured for Monte Carlo. As an example, let's look at the former.

`hitanalyzerrefit_A1_cfg.py` grabs a chunk of files from the filelist. Specify the total number of files in one filelist in line 12 (in this case, 80):
```myfilelist.extend( chunks(formattedfilespyA.fileListA, 80)[part - 1] )```
and specify the chunk number in line 6 (in this case, 1):
```part=1```

Thus, you can copy this cfg file and have different versions with different chunk numbers (change one value on line 6). This allows you to break a large filelist down into smaller chunks for better parallelization when using LSF batch submission. My convention is that `A1` specifies that this will process "part 1" (chunk 1) of Run 2012A.



