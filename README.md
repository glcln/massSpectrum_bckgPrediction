# massSpectrum_bckgPrediction

This repositery is dedicated to HSCP analysis, for the background estimate of the mass spectrum. 

## Setup working area

```bash
export SCRAM_ARCH=slc7_amd64_gcc700
cmsrel CMSSW_10_6_30
cd CMSSW_10_6_30/src/
cmsenv
```

For the following step you should have a ssh key associated to your GitHub account.
For more information, see [connecting-to-github-with-ssh-key](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent).

```bash
git clone -b master git@github.com:dapparu/massSpectrum_bckgPrediction.git massSpectrum_bckgPrediction 
``` 

## Run the background estimate code 

This part concerns the run of the background estimate method. 

If you want to run separately each configuration (for each systematic), you need to uncomment the chosen line and run :

```bash
python LaunchBckgOn1Syst.py
```

You can change on what you want to run directly in ```LaunchBckgOn1Syst.py``` and ```step2_backgroundPrediction.C``` files. 

In the first file (```LaunchBckgOn1Syst.py```): 
- You can set on which datasets you want to run, giving the path to the root file with all the needed histograms, in the ```datasetList``` array. 
- The ```config``` array is used to indicate which kind of estimates are ran: nominal or the different systematics. 
- The ```nPe```variable set the number of pseudo-experiments done during the background estimate. 
- The ```odir``` array gives the directory where you can find fast produced plots. 

The code runs on 25 cores in parallel (in local) and it can be changed in the file ```Regions.h```.

In the second file (```step2_backgroundPrediction.C```), you can set which regions you want for estimation.

After you ran the code, a new file is created with the histograms corresponding to the background estimate in the wanted regions, and labels are set to correspond to the different systematics cases. 

*If you want to run everything automatically at once (will take several hours, use a screen during the night), use :*
```bash
python LaunchBckgOnAllSyst.py
```


## Plotting part

This part concerns the run of mass spectrum plotter code, with nice style. 

The code runs with the code : 
```bash
python ShowPreds.py
```

which is calling :
```bash
python2.7 MyMacroMass.py --ifile inputfile --ofile mass_plot --region reg --odir outputdir
```

Concerniing the options:
- option ```--ifile``` is for the input file; obtained at the previous step. 
- option ```--ofile``` is for the output file label. 
- option ```--region``` is for the region on which one wants to run. 
- option ```--odir``` is for the output directory. 