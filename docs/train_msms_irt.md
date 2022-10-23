# DeepDIA Tutorial: Training New Models for MS/MS and iRT Prediction
Training new models for MS/MS and iRT prediction using data-dependent acquisition (DDA) data.

## 1. Prerequisites
See [**DeepDIA Tutorial: Spectral Library Generation From Peptide Lists**](predict_msms_irt.md).

NVIDIA graphics cards with CUDA are recommended.
- [CUDA Toolkit](https://developer.nvidia.com/cuda-toolkit-archive) (version 11.2)
- [cuDNN SDK](https://developer.nvidia.com/cudnn) (version 8.1.0)

The following software is not needed in the process of model training, but is used in the DDA analysis to get training data:
- [SpectroMine](https://biognosys.com/software/spectromine) (Biognosys AG, Schlieren, Switzerland)


## 2. Tutorial Data
### 2.1. Starting Materials for Model Training
Starting materials of this tutorial are available at ProteomeXchange and iProX with identifier [`PXD014108`](http://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD014108) or [`IPX0001628000`](https://www.iprox.cn/page/project.html?id=IPX0001628000). 

Fragment reports exported from SpectroMine.
- HeLa_DDA.csv.zip

### 2.2. Tutorial Data for DDA Analysis
These files are not needed in the process of model training, but are used in the DDA analysis to get training data.

LC-MS/MS DDA data of HeLa cells on Q Exactive HF are available at ProteomeXchange with the data set identifier [`PXD005573`](http://proteomecentral.proteomexchange.org/GetDataset?ID=PXD005573) (Bruderer, R. et al. Mol. Cell. Proteomics 2017, 16, 2296-2309, [doi:10.1074/mcp.RA117.000314](https://doi.org/10.1074/mcp.RA117.000314)).
- C_D160304_S251-Hela-2ug-2h_MSG_R01_T0.raw
- C_D160304_S251-Hela-2ug-2h_MSG_R02_T0.raw
- C_D160304_S251-Hela-2ug-2h_MSG_R03_T0.raw
- C_D160331_S209-HPRP-HeLa-05_MSG_R01_T0.raw
- C_D160331_S209-HPRP-HeLa-10_MSG_R01_T0.raw
- C_D160331_S209-HPRP-HeLa-15_MSG_R01_T0.raw
- C_D160331_S209-HPRP-HeLa-20_MSG_R01_T0.raw
- C_D160331_S209-HPRP-HeLa-25_MSG_R01_T0.raw
- C_D160331_S209-HPRP-HeLa-50_MSG_R01_T0.raw
- C_D160331_S209-HPRP-HeLa-FT_MSG_R01_T0.raw
- C_D160401_S209-HPRP-HeLa-05_MSG_R01_T0.raw
- C_D160401_S209-HPRP-HeLa-10_MSG_R01_T0.raw
- C_D160401_S209-HPRP-HeLa-15_MSG_R01_T0.raw
- C_D160401_S209-HPRP-HeLa-20_MSG_R01_T0.raw
- C_D160401_S209-HPRP-HeLa-25_MSG_R01_T0.raw
- C_D160401_S209-HPRP-HeLa-50_MSG_R01_T0.raw
- C_D160401_S209-HPRP-HeLa-FT_MSG_R01_T0.raw

The raw DDA data have been searched against the SwissProt *Homo sapiens* database (FASTA), which can be downloaded from [UniProt](https://www.uniprot.org/). The FASTA file (2018-04 version, 20,301 entries)
has been deposited to ProteomeXchange via the iProX partner repository with the data set identifier [`PXD014108`](http://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD014108) or [`IPX0001628000`](https://www.iprox.cn/page/project.html?id=IPX0001628000). 
- swissprot_human_201804_validated.fasta

The saved SpectroMine project is available at ProteomeXchange and iProX with identifier [`PXD014108`](http://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD014108) or [`IPX0001628000`](https://www.iprox.cn/page/project.html?id=IPX0001628000).
- HeLa_DDA.psar.zip

SpectroMine reports should be exported with the schema provided in the [`misc\SpectroMine_Report_Schema`](misc/SpectroMine_Report_Schema) folder.
- FragmentReport.rs

For detailed instructions, see [SpectroMine Manual](https://biognosys.com/resources/spectromine-manual).


## 3. Train a MS/MS Model
### 3.1. Prepare Training Data
Create a directory for the project. Rename and place the fragment report in the project folder:
- HeLa.FragmentReport.csv

Open Anaconda PowerShell Prompt, set the project folder as working directory, and set path to the scripts as global parameter.
``` powershell
cd "Path_to_project_data"
$script_path = "Path_to_DeepDIA\src"
```

Training data can be converted from a SpectroMine fragment report.
``` powershell
python $script_path\extract_ms2_from_SpectroMine.py `
--in HeLa.FragmentReport.csv `
--out HeLa.assay.pickle
```

Run `remove_redundant_assays.py` to get a unique MS/MS spectrum for each peptide precursor.
``` powershell
python $script_path\remove_redundant_assays.py `
--in HeLa.assay.pickle `
--out HeLa_nonredundant.assay.pickle `
--action best `
--score qvalue `
--score_ascending `
--across_run
```

By setting `--score qvalue` and `--score_ascending`, the one with minimum q-value is selected from replicate spectra of each peptide precursor.

Run `convert_assays_to_ions.py`.
``` powershell
python $script_path\convert_assays_to_ions.py `
--in HeLa_nonredundant.assay.pickle `
--out HeLa.ions.json
```

By default, We only use the peptide precursors of charge 2+ and 3+ since others do not contain sufficient numbers of spectra for model training.

The MS/MS ion intensities are saved in JSON files. 
- HeLa_charge2.ions.json
- HeLa_charge3.ions.json

Move them into saperate subfolders `charge2` and `charge3`.
``` powershell
mkdir charge2
mv HeLa_charge2.ions.json charge2
mkdir charge3
mv HeLa_charge3.ions.json charge3
```

### 3.2. Train a Model
Run `train_ms2.py` in the `charge2` directory.
``` powershell
cd charge2
python $script_path\train_ms2.py `
--in HeLa_charge2.ions.json
```

Expected run time depends on the number of peptide spectra and the performance of the computer. In this toturial, this command may take several hours.

In the `models` folder, we find the trained model (with checkpoints during training) for charge 2+ peptides.
- charge2\models\epoch_*.hdf5

### 3.3. Test the Trained Model
Predict MS/MS ion intensities using the trained model and calculate similarties between predicted and experimental spectra.
``` powershell
python $script_path\predict_ms2.py `
--in HeLa_charge2.ions.json `
--score `
--charge 2
```

It is recommanded to test the model with another dataset. Here the training data (`HeLa_charge2.ions.json`) is used only as demo.

The predicted MS/MS ion intensities are saved in a JSON file. The similarty scores are saved in a CSV file.
- HeLa_charge2.prediction.ions.json
- HeLa_charge2.prediction.ions_score.csv

The median and quantiles of the similarty scores are printed in the prompt.


Train and test the model for charge 3+ following the same steps.
``` powershell
cd ..\charge3
python $script_path\train_ms2.py `
--in HeLa_charge3.ions.json

python $script_path\predict_ms2.py `
--in HeLa_charge3.ions.json `
--score `
--charge 3
```


## 4. Train an iRT Model
### 4.1. Prepare Training Data
Set the project folder as working directory.
``` powershell
cd "Path_to_project_data"
```

Training data can be converted can be converted from a SpectroMine fragment report.
``` powershell
python $script_path\extract_rt_from_SpectroMine.py `
--in HeLa.FragmentReport.csv `
--out HeLa.irt.csv
```

The iRT are saved in a CSV file.
- HeLa.irt.csv

Move them into a subfolder `irt`.
``` powershell
mkdir irt
mv HeLa.irt.csv irt
```

### 4.2. Train a Model
Run `train_rt.py` in the `irt` directory.
``` powershell
cd irt
python $script_path\train_rt.py `
--in HeLa.irt.csv
```

In the `models` folder, we find the trained model (with checkpoints during training).
- irt\models\epoch_*.hdf5

### 4.3. Test the Trained Model
Predict peptide iRT using the trained model and calculate correlations and error between predicted and experimental values.
``` powershell
python $script_path\predict_rt.py `
--in HeLa.irt.csv `
--score
```

It is recommanded to test the model with another dataset. Here the training data (`HeLa.irt.csv`) is used only as demo.

The predicted iRT and the prediction errors (differences between predicted and experimental values) are saved in a CSV file.
- HeLa.prediction.irt.csv
- HeLa.prediction.irt_score.csv

The correlation between predicted and experimental values, as well as the ranges of the prediction errors are printed in the prompt.
