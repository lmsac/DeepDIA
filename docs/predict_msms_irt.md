# DeepDIA Tutorial: Spectral Library Generation From Peptide Lists
Using deep learning to generate in silico spectral libraries from peptide lists for data-independent acquisition (DIA) analysis. 

## 1. Prerequisites
### 1.1. System Requirements
In this tutorial, spectral library generation has been tested on a workstation with Intel Xeon E5-2690 v3 CPU, 16 GB RAM, and Microsoft Windows Server 2016 Version 1607 (OS Build 14393.2430) operating system.

After spectral library generation, DIA data analysis is performed on a workstation with Intel Core i9-7960X CPU, 128 GB RAM, and Microsoft Windows 10 Version 1809 (OS Build 17763.503) 64-bit operating system.

### 1.2. Software Dependency
The following software and packages are required:
- Python (version 3.7 or later, [Anaconda](https://www.anaconda.com/) distribution is recommended)
- [TensorFlow](https://www.tensorflow.org/) (version 2.0 or later)
- [Keras](https://keras.io/) (packaged with TensorFlow)

DeepDIA requires the following Python packages integrated in Anaconda:
- numpy (version 1.18.5)
- pandas (version 0.25.3)
- scipy (version 1.4.1)
- statsmodels (version 0.13.2)

Later versions may be compatible, but have not been tested.

The following software is not needed in the process of spectral library generation, but is used in the complete DIA data analysis workflow:
- [Spectronaut](https://biognosys.com/software/spectronaut) (Biognosys AG, Schlieren, Switzerland)


## 2. Tutorial Data
### 2.1. Starting Materials for Spectral Library Generation
Starting materials of this tutorial are available at ProteomeXchange and iProX with identifier [`PXD014108`](http://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD014108) or [`IPX0001628000`](https://www.iprox.cn/page/project.html?id=IPX0001628000). 

#### 2.1.1. Models for MS/MS and iRT Prediction
Pre-trained models: 
- HeLa.model.zip

They are also provided in [`data\models`](data/models) folder in this repositry.

The models have been trained with HeLa data on Q Exactive HF (Bruderer, R. et al. Mol. Cell. Proteomics 2017, 16, 2296-2309, [doi:10.1074/mcp.RA117.000314](https://doi.org/10.1074/mcp.RA117.000314)).

#### 2.1.2. Peptide Lists
A peptide list is stored in a comma-separated values (CSV) file including a columns named `protein` and `sequence`.  

Example peptide lists are provided in [`data\peptide`](data/peptide) folder:
- Pan_human.peptide.csv
- Pan_human_charge2.peptide.csv
- Pan_human_charge3.peptide.csv

The peptide lists have been collected from the Pan Human Library (Rosenberger, G. et al. Sci. Data 2014, 1, 140031, [doi:10.1038/sdata.2014.31](https://doi.org/10.1038/sdata.2014.31)).

Note that DeepDIA only supports peptide sequences with standard amino acids (ACDEFGHIKLMNPQRSTVWY) and length <= 50.

### 2.2. Tutorial Data for DIA Analysis
These files are not needed in the process of spectral library generation, but are used in the complete DIA data analysis workflow.

#### 2.2.1. LC-MS/MS Raw Data
LC-MS/MS data of 3 DIA technical replicates of 2 h gradient of HeLa cells on Q Exactive HF are available at ProteomeXchange with the data set [`PXD005573`](http://proteomecentral.proteomexchange.org/GetDataset?ID=PXD005573) (Bruderer, R. et al. Mol. Cell. Proteomics 2017, 16, 2296-2309, [doi:10.1074/mcp.RA117.000314](https://doi.org/10.1074/mcp.RA117.000314)).
- Fig1_MP-DIA-120min120kMS1-22W30k-8dppp_MHRM_R01.raw
- Fig1_MP-DIA-120min120kMS1-22W30k-8dppp_MHRM_R02.raw
- Fig1_MP-DIA-120min120kMS1-22W30k-8dppp_MHRM_R03.raw

#### 2.2.2. Protein Sequence Database
SwissProt *Homo sapiens* database (FASTA) can be downloaded from [UniProt](https://www.uniprot.org/). The FASTA file (2018-04 version, 20,301 entries)
has been deposited to ProteomeXchange via the iProX partner repository with the data set identifier [`PXD014108`](http://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD014108) or [`IPX0001628000`](https://www.iprox.cn/page/project.html?id=IPX0001628000). 
- swissprot_human_201804_validated.fasta


## 3. Spectral Library Generation
### 3.1. Prepare Starting Materials
Create a directory for the project. Place the peptide lists and models in the project folder in the following paths:
- Pan_human.peptide.csv
- Pan_human_charge2.peptide.csv
- Pan_human_charge3.peptide.csv
- models\charge2\epoch_035.hdf5
- models\charge3\epoch_034.hdf5
- models\irt\epoch_082.hdf5

Open Anaconda PowerShell Prompt, set the project folder as working directory, and set path to the scripts as global parameter.
``` powershell
cd "Path_to_project_data"
$script_path = "Path_to_DeepDIA\src"
```

### 3.2. Predict MS/MS Spectra
Run `predict_ms2.py` to predict MS/MS ion intensities for peptide precursors with charge 2+.
``` powershell
python $script_path\predict_ms2.py `
--in Pan_human_charge2.peptide.csv `
--model models\charge2\epoch_035.hdf5 `
--charge 2 `
--out Pan_human_charge2.prediction.ions.json
```

The predicted MS/MS ion intensities are saved in a JSON file.
- Pan_human_charge2.prediction.ions.json

Predict MS/MS spectra for charge 3+ following the same steps.
``` powershell
python $script_path\predict_ms2.py `
--in Pan_human_charge3.peptide.csv `
--model models\charge3\epoch_034.hdf5 `
--charge 3 `
--out Pan_human_charge3.prediction.ions.json
```

### 3.3. Predict iRT
Run `predict_rt.py`.
``` powershell
python $script_path\predict_rt.py `
--in Pan_human.peptide.csv `
--model models\irt\epoch_082.hdf5 `
--out Pan_human.prediction.irt.csv
```

The predicted iRT values are saved in a CSV file.
- Pan_human.prediction.irt.csv

### 3.4. Generate Spectral Library
Ensure that the predicted MS/MS and iRT files are present in the project folder.
- Pan_human.peptide.csv
- Pan_human_charge2.prediction.ions.json
- Pan_human_charge3.prediction.ions.json
- Pan_human.prediction.irt.csv

Run `build_assays_from_prediction.py`.
``` powershell
python $script_path\build_assays_from_prediction.py `
--peptide Pan_human.peptide.csv `
--ions Pan_human_charge2.prediction.ions.json `
       Pan_human_charge3.prediction.ions.json `
--rt Pan_human.prediction.irt.csv `
--out Pan_human.prediction.assay.pickle
```

The generated spectral library is saved in a Python binary file.
- Pan_human.prediction.assay.pickle

Run `convert_assays_to_Spectronaut_library.py`.
``` powershell
python $script_path\convert_assays_to_Spectronaut_library.py `
--in Pan_human.prediction.assay.pickle `
--out Pan_human.prediction.library.csv
```

The generated spectral library is converted to a CSV file that is compatible with Spectronaut.
- Pan_human.prediction.library.csv


## 4. DIA Data Analysis
The predicted spectral library (`Pan_human.prediction.library.csv`) can be imported into Spectronaut with the FASTA file (`swissprot_human_201804_validated.fasta`) to analyze the DIA raw data.

For detailed instructions, see [Spectronaut Manual](https://biognosys.com/resources/spectronaut-manual
).

