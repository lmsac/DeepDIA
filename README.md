# DeepDIA
Using deep learning to generate in silico spectral libraries for data-independent acquisition (DIA) analysis. 

## Updates
1.1.0
- Dependency of R removed
- FASTA digestion
- Ion mobility prediction (experimental)

For the version of the [*Nat Commun* 2020 publication](https://doi.org/10.1038/s41467-019-13866-z), please refer to the commit [#674e2fb](https://github.com/lmsac/DeepDIA/tree/674e2fb8f85542243eed6340939df54864178552).

## Dependency
The following software and packages are required:
- Python (version 3.7 or later, [Anaconda](https://www.anaconda.com/) distribution is recommended)
- [TensorFlow](https://www.tensorflow.org/) (version 2.0 or later)
- [Keras](https://keras.io/) (packaged with TensorFlow)

For spectral library generation from FASTA files and data preprocessing for training detectability models, the following package is required:
- [Biopython](https://biopython.org/) (version 1.70 or later)

DeepDIA requires the following Python packages integrated in Anaconda:
- numpy (version 1.18.5)
- pandas (version 0.25.3)
- scipy (version 1.4.1)
- statsmodels (version 0.13.2)

Later versions may be compatible, but have not been tested.

For model training, NVIDIA graphics cards with CUDA are recommended.
- [CUDA Toolkit](https://developer.nvidia.com/cuda-toolkit-archive) (version 11.2)
- [cuDNN SDK](https://developer.nvidia.com/cudnn) (version 8.1.0)

## Installation
### 1. Install Python (Anaconda)
Download and install [Anaconda](https://www.anaconda.com).

Check successful installation by in the Anaconda Prompt:
```
pip list
```

Ensure that the following Python packages are installed: numpy, pandas, scipy, and statsmodels. If not, install the missing packages using the following command (as an example for statsmodels):
```
pip install statsmodels
```

### 2. Install TensorFlow
Ensure that NVIDIA GPU driver has been installed. Install the CUDA and cuDNN with conda. This step can be skipped if you run TensorFlow on CPU only.
```
conda install -c conda-forge cudatoolkit=11.2 cudnn=8.1.0
```

Install TensorFlow using `pip`:
```
pip install tensorflow
```

### 3. Install Biopython
Install Biopython using `pip`:
```
pip install biopython
```
or `conda`:
```
conda install -c conda-forge biopython
```

## Getting Started
### 1. Prepare a Peptide List
A peptide list is stored in a comma-separated values (CSV) file including columns named `protein` and `sequence`.  
```
"protein","sequence"
"O43504","HDGITVAVHK"
"P56470","VGSSGDIALHINPR"
"Q9UHL4","LDHFNFER"
"P68371","IREEYPDR"
"P01024","AKDQLTCNK"
```
Peptides can be collected from public resources. 
From the Pan Human Library (Rosenberger, G. et al. Sci. Data 2014, 1, 140031, [doi:10.1038/sdata.2014.31](https://doi.org/10.1038/sdata.2014.31)), peptide lists have been collected and provided as an example in [`data\peptide`](data/peptide) folder:
- Pan_human.peptide.csv
- Pan_human_charge2.peptide.csv
- Pan_human_charge3.peptide.csv

DeepDIA only supports peptide sequences with standard amino acids (ACDEFGHIKLMNPQRSTVWY) and length <= 50.

### 2. Predict MS/MS Spectra
Prepare a model for MS/MS prediction.
You can use pre-trained models or train your own models. A model trained with HeLa data on Q Exactive HF (Bruderer, R. et al. Mol. Cell. Proteomics 2017, 16, 2296-2309, [doi:10.1074/mcp.RA117.000314](https://doi.org/10.1074/mcp.RA117.000314)) is provided as an example in [`data\models`](data/models) folder: 
- data\models\charge2\epoch_035.hdf5
- data\models\charge3\epoch_034.hdf5

Run `predict_ms2.py` to predict MS/MS ion intensities for peptide precursors with charge 2+.
``` powershell
python src\predict_ms2.py `
--in data\peptide\Pan_human_charge2.peptide.csv `
--model data\models\charge2\epoch_035.hdf5 `
--charge 2 `
--out data\Pan_human_charge2.prediction.ions.json
```

The predicted MS/MS ion intensities are saved in a JSON file (`*.prediction.ions.json`).

Predict MS/MS for charge 3+ following the same steps.
``` powershell
python src\predict_ms2.py `
--in data\peptide\Pan_human_charge3.peptide.csv `
--model data\models\charge3\epoch_034.hdf5 `
--charge 3 `
--out data\Pan_human_charge3.prediction.ions.json
```

### 3. Predict iRT
Prepare a model for iRT prediction.
You can use pre-trained models or train your own models. A pretrained model is provided as an example in [`data\models`](data/models) folder: 
- data\models\irt\epoch_082.hdf5

Run `predict_rt.py`.
``` powershell
python src\predict_rt.py `
--in data\peptide\Pan_human.peptide.csv `
--model data\models\irt\epoch_082.hdf5 `
--out data\Pan_human.prediction.irt.csv
```

The predicted iRT values are saved in a CSV file (`*.prediction.irt.csv`).

### 4. Generate Spectral Library
Ensure that the predicted MS/MS and iRT files are present in the `data` folder.

Run `build_assays_from_prediction.py`.
``` powershell
python src\build_assays_from_prediction.py `
--peptide data\peptide\Pan_human.peptide.csv `
--ions data\Pan_human_charge2.prediction.ions.json `
       data\Pan_human_charge3.prediction.ions.json `
--rt data\Pan_human.prediction.irt.csv `
--out data\Pan_human.prediction.assay.pickle
```

The generated spectral library is saved in a Python binary file (`*.assay.pickle`).

Run `convert_assays_to_Spectronaut_library.py`.
``` powershell
python src\convert_assays_to_Spectronaut_library.py `
--in data\Pan_human.prediction.assay.pickle `
--out data\Pan_human.prediction.library.xls
```

The generated spectral library is converted to a speadsheet file (`*.library.xls`) that is compatible with Spectronaut and DIA-NN.


## Tutorial
Tutorials are avaliable in the [`docs`](docs) folder.

### Spectral Library Pretiction
[**DeepDIA Tutorial: Spectral Library Generation From Peptide Lists**](docs/predict_msms_irt.md) describes the workflow to generate in silico spectral libraries from peptide lists.

### Detectability Prediction
[**DeepDIA Tutorial: Spectral Library Generation with Detectability Prediction**](docs/predict_detectability.md) describes the complete workflow to generate in silico spectral libraries from proteome databases with detectability filtering.

### Model Training
[**DeepDIA Tutorial: Training New Models for MS/MS and iRT Prediction**](docs/train_msms_irt.md) describes the workflow for training new models for MS/MS and iRT prediction using data-dependent acquisition (DDA) data.

[**DeepDIA Tutorial: Training a New Model for Detectability Prediction**](docs/train_detectability.md) describes the workflow for training a new model for MS detectability prediction using data-dependent acquisition DDA data.


## Publications
Yang, Y., Liu, X., Shen, C., Lin, Y., Yang, P., Qiao, L. In silico spectral libraries by deep learning facilitate data-independent acquisition proteomics. *Nat Commun* **11**, 146 (2020). https://doi.org/10.1038/s41467-019-13866-z.

## License
DeepDIA is distributed under a BSD license. See the LICENSE file for details.

## Contacts
Please report any problems directly to the github issue tracker. Also, you can send feedback to liang_qiao@fudan.edu.cn.
