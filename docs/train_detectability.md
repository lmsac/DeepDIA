# DeepDIA Tutorial: Training a New Model for Detectability Prediction
Training a new model for MS detectability prediction using data-dependent acquisition (DDA) data.

## 1. Prerequisites
See [**DeepDIA Tutorial: Spectral Library Generation From Peptide Lists**](predict_msms_irt.md) and [**DeepDIA Tutorial: Training New Models for MS/MS and iRT Prediction**](train_msms_irt.md).

For data preprocessing for training detectability models, the following package is required:
- [Biopython](https://biopython.org/) (version 1.70 or later)


## 2. Tutorial Data
### 2.1. Starting Materials for Model Training
Starting materials of this tutorial are available at ProteomeXchange and iProX with identifier [`PXD014108`](http://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD014108) or [`IPX0001628000`](https://www.iprox.cn/page/project.html?id=IPX0001628000). 

Peptide and protein reports exported from SpectroMine.
- HEK293_DDA.csv.zip

SwissProt *Homo sapiens* database (FASTA), which can be downloaded from [UniProt](https://www.uniprot.org/). The FASTA file (2018-04 version, 20,301 entries)
has been deposited to ProteomeXchange via the iProX partner repository with the data set identifier [`PXD014108`](http://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD014108) or [`IPX0001628000`](https://www.iprox.cn/page/project.html?id=IPX0001628000). 
- swissprot_human_201804_validated.fasta

### 2.2. Tutorial Data for DDA Analysis
These files are not needed in the process of model training, but are used in the DDA analysis to get training data.

LC-MS/MS DDA data of HEK-293 cells on Q Exactive HF are available at ProteomeXchange (http://proteomecentral.proteomexchange.org/) with the data set identifier [`PXD005573`](http://proteomecentral.proteomexchange.org/GetDataset?ID=PXD005573) (Bruderer, R. et al. Mol. Cell. Proteomics 2017, 16, 2296-2309, [doi:10.1074/mcp.RA117.000314](https://doi.org/10.1074/mcp.RA117.000314)).
- Fig4_HEK293-1m-HPRP-10perc_DDA_R01_T0.raw
- Fig4_HEK293-1m-HPRP-15perc_DDA_R01_T0.raw
- Fig4_HEK293-1m-HPRP-20perc_DDA_R01_T0.raw
- Fig4_HEK293-1m-HPRP-25perc_DDA_R01_T0.raw
- Fig4_HEK293-1m-HPRP-50perc_DDA_R01_T0.raw
- Fig4_HEK293-1m-HPRP-5perc_DDA_R01_T0.raw
- Fig4_HEK293-1m-HPRP-FT_DDA_R01_T0.raw
- Fig4_HEK293-1m_DDA_R01_T0.raw
- Fig4_HEK293-1m_DDA_R02_T0.raw
- Fig4_HEK293-1m_DDA_R03_T0.raw

The raw DDA data have been searched against the SwissProt *Homo sapiens* database. The saved SpectroMine projects are available at ProteomeXchange and iProX with identifier [`PXD014108`](http://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD014108) or [`IPX0001628000`](https://www.iprox.cn/page/project.html?id=IPX0001628000).
- HEK293_DDA.psar.zip

SpectroMine reports should be exported with the schema provided in the [`misc\SpectroMine_Report_Schema`](misc/SpectroMine_Report_Schema) folder.
- PeptideReport.rs
- ProteinReport.rs

For detailed instructions, see [SpectroMine Manual](https://biognosys.com/resources/spectromine-manual).


## 3. Train a Detectability Model
### 3.1. Prepare Training Data
Create a directory for the project. Rename and place the peptide and protein reports in the project folder:
- HEK293.PeptideReport.csv
- HEK293.ProteinReport.csv

Open Anaconda PowerShell Prompt, set the project folder as working directory, and set path to the scripts as global parameter.
``` powershell
cd "Path_to_project_data"
$script_path = "Path_to_DeepDIA\src"
```

Training data can be converted from SpectroMine peptide and protein reports.
``` powershell
python $script_path\extract_detect_from_SpectroMine.py `
--protein HEK293.ProteinReport.csv `
--peptide HEK293.PeptideReport.csv `
--fasta swissprot_human_201804_validated.fasta `
--out HEK293_excludeSingleHit_coverage25.detectability.csv `
--supplement_run_pattern "HPRP"
```

By setting the `--supplement_run_pattern` parameter, LC-MS/MS runs with fractionation (with file name containing `HPRP`) are used as supplement runs to those without fractionation. This setting is based on
the assumption that peptides observed in DDA without fractionation are more likely to be detected than those only observed in DDA with fractionation, and detectability scores of peptides only observed with fractionation are set to the average of 0.5 and the minimum of detectability scores of those detected without fractionation.

By default, `Trypsin/P` is selected as enzyme and digestion is performed with the following parameters:
- Maximum missed cleavages: 2
- Minimum peptide length: 7
- Maximum peptide length: 50
- Minimum peptide mass: 0
- Maximum peptide mass: 4000
- Remove N-terminal methionine: True

Ensure that these parameters are consistent with those used in database searching by SpectroMine.

Single hits are excluded and only proteins with sequence coverage >= 25% are taken into considerasion.

You can the view the parameters using the following command:
``` powershell
python $script_path\extract_detect_from_SpectroMine.py --help
```

The detectability scores are saved in CSV files.
- HEK293_excludeSingleHit_coverage25.detectability.csv
- HEK293_excludeSingleHit_coverage25_negative.detectability.csv

The negative file (`*_negative.detectability.csv`) contains the theoretical peptides not found in the experimental data.

Move them into a subfolder `detectability`.
``` powershell
mkdir detectability
mv *.detectability.csv detectability
```

### 3.2. Train a Model
Run `train_detectability_hard_negative.py` in the `detectability` directory.
``` powershell
cd detectability
python $script_path\train_detectability_hard_negative.py
```

Expected run time depends on the number of peptides and the performance of the computer. In this toturial, this command may take several hours to a day.

In the `training_*\models` folder, we find the trained model (with checkpoints during training).
- detectability\models\training_*\epoch_*.hdf5
