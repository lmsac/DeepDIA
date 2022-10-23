import argparse

parser = argparse.ArgumentParser(
    description='Predict peptide retention time.'
)
parser.add_argument(
    '--in', nargs='+',
    help='input peptide list files'
)
parser.add_argument(
    '--model',
    help='model file'
)
parser.add_argument(
    '--reference',
    help='reference file'
)
parser.add_argument(
    '--out', nargs='+',
    help='output retention time files'
)
parser.add_argument(
    '--score', action='store_true', default=False,
    help='calculate correlations and error between predicted and experimental values (input files must contain experimental retention time values)'
)

args = parser.parse_args()
peptide_files = getattr(args, 'in')
model_file = args.model
reference_file = args.reference
out_files = args.out
score = args.score

# %%
import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(filename)s: [%(levelname)s] %(message)s'
)

# %%
from util import list_files

if globals().get('peptide_files', None) is None:
    peptide_files = list_files(
        path='.',
        pattern='\\.peptide\\.csv$' if not globals().get('score', False) \
            else '(?<!prediction)\\.irt\\.csv$'
    )

if len(peptide_files) == 0:
    raise ValueError('no peptide list files')

# %%
import os

if globals().get('model_file', None) is None:
    model_files = list_files(
        path='models' if os.path.isdir('models') else '.',
        pattern='^epoch_[0-9]+\\.hdf5$',
        recursive=True
    )

    if len(model_files) == 0:
        raise ValueError('no model file')

    model_file = model_files[-1]

# %%
if globals().get('out_files', None) is None:
    out_files = []
    for peptide_file in peptide_files:
        out_file = os.path.splitext(peptide_file)[0]
        if out_file.endswith('.peptide'):
            out_file = out_file[:-len('.peptide')]
        elif out_file.endswith('.irt'):
            out_file = out_file[:-len('.irt')]
        out_file += '.prediction.irt.csv'
        out_files.append(out_file)

if len(out_files) != len(peptide_files):
    raise ValueError('numbers of peptide list files and output files not match')


# %%
import pandas as pd

from peprt import PeptideRTPredictor, PeptideRTOptions


options = PeptideRTOptions.default()


logging.info('use model: ' + model_file)

predictor = PeptideRTPredictor(
    options=options,
    model_path=model_file
)

# %%
if globals().get('reference_file', None) is not None:
    from assay.rtcalibration import RetentionTimeCalibrator

    rt_calibrator = RetentionTimeCalibrator(
        model='interpolate', 
        smooth='lowess',
        smooth_args={'frac': 0.667, 'it': 0}
    )
    
    logging.info('load references: ' + reference_file)

    reference = pd.read_csv(reference_file)

    if 'rt' not in reference.columns and 'irt' in reference.columns:
        reference.rename(columns={'irt': 'rt'}, inplace=True)

    rt_calibrator.load_reference_data(reference)
    
    logging.info('references loaded: {0} peptides' \
                 .format(len(rt_calibrator.reference_data)))


# %%
for peptide_file, out_file in zip(peptide_files, out_files):
    logging.info('load peptides: ' + peptide_file)

    peptides = pd.read_csv(peptide_file)

    peptides = peptides.loc[
        (peptides['sequence'].str.len() <= options.max_sequence_length) & \
        peptides['sequence'].map(lambda s: \
            all(map(lambda a: a in options.amino_acids, s))) \
    , :]

    logging.info('peptides loaded: {0} valid peptides' \
                 .format(len(peptides)))

    logging.info('predict peptide iRT: ' + peptide_file)

    sequences = peptides['sequence']
    modifications = peptides.get('modification', None)
    prediction = predictor.predict(sequences, modifications)

    logging.info('peptide iRT predicted: {0} peptides' \
                 .format(len(prediction)))

    if globals().get('reference_file', None) is not None:
        logging.info('calibrate peptide iRT: ' + peptide_file) 

        prediction.rename(columns={'irt': 'rt'}, inplace=True)

        prediction = rt_calibrator.calibrate_rt_data(prediction)

        prediction.rename(columns={
           'rt_old': 'irt_uncalibrated', 
           'rt_new': 'irt'
        }, inplace=True)

        logging.info('peptide iRT calibrated: ' + peptide_file) 

    logging.info('saving peptide iRT: {0}' \
        .format(out_file))

    prediction.to_csv(out_file, index=False)

    logging.info('peptide iRT saved: {0}, {1} peptides' \
        .format(out_file, len(prediction)))

    if globals().get('score', False):
        logging.info('calculate RT prediction accuracy: ' + out_file)
    
        scores = pd.DataFrame.from_dict({
            'sequence': sequences,
            'modification' : modifications,
            'irt_experimental': peptides['irt'],
            'irt_predicted': prediction['irt'],
        })

        scores['delta_irt'] = scores['irt_predicted'] - scores['irt_experimental']

        out_score_file = os.path.splitext(out_file)[0]
        if out_score_file.endswith('.irt'):
            out_score_file = out_score_file[:-len('.irt')]
        out_score_file += '.irt_score.csv' 
        
        scores.to_csv(out_score_file, index=False)
        
        logging.info('scores saved: {0}, {1} peptides' \
            .format(out_score_file, len(scores)))

        logging.info('iRT correlation: pearson={0}'.format(
            scores['irt_predicted'].corr(scores['irt_experimental'])
        ))
        
        logging.info('iRT difference: median={0}, IQR={1}, range(95%)={2}'.format(
            scores['delta_irt'].median(),
            scores['delta_irt'].quantile(0.75) - scores['delta_irt'].quantile(0.25),
            scores['delta_irt'].quantile(0.975) - scores['delta_irt'].quantile(0.025),
        ))

