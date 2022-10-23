import os
import numpy as np

from pepdetect import PeptideDetectabilityTrainer, \
    PeptideDetectabilityPredictor
from util import list_files

def find_initial_negative_subset(data1, data2, ratio=1, seed=None):
    if seed is None:
        random = np.random
    else:
        random = np.random.RandomState(seed=seed)

    if (data1.shape[0] * ratio >= data2.shape[0]):
        return np.arange(data2.shape[0])
    else:
        indexes = random.permutation(data2.shape[0])
        return indexes[:(data1.shape[0] * ratio)]


def evaluate_model(data1, data2, options,
                   model_path=None, model=None,
                   positive_threshold=0.5):
    predictor = PeptideDetectabilityPredictor(
        options=options,
        model_path=model_path, model=model
    )
    prediction1 = predictor.predict(data1)
    prediction2 = predictor.predict(data2)
    tp = np.where(prediction1["detectability"].values >= positive_threshold)[0]
    fn = np.where(prediction1["detectability"].values < positive_threshold)[0]
    fp = np.where(prediction2["detectability"].values >= positive_threshold)[0]
    tn = np.where(prediction2["detectability"].values < positive_threshold)[0]
    tpr = len(tp) / (len(tp) + len(fn))
    fpr = len(fp) / (len(fp) + len(tn))
    precision = len(tp) / (len(tp) + len(fp))
    return tpr, fpr, precision, \
        tp, fn, fp, tn, \
        prediction1[["detectability"]].values, \
        prediction2[["detectability"]].values


def find_hard_negative_subset(data1, data2, options,
                              model_path=None, model=None,
                              positive_threshold=0.5,
                              min_ratio=None, max_ratio=None):
    tpr, fpr, precision, \
        tp, fn, fp, tn, \
        prediction1, prediction2 = \
        evaluate_model(
            data1, data2, options=options,
            model_path=model_path, model=model,
            positive_threshold=positive_threshold
        )

    indexes = fp

    if min_ratio is not None and min_ratio > 0:
        min_count = data1.shape[0] * min_ratio
    else:
        min_count = None
    if max_ratio is not None and max_ratio > 0:
        max_count = data1.shape[0] * max_ratio
    else:
        max_count = None
    if min_count is not None and \
        max_count is not None and \
        min_count > max_count:
            raise ValueError('min_count > max_count')

    if min_count is not None and min_count > len(indexes):
        indexes = np.concatenate(
            (indexes, tn[np.argsort(-prediction2[tn], axis=None)])
        )[:min_count]
    elif max_count is not None and max_count < len(indexes):
        indexes = indexes[np.argsort(-prediction2[indexes], axis=None) \
                          [:max_count]]

    return indexes, tpr, fpr, precision


def train_hard_negative_iter(data_positive, data_negative, options,
                             model_path=None,
                             working_dir='.',
                             max_rounds=10,
                             positive_threshold=0.5,
                             fpr_threshold=0.01,
                             seed=None,
                             lr=0.001,
                             reduce_lr_patience=5,
                             reduce_lr_factor=0.1,
                             log=print,
                             **kwargs):
    result = None

    for i in range(0, max_rounds + 1):
        if i == 0:
            indexes_negative = find_initial_negative_subset(
                data_positive, data_negative, seed=seed
            )

        else:
            model_path = list_files(
                path=os.path.join(
                    working_dir,
                    'training_{i}'.format(i=i - 1),
                    'models'
                ),
                pattern=r'^epoch_[0-9]+\.hdf5$',
                recursive=True
            )[-1]
            indexes_negative, tpr, fpr, precision = find_hard_negative_subset(
                data_positive, data_negative, options=options,
                model_path=model_path,
                positive_threshold=positive_threshold,
                min_ratio=1, max_ratio=1
            )

            if log is not None:
                log('TPR: {tpr}, FPR: {fpr}, Precision: {precision}'.format(
                    tpr=tpr, fpr=fpr, precision=precision
                ))

            if result is not None:
                result['evaluate'] = {
                    'tpr': tpr,
                    'fpr': fpr,
                    'precision': precision
                }

                yield result

            if (fpr <= fpr_threshold):
                if log is not None:
                    log('Early stopping')
                break

        if (i >= max_rounds):
            break

        if log is not None:
            log('Round {i}: train on {n} positive samples, {m} negative samples' \
                .format(
                    i=i, n=len(data_positive), m=len(indexes_negative)
                )
            )

        train_dir = os.path.join(
            working_dir,
            'training_{i}'.format(i=i)
        )
        os.makedirs(train_dir, exist_ok=True)
        os.makedirs(os.path.join(train_dir, 'models'), exist_ok=True)

        data1 = data_positive
        data2 = data_negative.iloc[indexes_negative]

        trainer = PeptideDetectabilityTrainer(
            options=options,
            save_path=os.path.join(
                train_dir, 'models', 'epoch_{epoch:03d}.hdf5'
            ),
            log_path=os.path.join(train_dir, 'training.log')
        )
        if model_path is not None:
            trainer.load_model(model_path)

        if reduce_lr_patience and reduce_lr_factor:
            trainer.use_reduced_lr(
                patience=reduce_lr_patience,
                factor=reduce_lr_factor
            )
        
        if lr is not None:
            log('use lr: ' + str(lr))

        result = trainer.train(
            data1, data2, seed=seed, 
            lr=lr,
            **kwargs
        )
        trainer.save_model(os.path.join(
            train_dir, 'models', 'last_epoch.hdf5'
        ))

        result['split'] = {
            'seed': seed,
            'round': i,
            'train': {
                'positive': result['split']['train']['positive'],
                'negative': indexes_negative \
                    [np.array(result['split']['train']['negative'])].tolist()
            },
            'validate': {
                'positive': result['split']['validate']['positive'],
                'negative': indexes_negative \
                    [np.array(result['split']['validate']['negative'])] \
                    .tolist()
            }
        }

