from assay.assay2table import AssayToDataFrameConverter
from assay.modseq import stringify_modification

import pandas as pd
import numpy as np
import copy


class RetentionTimeCalibrator():
    def __init__(self, model='interpolate', 
                 smooth='savgol', smooth_args=None):
        columns = [
            {
                'name': 'run',
                'path': ['metadata', 'file'],
                'convert': str,
                'default': ''
            },
            {
                'name': 'peptideSequence',
                'path': 'peptideSequence'
            },
            {
                'name': 'modification',
                'path': 'modification',
                'convert': stringify_modification,
            },
            {
                'name': 'rt',
                'path': 'rt'
            }
        ]
        self.data_converter = AssayToDataFrameConverter(columns=columns)
        
        if model == 'linear':
            def linear(x, y, new_x):
                coef = np.polyfit(x, y, 1)
                return coef[0] * new_x + coef[1]
            
            self.model_func = linear
        elif model == 'interpolate':
            from scipy.interpolate import interp1d
            
            def interpolate(x, y, new_x):  
                x, index = np.unique(x, return_index=True)
                y = y[index]
                interp = interp1d(x, y, fill_value='extrapolate')
                return interp(new_x)
            
            self.model_func = interpolate
        else:
            raise ValueError('invalid model: ' + str(model))
        
        if smooth == 'savgol':
            from scipy.signal import savgol_filter
            
            def savgol(x, y, **kwargs):
                smooth_args = {
                    'window_length': 7, 
                    'polyorder': 1
                }
                smooth_args.update(kwargs)
                
                x = savgol_filter(
                    x, **smooth_args
                )
                return x, y
            
            self.smooth_func = savgol
        elif smooth == 'lowess':
            import statsmodels.api as sm
            
            def lowess(x, y, **kwargs):
                r = sm.nonparametric.lowess(y, x, **kwargs)
                
                # https://github.com/statsmodels/statsmodels/issues/2449
                if any(np.isnan(r[:, 1])):
                    data = pd.DataFrame.from_dict({'x': x, 'y': y}) \
                        .groupby(x).mean()
                    x = data['x']
                    y = data['y']
                    r = sm.nonparametric.lowess(y, x, **kwargs)
                    
                return r[:, 0], r[:, 1]
            
            self.smooth_func = lowess
        elif smooth is None or smooth == 'none' or smooth == 'None':
            self.smooth_func = None
        else:
            raise ValueError('invalid smooth: ' + str(smooth))
        
        self.smooth_args = smooth_args
            
        
    def load_reference(self, reference_assays):
        reference_data = self.data_converter \
            .assays_to_dataframe(reference_assays)
        self.load_reference_data(reference_data)
    
    def load_reference_data(self, reference_data):
        reference_data.insert(0, 'index', reference_data.index)
        self.reference_data = reference_data
        
    
    def calculate_rt(self, data):        
        merged_data = pd.merge(
            self.reference_data.drop(columns=['index', 'run']), data, 
            on=self.reference_data.columns.drop(['index', 'run', 'rt']).tolist(),
            suffixes=['_reference', '']
        )   
        
        index = np.argsort(merged_data['rt_reference'])
        y = merged_data['rt_reference'][index].values
        x = merged_data['rt'][index].values
           
        if self.smooth_func is not None:
            x, y = self.smooth_func(x, y, **(self.smooth_args or {}))
        
        if any(map(lambda x: not x > 0, x)):
            raise ValueError(x)
        if any(map(lambda x: not x > 0, y)):
            raise ValueError(y)
            
        y_new = self.model_func(x, y, data['rt'].values)
        
#        if any(map(lambda x: not x > 0, y_new)):
#            raise ValueError(y_new)

        return y_new
        
            
    def calibrate_rt(self, assays, multiple_runs=False, inplace=False, return_data=False):
        assay_data = self.data_converter \
            .assays_to_dataframe(assays)
            
        assay_data = self.calibrate_rt_data(assay_data, multiple_runs=multiple_runs)               
                
        if not inplace:
            assays = copy.deepcopy(assays)
            
        for i, assay in enumerate(assays):
            assay['rt'] = float(assay_data['rt_new'][i])
        
        if return_data:    
            assay_data = pd.merge(
                assay_data, 
                self.reference_data.drop(columns=['run']) \
                    .rename(columns={
                        'rt': 'rt_reference',
                        'index': 'index_reference'
                    }), 
                on=self.reference_data.columns. \
                    drop(['index', 'run', 'rt']).tolist(),
                how='left',
                suffixes=['', '_reference']
            )
            assay_data['index_reference'] = assay_data['index_reference'] \
                .fillna(-1).astype(int)
            
            if not inplace:
                return assays, assay_data
            else:
                return assay_data
        
        return assays

    def calibrate_rt_data(self, assay_data, multiple_runs=False):           
        if not multiple_runs:                        
            rt_new = self.calculate_rt(assay_data)               
            
        else:
            rt_new = pd.Series(None, index=assay_data.index)
            for run, data in assay_data.groupby(by=['run'], group_keys=False):
                rt_new[data.index] = self.calculate_rt(data) 
           
        assay_data.rename(columns={'rt': 'rt_old'}, inplace=True)
        assay_data['rt_new'] = rt_new
            
        return assay_data
            
        