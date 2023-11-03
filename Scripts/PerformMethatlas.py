import numpy as np
import pandas as pd
from scipy import optimize
import os.path as op
import sys
from multiprocessing import Pool

class Deconvolve:
    def __init__(self, atlas_path, samp_path, out_dir, resid, slim=False, plot=False):
        self.out_dir = out_dir  # Output dir to save mixture results and plot
        self.slim = slim  # Write results table w\o indexes and header (bool)
        self.plot = plot  # Plot results (bool)
        self.resid = resid  # Output residuals as well
        self.out_bname = self.get_bname(samp_path)  # output files path w/o extension

        # Load input files:
        self.atlas = self.load_atlas(atlas_path)  # Atlas
        self.samples = self.load_sample(samp_path)  # Samples to deconvolve

    def get_bname(self, samp_path):
        """
        Compose output files path:
        join the out_dir path with the basename of the samples file
        remove csv and gz extensions.
        """
        base_fname = op.basename(samp_path)

        if base_fname.endswith('.gz'):
            base_fname = op.splitext(base_fname)[0]
        base_fname = op.splitext(base_fname)[0]
        return op.join(self.out_dir, base_fname)

    @staticmethod
    def load_atlas(atlas_path):
        """
        Read the atlas csv file, save data in self.atlas
        :param atlas_path: Path to the atlas csv file
        """
        # validate path:
        # Deconvolve._validate_csv_file(atlas_path)

        # Read atlas, sort it and drop duplicates
        # print('atlas_path', atlas_path)
        df = pd.read_csv(atlas_path)
        df.rename(columns={list(df)[0]: 'acc'}, inplace=True)
        df = df.sort_values(by='acc').drop_duplicates(subset='acc').reset_index(drop=True)
        return df

    @staticmethod
    def decon_single_samp(samp, atlas):
        """
        Deconv ng NNLS, to get the mixture coefficients.
        :param samp: a vector of a single sample
        :param atlas: the atlas DadtaFrame
        :return: the mixture coefficients (of size 25)
        """

        name = samp.columns[1]

        # remove missing sites from both sample and atlas:
        data = samp.merge(atlas, on='acc', how='inner').copy().dropna(axis=0)
        if data.empty:
            print('Warning: skipping an empty sample {}'.format(name), file=sys.stderr)
            # print('Dropped {} missing sites'.format(self.atlas.shape[0] - red_atlas.shape[0]))
            return np.nan
        print('{}: {} sites'.format(name, data.shape[0]), file=sys.stderr)
        del data['acc']

        samp = data.iloc[:, 0]
        red_atlas = data.iloc[:, 1:]

        # get the mixture coefficients by deconvolution (non-negative least squares)
        mixture, residual = optimize.nnls(red_atlas, samp)
        mixture /= np.sum(mixture)
        return mixture, residual

    def load_sample(self, samp_path):
        """
        Read samples csv file. Reduce it to the atlas sites, and save data in self.samples
        Note: samples file must contain a header line.
        """

        # validate path:

        samples = pd.read_csv(samp_path)
        samples.rename(columns={list(samples)[0]: 'acc'}, inplace=True)
        samples = samples.sort_values(by='acc').drop_duplicates(subset='acc').reset_index(drop=True)
        samples = samples.merge(self.atlas['acc'].to_frame(), how='inner', on='acc')
        return samples

    def run(self):

        # run deconvolution on all samples in parallel
        processes = []
        with Pool() as p:
            for i, smp_name in enumerate(list(self.samples)[1:]):
                params = (self.samples[['acc', smp_name]], self.atlas)
                processes.append(p.apply_async(Deconvolve.decon_single_samp, params))
            p.close()
            p.join()

        self.samples = self.samples.iloc[:, 1:]

        # collect the results to 'res_table':
        arr = [pr.get() for pr in processes]
        res_table = np.empty((self.atlas.shape[1] - 1, self.samples.shape[1]))
        resids_table = np.empty((self.samples.shape[1], 1))
        for i in range(len(arr)):
            res_table[:, i], resids_table[i] = arr[i]
        df = pd.DataFrame(res_table, columns=self.samples.columns, index=list(self.atlas.columns)[1:])
        rf = pd.DataFrame(resids_table, columns=['Residuals'], index=self.samples.columns)
        return df, rf



real = pd.read_csv([sys.argv][0], index_col = 0, sep = '\t')
methods = ['none', 'z_score', 'min_max', 'col_z_score', 'col_min_max', 'QN', 'lognorm']
for method in methods:
    path = f'./Cache/{method}_reference.csv'
    ref = pd.read_csv(path, index_col=0)
    df = pd.read_csv(f'./Cache/{method}_mixture.csv', index_col=0)



    deco = Deconvolve(atlas_path=path, samp_path=f'./Cache/{method}_mixture.csv',
                        out_dir='.', resid=True)
    df2 = pd.DataFrame(np.empty((len(df.columns), len(ref.columns),)))
    df2[:] = np.nan
    df2.index = df.columns
    df2.columns = ref.columns
    for i in deco.samples.columns[1:]:
        df2.loc[i] = deco.decon_single_samp(samp=pd.DataFrame(deco.samples[['acc', i]]).replace([np.inf, -np.inf], np.nan).dropna(), atlas=deco.atlas.replace([np.inf, -np.inf], np.nan).dropna())[0]
    df2.to_csv(f'./Cache/{method}_results_meth_atl.csv')