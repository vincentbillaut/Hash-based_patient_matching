import numpy as np
import pandas as pd
import gzip
import pickle
import os
import distance

import common.preprocessing as pp

from multiprocessing import Pool
from tqdm import tqdm


def file_len(path):
    with gzip.open(path, 'rb') as f:
        for i, _ in enumerate(f):
            pass
    return i + 1

def get_initial_ordering(path="Data/ALL.chr22.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf.gz"):
    """
    Returns initial ordering of patients in the database
    """
    with gzip.open(path, 'rb') as f:
        for i,l in enumerate(f):
            s = l.decode('utf-8')
            if i == 29:
                initial_order = pd.DataFrame(s.split()[9:])
                break
    return initial_order


class GenomeData:
    def __init__(self, name="chrom22", path="Data/22hap0.gz", verbose = True, test = False):
        self.name = name
        self.source_path = path if not test else "Data/22hap0_test.gz"
        self.verbose = verbose
        if self.verbose:
            print("Retrieving file length...")
        self.n_indiv = file_len(self.source_path)
        self.haploid0, self.haploid1, self.diploid = None, None, None
        self.ethnicity, self.initial_order, self.new_order = None, None, None

    def extract_data(self, path=None, verbose=None):
        if path is None:
            path=self.source_path
        if verbose is None:
            verbose=self.verbose
        with gzip.open(path, 'rb') as f:
            for l in f:
                break
            s = l.decode('utf8')
            n = len(s.split())
            data = np.zeros((self.n_indiv,n), dtype=int)
            f.seek(0)
            if verbose:
                for i,l in tqdm(enumerate(f), total=self.n_indiv,
                                desc="data extr. ({})".format(path)):
                    s = l.decode('utf8')
                    data[i,] = np.array(list(map(float, s.split())))
            else:
                for i,l in enumerate(f):
                    s = l.decode('utf8')
                    data[i,] = np.array(list(map(float, s.split())))
        if verbose:
            print("\tdone.")
        return data

    def extract_labels(self, path="Data/ALL.chr22.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf.gz"):
        if self.verbose:
            print("Extracting patient IDs...")
        with gzip.open(path, 'rb') as f:
            for i,l in enumerate(f):
                if i == 29:
                    s = l.decode('utf-8')
                    idList = s.split()[9:]
                    break
        self.index_to_id = idList
        self.id_to_index = {idList[i]:i for i in range(len(idList))}

    def extract_haploid1(self):
        if self.verbose:
            print("Extracting second haploid data...")
        self.haploid1 = self.extract_data(self.source_path.replace('0', '1'))

    def build_diploid(self):
        if self.verbose:
            print("Building diploid data...")
        self.diploid = self_get_haploid0() + self_get_haploid1()

    def extract_ethnicity(self, path='Data/integrated_call_samples.20101123.ALL.panel'):
        data = pd.DataFrame(columns=['patient', 'country', 'continent'])
        with open(path, 'r') as f:
            for i,l in enumerate(f):
                info = l.strip().split()
                data.loc[i] = info[:3]
        self.ethnicity = data.set_index('patient').sort_values(by=['continent', 'country'], ascending=[False, True])

    def extract_ethnicity_index(self):
        if self.ethnicity is None:
            self.extract_ethnicity()
        ethnicity = self.ethnicity
        if self.initial_order is None:
            self.initial_order = get_initial_ordering()
        self.initial_order['index1'] = self.initial_order.index
        self.initial_order.set_index(0)
        self.new_order = ethnicity.join(self.initial_order.set_index(0), how='left').index1

    def apply(self, method, i, j):
        pass

    def get_name(self):
        return self.name

    def get_n_indiv(self):
        return self.n_indiv

    def get_n_positions(self):
        if self.haploid0 is None:
            if self.haploid1 is None:
                return self.get_haploid0().shape[1]
            else:
                return self.haploid1.shape[1]
        else:
            return self.haploid0.shape[1]

    def get_haploid0(self):
        if self.haploid0 is None:
            self.haploid0 = self.extract_data()
        if self.new_order is None:
            self.extract_ethnicity_index()
            return self.haploid0[self.new_order]
        else:
            return self.haploid0

    def get_haploid1(self):
        if self.haploid1 is None:
            self.extract_haploid1()
        if self.new_order is None:
            self.extract_ethnicity_index()
            return self.haploid1[self.new_order]
        else:
            return self.haploid1

    def get_diploid(self):
        if self.diploid is None:
            self.build_diploid()
        if self.new_order is None:
            self.extract_ethnicity_index()
            return self.diploid[self.new_order]
        else:
            return self.diploid



class ComparisonEngine:
    """

    """
    def __init__(self, data, preprocessing, channel = 'hap0', metric='similarity',
                 option = '', ethnicity_reorder = True, verbose = True):
        self.genome_data = data
        if preprocessing is not None:
            self.preprocessing = preprocessing
        else:
            self.preprocessing = pp.Preprocessing()
        self.channel = channel
        self.option = option
        self.verbose = verbose
        self.ethnicity_reorder = ethnicity_reorder
        # retrieving the data
        self.get_data()
        # setting the metric
        self.metric_name = metric
        self.set_metric(self.metric_name)
        # setting cache_filename
        self.cache_filename = None

    def get_data(self):
        if self.channel == 'hap0':
            data = self.genome_data.get_haploid0()
        elif self.channel == 'hap1':
            data = self.genome_data.get_haploid1()
        elif self.channel == 'dip':
            data = self.genome_data.get_diploid()
        self.data = self.preprocessing.apply(data)

    def set_metric(self, metric):
        def hamming(x,y):
            return np.sum(np.abs(x-y))
        def similarity(x,y):
            d = hamming(x,y)
            return (len(x) + len(y) - d) * 1.0 / (len(x) + len(y) + d)
        if metric == 'hamming':
            self.metric = hamming
        if metric == 'similarity':
            self.metric = similarity

    def get_metric(self):
        return self.metric

    def init_filename(self):
        if not os.path.exists('Data/cache/'):
            os.mkdir('Data/cache/')
        self.cache_filename = 'Data/cache/{}_{}_{}{}{}{}.pkl'.format(self.genome_data.get_name(),
                                                                  self.channel,
                                                                  self.metric_name,
                                                                  self.preprocessing.get_name(),
                                                                  "_originalorder" if not self.ethnicity_reorder else "",
                                                                  ("_"+self.option) if self.option != "" else "")

    def load_from_cache(self):
        if self.cache_filename is None:
            self.init_filename()
        result = pickle.load(open(self.cache_filename, 'rb'))
        return result

    def dump_to_cache(self, data):
        if self.cache_filename is None:
            self.init_filename()
        pickle.dump(data, open(self.cache_filename, 'wb'))

    def compute(self):
        n = self.genome_data.get_n_indiv()
#         n = self.genome_data.n_indiv
        matrix = np.zeros((n, n))
        test = np.array([0])
        for i in range(n):
            matrix[i,i] = self.metric(test, test)

        def compute_row(i):
            res = []
            for j in range(i+1,n):
                res.append(self.metric(self.data[i,], self.data[j,]))
            return res

#         pool = Pool(8)
#         results = pool.map(compute_row, list(range(n-1)))
#         pool.close()
#         pool.join()

        results = []
        if self.verbose:
            for i in tqdm(range(n-1), desc="pairwise metric"):
                results.append(compute_row(i))
        else:
            for i in range(n-1):
                results.append(compute_row(i))

        for i in range(n-1):
            matrix[i+1:,i] = np.array(results[i])
            matrix[i,i+1:] = np.array(results[i])

        return matrix

    def compute_individual(self, i):
        pass

    def apply(self, force_recompute = False, force_dump = True):
        if not force_recompute:
            try:
                return self.load_from_cache()
            except:
                result = self.compute()
                if force_dump:
                    self.dump_to_cache(result)
                return result
        else:
            result = self.compute()
            if force_dump:
                self.dump_to_cache(result)
            return result

    def heatmap(self, **args):
        """
        Only on notebook ; requires graphic output
        """
        try:
            sns.heatmap(self.load_from_cache(), **args)
            plt.show()
        except:
            print("Exception occured when running heatmap.\nProbably no graphic handler, or you haven't computed the data yet.")
