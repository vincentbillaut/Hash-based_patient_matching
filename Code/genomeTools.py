import numpy as np
import pandas as pd
import gzip
import pickle
import os
import distance

from tqdm import tqdm


def file_len(path):
    with gzip.open(path, 'rb') as f:
        for i, _ in enumerate(f):
            pass
    return i + 1


class GenomeData:
    def __init__(self, name="chrom22", path="Data/22hap0.gz", verbose = True):
        self.name = name
        self.source_path = path
        self.verbose = verbose
        if self.verbose:
            print("Retrieving file length...")
        self.n_indiv = file_len(self.source_path)
        self.haploid0 = self.extract_data()
        self.haploid1, self.diploid = None, None

    def extract_data(self, path=None, verbose=None):
        if path is None:
            path=self.source_path
        if verbose is None:
            verbose=self.verbose
        if verbose:
            print("Extracting data from {}...".format(path))
        with gzip.open(path, 'rb') as f:
            for l in f:
                break
            s = l.decode('utf8')
            n = len(s.split())
            data = np.zeros((self.n_indiv,n), dtype=int)
            f.seek(0)
            for i,l in tqdm(enumerate(f), total=self.n_indiv):
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
        self.diploid = self.haploid0 + self.haploid1
    
    def extract_ethnicity(self, path='Data/integrated_call_samples.20101123.ALL.panel'):
        data = pd.DataFrame(columns=['patient', 'country', 'continent'])
        with open(path, 'r') as f:
            for i,l in enumerate(f):
                info = l.strip().split()
                data.loc[i] = info[:3]
        self.ethnicity = data
        
    def apply(self, method, i, j):
        pass
        
    def get_name(self):
        return self.name
        
    def get_haploid0(self):
        return self.haploid0
    
    def get_haploid1(self):
        if self.haploid1 is None:
            self.extract_haploid1()
        return self.haploid1
    
    def get_diploid(self):
        if self.diploid is None:
            self.build_diploid()
        return self.diploid
    
    def get_n_indiv(self):
        return self.n_indiv
    
    
    
class ComparisonEngine:
    def __init__(self, data=None, channel = 'hap0', window_size = 1, metric='similarity', option = ''):
        self.genome_data = data
        self.channel = channel
        if self.genome_data is not None:
            self.get_data()
        self.window_size = window_size
        self.set_metric(metric)
        self.option = option
        self.cache_filename = None
        
    def get_data(self):
        if self.channel == 'hap0':
            self.data = self.genome_data.getHaploid0()
        elif self.channel == 'hap1':
            self.data = self.genome_data.getHaploid1()
        elif self.channel == 'dip':
            self.data = self.genome_data.getDiploid()
            
    def set_metric(self, metric):
        def similarity(x,y):
            d = distance.hamming(x,y)
            return (len(x) + len(y) - d) * 1.0 / (len(x) + len(y) + d)
        if metric == 'hamming':
            self.metric = lambda x,y: distance.hamming(x,y,normalized=True)
        if metric == 'similarity':
            self.metric = similarity
            
    def get_metric(self):
        return self.metric
            
    def init_filename(self):
        if not os.path.exists('data/cache/'):
            os.mkdir('data/cache/')
        self.cache_filename = 'data/cache/{}_{}_{}_{}.pkl'.format(self.genome_data.get_name(),
                                                                  self.metric,
                                                                  self.window_size,
                                                                  self.option)
        
    def load_from_cache(self):
        if self.cache_filename is None:
            self.init_filename()
        result = pickle.load(open(self.cache_filename, 'rb'))
        return result
    
    def dump_to_cache(self, data):
        if self.cache_filename is None:
            self.init_filename()
        pickle.dump(result, open(self.cache_filename, 'wb'))
    
    def compute(self):
        n = self.genome_data.get_n_indiv()
        matrix = np.zeros((n, n))
        test = [0]
        for i in range(n):
            matrix[i,i] = self.metric(test, test)
        
        def compute_row(i):
            res = []
            for j in range(i+1,n):
                res.append()
            return res
        
        
        
    
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