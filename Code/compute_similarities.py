import sys
import gzip
import re
import pandas as pd
import numpy as np
import distance
import random

from multiprocessing.dummy import Pool as ThreadPool
from tqdm import tqdm
from pprint import pprint


def file_len(path):
    with gzip.open(path, 'rb') as f:
        for i, _ in enumerate(f):
            pass
    return i + 1


def computeSimilarity(dataset, i, j):
    l = dataset.shape[1]
    d = distance.hamming(dataset[i], dataset[j])
    return (2*l - d) * 1.0 / (2*l + d)


def similarityWithFixed(dataset, threads = 8):
    n, l = dataset.shape[0], dataset.shape[1]
    i = random.randint(0,n-1)
#     i = 1
    ref = dataset[i]
    
    lst = list(range(n))

    pool = ThreadPool(threads)
    results = pool.map(lambda j: computeSimilarity(dataset, i, j), lst)
    pool.close()
    pool.join()

    return np.array(results), i


def extractData(path, m=None):
    if m is None:
        print("\tretrieving file size...")
        m = file_len(path)
    print("\textracting data...")
    with gzip.open(path, 'rb') as f:
        for l in f:
            break
        s = l.decode('utf8')
        n = len(s.split())
        data = np.zeros((m,n), dtype=int)
        f.seek(0)
        for i,l in tqdm(enumerate(f), total=m):
            s = l.decode('utf8')
            data[i,] = np.array(list(map(float, s.split())))
    print("\tdone.")
    return data


if __name__ == "__main__":
    # data extraction
    print("extracting data...")
    try:
        data = extractData(sys.argv[1], 1092)
    except Exception as inst:
        print("Error in opening data file")
        print(inst)
        sys.exit(1)

    # similarily
    print("computing similarities...")
    similarities, i = similarityWithFixed(data, 12)
    
    # data writing
    print("outputting similarities ({} reference)...".format(i))
    newFileName = sys.argv[1].split('.')[0] + "_sim_{}.gz".format(i)
    np.savetxt(fname=newFileName, X=similarities)
    
    print("done.")