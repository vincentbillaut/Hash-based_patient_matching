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

# def buildDiploid(path0, path1, m=1092):
#     with gzip.open(path0, 'rb') as f0:
#         with gzip.open(path1, 'rb') as f1:
#             for l0 in f0:
#                 break
#             s0 = l0.decode('utf8')
#             n = len(s0.split())
#             data = np.zeros((m,n), dtype=int)
#             f0.seek(0)
#             for i,(l0,l1) in tqdm(enumerate(zip(f0,f1)), total=m):
#                 s0,s1 = l0.decode('utf8'), l1.decode('utf8')
#                 data[i,] = np.array(list(map(float, s0.split()))) + np.array(list(map(float, s1.split())))
#     print("done.")
#     return data


if __name__ == "__main__":
    # data extraction
    print("extracting data...")
    try:
        data0 = extractData(sys.argv[1], 1092)
        data1 = extractData(sys.argv[2], 1092)
    except Exception as inst:
        print("Error while opening data file.")
        print(inst)
        sys.exit(1)
    
    # data writing
    print("outputting diploid data...")
    newFileName = sys.argv[1].replace("hap0", "dip")
    np.savetxt(fname=newFileName, X=(data0+data1))
    
    print("done")