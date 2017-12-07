import sys
import genomeTools as gt

from tqdm import tqdm

w = int(sys.argv[1])
n = int(sys.argv[2])

print("retrieve data")
genomedata = gt.GenomeData()

for i in tqdm(range(n)):
    ce = gt.ComparisonEngine(genomedata, window_size=w, verbose=False)
    ce.apply()
    
print("done")