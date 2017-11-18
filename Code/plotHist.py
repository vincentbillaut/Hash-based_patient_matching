import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import sys

data = np.genfromtxt(sys.argv[1])
sns.hist(data, bins=200)
plt.title("Similarity Distribution relative to sequence #121")
plt.show()