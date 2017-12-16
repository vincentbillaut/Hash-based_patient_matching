import numpy as np
import pandas as pd
import hashlib

from tqdm import tqdm


# verbose function
def tqdm_v(gen, verbose=True, **kwargs):
    if verbose:
        return tqdm(gen, **kwargs)
    else:
        return gen


class Preprocessing:
    def __init__(self):
        pass

    def get_name(self):
        return ""

    def apply(self, data):
        return data


class RandomSlice(Preprocessing):
    """
    Preprocessing object that extracts a random slice from the data.
    """
    def __init__(self, slice_size = -1, slice_start = None):
        self.slice_param = slice_size
        self.start_param = slice_start

    def get_name(self):
        """Returns the name of the object.

        Returns
        -------
        str
            Name of the instanciated RandomSlice object,
            based on its parameters.

        """
        return "_rand-slice-{}{}".format(self.slice_param,
        ("-"+str(self.start_param)) if self.start_param is not None else "")

    def apply(self, data):
        """Applies slicing to given data.

        Parameters
        ----------
        data : numpy.ndarray
            The data to apply the slicing on.

        Returns
        -------
        numpy.ndarray
            The data, sliced accordingly.

        """
        n = data.shape[1]
        # slice_size
        if self.slice_param == -1 or self.slice_param > n:
            slice_size = n
        elif type(self.slice_param) == float:
            slice_size = int(self.slice_param * n)
        else:
            slice_size = self.slice_param
        # slice_start
        if self.start_param is None:
            slice_start = np.random.randint(0,n - slice_size + 1)
        else:
            slice_start = self.start_param
        # slicing
        return data[:, slice_start:(slice_start + slice_size)]


class HashingWindow(Preprocessing):
    def __init__(self, granularity=1, hashing=None):
        """Initializes HashingWindow object.

        Parameters
        ----------
        granularity : int
            Size of window to consider as atomic element of a sequence.
            Default, 1, means that we take the data as is.
            10 means that we slice the data into slices of width 10, and work on hashes of those length-10 slices.
        hashing : type
            Hashing method. Must be a hashlib method.

        """
        self.granularity = granularity
        self.hasher = hashing if hashing is not None else hashlib.md5

    def get_name(self):
        """Returns the name of the HashingWindow object.

        Returns
        -------
        str
            Name of the object, taking into account granularity and hashing method.

        """
        return "hash-{}-{}".format(self.hasher().name, self.granularity)

    def apply(self, data, verbose=False):
        # slicing data into evenly sized chunks
        def chunk_sizes():
            l = list()
            curr = 0
            for i in range(0, data.shape[1], self.granularity):
                curr += min(self.granularity, data.shape[1]-i)
                l.append(curr)
            return np.array(l[:-1])

        chunks = chunk_sizes()
        result = np.zeros((data.shape[0], len(chunks)+1))
        for j,subarray in tqdm_v(enumerate(np.split(data, chunks, axis=1)),
                                    verbose, desc="hashing"):
            for i in range(subarray.shape[0]):
                result[i,j] = int(self.hasher((''.join(map(str,subarray[i,:]))).encode('UTF-8')).hexdigest(),16)
        return result
