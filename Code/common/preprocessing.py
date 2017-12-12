import numpy as np
import pandas as pd

class Preprocessing:
    def get_name(self):
        pass

    def apply(self, data):
        pass

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
        string
            Name of the instanciated RandomSlice object,
            based on its parameters.

        """
        return "rand-slice_{}{}".format(self.slice_param,
        "_"+str(self.start_param) if self.start_param is not None else "")

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
        n = data.get_n_positions()
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
    def __init__(self):
        pass

    def get_name(self):
        pass

    def apply(self, data):
        pass
    
