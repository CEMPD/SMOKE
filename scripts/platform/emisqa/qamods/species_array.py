from __future__ import division
from __future__ import print_function
from builtins import object
from past.utils import old_div
import numpy as np
from qamods.chem_mechs import *

class SpeciesArray(object):
    """
    A class container for an array of a species.
    Takes a species name and initial array upon creation.

    obj.add_array(inarray): Adds another array of data to the current array
    obj.pct_err(inarray): Returns the percent error with another array as an array.
    obj.maxMin(): Returns the maximum and minimum values for the current array.
    obj: returns current species name when called as a string 
    obj(): returns current array
    """

    def __init__(self, init_array, species_name):
        self.species_name = species_name 
        self.array = init_array

    def __str__(self):
        return self.species_name

    def __call__(self):
        return self.array

    def add_array(self, inarray):
        """
        Adds another array of data to the current array.
        """
        self.array = self.array + inarray

    def diff_arr(self, inarray):
        """
        Gets the difference between the current array and another array.
        """
        if inarray.shape != self.array.shape:
            raise IndexError('Array size mismatch in percent error calculation')

        outarray = (inarray - self.array)
        outarray = np.where(outarray != outarray, 0.0, outarray)
        outarray = np.where(outarray > 1e10, 0.0, outarray)
        return outarray
    
    def pct_err(self, inarray):
        """
        Gets the percent error between the current array and another array.  Outputs to an array of
        the same size.
        """
        if inarray.shape != self.array.shape:
            raise IndexError('ERROR: Array size mismatch in percent error calculation')

        outarray = (old_div((inarray - self.array), self.array)) * float(100)
        outarray = np.where(outarray != outarray, 0.0, outarray)
        outarray = np.where(outarray > 1e10, 0.0, outarray)
        return outarray

    def maxMin(self):
        """
        Gives the maximum and minimum values in the array.
        """
        minVal = argmin(self.array.flat)
        maxVal = argmax(self.array.flat)
        return minVal, maxVal

    def moles2tons(self, informat, mech = 'cmaq_cb6'):
        """
        Converts a value or array of values from moles/s to tons/hr
        """
        mech_dct = molecDct[mech]

        if self.species_name in mech_dct: 
            factor = mech_dct[self.species_name]
        else: 
            print('WARNING: No match found for %s in mech table' %self.species_name)
            factor = 1

        self.array = self.array * float(factor)   # Convert moles per second to grams per second

        if informat != 'UAM':
            self.array = self.array * 3600.   # Convert tons per second to tons per hour

        self.array = self.array * (0.00000110231131)  # Convert grams per hour to tons per hour 

    def sum_dims(self):
        self.array = np.sum(self.array, axis=(1,2), keepdims=True)

