from __future__ import print_function
from builtins import range
from builtins import object
import sys
import netCDF4 as ncf

class StkGrp(object):
    """
    Develop the stack group col/row x-ref information
    """

    def __init__(self, infile_name, verbosity=False):
        """
        """
        self.infile_name = infile_name
        self.verbosity = verbosity
        self.infile = self._load_infile()
        self.pt_xref = {}
        self.stk_num = 0
        self._get_xref()

    def _open_NCF(self):
        '''
        Opens the netCDF input file and returns an open file object.
        '''
        try: 
            infile = ncf.Dataset(self.infile_name)
        except TypeError: 
            sys.exit('ERROR: %s not available for access or not a netCDF.' %self.infile_name)
        else: 
            return infile 

    def _load_infile(self):
        """
        Set the infile name based on the SMOKE conventions.
        """
        if self.verbosity:
            print("Stack groups: " + self.infile_name)
        return self._open_NCF()

    def _get_xref(self):
        """
        Process the col and row to create a x ref to stack
        """
        # Fetch a variable from the in file
        row_spec = self.infile.variables['ROW']
        row_in = row_spec[:]
        col_spec = self.infile.variables['COL']
        col_in = col_spec[:]
        self.stk_num = row_spec.shape[2]
        for stack in range(row_spec.shape[2]):
            row = int(row_in[0][0][stack][0])
            col = int(col_in[0][0][stack][0])
            self.pt_xref[stack] = {'col': col, 'row': row}

