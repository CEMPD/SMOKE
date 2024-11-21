#!/work/EMIS/python/miniforge3/envs/psemplot/bin/python

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import emisqa
from emisqa.formulas import calc_form
from emisqa.helpers import clean_temp
from emisqa.dataout.data_out import *
from emisqa.run_parse import *

# Get the command line options in an object 
opts = emisqa.run_parse.RunOpts()

# Get the output dictionary of species arrays
out_dict = emisqa.run_select.runQA(opts)

# Caculate any output species from formulas
out_dict = calc_form(out_dict, opts.formK, opts.formNK, opts.ignore_spec, opts.verbosity)

# Output the dictionary into a csv or ncf file
outfile = DataOut(opts.outfile_name, opts.out_type, opts.gsdate, opts.verbosity)
outfile.write_outfile(out_dict, opts.grid, opts.region, opts.tons, opts.units, opts.srgfile)

#Clean up temporary zip files
clean_temp(opts.zip_dict)
