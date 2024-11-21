#!/usr/bin/env python

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import qamods
from qamods.formulas import calc_form
from qamods.helpers import clean_temp
from qamods.dataout.data_out import *
from qamods.run_parse import *

# Get the command line options in a dictionary
opt_dict = qamods.run_parse.get_opts()

# Get the output dictionary of species arrays
out_dict = qamods.run_select.runQA(opt_dict['run_type'], opt_dict['species_list'], opt_dict)

# Caculate any output species from formulas
out_dict = calc_form(out_dict, opt_dict['formK'], opt_dict['formNK'], opt_dict['ignore_spec'], opt_dict['verbosity'])

# Output the dictionary into a csv or ncf file
outfile = data_out(opt_dict['outfile_name'], opt_dict['out_type'], opt_dict['gsdate'], opt_dict['verbosity'])
outfile.write_outfile(out_dict, opt_dict['grid'], opt_dict['grid_desc'], opt_dict['region'], 
  opt_dict['tons'], opt_dict['units'], opt_dict['srgfile'])

#Clean up temporary zip files
clean_temp(opt_dict['zip_dict'])
