This afdust_xportfrac tool includes a sample script for creating the transportable fraction
file (XPORTFRAC) used to adjust afdust emissions after SMOKE is run.

This tool has been updated from original version that was provided with 2016 beta platform to support 
outputs from MCIP pre- and post-version 5.0. Specifically, MCIP version 4.5 and earlier generate 
GRIDCRO2D file that contains LUFRAC variable; MCIP version 5.0 and later writes LUFRAC variables to 
separate file LUFRAC_CRO, therefore both GRIDCRO2D and LUFRAC_CRO are required input files.

An sample script gen_afdust_tfrac.12US1.GRIDCRO2D.csh and is set up
for the 12US1 grid. This script can be run for any grid for which you have
a GRIDCRO2D file. You may use a GRIDCRO2D file for any day of the year.
The script calls a precompiled program, gen_afdust_tfrac.

In addition to the GRIDCRO2D (and LUFRAC_CRO), the program uses two ancillary files:
- gridcro2d_lu_types.12US2.csv: This provides descriptions for each land use
  variable in the GRIDCRO2D. **IMPORTANT** This must match your GRIDCRO2D file!
  You may verify this by looking in the NetCDF header of your GRIDCRO2D file,
  which should include descriptions for each LUFRAC variable. This particular
  file says "12US2", but 12US1 uses the same variables.
  An alternate file, gridcro2d_lu_types.USGS24.csv, is provided in this package
  for GRIDCRO2D files based on USGS land types (24 LUFRAC variables instead of 40).

- captureclass_fractions.txt: This specifies a transportable fraction for each
  land use type, as defined by the gridcro2d_lu_types file.
