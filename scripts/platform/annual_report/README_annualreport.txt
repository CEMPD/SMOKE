README annual_report.py


1. Options

usage: annual_report.py [options] parameterFile

options:
  -h, --help            show this help message and exit
  -s  sector_name       Sector name.  Overrides environmental variable SECTOR
  -c  case_abbrev       Case abbreviation.  Overrides environmental variable
                        CASE
  -y  year              Year.  Overrides environmental variable BASE_YEAR.
                        This is used for constructing the smkmerge dates
                        filename.
  -g  grid              Grid abbreviation.  Overrides environmental variable
                        GRID
  -f  output_format     Output format. The 'state' format is 1 state per row
                        and species as columns. The 'species' format is 1
                        species per row. Default is 'both', create both
                        formatted outputs.
  -m  mole_o_mass       Mass or mole. This is used to construct the SMOKE
                        report names.  The default values is 'mole'.
                        Alternatively, could be set to 'mass'
  -S  sector_list       Sector list file.  Overrides environmental variable
                        SECTORLIST
  -R  report_dir        Report directory.  Where the SMOKE reports reside.
                        Below the report_dir, the script expects sector
                        specific directories.
  -D  dates_dir         Smkmerge date files directory.  Where the smkmerge
                        dates file reside.
  -o  out_name          Output name for the annual summary reports.  For the
                        'state' format it will append a '.txt' and for the
                        'species' format it will append '_emf.csv'.  For
                        example, out_name = 'onroad_annual_2002' will become
                        onroad_annual_2002.txt and onroad_annual_2002_emf.csv.
                        Default value is 'annual_summary'.
  -O  out_dir           Output directory for the annual reports.
  --chem= chem_mech     Chemical mechanism.  Overrides chem_mechanism in the
                        parameter file
  --start_month= month  Start month. First month to aggregate over. Default is
                        1, January.
  --end_month= month    End month. Last month to aggregate over. Default is
                        12, December.
  --debug               Turn debugging on.
  --debug_extra         Turn debugging on. Prints which smkmerge reports were
                        processed.

2. Example usage

     $ $script_dir/annual_report.py \
     -R ~/emf/smoke_optim/reports/smkmerge \
     -D ~/emf/EPA_scripts/inputs/2002ae_02b/mrggrid -O ~/tmp \
     -o {$SECTOR}_annual_summary_$CASE $script_dir/parameter_file.txt

   The script will look for the smkmerge daily reports under a sector
   directory below "~/emf/smoke_optim/reports/smkmerge", and the
   mrggrid dates files in
   "~/emf/EPA_scripts/inputs/2002ae_02b/mrggrid".

   It will write the annual summary report(s) to "~/tmp" directory and
   will create two outputs: $SECTOR_annual_summary_$CASE.txt and
   $SECTOR_annual_summary_$CASE_emf.csv


   Here, script_dir has been set to the location of the
   annual_report.py and the parameter_file.txt.  In addition, the
   following environmental variables were set: CASE, BASE_YEAR, GRID,
   SECTOR, SECTORLIST.  Note, that the environmental variables do not
   need to be set, their values can be passed as options.

   example of running without environmental variables set:
      $ ./annual_report.py -s othpt -c 2002ae_02b  -g 36US1 \
      -y 2002 -S ~/emf/EPA_scripts/inputs/2002ae_02b/mrggrid/work/sectorlist_2002ae_02b.txt \
      -R ~/emf/smoke_optim/reports/smkmerge \
      -D ~/emf/EPA_scripts/inputs/2002ae_02b/mrggrid \
      -o {$SECTOR}_annual_summary_$CASE -O ~/tmp/2002ae parameter_file.txt       

3. Output format

   The default is to create two outputs, state and species.  The
   'state' format is the original format from CSC's python scripts.
   The only change is that it has aggregate species included as
   additional columns, e.g. "NOX".

   The 'species' format is in CSV format and has 4 columns:
      State,Sector,Species,annu_emis

   NOTE: If you run the annual_report.py multiple times with the same
   output name, it will create new 'state' reports but will
   continuously append to the 'species' report.

4. Parameter file

   This contains the default chemical mechanism, eg. 'cmaq_cb05'.  It
   also contains 2 dictionaries.  The first dictionary is the
   molecular weights dictionary for this chemical mechanism.  The
   second is an aggregate species dictionary that defines the
   composition of the aggregate species. If you add a new aggregate
   species, it must also be added to the molecular weights
   dictionary. Note: the aggregate species dictionary elements are
   sets.


