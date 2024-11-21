from __future__ import print_function
from qamods.helpers import checkEV, conv2greg
from qamods.default_paths import *
from optparse import OptionParser,OptionGroup
from qamods.data_file import DataFile
from qamods.dateloop import inday as dl
from qamods.inline.stack_group import StkGrp
import os.path
from sys import exit

### Main

def listRunTypes():
    '''
    Lists the available run types and syntax.
    '''
    print('''
    \n\t\tRun types for the pyQA script.\n\t\t      ------------  \n
    pe: percent error run between two files and the specified species.  Requires two files listed with the -f option.
    rd: raw difference run between two files and the specified species.  Requires two files listed with the -f option.
    add: adds together two files and the specified species. Returns CSV.  Requires two files listed with the -f option.
    dv: dumps the daily value or timestep sum of a specified file.  Requires one file listed with the -f option.
    sum: sums up a group of files. Requires case, sector, speciation, input path, grid, GSDATE, and rundays.
    avg: averages a group of files over the number of days. Requires case, sector, speciation, input path, grid, GSDATE, and rundays.
    domain: writes the daily domain of a file.  Returns as one line CSV.  Requires one file listed with the -f option.
    mm: gives the maximum and minimum species values for a file.  Requires one file listed with the -f option.
    dump: raw dump of a gridded file.  Requires one file with the -f option.
    dh: Dump hourly data from a gridded file. Useful for type conversions.\n
    ''')
    exit(1)

def get_opts():
    ### Test for environment variables and fill in the variables locally.
    grid = checkEV('GRID')
    gsdate = checkEV('GSDATE')
    case = checkEV('CASE')
    sector = checkEV('SECTOR')
    inPath = checkEV('IMD_ROOT')
    spec = checkEV('SPEC')
    grid_desc = checkEV('GRIDDESC')

    # Handle command line arguments and options.
    parser = OptionParser(usage = 'usage: %prog [options] OUTFILE RUNTYPE')

    outputGroup = OptionGroup(parser, "Output File Configuration Options")
    loopGroup = OptionGroup(parser, "Options Required for sum and avg Methods")
    regionGroup = OptionGroup(parser, "Options Required for Region Output","Output type must be set to CSV or left default")
    parser.add_option('-s', '--speciesname', dest='species_name', help='List of species to process.  Ex. -s SPEC1,SPEC2,...', default='')
    parser.add_option('-a', '--allspecies', action='store_true', dest='all_species', help='Run for all species.', default=False)
    parser.add_option('-l', '--listruns', action='store_true', dest='listRuns', help='List run types.', default=False) 
    parser.add_option('-f', '--file', dest='file_name', help='Specify a file or a list of up to two files for access.', metavar='FILE', default='')
    parser.add_option('-v', '--verbose', action='store_true', dest='verbosity', help='Give more information during run.', default=False)
    parser.add_option('-i', '--inline', action='store_true', dest='inln', help='Use inline input files with stack groups rather than 2D', default=False)
    parser.add_option('-g', '--stack-groups', dest='stackFile', help='Explicitly set stack groups file location.  Set to default path when looping with -i option.', default = '')
    parser.add_option('-c', '--informat', dest='informat', help='Input file type. Specify NCF, CSV, or UAM.', default='NCF', metavar='INTYPE')
    parser.add_option('--camx-ptsr', action='store_true', dest='ptsr', help='Use PTSOURCE CAMx emissions format as input.', default=False)
    parser.add_option('--layer', dest='layer', help='Specify an individual layer.  Defaults to flatten all layers.', metavar='#', default='')
    parser.add_option('--ignore-missing', action='store_true', dest='ignore_spec', help='Ignore missing input species instead of producing fatal error.  May produce inaccurate results.', default=False)
    parser.add_option('--interpolate', action='store_true', dest='interpolate', help='Average between hours to get interpolated results for CAMx comparison.', default=False)
    regionGroup.add_option('-e', '--region', dest='region', help='Specify CSV output type as county or state.  Defaults to gridded.\nRegion requires -G [grid] command line option.', metavar='region', default='')
    regionGroup.add_option('--srg_file', dest='srgfile', help='Specify path to surrogate file used for getting regions. Defaults to 12km land area.', default=defSrgFile)
    loopGroup.add_option('-r', '--repdays', dest='rep_days', help='Use representative days of specified type: aveday_N, aveday_Y, mwdss_N, mwdss_Y, week_N, week_Y, all\n  Defaults to all', metavar='TYPE', default='')
    loopGroup.add_option('-G', '--grid', dest='grid', help='Set the grid name when looping or when writing state/county based reports.', metavar='GRID', default='')
    loopGroup.add_option('-D', '--gsdate', dest='gsdate', help='Set the starting GSDATE when looping.', metavar='YYYYMMDD', default='')
    loopGroup.add_option('-R', '--rundays', dest='run_days', help='Number of days from the starting GSDATE to run for looping.', metavar='1-366', default='')
    loopGroup.add_option('-C', '--case', dest='case', help='Name of case for the input data when looping.', metavar='CASE', default='')
    loopGroup.add_option('-S', '--sector', dest='sector', help='Sector of the input data when looping.', metavar='SECTOR', default='')
    loopGroup.add_option('-P', '--speciation', dest='spec', help='Speciation of the input data for mole->ton conversion and filenames when looping.', metavar='SPEC', default='cmaq_cb6')
    loopGroup.add_option('-I', '--inputpath', dest='inPath', help='Path to the input data.  Used when looping.', metavar='PATH', default='')
    outputGroup.add_option('-o', '--outtype', dest='out_type', help='Output to NCF or CSV.  Default is CSV.', default='CSV')
    outputGroup.add_option('-t', '--tons', action='store_true', dest='tons', help='Convert units from moles/sec to tons/hr.  Defaults to no conversion.', default=False)
    outputGroup.add_option('-F', '--formula', dest='formula', help='Adds existing species into an output species. Format: CALCSPEC1=SPECa+SPECb,CALCSPEC2=SPECc+SPECd,...', default='')
    outputGroup.add_option('-K', '--formula-no-keep', dest='formulaNK', help='Same as -F, except that all species used in a calculation will not be in output.', default='')
    outputGroup.add_option('-u', '--units', dest='units', help='Override units name.  Default is moles/s or tons/day if -t flag is used.', default='')
    outputGroup.add_option('--all-hours', action='store_true', dest='all_hours', help='Sum up species over all time steps/hours in a file.  Default is to sum first 24 hours.', default=False)
    parser.add_option('--griddesc', dest='grid_desc', help='Specify the path to the grid description file', metavar='#', default='')
    parser.add_option_group(regionGroup)
    parser.add_option_group(loopGroup)
    parser.add_option_group(outputGroup)

    (options, args) = parser.parse_args()

    if options.listRuns: 
        listRunTypes() 

    if len(args) != 2: 
        print('Must specify an outfile and a run type.')
        print('Use the -l option to list run types or -h for futher options.')
        exit()

    # Override environment variables
    if options.grid: 
        grid = options.grid
    if options.gsdate: 
        gsdate = options.gsdate
    if options.case: 
        case = options.case
    if options.sector: 
        sector = options.sector
    if options.inPath: 
        inPath = options.inPath
    if options.spec: 
        spec = options.spec
    if options.run_days: 
        run_days = int(options.run_days)
    else:
        run_days = ''
    if options.grid_desc: 
        grid_desc = options.grid_desc
    else:
        grid_desc = defGridDesc
    # Set command line options
    layer = options.layer
    verbosity = options.verbosity
    tons = options.tons
    region = options.region.strip().lower()
    rep_days = options.rep_days.strip().upper()
#    formK = options.formula.strip().upper()
#    formNK = options.formulaNK.strip().upper()
    formK = options.formula.strip()
    formNK = options.formulaNK.strip()
    inln = options.inln
    ignore_spec = options.ignore_spec
    units = options.units
    stackFile = options.stackFile
    all_hours = options.all_hours
    interpolate = options.interpolate
    species_name = options.species_name
    all_species = options.all_species
    file_name = options.file_name
    informat = options.informat.upper()
    ptsr = options.ptsr

    #  Load the molecular weight conversion dictionary from the external parameter file.
    if spec in mwDict:
        mwSpec = spec.lower()
    else:
        mwSpec = mwDefault
    if verbosity: 
        print('Loading molecular weight file for speciation: %s' %mwSpec)
#    execfile(mwDict[mwSpec])

    outfile_name = args[0]
    run_type = args[1].strip().lower()
    out_type = options.out_type.upper()
    zip_dict = {} # Dictionary pointing to unzipped input files

    if out_type == 'NCF':
        if not grid or not grid_desc:
            parser.error('Must specify a grid name and a grid description file when outputting to NCF.')

    # Set the input file name prefix for inline versus 2D
    if inln:
        inprefix = 'inln'

        if not grid:
            parser.error('Inline to 2D conversion requires setting a grid.  Please specify a grid with -G')
        
        if not stackFile:
            if inPath and sector and grid and case:
                stackFile = os.path.join(inPath, 'stack_groups_%s_%s_%s.ncf' %(sector, grid, case))
            else:
                parser.error('Stack groups file not set.  Please set path using -g.')

        stacks = StkGrp(stackFile)
    else:
        inprefix = 'emis'
        stacks = ''

    # Handle species list options
    if not species_name and not all_species: 
        parser.error('No species specified.  Must either specify the -s or -a option.')
    if species_name and all_species: 
        parser.error('You must only specify either the -s or the -a option.')
#    species_list = species_name.strip().upper().split(',')
    species_list = species_name.strip().split(',')

    # Load first file for obtaining species list
    if len(file_name.split(',')) > 1: 
        speciesfile_name = file_name.split(',')[0]
    else: 
        speciesfile_name = file_name

    if speciesfile_name:
        infile = DataFile(speciesfile_name, verbosity, informat, ptsr, zip_dict)

        if not gsdate.strip():
            esdate = infile.sdate
            if int(esdate) < 190000:
                esdate = 2011001
            gsdate = conv2greg(esdate)
    else:
        if not grid or not gsdate or not case or not sector or not inPath or not spec or not run_days:
            parser.error('Must set an input path, case, sector, gsdate, grid, rundays, and speciation OR an input file name (-f).')

        # Handle representative days    
        if run_days:
            sdate = dl.InDay(gsdate, rep_days, run_days, smkDatesPath) 
        else:
            sdate = gsdate

        if sector.lower() == 'mrggrid':
            infile_name = os.path.join(inPath, 'emis_mole_all_%s_%s_%s_%s.ncf' %(sdate, grid, spec, case))
        elif sector.lower() == 'mrggrid_withbeis':
            infile_name = os.path.join(inPath, 'emis_mole_all_%s_%s_withbeis_%s.ncf' %(sdate, grid, case))
        elif sector.lower() == 'mrggrid_nobeis':
            infile_name = os.path.join(inPath, 'emis_mole_all_%s_%s_nobeis_%s.ncf' %(sdate, grid, case))
        else:
            infile_name = os.path.join(inPath, sector, '%s_mole_%s_%s_%s_%s_%s.ncf' %(inprefix, sector, sdate, grid, spec, case))  # v5 directory structure 

        infile = DataFile(infile_name, verbosity, informat, ptsr, zip_dict)

    # Get species list from open file
    if all_species:
        species_list = list( species_name for species_name in infile.species_list if species_name != 'TFLAG' )
#    infile.closeFile()

    if verbosity and layer: 
        print('Running for layer %s' %layer)

    opt_dict = {'species_list': species_list, 'file_name': file_name, 'verbosity': verbosity, 'tons': tons, 'region': region, 
        'rep_days': rep_days, 'formK': formK, 'formNK': formNK, 'inln': inln, 'stacks': stacks, 'units': units, 
        'ignore_spec': ignore_spec, 'all_hours': all_hours, 'interpolate': interpolate, 'informat': informat, 'ptsr': ptsr, 
        'grid': grid, 'gsdate': gsdate, 'case': case, 'sector': sector, 'inPath': inPath, 'spec': spec, 'run_days': run_days, 
        'layer': layer, 'grid_desc': grid_desc, 'out_type': out_type, 'run_type': run_type, 'zip_dict': zip_dict, 'outfile_name': outfile_name,
        'srgfile': options.srgfile}

    return opt_dict

