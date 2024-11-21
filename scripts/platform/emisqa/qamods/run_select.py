from qamods.runtypes import *

def runQA(run_type, species_list, run_opts):
    # Select run type and run the script.
    # Run types dictionary contains x-reference from command line run type name to run type function name 
    run_types = {'pe': percentError, 'add': addFiles, 'dv': dumpDV, 'sum': sumDV, 'avg': avgDV, 'domain': singleDomain, 'mm': mMDomain, 'rd': rawDiff, 'dump': rawDump, 'hd': hourDump, 'yd': hrlyDomain}
    if run_type not in run_types: 
        raise ValueError('Specified run type not available.  Please refer to the list of run types using the -l argument.')

    out_dict = run_types[run_type](species_list, run_opts)

    # Convert to tons by species if conversion is turned on
    if run_opts['tons']: 
        for speciesName in list(out_dict.keys()):
            out_dict[speciesName].moles2tons(run_opts['informat'], run_opts['spec'])

    return out_dict

'''
Individual run types by function name
Extending:
Function arguments must always be the "species_list" list and the "run_opts" dictionary
Always returns dictionary of 3d species arrays [hrs,cols,rows]
'''

def percentError(species_list, run_opts):
    '''
    Calculate the percentError between two NCF files.
    '''
    out_dict = pe.get_dict(species_list, run_opts['file_name'], run_opts['all_hours'], run_opts['grid'], run_opts['grid_desc'], run_opts['ignore_spec'], run_opts['inln'], run_opts['interpolate'], 
        run_opts['layer'], run_opts['region'], run_opts['stacks'], run_opts['ptsr'], run_opts['informat'], run_opts['verbosity'], run_opts['zip_dict'])
    return out_dict

def rawDiff(species_list, run_opts):
    '''
    Calculate the raw difference between two NCF files.
    '''
    out_dict = raw_diff.get_dict(species_list, run_opts['file_name'], run_opts['all_hours'], run_opts['grid'], run_opts['grid_desc'], run_opts['ignore_spec'], run_opts['inln'], run_opts['interpolate'], 
        run_opts['layer'], run_opts['region'], run_opts['stacks'], run_opts['ptsr'], run_opts['informat'], run_opts['verbosity'], run_opts['zip_dict'])
    return out_dict

def addFiles(species_list, run_opts):
    '''
    Adds together two NCF files
    '''
    out_dict = add_files.get_dict(species_list, run_opts['file_name'], run_opts['all_hours'], run_opts['grid'], run_opts['grid_desc'], run_opts['ignore_spec'], run_opts['inln'], run_opts['interpolate'], 
        run_opts['layer'], run_opts['region'], run_opts['stacks'], run_opts['ptsr'], run_opts['informat'], run_opts['verbosity'], run_opts['zip_dict'])
    return out_dict

def sumDV(species_list, run_opts):
    '''
    Sums the daily values of NCF files from a start date through the number of run dates.
    '''
    out_dict = sumdv.get_dict(species_list, run_opts['all_hours'], run_opts['grid'], run_opts['gsdate'], run_opts['case'], run_opts['sector'], run_opts['inPath'], run_opts['spec'], run_opts['rep_days'], run_opts['run_days'], 
        run_opts['grid_desc'], run_opts['ignore_spec'], run_opts['inln'], run_opts['interpolate'], run_opts['layer'], run_opts['region'], run_opts['stacks'], run_opts['ptsr'], run_opts['informat'], 
        run_opts['verbosity'], run_opts['zip_dict'])
    return out_dict

def avgDV(species_list, run_opts):
    '''
    Averages the daily values of NCF files from a start date through the number of run dates by the number of run dates.
    '''
    out_dict = avgdv.get_dict(species_list, run_opts['all_hours'], run_opts['grid'], run_opts['gsdate'], run_opts['case'], run_opts['sector'], run_opts['inPath'], run_opts['spec'], run_opts['run_days'], 
        run_opts['grid_desc'], run_opts['ignore_spec'], run_opts['inln'], run_opts['interpolate'], run_opts['layer'], run_opts['region'], run_opts['stacks'], run_opts['ptsr'], run_opts['informat'], 
         run_opts['verbosity'], run_opts['zip_dict'])
    return out_dict

def dumpDV(species_list, run_opts):
    '''
        Dumps the daily value data of a single NCF file to a file.
        '''
    out_dict = dump_dv.get_dict(species_list, run_opts['file_name'], run_opts['all_hours'], run_opts['grid'], run_opts['grid_desc'], run_opts['ignore_spec'], run_opts['inln'], run_opts['interpolate'], 
        run_opts['layer'], run_opts['region'], run_opts['stacks'], run_opts['ptsr'], run_opts['informat'], run_opts['verbosity'], run_opts['zip_dict'])
    return out_dict

def singleDomain(species_list, run_opts):
    '''
    Sums up every grid cell for every hour for each species.
    '''
    out_dict = single_domain.get_dict(species_list, run_opts['file_name'], run_opts['all_hours'], run_opts['grid'], run_opts['grid_desc'], run_opts['ignore_spec'], run_opts['inln'], run_opts['interpolate'], 
        run_opts['layer'], run_opts['region'], run_opts['stacks'], run_opts['ptsr'], run_opts['informat'], run_opts['verbosity'], run_opts['zip_dict'])
    return out_dict

def hrlyDomain(species_list, run_opts):
    '''
    Sums up every grid cell for every hour by hour for each species.
    '''
    out_dict = hourly_domain.get_dict(species_list, run_opts['file_name'], run_opts['all_hours'], run_opts['grid'], run_opts['grid_desc'], run_opts['ignore_spec'], run_opts['inln'], run_opts['interpolate'], 
        run_opts['layer'], run_opts['region'], run_opts['stacks'], run_opts['ptsr'], run_opts['informat'], run_opts['verbosity'], run_opts['zip_dict'])
    return out_dict

def mMDomain(species_list, run_opts):
    '''
    Sums up every grid cell for a day for each species.
    Returns the min and max
    '''
    out_dict = mm_domain.get_dict(species_list, run_opts['file_name'], run_opts['all_hours'], run_opts['grid'], run_opts['grid_desc'], run_opts['ignore_spec'], run_opts['inln'], run_opts['interpolate'], 
        run_opts['layer'], run_opts['region'], run_opts['stacks'], run_opts['ptsr'], run_opts['informat'], run_opts['verbosity'], run_opts['zip_dict'])
    return out_dict

def rawDump(species_list, run_opts):
    '''
    Does a raw dump of a single NCF file
    '''
    out_dict = raw_dump.get_dict(species_list, run_opts['file_name'], run_opts['all_hours'], run_opts['grid'], run_opts['grid_desc'], run_opts['ignore_spec'], run_opts['inln'], run_opts['interpolate'], 
        run_opts['layer'], run_opts['region'], run_opts['stacks'], run_opts['ptsr'], run_opts['informat'], run_opts['verbosity'], run_opts['zip_dict'])
    return out_dict

def hourDump(species_list, run_opts):
    '''
    Does a raw dump of a sum of all hours in a single NCF file to a file.
    '''
    out_dict = hour_dump.get_dict(species_list, run_opts['file_name'], run_opts['all_hours'], run_opts['grid'], run_opts['grid_desc'], run_opts['ignore_spec'], run_opts['inln'], run_opts['interpolate'], 
        run_opts['layer'], run_opts['region'], run_opts['stacks'], run_opts['ptsr'], run_opts['informat'], run_opts['verbosity'], run_opts['zip_dict'])
    return out_dict

