###
# Path and file defaults
# Path and filename of molecular weight conversion files
mwDict = {'cmaq_cb05_soa': '/work/EMIS/users/bte/apps/old/pyQA/parameter_file_cmaq_cb05.txt', 'cmaq_cb05': '/work/EMIS/users/bte/apps/old/pyQA/parameter_file_cmaq_cb05.txt', 'cmaq_cb05_soa': '/work/EMIS/users/bte/apps/old/pyQA/parameter_file_cmaq_cb05.txt', 'saprc07t': '/work/EMIS/users/bte/apps/old/pyQA/parameter_file_saprc07t.txt', 'cmaq_cb6': '/work/EMIS/users/bte/apps/old/pyQA/parameter_file_cmaq_cb6.txt'} 
mwDefault = 'cmaq_cb6'  #  Default molecular weight file 

#Path to the SMOKE dates files
smkDatesPath = '/work/EMIS/em_v6.3/ge_dat/smk_dates/'

#Path to temporary directory
tmp_dir = '/work/EMIS/users/bte/tmp'

# Set the default grid description file.  This can be overridden with a GRIDDESC environment variable.
defGridDesc = '/work/EMIS/em_v7/ge_dat/gridding/griddesc_lambertonly_26may2017_nf_v78.txt'
#defGridDesc = '/work/EMIS/em_v6.3/ge_dat/gridding/griddesc_hemispheric_17mar2016_v0.txt'
#defGridDesc = '/home/callen05/tmp/griddesc_hemispheric_17mar2016_v0.txt'

# Default surrogate gridding file
defSrgFile = '/work/EMIS/em_v6.3/ge_dat/gridding/surrogates/New_Surrogates/CONUS12_2010_v5_20141015/USA_340_NOFILL.txt'
