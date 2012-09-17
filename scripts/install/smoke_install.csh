#!/bin/csh -f

set exitstat = 0
set here = `pwd`
set tmpfile = .install

/bin/rm -rf $tmpfile

# Check that downloaded files are available
set file = smoke_v31.nctox.data.tar.gz
echo "Checking for $file file..."

ls $file > $tmpfile
if ( $status > 0 ) then
    echo "ERROR: Could not find $file"
    echo "       Please download the file $file and try again"
    set exitstat = 1
endif

set file = smoke_v31.Linux2_x86_64pg.tar.gz
echo "Checking for $file file..."

ls $file > $tmpfile
if ( $status > 0 ) then
    echo "ERROR: Could not find $file"
    echo "       Please download the file $file and try again"
    set exitstat = 1
endif

if ( $exitstat > 0 ) then
   exit( $exitstat )
endif

# Check that SMK_HOME is set
if ( ! $?SMK_HOME ) then
    echo "ERROR: You must define the SMK_HOME environment variable"
    echo "       before running this script.  Use the following command"
    echo "       to set your SMK_HOME variable:"
    echo "          setenv SMK_HOME <your chosen SMOKE installation dir>"
    exit( 1 )
endif

echo "SMOKE v3.1 will be installed in the following directory:"
echo "      $SMK_HOME"
echo " "

# Install files
cd $SMK_HOME

set file = smoke_v31.nctox.data.tar.gz
echo "Unpacking file $file..."
tar xzvf $here/$file
if ( $status > 0 ) then
   echo "ERROR: Could not unpack the file $file"
   echo "       Confirm that the 'tar' command is available and"
   echo "       that you have the correct permissions"
   set exitstat = 1
endif

set file = smoke_v31.Linux2_x86_64pg.tar.gz
echo "Unpacking file $file..."
tar xzvf $here/$file
if ( $status > 0 ) then
   echo "ERROR: Could not unpack the file $file"
   echo "       Confirm that the 'tar' command is available and"
   echo "       that you have the correct permissions"
   set exitstat = 1
endif

if ( $exitstat > 0 ) then
   exit( $exitstat )
endif

# Source the assigns file
cd subsys/smoke/assigns
set file = ASSIGNS.nctox.cmaq.cb05_soa.us12-nc 
source $file 
if ( $status > 0 ) then
   echo "ERROR: Could not source the Assigns file $file"
   echo "       Please contact the CMAS Help Desk at http://www.cmascenter.org"
   set exitstat = 1
endif

if ( $exitstat > 0 ) then
   exit( $exitstat )
endif

# Make necessary symbolic links
echo "Creating symbolic links..."
cd $SMKROOT/src
foreach dir ( biog cntlmat emmod emqa emutil grdmat inc lib \
              point smkinven smkmerge movesmrg spcmat temporal )
    cd $dir
    ln -s ../../scripts/make/Makeit ./
    if ( $status > 0 ) then
       echo "ERROR: Could not create a symbolic link in directory $dir"
       echo "       Please contact the CMAS Help Desk at http://www.cmascenter.org"
       set exitstat = 1
    endif
    cd ..
end

if ( $exitstat > 0 ) then
   exit( $exitstat )
endif

# Create the inventory list files
echo "Creating inventory list files..."

cd $INVDIR/area 
echo "#LIST" > arinv.area.lst
ls $INVDIR/area/arinv_avefire_2002ce_ida_nc.txt >> arinv.area.lst
ls $INVDIR/area/arinv_lm_no_c3_cap2002v3_orl_nc.txt >> arinv.area.lst
ls $INVDIR/area/arinv_nonpt_pf4_cap_nopfc_orl_nc.txt >> arinv.area.lst

cd $INVDIR/nonroad
echo "#LIST" > arinv.nonroad.lst
ls $INVDIR/nonroad/arinv_nonroad_caps_2005v2_jul_orl_nc.txt >> arinv.nonroad.lst

cd $INVDIR/point
echo "#LIST" > ptinv.point.lst
ls $INVDIR/point/ptinv_ptipm_cap2005v2_orl_nc.txt >> ptinv.point.lst
ls $INVDIR/point/ptinv_ptnonipm_xportfrac_cap2005v2_orl_nc.txt >> ptinv.point.lst

cd $INVDIR/rateperdistance_noRFL
echo "#LIST" > mbinv.rateperdistance_noRFL.lst
ls $INVDIR/rateperdistance_noRFL/mbinv.VMT.nc.Apr2011.txt >> mbinv.rateperdistance_noRFL.lst 
ls $INVDIR/rateperdistance_noRFL/mbinv.SPEED.nc.Apr2011.txt >> mbinv.rateperdistance_noRFL.lst 

mkdir $INVDIR/rateperdistance_RFLonly
cd $INVDIR/rateperdistance_RFLonly
ln -s ../rateperdistance_noRFL/mbinv.rateperdistance_noRFL.lst mbinv.rateperdistance_RFLonly.lst

cd $INVDIR/ratepervehicle_noRFL
echo "#LIST" > mbinv.ratepervehicle_noRFL.lst
ls $INVDIR/ratepervehicle_noRFL/mbinv.VPOP.nc.Apr2011.txt >> mbinv.ratepervehicle_noRFL.lst 

mkdir $INVDIR/ratepervehicle_RFLonly
cd $INVDIR/ratepervehicle_RFLonly
ln -s ../ratepervehicle_noRFL/mbinv.ratepervehicle_noRFL.lst mbinv.ratepervehicle_RFLonly.lst

mkdir $INVDIR/rateperprofile
cd $INVDIR/rateperprofile
ln -s ../ratepervehicle_noRFL/mbinv.ratepervehicle_noRFL.lst mbinv.rateperprofile.lst

# Return to original directory
cd $here
rm -rf $tmpfile

echo "Installation completed successfully."
echo " "
echo "Please follow the instructions in Section 4.3 of the SMOKE User's Manual"
echo "   to run the nctox test case."
echo "http://www.smoke-model.org/version3.0/html/ch04s03.html"
echo " "

exit( 0 )

