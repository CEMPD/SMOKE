#!/bin/csh -f

set exitstat = 0
set here = `pwd`
set tmpfile = .install

/bin/rm -rf $tmpfile

# Check that downloaded files are available
set file = smoke.nctox.data.tar.gz
echo "Checking for $file file..."

ls $file > $tmpfile
if ( $status > 0 ) then
    echo "ERROR: Could not find $file"
    echo "       Please download the file $file and try again"
    set exitstat = 1
endif

set file = smoke.Linux2_x86pg.tar.gz
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

# Check that EDSS_ROOT is set
if ( ! $?EDSS_ROOT ) then
    echo "ERROR: You must define the EDSS_ROOT environment variable"
    echo "       before running this script.  Use the following command"
    echo "       to set your EDSS_ROOT variable:"
    echo "          setenv EDSS_ROOT <your chosen SMOKE installation dir>"
    exit( 1 )
endif

echo "SMOKE v2.2 will be installed in the following directory:"
echo "      $EDSS_ROOT"
echo " "

# Install files
cd $EDSS_ROOT

set file = smoke.nctox.data.tar.gz
echo "Unpacking file $file..."
tar xzvf $here/$file
if ( $status > 0 ) then
   echo "ERROR: Could not unpack the file $file"
   echo "       Confirm that the 'tar' command is available and"
   echo "       that you have the correct permissions"
   set exitstat = 1
endif

set file = smoke.Linux2_x86pg.tar.gz
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
set file = ASSIGNS.nctox.cmaq.cb4p25_wtox.us12-nc
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
              mo6 mobile point smkinven smkmerge spcmat temporal )
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

cd $ARDAT
echo "#LIST" > arinv.stationary.lst
ls $ARDAT/arinv.nonpoint.nti99_NC.new.txt >> arinv.stationary.lst
ls $ARDAT/arinv.stationary.nei96_NC.ida.txt >> arinv.stationary.lst

cd $MBDAT
echo "#LIST" > mbinv.lst
ls $MBDAT/mbinv.nei99_NC.ida.txt >> mbinv.lst
ln -s mcref.nctox.txt mcref.nctox_18.txt
ln -s mvref.nctox.txt mvref.nctox_18.txt

cd $INVDIR/nonroad
echo "#LIST" > arinv.nonroad.lst
ls $INVDIR/nonroad/arinv.nonroad.n* >> arinv.nonroad.lst

cd $PTDAT
echo "#LIST" > ptinv.lst
ls $PTDAT/ptinv.n* >> ptinv.lst

# Return to original directory
cd $here
rm -rf $tmpfile

echo "Installation completed successfully."
echo " "
echo "Please follow the instructions in Section 4.4 of the SMOKE User's Manual"
echo "   to run the nctox test case."
echo "http://www.cep.unc.edu/empd/products/smoke/version2.1/html/ch04s04.html"
echo " "

exit( 0 )

