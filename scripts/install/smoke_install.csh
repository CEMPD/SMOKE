#!/bin/csh -f

set exitstat = 0
set blah_foo_arg = `pwd`

/bin/rm -rf .install_log

ls smoke.nctox* > .install_log
if ( $status > 0 ) then
    echo "ERROR: smoke.nctox.data.tar file is not available"
    echo "       Please download that file and retry."
    set exitstat = 1
endif

set found1 = n
if ( -e smoke.IRIXn32f90.tar.gz || -e smoke.IRIXn32f90.tar ) then
    echo "Installing smoke.IRIXn32f90.tar..."
    set found1 = y
endif
if ( -e smoke.Linux2_x86ifc.tar.gz || -e smoke.Linux2_x86ifc.tar ) then
    echo "Installing smoke.Linux2_x86ifc.tar..."
    set found1 = y
endif
if ( -e smoke.SunOS5f90.tar.gz || -e smoke.SunOS5f90.tar ) then
    echo "Installing smoke.SunOS5f90.tar..."
    set found1 = y
endif

if ( $found1 == n ) then
   echo "ERROR: No smoke installation file for a specific platform"
   echo "       is available. Please download that file and retry."
   echo "       
   set exitstat = 1
endif

if ( $exitstat > 0 ) then
   exit( $exitstat )
endif

set edssstat = 0
if ( $?EDSS_ROOT ) then
      set cnt = `echo $EDSS_ROOT | wc -c`
      @ cnt = $cnt - 1
      if ( $cnt > 24 ) then
         set edssstat = 1
      endif
else
    echo "ERROR: You must define the EDSS_ROOT environment variable"
    echo "       prior to running this script.  Use the command line"
    echo "       instruction as follows:"
    echo "          setenv EDSS_ROOT <your chosen SMOKE installation dir>"
    exit( 1 )
endif

# Check if gunzip is needed
ls smoke*tar* | grep gz > .install_log
set nline = `cat .install_log | wc -l`

# If needed, unzip files
if ( $nline > 0 ) then
   which gunzip
   if ( $status > 0 ) then
       echo "ERROR: Your computer does not have the gunzip utility for"
       echo "       unzipping the SMOKE tar files. You can mannually unzip"
       echo "       these files using the pkzip or zip utilities, if these"
       echo "       are available instead.  Then rerun this script."
       echo " "
   endif

   gunzip smoke*tar.gz

endif

# Install files
cd $EDSS_ROOT

tar xvf $blah_foo_arg/smoke.nctox*tar
if ( $status > 0 ) then
   set exitstat = 1
endif

set found1 = n
if ( -e $blah_foo_arg/smoke.IRIXn32f90.tar ) then
   tar xvf $blah_foo_arg/smoke.IRIXn32f90.tar
   set found1 = y
endif

if ( -e $blah_foo_arg/smoke.Linux2_x86ifc.tar ) then
   tar xvf $blah_foo_arg/smoke.Linux2_x86ifc.tar
   set found1 = y
endif

if ( -e $blah_foo_arg/smoke.SunOS5f90.tar ) then
   tar xvf $blah_foo_arg/smoke.SunOS5f90.tar
   set found1 = y
endif

if ( $found1 == n ) then
   echo "ERROR: Cannot find any executable tar files in directory"
   echo "       "$blah_foo_arg"
   echo "       "
   set exitstat = 1
endif

cd subsys/smokev1/assigns
source ASSIGNS.nctox.cmaq.cb4p25.us36-nc
if ( $status > 0 ) then
   set exitstat = 1
endif

if ( $exitstat > 0 ) then
    echo "ERROR: installation did not complete successfully."
    exit ( $exitstat )
else
    echo "CONGRATULATIONS: installation completed successfully."
endif

cd $ARDAT
ls $ARDAT/arinv.stationary.nei96_NC.ida.txt > arinv.stationary.lst
cd $MBDAT
ls $MBDAT/mbinv99v1_vmt_NC.ida > mbinv.lst
cd $INVDIR/nonroad
ls $INVDIR/nonroad/arinv.nonroad.n* > arinv.nonroad.lst
cd $PTDAT
ls $PTDAT/ptinv.nei96_NC.ida.txt > ptinv.lst
echo "DATERANGE 0709 0710" > pthour.lst
ls $SMKDAT/cem/1996/q3/* >> pthour.lst

cd $blah_foo_arg

if ( $edssstat > 0 ) then
     echo "NOTE: EDSS_ROOT may be too long to be able to run"
     echo "      MOBILE6.  The SMK_M6PATH will need to be"
     echo "      50 characters or less, and this will happen"
     echo "      when EDSS_ROOT > 24 characters and ESCEN >=5"
     echo "      characters. If you plan to use MOBILE6, reinstall"
     echo "      with a shorter EDSS_ROOT path. This limitation"
     echo "      will be removed in future releases."
endif

exit( 0 )

