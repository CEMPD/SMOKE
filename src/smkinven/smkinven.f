
         PROGRAM SMKINVEN

C***********************************************************************
C  program body starts at line 171
C
C  DESCRIPTION:
C    The smkinven program reads the any source inventory in one of four
C    formats: EMS-95, EPS, IDA, and SMOKE list format.  It permits a flexible
C    definition of a point source, which depends on the inventory input
C    formats.  It allows any number of inventory pollutants, within the limit
C    of I/O API (this works out to only 15 pollutants per file).
C
C  PRECONDITIONS REQUIRED:
C    Set environment variables:
C      PROMPTFLAG:    If N, default inputs are used
C      RAW_DUP_CHECK: If Y, duplicate sources disallowed
C      VELOC_RECALC:  If Y, recalculates velocity based on flow
C      WEST_HSPHERE:  If N, does not reset positive longitudes to negative
C    Input files:  
C      PTINV: ASCII point sources inventory
C      PSTK: Replacement stack parameters file
C      ZONES: Time zones files
C      SIPOLS: Master list of pollutant codes and names (in output order)
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C    Subroutines: I/O API subroutines, INITEM, CHECKMEM, RDTZONE, RDSIPOLS, 
C       RDPTINV, FMTCSRC, FIXSTK, WRPTSCC, WRPTREF, OPENPNTS, WPNTSCHR, 
C       WPNTSPOL
C    Functions: I/O API functions, GETFLINE, GETTZONE 
C
C  REVISION  HISTORY:
C    started 10/98 by M Houyoux as rawpoint.f from emspoint.F 4.3
C    smkinven changes started 4/98
C    toxics changes 11/2002  A. Holland
C
C***************************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2004, Environmental Modeling for Policy Development
C All Rights Reserved
C 
C Carolina Environmental Program
C University of North Carolina at Chapel Hill
C 137 E. Franklin St., CB# 6116
C Chapel Hill, NC 27599-6116
C 
C smoke@unc.edu
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C***************************************************************************

C...........   MODULES for public variables
C...........   This module is the inventory arrays
        USE MODSOURC, ONLY: CIFIP, CISIC, CSRCTYP, TZONES, CSCC, IDIU, IWEK,
     &                      CINTGR, CEXTORL
C.........  This module contains the lists of unique inventory information
        USE MODLISTS, ONLY: MXIDAT, INVSTAT, INVDNAM, FIREFLAG, INVDCOD, 
     &                      FF10FLAG, MEDSFLAG, NCDFLAG, APIFLAG, FIREFF10
C.........  This module contains the information about the source category
        USE MODINFO, ONLY: CATEGORY, NIPOL, NIACT, NIPPA, EIIDX, INV_MON,
     &                     EINAM, AVIDX, ACTVTY, EANAM, NSRC

C.........  This module contains data for day- and hour-specific data
        USE MODDAYHR, ONLY: DAYINVFLAG, HRLINVFLAG

        IMPLICIT NONE

C...........   INCLUDES:
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.
        INCLUDE 'SETDECL.EXT'   !  FileSetAPI variables and functions

C...........   EXTERNAL FUNCTIONS and their descriptions:
        
        CHARACTER(2)  CRLF
        INTEGER       ENVINT
        LOGICAL       ENVYN
        INTEGER       INDEX1
        INTEGER       GETFLINE
        INTEGER       GETTZONE
        INTEGER       STR2INT
        LOGICAL       USEEXPGEO

        EXTERNAL      CRLF, ENVINT, ENVYN, INDEX1, GETFLINE, GETTZONE, 
     &                STR2INT, USEEXPGEO

C...........  LOCAL PARAMETERS and their descriptions:

        CHARACTER(50), PARAMETER :: CVSW = '$Name SMOKEv4.8.1_Jan2021$' ! CVS release tag

C.........  LOCAL VARIABLES and their descriptions:

C.........  Day-specific and hour-specific variable indices
        INTEGER     ::  DEAIDX( MXVARS3 ) = 0
        INTEGER     ::  DSPIDX( MXSPDAT ) = 0
        INTEGER     ::  HEAIDX( MXVARS3 ) = 0
        INTEGER     ::  HSPIDX( MXSPDAT ) = 0

C.........  Array that contains the names of the inventory variables needed for
C           this program
        CHARACTER(IOVLEN3) IVARNAMS( MXINVARR )

C.........  File units and logical/physical names

        INTEGER    :: ADEV = 0  !  unit no. for REPINVEN file
        INTEGER    :: CDEV = 0  !  unit no. for SCCs description
        INTEGER    :: DDEV = 0  !  unit no. for day-specific input file 
        INTEGER    :: HDEV = 0  !  unit no. for hour-specific input file 
        INTEGER    :: IDEV = 0  !  unit no. for inventory file (various formats)
        INTEGER    :: LDEV = 0  !  unit no. for log file
        INTEGER    :: MDEV = 0  !  unit no. for map inventory file
        INTEGER    :: ODEV = 0  !  unit number for ORIS description
        INTEGER    :: PDEV = 0  !  unit number for inventory data table
        INTEGER    :: RDEV = 0  !  unit no. for def stack pars or mobile codes
        INTEGER    :: UDEV = 0  !  unit no. for non-HAP inclusions/exclusions file
        INTEGER    :: SDEV = 0  !  unit no. for ASCII output inventory file
        INTEGER    :: TDEV = 0  !  unit no. for MEDS daily/hourly inventory file
        INTEGER    :: YDEV = 0  !  unit no. for area-to-point factors file
        INTEGER    :: ZDEV = 0  !  unit no. for time zone file

        CHARACTER(NAMLEN3) :: ANAME = ' '! inven ASCII output logical name
        CHARACTER(NAMLEN3) :: DNAME = ' '! day-specific input logical name
        CHARACTER(NAMLEN3) :: ENAME = ' '! inven I/O API output logical name
        CHARACTER(NAMLEN3) :: GNAME = ' '! gridded I/O API input logical
        CHARACTER(NAMLEN3) :: HNAME = ' '! hour-specific input logical name
        CHARACTER(NAMLEN3) :: INAME = ' '! inven input logical name
        CHARACTER(NAMLEN3) :: TNAME = ' '! tmp day- or hour-specific input logical name

C...........   Other local variables
                                
        INTEGER         S, I, J, J1, J2, J3, K, L, L2, V !  counters and indices

        INTEGER      :: DNSTEP = 0 ! day-specific data time step number
        INTEGER      :: DSDATE = 0 ! day-specific data start date
        INTEGER      :: DSTIME = 0 ! day-specific data start time

        INTEGER      :: HNSTEP = 0 ! day-specific data time step number
        INTEGER      :: HSDATE = 0 ! day-specific data start date
        INTEGER      :: HSTIME = 0 ! day-specific data start time
        INTEGER         INSTEP     ! expected input time step HHMMSS
        INTEGER         IOS        ! I/O status
        INTEGER         MAXK       ! test for maximum value of K in output loop
        INTEGER      :: MXSRCDY= 0 ! max no. day-specific sources
        INTEGER      :: MXSRCHR= 0 ! max no. hour-specific sources
        INTEGER      :: NDAT = 0   ! tmp no. actual pols & activities
        INTEGER      :: NINVARR = 0! no. inventory variables to read
        INTEGER      :: NRAWBP = 0 ! number of sources with pollutants
        INTEGER      :: NRAWSRCS= 0! number of unique sources
        INTEGER      :: NVARDY = 0 ! no. day-specific variables
        INTEGER      :: NVSPDY = 0 ! no. day-specific special variables
        INTEGER      :: NVARHR = 0 ! no. hour-specific variables
        INTEGER      :: NVSPHR = 0 ! no. hour-specific special variables
        INTEGER         OUTSTEP    ! output time step HHMMSS for day/hour data
        INTEGER         TZONE      ! output time zone for day- & hour-specific

        LOGICAL         A2PFLAG          ! true: using area-to-point processing
        LOGICAL      :: CFLAG = .FALSE.  ! true: CEM processing
        LOGICAL         IFLAG            ! true: average inventory inputs used
        LOGICAL         NONPOINT         ! true: importing nonpoint inventory
        LOGICAL         ORLFLG           ! true: ORL format inventory
        LOGICAL         STKFLG           ! true: check stack parameters

        CHARACTER(5)          TYPNAM      !  'day' or 'hour' for import
        CHARACTER(256)        MESG        !  message buffer
        CHARACTER(IOVLEN3) :: GRDNM = ' ' !  I/O API input file grid name
        CHARACTER(PHYLEN3) :: VARPATH = './' ! path for pol/act files
        CHARACTER(FIPLEN3)    CFIP        ! temporary FIPS code
        CHARACTER(FIPLEN3)    PCFIP       ! previous FIPS code

        CHARACTER(16) :: PROGNAME = 'SMKINVEN'   !  program name

C***********************************************************************
C   begin body of program SMKINVEN

        LDEV = INIT3()

C.........  Write out copyright, version, web address, header info, and prompt
C           to continue running the program.
        CALL INITEM( LDEV, CVSW, PROGNAME )

C.........  Set source category based on environment variable setting
        CALL GETCTGRY

C.........  Output time zone
        TZONE = ENVINT( 'OUTZONE', 'Output time zone', 0, IOS )

C.........  Get names of input files
        CALL OPENINVIN( CATEGORY, IDEV, DDEV, HDEV, RDEV, SDEV,
     &                  PDEV, ZDEV, CDEV, ODEV, UDEV, YDEV,
     &                  ENAME, INAME, DNAME, HNAME )

C.........  Set controller flags depending on unit numbers
        IFLAG   = ( IDEV .NE. 0 )
        A2PFLAG = ( YDEV .NE. 0 )

C.........  Set IOAPI or NetCDF gridded input file name, if available
        IF( APIFLAG .OR. NCDFLAG ) GNAME = ENAME

        MESG = 'Setting up to read inventory data...'
        CALL M3MSG2( MESG )

C.........  Read country, state, and county file for time zones
        IF( USEEXPGEO() ) THEN
            CALL RDGEOCODES( 1, I )
        ELSE IF( ZDEV .GT. 0 ) THEN
            CALL RDSTCY( ZDEV, 1, I )   !  "I" used as a dummy
        END IF

C.........  Read, sort, and store inventory data table file
        CALL RDCODNAM( PDEV )

C.........  Read/store informatino for MEDS inv processing
        IF( MEDSFLAG ) CALL RDMEDSINFO  ! read GAI_LOOKUP_TABLE (col/row to lat/lon) for MEDS format inv

C.........  Define the month of inventory processing
        MESG = 'Define the processing inventory month'
        INV_MON = ENVINT( 'SMKINVEN_MONTH', MESG, 0, IOS )

C.........  Process for ASCII average day or annual inventory
        IF( IFLAG ) THEN

C.............  Read the source information from the raw inventory files, 
C               store in unsorted order, and determine source IDs
C.............  The arrays that are populated by this subroutine call
C               are contained in the module MODSOURC
            IF( INV_MON == 0 ) THEN
                MESG = 'Processing Annual inventory....'
                CALL M3MSG2( MESG )
            ELSE
                MESG = 'Processing ' // TRIM( MON_NAME( INV_MON ) ) // 
     &                 ' inventory'
                CALL M3MSG2( MESG )
            END IF

            CALL M3MSG2( 'Reading inventory sources...' )

            CALL RDINVSRCS( IDEV, INAME, NRAWBP, NRAWSRCS, ORLFLG )

C.............  Process source information and store in sorted order
            CALL M3MSG2( 'Processing inventory sources...' )
            
            CALL PROCINVSRCS( NRAWSRCS )

C.............  Read the data from the raw inventory files and store in 
C               sorted order
            CALL M3MSG2( 'Reading inventory data...' )
            
            CALL RDINVDATA( IDEV, INAME, NRAWBP, NONPOINT )

C.............  Check if ORL inventory and reset NHAPEXCLUDE if needed
            IF( .NOT. ORLFLG .AND. UDEV > 0 ) THEN
                MESG = 'NOTE: Ignoring SMK_PROCESS_HAPS setting ' //
     &                 'since an ORL or FF10 inventory is not being processed'
                CALL M3MSG2( MESG )
                UDEV = 0
            END IF

C.............  Process inventory records and store in sorted order
            CALL M3MSG2( 'Processing inventory data...' )

            CALL PROCINVEN( NRAWBP, UDEV, YDEV, CDEV, LDEV ) 

C.............  Integrate criteria and toxic pollutants
            IF( ORLFLG ) THEN
                CALL SETNONHAP( NRAWBP )
            END IF

C.............  Determine memory needed for actual pollutants list and actual
C               activities list and allocate them. Invstat has been updated
C               to be +/- 2 depending on whether the pollutant or activity was
C               present in the inventory.
            NIPOL = 0
            NIACT = 0
            DO I = 1, MXIDAT
                IF( INVSTAT( I ) .GT.  1 ) NIPOL = NIPOL + 1
                IF( INVSTAT( I ) .LT. -1 ) NIACT = NIACT + 1
            ENDDO

            NIPPA = NIPOL + NIACT

            ALLOCATE( EIIDX( NIPOL ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EIIDX', PROGNAME )
            ALLOCATE( EINAM( NIPOL ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EINAM', PROGNAME )
            ALLOCATE( AVIDX( NIACT ), STAT=IOS )
            CALL CHECKMEM( IOS, 'AVIDX', PROGNAME )
            ALLOCATE( ACTVTY( NIACT ), STAT=IOS )
            CALL CHECKMEM( IOS, 'ACTVTY', PROGNAME )
            ALLOCATE( EANAM( NIPPA ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EANAM', PROGNAME )

C.............  Create list of actual pollutants and activities and indexes to 
C               the master list. The order in EINAM and ACTVTY will be the  
C               output order. The indexes are for accessing INVDCOD, if needed.
C.............  These are for opening output file and processing output data
            J1 = 0
            J2 = 0
            J3 = 0
            DO I = 1, MXIDAT

               IF( INVSTAT( I ) .GT. 0 ) THEN
                   J1 = J1 + 1
                   J3 = J3 + 1
                   EIIDX( J3 ) = I
                   EINAM( J3 ) = INVDNAM( I )
                   EANAM( J1 ) = INVDNAM( I )
               END IF

               IF( INVSTAT( I ) .LT. 0 ) THEN
                   J1 = J1 + 1
                   J2 = J2 + 1
                   AVIDX ( J2 ) = I
                   ACTVTY( J2 ) = INVDNAM( I )
                   EANAM ( J1 ) = INVDNAM( I )
               END IF
            
            END DO

C.............   Fix stack parameters for point sources
C.............   Some of these arguments are variables that are defined in the
C                module MODSOURC
            IF( CATEGORY .EQ. 'POINT' ) THEN
                MESG = 'Check stack parameter values'
                STKFLG = ENVYN( 'CHECK_STACKS_YN', MESG, .TRUE., IOS )
                IF( FIREFLAG .OR. FIREFF10 ) STKFLG = .FALSE.    ! skip check stack para when wildfire
                IF( STKFLG ) CALL FIXSTK( RDEV, NSRC )
            END IF

C.............  Set time zones based on country/state/county code. Note that a
C               few counties in the Western U.S. are divided by a time zone, so 
C               this is not perfectly accurate for all counties.
            ALLOCATE( TZONES( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'TZONES', PROGNAME )

            PCFIP = ''
            DO S = 1, NSRC
                CFIP = CIFIP( S )
 
                IF( CFIP /= PCFIP ) THEN
                    TZONES( S ) = GETTZONE( CFIP )
                    PCFIP = CFIP
                ELSE
                    TZONES( S ) = TZONES( S - 1 )
                END IF
            END DO

C.............  Write out primary inventory files. Do this before the day- or 
C               hour-specific processing so that if there is a problem, the
C               lengthy inventory import does not need to be redone...

C.............  Write out SCC file
            CALL WRCHRSCC( CSCC )

        END IF  ! For ASCII annual/ave-day inputs

C.........  Input gridded I/O API inventory data
        IF( APIFLAG .OR. NCDFLAG ) THEN

            IF( APIFLAG ) CALL RDGRDAPI( GNAME, GRDNM ) 
            IF( NCDFLAG ) CALL RDGRDNCF( INAME, ENAME ) 

C.............  initialize arrays for later 
            ALLOCATE( CISIC ( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CISIC', PROGNAME )
            ALLOCATE( CSRCTYP( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CSRCTYP', PROGNAME )
            ALLOCATE( CINTGR( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CINTGR', PROGNAME )
            ALLOCATE( CEXTORL( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CEXTORL', PROGNAME )

            CISIC   = ' '       ! array
            CSRCTYP = ' '       ! array
            CINTGR  = ' '       ! array
            CEXTORL = ' '       ! array

            IFLAG = .TRUE.

        END IF  ! For gridded I/O API NetCDF inventory

C.........  Output SMOKE inventory files
        IF( IFLAG ) THEN

C.............  Generate message to use just before writing out inventory files
C.............  Open output I/O API and ASCII files 

            CALL OPENINVOUT( A2PFLAG, GRDNM, ENAME, ANAME, MDEV, SDEV,
     &                       ADEV, VARPATH )

            MESG = 'Writing SMOKE ' // TRIM( CATEGORY ) // 
     &             ' SOURCE INVENTORY file...'

            CALL M3MSG2( MESG )

C.............  Write source characteristics to inventory files (I/O API and
C               ASCII)
            CALL WRINVCHR( ENAME, SDEV, A2PFLAG, NONPOINT )

C.............  Deallocate sorted inventory info arrays, except CSOURC
            CALL SRCMEM( CATEGORY, 'SORTED', .FALSE., .FALSE., 1, 1, 1 )

C.............  Write out average inventory data values
C.............  Compute inventory data values, if needed
            CALL WRINVEMIS( MDEV, VARPATH )

C.............  Deallocate sorted inventory info arrays
            CALL SRCMEM( CATEGORY, 'SORTED', .FALSE., .TRUE., 1, 1, 1 )

C.........  If the inventory is not being created, then read necessary
C           information from existing inventory files, which will be used
C           for day- and hour-specific data import.
        ELSE

C.............  Store source-category-specific header information, 
C           including the inventory pollutants in the file (if any).  Note that 
C           the I/O API head info is passed by include file and the
C           results are stored in module MODINFO.
            CALL GETSINFO( ENAME )

C.............  Since not reading inventory files yet, but need to know
C               if this is a fire inventory, scan inventory data names
C               to see if HFLUX is present. If so, FIREFLAG = .true.
            IF ( INDEX1( 'HFLUX',NIPPA,EANAM ) .GT. 0 ) FIREFLAG= .TRUE.

            NINVARR = 5
            IVARNAMS( 1 ) = 'CIFIP'   ! In case CEM input
            IVARNAMS( 2 ) = 'CSOURC'  ! In case non-CEM input
            IVARNAMS( 3 ) = 'CSCC'    ! In case CEM input (for reporting)
            IVARNAMS( 4 ) = 'CPDESC'  ! In case CEM input
            IVARNAMS( 5 ) = 'CINTGR'  ! for VOC + HAPs integration

            IF ( .NOT. FIREFLAG ) THEN
                NINVARR = 7
                IVARNAMS( 6 ) = 'CORIS'   ! In case CEM input
                IVARNAMS( 7 ) = 'CBLRID'  ! In case CEM input
            END IF

            CALL RDINVCHR( CATEGORY, ENAME, SDEV, NSRC, NINVARR,
     &                     IVARNAMS )

        END IF !   End processing of average annual import or not

C.........  Read in daily or hourly MEDS emission values and output to a SMOKE inter output file
        IF( MEDSFLAG ) THEN
C.............  Read and output day-specific data
            IF( DAYINVFLAG ) THEN
                TYPNAM = 'day'
                TDEV   = DDEV
                CALL GENMEDSOUT( TDEV, TNAME, TZONE, TYPNAM ) 

            ELSE IF( HRLINVFLAG ) THEN
                TYPNAM = 'hour'
                TDEV   = HDEV
                CALL GENMEDSOUT( TDEV, TNAME, TZONE, TYPNAM ) 

            END IF

        END IF

C.........  Read in daily emission values and output to a SMOKE file
        IF( DAYINVFLAG .AND. .NOT. MEDSFLAG ) THEN

            INSTEP  = 240000
            OUTSTEP = 10000
            TYPNAM  = 'day'

C.............  Preprocess day-specific file(s) to determine memory needs.
C               Also determine maximum and minimum dates for output file.
            CALL GETPDINFO( DDEV, TZONE, INSTEP, OUTSTEP, TYPNAM, DNAME, 
     &                      DSDATE, DSTIME, DNSTEP, NVARDY, NVSPDY, 
     &                      MXSRCDY, DEAIDX, DSPIDX, CFLAG )

C.............  Read and output day-specific data
            CALL GENPDOUT( DDEV, CDEV, ODEV, ADEV, TZONE, DSDATE,DSTIME, 
     &                     DNSTEP, INSTEP, OUTSTEP, NVARDY, NVSPDY, 
     &                     MXSRCDY, TYPNAM, DNAME, DEAIDX, DSPIDX, CFLAG )

        END IF

C.........  Read in hourly emission values and output to a SMOKE file
        IF( HRLINVFLAG .AND. .NOT. MEDSFLAG ) THEN

            INSTEP  = 10000
            OUTSTEP = 10000
            TYPNAM  = 'hour'

C.............  Preprocess hour-specific file(s) to determine memory needs.
C               Also determine maximum and minimum dates for output file.
            CALL GETPDINFO( HDEV, TZONE, INSTEP, OUTSTEP, TYPNAM, HNAME, 
     &                      HSDATE, HSTIME, HNSTEP, NVARHR, NVSPHR, 
     &                      MXSRCHR, HEAIDX, HSPIDX, CFLAG )

C.............  Read and output hour-specific data
            CALL GENPDOUT( HDEV, CDEV, ODEV, ADEV,TZONE, HSDATE,HSTIME,
     &                     HNSTEP, INSTEP, OUTSTEP, NVARHR, NVSPHR,  
     &                     MXSRCHR, TYPNAM, HNAME, HEAIDX, HSPIDX, CFLAG )

        END IF

C.............  Write inventory report file
C        CALL M3MSG2( 'Writing inventory report file...' )
C        IF( .NOT. MEDSFLAG ) CALL WREPINVEN( ADEV, CDEV )

C.........  End program successfully
999     MESG = ' '
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 0 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )

92010   FORMAT( 5X, A, :, I10 )


C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C93041   FORMAT( I5, X, I5, X, I3, X, I8, X, I5.5, 3( X, I3 ) )

93060   FORMAT( 10( A, :, E10.3, :, 1X ) )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94020   FORMAT( A, 1X, I5.5, 1X, A, 1X, I8.8, 1X,
     &          A, I6, 1X, A, I6, 1X, A, :, I6 )

94040   FORMAT( A, I2.2 )

94060   FORMAT( 10( A, :, E10.3, :, 1X ) )

94080   FORMAT( '************  ', A, I7, ' ,  ' , A, I12 )

        END PROGRAM SMKINVEN
