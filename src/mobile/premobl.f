 
        PROGRAM PREMOBL
 
C***********************************************************************
C  program body starts at line 
C
C  DESCRIPTION:
C     This program inputs gridded, time-dependent temperature data, a mobile 
C     list file, and an ungridding matrix and
C     determines the min/max temperature combinations for each mobile source 
C     and for each emission factor set.
C
C     This version of the program assumes that MOBILE5 is the emission
C     factor model of interest, and therefore is focused on temperature
C     as the met data needing preprocessing, and sets the permitted min/max
C     temperature ranges accordingly.  Ultimately, this program could become
C     PREEMFAC, and do more general met pre-processing.  One input would then
C     need to be the activity data names associated with a method, such as
C     MOBILE5, for computing emissions.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Copied from prediur.F 3.4 by M. Houyoux
C
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1999, MCNC--North Carolina Supercomputing Center
C All Rights Reserved
C
C See file COPYRIGHT for conditions of use.
C
C Environmental Programs Group
C MCNC--North Carolina Supercomputing Center
C P.O. Box 12889
C Research Triangle Park, NC  27709-2889
C
C env_progs@mcnc.org
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C****************************************************************************

C...........   MODULES for public variables
C...........   This module is the source inventory arrays
        USE MODSOURC

C...........   This module contains the cross-reference tables
        USE MODXREF

C.........  This module contains emission factor tables and related
        USE MODEMFAC

C...........   This module contains the information about the source category
        USE MODINFO

C...........   This module is the derived meteorology data for emission factors
        USE MODMET

        IMPLICIT NONE
 
C...........   INCLUDES:
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS and their descriptions:

        CHARACTER*2     CRLF
        INTEGER         GETIFDSC
        CHARACTER*14    MMDDYY
        INTEGER         PROMPTFFILE
        CHARACTER*16    PROMPTMFILE
        INTEGER         WKDAY

        EXTERNAL     CRLF, GETIFDSC, MMDDYY, PROMPTFFILE, PROMPTMFILE,
     &               WKDAY

C...........   LOCAL PARAMETERS
        CHARACTER*50, PARAMETER :: SCCSW = ! SCCS string with version no. at end
     &               '@(#)$Id$'

C...........   LOCAL VARIABLES and their descriptions:

C...........   Gridded meteorology data (dim: NGRID)
        REAL   , ALLOCATABLE :: TA( : )   !  one layer of temperature

C...........   Ungridding Matrix
        INTEGER, ALLOCATABLE :: UMAT( : ) ! Contiguous ungridding matrix

C...........  Allocatable per-source arrays
        INTEGER, ALLOCATABLE :: DAYBEGT ( : ) ! daily start time HHMMSS
        INTEGER, ALLOCATABLE :: DAYENDT ( : ) ! daily end time HHMMSS
        LOGICAL, ALLOCATABLE :: LDAYSAV ( : ) ! true: src uses daylight time

C...........  Alloctable array for flagging needed PSI/min-max-combos
C...........  To conserve memory, this could be allocated and reallocated
C             at each time step for only the used min-max-combos
C...........  This is used only by the genefmet routine
        INTEGER, ALLOCATABLE :: TIPSI( :, :, : )

C...........  Array for the variable names of the SRC-PSI file
        CHARACTER(LEN=IOVLEN3), ALLOCATABLE :: PSVNAME( : )

C...........  Array that contains the names of the inventory variables needed 
C             for this program
        CHARACTER(LEN=IOVLEN3) IVARNAMS( MXINVARR )

C...........   File units and logical names:

        INTEGER      EDEV  ! unit number for output METIDX per PSI
        INTEGER      IDEV  ! unit number for emissions factors xref file
        INTEGER      LDEV  ! unit number for log file
        INTEGER      RDEV  ! unit number for mobile codes conversions file
        INTEGER      SDEV  ! unit number for ASCII inventory file

        CHARACTER*16 ANAME ! logical name for mobile ASCII inventory file
        CHARACTER*16 ENAME ! logical name for mobile I/O API inventory file
        CHARACTER*16 FNAME ! logical name for output METIDX per PSI
        CHARACTER*16 MNAME ! logical name for output ungridded min/max temp
        CHARACTER*16 PNAME ! logical nm for output PSIs by src, actvty, & 24 hrs
        CHARACTER*16 TNAME ! logical name for surface temp input file
        CHARACTER*16 UNAME ! logical name for ungridding-matrix input file

C...........   Other local variables:

        INTEGER    H, I, J, K, L, S, T, V  ! Counters and pointers

        INTEGER    DAY     !  tmp day of week number
        INTEGER    EDATE   !  ending input date counter (YYYYDDD) in GMT
        INTEGER    ENLEN   !  length of the emissions inven name
        INTEGER    ETIME   !  ending input time counter (HHMMSS)  in GMT
        INTEGER    IDATE   !  output date for min/max
        INTEGER    IOS     !  temporary I/O status
        INTEGER    ITIME   !  output time for min/max
        INTEGER    JDATE   !  input date counter (YYYYDDD) in GMT
        INTEGER    JTIME   !  input time counter (HHMMSS)  in GMT
        INTEGER    LDATE   !  date from previous loop iteration
        INTEGER    NCOLS   !  no. grid columns
        INTEGER    NCOLSU  !  no. grid columns in ungridding matrix
        INTEGER    NGRID   !  no. grid cells
        INTEGER    NINVARR !  no. inventory variables to read
        INTEGER    NMATX   !  size of ungridding matrix
        INTEGER    NROWS   !  no. grid rows
        INTEGER    NROWSU  !  no. grid rows in ungridding matrix
        INTEGER    NSTEPS  !  number of time steps to process temperature data
        INTEGER    ODATE   !  output date
        INTEGER    OTIME   !  time in GMT for determining when to output
        INTEGER    OSRC    !  number of sources outside grid
        INTEGER    PSI     !  tmp parameter scheme index
        INTEGER    SDATE   !  output start date
        INTEGER    SDATE_MET ! met file start date
        INTEGER    STIME   !  output start time
        INTEGER    STIME_MET ! met file start time
        INTEGER    TSTEP   !  time step of input temperature data (HHMMSS)
        INTEGER    TZONE   !  zone to determine output days

        REAL       TMAX   !  deg F temporary source max temperature on interval
        REAL       TMIN   !  deg F temporary source min temperature on interval
        REAL       MSAV   !  temporary saved minimum temperature

        LOGICAL :: EFLAG    = .FALSE.  !  true: error found
        LOGICAL :: LASTTIME = .FALSE.  !  true: final time step
        LOGICAL :: OFLAG    = .FALSE.  !  true: ungridding is 0 for some srcs

        CHARACTER(LEN=IOVLEN3) :: TVARNAME    !  temperature variable name
        CHARACTER*300             MESG        !  message buffer

        CHARACTER*16 :: PROGNAME = 'PREMOBL'   !  program name

C***********************************************************************
C   begin body of program PREMOBL

        LDEV = INIT3()

C.........  Write out copywrite, version, web address, header info, and prompt
C           to continue running the program.
        CALL INITEM( LDEV, SCCSW, PROGNAME )

C.........  Set source category based on environment variable setting
        CALL GETCTGRY

C.........  End program if source category is not mobile sources
        IF( CATEGORY .NE. 'MOBILE' ) THEN
            L = LEN_TRIM( PROGNAME )
            MESG = 'Program ' // PROGNAME( 1:L ) // ' does not ' //
     &             'support ' // CATEGORY( 1:CATLEN ) // ' sources.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF

C.........  Obtain settings from the environment...
C.........  Get the name of the activity to use for one run
        MESG = 'Temperature variable name'
        CALL ENVSTR( 'TVARNAME', MESG, 'TEMP1P5', TVARNAME, IOS )

C.........  Set default name of meterology file, depending on the name of the
C           temperature variable
        TNAME = 'MET_CRO_2D'
        IF( TVARNAME .EQ. 'TA' ) TNAME = 'MET_CRO_3D'

C.........  Get inventory file names given source category
        CALL GETINAME( CATEGORY, ENAME, ANAME )

C.......   Get file names and units; open input files

        ENAME = PROMPTMFILE( 
     &          'Enter logical name for I/O API INVENTORY file',
     &          FSREAD3, ENAME, PROGNAME )
        ENLEN = LEN_TRIM( ENAME )

        SDEV = PROMPTFFILE( 
     &           'Enter logical name for ASCII INVENTORY file',
     &           .TRUE., .TRUE., ANAME, PROGNAME )

        UNAME = PROMPTMFILE(
     &          'Enter logical name for UNGRIDDING MATRIX file',
     &          FSREAD3, CRL // 'UMAT', PROGNAME )

        TNAME = PROMPTMFILE(
     &         'Enter logical name for SURFACE TEMPERATURE file',
     &          FSREAD3, TNAME, PROGNAME )

        IDEV = PROMPTFFILE(
     &         'Enter logical name for EMISSION FACTOR INDEX LIST file',
     &         .TRUE., .TRUE., CRL // 'PLIST', PROGNAME )

C.........  Get file name for converting road-class to road type &
C           vehicle type name to vehicle type number.
        MESG = 'Enter logical name for MOBILE CODES file'
        RDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'MCODES', PROGNAME )

C.........  Get header description of inventory file, error if problem
        IF( .NOT. DESC3( ENAME ) ) THEN
            MESG = 'Could not get description of file "' //
     &             ENAME( 1:ENLEN ) // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C.........  Otherwise, store source-category-specific header information, 
C           including the inventory pollutants in the file (if any).  Note that 
C           the I/O API header info is passed by include file and the
C           results are stored in module MODINFO.
        ELSE

            CALL GETSINFO

        END IF

C.........  Create note about time zone expected in meteorology file
        WRITE( MESG, 94010 )
     &     'NOTE: Time stamps of input meteorology file are assumed ' //
     &     CRLF() // BLANK5 // '      to be in GMT'
        CALL M3MSG2( MESG )

C.........  Read header of temperature file
        IF ( .NOT. DESC3( TNAME ) ) THEN

            MESG = 'Could not get description of file ' // TNAME
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C.........  Save header information that will be needed later
        ELSE
            SDATE_MET = SDATE3D
            STIME_MET = STIME3D
            TSTEP = TSTEP3D
            NSTEPS= MXREC3D
            NROWS = NROWS3D
            NCOLS = NCOLS3D
            NGRID = NROWS * NCOLS

        END IF

C.........  Read header of ungridding matrix...
        IF( .NOT. DESC3( UNAME ) ) THEN
            MESG = 'Could not get description for file ' // UNAME
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Store number of ungridding factors
        NMATX = NCOLS3D

C.........  Check dimensions of ungridding matrix...
C.........  Check the number of sources 
        CALL CHKSRCNO( CATDESC, UNAME, NROWS3D, NSRC, EFLAG )

C.........  Compare the gridded file settings with the ungridded file settings
        NCOLSU = GETIFDSC( FDESC3D, '/NCOLS3D/', .TRUE. )
        NROWSU = GETIFDSC( FDESC3D, '/NROWS3D/', .TRUE. )
        CALL CHECK_GRID_DIMS( 'NCOLS', 'columns', NCOLSU, NCOLS )
        CALL CHECK_GRID_DIMS( 'NROWS', 'rows'   , NROWSU, NROWS )

C......... If the dimensions were in error, abort
        IF( EFLAG ) THEN
            MESG = 'Ungridding matrix is inconsistent with inventory '//
     &             'and/or gridded ' // CRLF() // BLANK5 //
     &             'meteorology data.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Set inventory variables to read
        NINVARR = 6
        IVARNAMS( 1 ) = 'IFIP'
        IVARNAMS( 2 ) = 'IRCLAS'
        IVARNAMS( 3 ) = 'CSOURC'
        IVARNAMS( 4 ) = 'CSCC'
        IVARNAMS( 5 ) = 'CLINK'
        IVARNAMS( 6 ) = 'TZONES'

C.........  Allocate memory for and read in required inventory characteristics
        CALL RDINVCHR( CATEGORY, ENAME, SDEV, NSRC, NINVARR, IVARNAMS )

C.........  Build unique lists of SCCs and country/state/county codes
C           from the inventory arrays
        CALL GENUSLST

C.........  Retrieve environment variable settings for temperature ranges
C.........  Populate table of valid min/max temperatures in MODMET
        CALL TMPRINFO( .TRUE., 'BOTH' )

C.........  Allocate memory for other arrays in the program
        ALLOCATE( UMAT( NSRC + 2*NMATX ), STAT=IOS )
        CALL CHECKMEM( IOS, 'UMAT', PROGNAME )
        ALLOCATE( TA( NGRID ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TA', PROGNAME )
        ALLOCATE( DAYBEGT( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'DAYBEGT', PROGNAME )
        ALLOCATE( DAYENDT( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'DAYENDT', PROGNAME )
        ALLOCATE( LDAYSAV( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LDAYSAV', PROGNAME )
        ALLOCATE( TASRC( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TASRC', PROGNAME )
        ALLOCATE( TKMIN( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TKMIN', PROGNAME )
        ALLOCATE( TKMAX( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TKMAX', PROGNAME )
        ALLOCATE( TKMINOUT( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TKMINOUT', PROGNAME )
        ALLOCATE( TKMAXOUT( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TKMAXOUT', PROGNAME )
        ALLOCATE( METIDX( NSRC,4 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'METIDX', PROGNAME )
        ALLOCATE( EFSIDX( NSRC,NIACT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EFSIDX', PROGNAME )
        ALLOCATE( PSVNAME( NIACT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'PSVNAME', PROGNAME )

C.........  Create array of which sources are affected by daylight savings
        CALL GETDYSAV( NSRC, IFIP, LDAYSAV )

C.........  Read ungridding matrix 

        CALL RDUMAT( UNAME, NSRC, NMATX, NMATX, 
     &               UMAT( 1 ), UMAT( NSRC+1 ), UMAT( NSRC+NMATX+1 ) )

C.........  Get default episode information for processing met data
        SDATE = SDATE_MET
        STIME = STIME_MET
        TZONE = 0
        CALL GETM3EPI( TZONE, SDATE, STIME, NSTEPS )

C.........  Set end date and time for run
        EDATE = SDATE
        ETIME = STIME
        CALL NEXTIME( EDATE, ETIME, NSTEPS*10000 )

C.........  Preprocess dates, times, and time zones.  Write report for
C           time zones here for zones that do no have a complete day of
C           data in beginning or end.  This is only an issue because of
C           computing min/max data per day.
        CALL CHKFULLDY( NSRC, SDATE, STIME, EDATE, ETIME, 
     &                  TZONES, LDAYSAV )

C.........  Fill tables for translating mobile road classes and vehicle types
C.........  The tables are passed through MODINFO
        CALL RDMVINFO( RDEV )

C...........   Read list of parameter scheme indices and emission factor
C...........   table names to use for each source
 
C.........  Read the cross-reference for assigning the parameter scheme indices
C           to the sources. This will populate parts of the MODXREF module.
C.........  Also create the list of unique parameter scheme indices.
        CALL RDEFXREF( IDEV, .TRUE. )

C...........   Map parameter scheme indexes (PSIs) onto sources for all hours
        CALL ASGNPSI( NIACT, ACTVTY )

C.........  Allocate memory for flag of which PSI/temp-combos are used. Have
C           waited until now because need MXNPSI from RDEFXREF
        ALLOCATE( TIPSI( MXXNPSI, NVLDTMM, NIACT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TIPSI', PROGNAME )

C.........  Open output files...

C.........  Open file(s) for per-source meteorology
        MNAME = 'MINMAXT'
        CALL OPENSMET( ENAME, SDATE, STIME, TVARNAME, MNAME )

C.........  Open file(s) for parameter scheme index outputs
        FNAME = 'MEFTEMP'
        PNAME = CRL // 'SRCPSI'
        CALL OPENPSIOUT( ENAME, FNAME, PNAME, EDEV, PSVNAME )
 
C.........  Process temperature information...

        L = LEN_TRIM( TVARNAME )
        MESG = 'Processing temperature data for using variable "' //
     &         TVARNAME( 1:L ) // '" ...'
        CALL M3MSG2( MESG )

C.........  Loop through days/hours of temperature files
        IDATE = SDATE
        ITIME = 0
        JDATE = SDATE
        JTIME = STIME
        LDATE = -9
        DO T = 1, NSTEPS

C.................  Read current temperature file
            IF ( .NOT. READ3( TNAME, TVARNAME, 1, 
     &                        JDATE, JTIME, TA ) ) THEN
                L = LEN_TRIM( TVARNAME )
                MESG = 'Could not read ' // TVARNAME( 1:L ) //
     &                 ' from ' // TNAME 
               CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )

            END IF

C.............  Apply ungridding matrix 
            CALL APPLUMAT( NSRC, NMATX, TA, UMAT(1), UMAT( NSRC+1 ),
     &                     UMAT( NSRC+NMATX+1 ), TASRC )

C.............  When new day...
            IF ( JDATE .NE. LDATE ) THEN

C.................  Write message for day of week and date
                DAY = WKDAY( JDATE )
                MESG = 'Processing ' // DAYS( DAY ) // MMDDYY( JDATE )
                CALL M3MSG2( MESG )

C.................  Set start and end hours of day for all sources
C.................  The first time this routine is called, ODATE and OTIME
C                   are set as well
                CALL SETSRCDY( NSRC, JDATE, TZONES, LDAYSAV, 
     &                         DAYBEGT, DAYENDT )

            END IF

C.............  First iteration in loop, set the output date/time
            IF( T .EQ. 1 ) THEN
                CALL SETOUTDT( NSRC, JDATE, TZONES, DAYBEGT, DAYENDT, 
     &                         ODATE, OTIME )
            END IF

C.............  Update the min/max temperatures based on this hour's data
            LASTTIME = ( T .EQ. NSTEPS )
            CALL DYMINMAX( NSRC, JTIME, LASTTIME, DAYBEGT, DAYENDT, 
     &                     TASRC, TKMIN, TKMAX, TKMINOUT, TKMAXOUT )
            
C.............  Adjust and output min/max data
            IF( LASTTIME .OR.
     &        ( JDATE .GE. ODATE .AND.
     &          JTIME .EQ. OTIME       ) ) THEN

C.................  Adjust the by-source meteorology data before output
                CALL ADJSMET( NSRC, NTMPR, NVLDTMM, MINT_MIN, MINT_MAX, 
     &                        MAXT_MIN, MAXT_MAX, TMMINVL, TMXINVL,
     &                        'temperature', VLDTMPR, VLDTMIN, VLDTMAX, 
     &                        TKMINOUT, TKMAXOUT, METIDX )

C.................  Count these sources (easier, but unecessary, to repeat)
                DO S = 1, NSRC

                    IF( UMAT( S ) .EQ. 0 ) THEN
                        OFLAG = .TRUE.
                        OSRC = OSRC + 1
                    END IF

                END DO

C.................  Write the by-source meteorology data
                CALL WRSMET( NSRC, IDATE, ITIME, 'MINMAXT', 
     &                       TKMINOUT, TKMAXOUT, METIDX    )

C.................  Update meteorology information for each emission factor
                CALL GENEFMET( NSRC, MXXNPSI, NVLDTMM, NIACT, 
     &                         METIDX, TIPSI )

C.................  Increment output time for per-day file
                CALL NEXTIME( IDATE, ITIME, 240000 )

            END IF   ! End of section to output and store for emission factors

            LDATE = JDATE
            CALL NEXTIME( JDATE, JTIME, TSTEP )

        END DO   !  End loop on hours of temperature files
 
C......... Write temperature combinations for each PSI into an ASCII file

        CALL M3MSG2( 'Writing out EF-REF/TEMPERATURE file...' )

        K = 0
        DO V = 1, NIACT

            DO I = 1, NPSI( V )

                DO J = 1, NVLDTMM

                    K = TIPSI( I, J, V )

                    IF( K .NE. 0 ) THEN

                        PSI  = PSILIST( I,V )
                        TMIN = VLDTMIN( J )
                        TMAX = VLDTMAX( J )

                        WRITE( EDEV,93010 ) 
     &                         PSI, TMIN, TMAX, J, ACTVTY( V )

                    END IF

                END DO
            END DO
        END DO

C.........  Allocate memory for storing PSIs for 24 hours.
        ALLOCATE( SRCPSI( NSRC,24 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SRCPSI', PROGNAME )

        CALL M3MSG2( 'Writing out SOURCE PSIs file...' )

C.........  Transfer PSIs to source-based PSIs list for each activity
        DO V = 1, NIACT

            DO S = 1, NSRC

                K = EFSIDX( S,V )

                DO H = 1, 24
                    SRCPSI( S,H ) = IPSIA( K,H )
                END DO

            END DO             ! End loop on sources

C.............  Write source-PSI file for each of 24 hours for current activity
            IF( .NOT. WRITE3( PNAME, PSVNAME( V ), 0, 0, SRCPSI ) ) THEN

        	MESG = 'Could not write PSIs per source for activity '//
     &                  ACTVTY( V ) // ' to "' //
     &                  PNAME( 1:LEN_TRIM( PNAME ) ) //  '".'
        	CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )

            END IF

        END DO                 ! End loop on activities

C.........  Write message when sources were excluded during ungridding
        IF( OFLAG ) THEN

            WRITE( MESG, 94010 )
     &             'NOTE: During ungridding, ', OSRC,
     &             'sources excluded from grid.' // CRLF() // BLANK10 //
     &             'These sources may be outside the grid.'

            CALL M3MESG( MESG )

        END IF

C......... End program sucessfully

        CALL M3EXIT( PROGNAME, 0, 0, ' ', 0 )
 
C******************  FORMAT  STATEMENTS   ******************************
 
C...........   Formatted file I/O formats............ 93xxx
 
93010   FORMAT( I8, 1X, F13.5, 1X, F13.5, 1X, I8, 1X, A )
 
C...........   Internal buffering formats............ 94xxx
 
94010   FORMAT( 10( A, :, I8, :, 1X ) )

C******************  INTERNAL SUBPROGRAMS  *****************************

        CONTAINS

C.............  This internal subprogram checks the dimension of the grid

            SUBROUTINE CHECK_GRID_DIMS( VNAME, VDESC, UVAL, TVAL )

C.............  Subprogram arguments
            CHARACTER(*), INTENT (IN) :: VNAME   ! name of value
            CHARACTER(*), INTENT (IN) :: VDESC   ! description of value
            INTEGER     , INTENT (IN) :: UVAL    ! ungridding matrix value
            INTEGER     , INTENT (IN) :: TVAL    ! met data file value

C.............  Local variables
            INTEGER       L1, L2

C----------------------------------------------------------------------

            IF( UVAL .NE. TVAL ) THEN
                EFLAG = .TRUE.
                L1 = LEN_TRIM( VNAME )
                L2 = LEN_TRIM( VDESC )
                WRITE( MESG,94010 ) 'ERROR: Inconsistent number of ' //
     &                 VDESC( 1:L2 ) // '. In ungridding matrix, ' // 
     &                 VNAME( 1:L1 ) // '= ', UVAL,
     &                 CRLF() // BLANK5 // 
     &                 'but in meteorology file, ' // VNAME( 1:L1 ) // 
     &                 '= ', TVAL
                CALL M3MSG2( MESG )

            END IF

            RETURN

C------------------- SUBPROGRAM FORMAT STATEMENTS ----------------------

C...........   Internal buffering formats............ 94xxx

94010       FORMAT( 10( A, :, I8, :, 1X ) )

            END SUBROUTINE CHECK_GRID_DIMS

        END PROGRAM PREMOBL



