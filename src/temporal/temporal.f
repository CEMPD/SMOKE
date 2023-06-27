
        PROGRAM TEMPORAL

C***********************************************************************
C
C  DESCRIPTION:
C    This program computes the hourly emissions data from inventoryg emissions
C    and/or activity and emission factor data. It can read average-inventory,
C    day-specific and hour-specific emissions and activity data.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C    copied by: M. Houyoux 01/99
C    origin: tmppoint.F 4.3
C
C     Created by Marc Houyoux
C
C     Revised 07/2014 by C.Coats for  new GENTPRO CSV profiles and cross-references
C
C***********************************************************************
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
C*************************************************************************

C.........  MODULES for public variables
C.........  This module contains the inventory arrays
        USE MODSOURC, ONLY: TZONES, TPFLAG, FLTRDAYL

C.........  This module contains the temporal profile tables
        USE MODTMPRL, ONLY: NMON, NWEK, NHRL, STDATE, STTIME, RUNLEN,
     &                      ITDATE, METPROFLAG, METPROTYPE, IPOL2D, LTFLAG

C.........  This module contains emission factor tables and related
        USE MODEMFAC, ONLY: EMTPOL, NEPOL, TEMPEF, USETIME

C.........  This module contains data for day- and hour-specific data
        USE MODDAYHR, ONLY: DYPNAM, DYPDSC, NDYPOA, NDYSRC,
     &                      HRPNAM, HRPDSC, NHRPOA, NHRSRC,
     &                      LDSPOA, LHSPOA, LHPROF,
     &                      INDXD, EMACD, INDXH, EMACH,
     &                      EMAC, EMACV, EMIST, EMFAC, TMAT

C.........  This module contains the lists of unique source characteristics
        USE MODLISTS, ONLY: NINVIFIP, INVCFIP, MXIDAT, INVDNAM, INVDVTS

C.........  This module contains the information about the source category
        USE MODINFO, ONLY: CATEGORY, BYEAR, NIPPA, EANAM, NSRC,
     &                     INVPIDX, NIPOL, EAREAD, EINAM, ACTVTY

C.........  This module is used for MOBILE setup information
        USE MODMBSET, ONLY: DAILY, WEEKLY, MONTHLY, EPISLEN

        IMPLICIT NONE

C.........  INCLUDES:
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  i/o api parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.
        INCLUDE 'SETDECL.EXT'   !  FileSetAPI variables and functions

C..........  EXTERNAL FUNCTIONS and their descriptions:

        LOGICAL         CHKINT
        CHARACTER(2)    CRLF
        INTEGER         ENVINT
        LOGICAL         ENVYN
        INTEGER         FINDC
        INTEGER         JULIAN
        INTEGER         GETDATE
        INTEGER         GETFLINE
        INTEGER         GETNUM
        INTEGER         INDEX1
        LOGICAL         ISDSTIME
        CHARACTER(14)   MMDDYY
        INTEGER         PROMPTFFILE
        INTEGER         RDTPROF
        INTEGER         SECSDIFF
        INTEGER         STR2INT
        LOGICAL         USEEXPGEO

        EXTERNAL    CHKINT, CRLF, ENVINT, ENVYN, FINDC, JULIAN,
     &              GETDATE, GETFLINE, GETNUM, INDEX1, ISDSTIME, MMDDYY,
     &              PROMPTFFILE, RDTPROF, SECSDIFF, STR2INT, USEEXPGEO

C.........  LOCAL PARAMETERS and their descriptions:

        CHARACTER(50), PARAMETER :: CVSW = '$Name SMOKEv5.0_Jun2023$'  ! CVS revision tag

C.........  Array that contains the names of the inventory variables needed for
C           this program
        CHARACTER(IOVLEN3) IVARNAMS( MXINVARR )

C.........  Day-specific, hour-specific data, and elevated sources data.
C.........  These need only to allow enough dimensions for one read per
C           pollutant per time step.

        INTEGER                 NPELV        ! optional elevated source-count
        INTEGER, ALLOCATABLE :: INDXE( : )   ! SMOKE source IDs
        REAL   , ALLOCATABLE :: EMISE( : )   ! elevated source emissions

C.........  Names of pollutants and activities associated with output variables
        CHARACTER(IOVLEN3), ALLOCATABLE:: ALLIN( : )

C.........  Reshaped input variables and output variables
        INTEGER         NGRP                ! no. of pol/emis-types groups
        INTEGER         NGSZ                ! no. of pols/emis-types per group
        CHARACTER(IOVLEN3), ALLOCATABLE:: ALLIN2D( :,: )
        CHARACTER(IOVLEN3), ALLOCATABLE:: EANAM2D( :,: )
        CHARACTER(IOVLEN3), ALLOCATABLE:: EAREAD2D( : )

C...........   Logical names and unit numbers

        INTEGER      :: CDEV = 0!  unit number for region codes file
        INTEGER      :: HDEV = 0!  unit number for holidays file
        INTEGER         KDEV    !  unit number for time periods file
        INTEGER         LDEV    !  unit number for log file
        INTEGER         PDEV    !  unit number for supplemental tmprl file
        INTEGER         SDEV    !  unit number for ASCII inventory file

        CHARACTER(16) :: ANAME = ' '    !  logical name for ASCII inven input
        CHARACTER(16) :: DNAME = 'NONE' !  day-specific  input file, or "NONE"
        CHARACTER(16) :: ENAME = ' '    !  logical name for I/O API inven input
        CHARACTER(16) :: FNAME = ' '    !  emission factors file
        CHARACTER(16) :: HNAME = 'NONE' !  hour-specific input file, or "NONE"
        CHARACTER(16) :: PNAME = 'NONE' !  hour-specific profile input file, or "NONE"
        CHARACTER(16) :: TNAME = ' '    !  timestepped (low-level) output file

C...........   Other local variables

        INTEGER         I, II, IS, J, K, L, L1, L2, N, S, T

        INTEGER         IOS, IOS1, IOS2, IOS3, IOS4 ! i/o status
        INTEGER         IOS6, IOS7, IOS8, IOS9      ! i/o status
        INTEGER         AVERTYPE            ! time period averaging type
        INTEGER         DYSTPOS, DYENDPOS   ! start and end position in file name string
        INTEGER         EARLYDATE           ! earliest starting date based on time zones
        INTEGER         EARLYTIME           ! earliest starting time based on time zones
        INTEGER         EARLST              ! earliest starting date within entire episode periods
        INTEGER         EDATE, ETIME        ! ending Julian date and time
        INTEGER         EFSDATE, EFEDATE    ! start and end date of current ef file
        INTEGER         ENLEN               ! length of ENAME string
        INTEGER         ENDPOS              ! ending position in ef day array
        INTEGER         FIRSTPOS            ! temporary position in file name string
        INTEGER         FDATE, FTIME        ! emission factor date and time
        INTEGER         HYPPOS              ! position of hyphen in file name string
        INTEGER         JDATE, JTIME        ! Julian date and time
        INTEGER         NDAYS               ! no. days in episode
        INTEGER         JYEAR, JMNTH, JDAYS ! tmp year, month, and date
        INTEGER         LATEDATE            ! latest ending date based on time zones
        INTEGER         LATETIME            ! latest ending time
        INTEGER         LATEST              ! latest starting date within entire episode periods
        INTEGER         NINVARR             ! no. inventory variables to read
        INTEGER         NLINES              ! no. lines in ef list file
        INTEGER         NMATX               ! size of ungridding matrix
        INTEGER         NMAJOR              ! no. major sources
        INTEGER         NPING               ! no. ping sources
        INTEGER         NSTEPS              ! number of output time steps
        INTEGER         NTPERIOD            ! No of time periods
        INTEGER      :: PYEAR = 0           ! projected year
        INTEGER         SDATE, STIME        ! starting Julian date and time
        INTEGER         STPOS               ! starting position in ef day array
        INTEGER         TSTEP               ! output time step
        INTEGER         TZONE               ! output-file time zone
        INTEGER         TZMIN               ! minimum time zone in inventory
        INTEGER         TZMAX               ! maximum time zone in inventory
        INTEGER      :: TDMAX = 0           ! maximum episode days
        INTEGER         LDATE               ! date used in previous subroutine call

        REAL            RTMP                ! tmp float

        LOGICAL      :: DAYLIT    = .FALSE. ! true: TZONES are in daylight time
        LOGICAL         DFLAG               !  true: day-specific  file available
        LOGICAL      :: EFLAG = .FALSE.     !  error-flag
        LOGICAL      :: EFLAG2= .FALSE.     !  error-flag (2)
        LOGICAL         ENDFLAG             !  true: couldn't find file end date
        LOGICAL      :: FNDOUTPUT = .FALSE. ! true: found output hydrocarbon
        LOGICAL         HFLAG               !  true: hour-specific file available
        LOGICAL         NFLAG               !  true: use all uniform temporal profiles
        LOGICAL         PFLAG               !  true: episode time periods needed
        LOGICAL         WFLAG               !  true: write QA on current time step

        CHARACTER(8)         TREFFMT   ! tmprl x-ref format (SOURCE|STANDARD)
        CHARACTER(8)         SCDATE    ! tmprl date
        CHARACTER(14)        DTBUF     ! buffer for MMDDYY
        CHARACTER(3)         INTBUF    ! buffer for integer
        CHARACTER(20)        MODELNAM  ! emission factor model name
        CHARACTER(256)       CURFNM    ! current emission factor file name
        CHARACTER(16)        CURLNM    ! current ef logical file name
        CHARACTER(IOVLEN3)   VOLNAM    ! volatile pollutant name
        CHARACTER(300)       MESG      ! buffer for M3EXIT() messages
        CHARACTER(IOVLEN3)   CBUF      ! pollutant name temporary buffer
        CHARACTER(IOVLEN3)   EBUF      ! pollutant name temporary buffer
        CHARACTER(20)        SEARCHSTR ! string used in search
        CHARACTER(MXDLEN3)   TEMPLINE  ! line from file description

        CHARACTER(16) :: PROGNAME = 'TEMPORAL' ! program name

C***********************************************************************
C   begin body of program TEMPORAL

        LDEV = INIT3()

C.........  Write out copyright, version, web address, header info, and prompt
C           to continue running the program.
        CALL INITEM( LDEV, CVSW, PROGNAME )

C.........  Obtain settings from the environment...

C.........  Get the time zone for output of the emissions
        TZONE = ENVINT( 'OUTZONE', 'Output time zone', 0, IOS )

C.........  Output hourly emission in Local Time.
        MESG = 'Outputs hourly emissions in local time  (No time shift)'
        LTFLAG = ENVYN( 'OUTPUT_LOCAL_TIME', MESG, .FALSE., IOS )

C.........  Get environment variable that overrides temporal profiles and
C               uses only uniform profiles.
        NFLAG = ENVYN( 'UNIFORM_TPROF_YN', MESG, .FALSE., IOS )

C.........  Set source category based on environment variable setting
        CALL GETCTGRY

C.........  Get the name of the emission factor model to use for one run
C        IF ( CATEGORY .EQ. 'MOBILE' ) THEN
C            MESG = 'Emission factor model'
C            CALL ENVSTR( 'SMK_EF_MODEL', MESG, '', MODELNAM, IOS)
C        ELSE
C            MODELNAM = ' '
C        END IF

C.........  Get inventory file names given source category
        CALL GETINAME( CATEGORY, ENAME, ANAME )

C.........  Prompt for and open input files
C.........  Also, store source-category specific information in the MODINFO
C           module.
        CALL OPENTMPIN( NFLAG, PFLAG, ENAME, ANAME, DNAME,
     &                  HNAME, SDEV, CDEV, HDEV, KDEV, PYEAR )

C.........  Determine status of some files for program control purposes
        DFLAG = ( DNAME .NE. 'NONE' )  ! Day-specific emissions
        HFLAG = ( HNAME .NE. 'NONE' )  ! Hour-specific emissions

C.........  Get length of inventory file name
        ENLEN = LEN_TRIM( ENAME )

C.........  Set inventory variables to read for all source categories
        IVARNAMS( 1 ) = 'CIFIP'
        IVARNAMS( 2 ) = 'TZONES'
        IVARNAMS( 3 ) = 'TPFLAG'
        IVARNAMS( 4 ) = 'CSCC'
        IVARNAMS( 5 ) = 'CSOURC'

C.........  Set inventory variables to read for specific source categories
        IF( CATEGORY .EQ. 'AREA' ) THEN
            NINVARR = 5

        ELSE IF( CATEGORY .EQ. 'MOBILE' ) THEN
            NINVARR = 9
            IVARNAMS( 6 ) = 'IRCLAS'
            IVARNAMS( 7 ) = 'IVTYPE'
            IVARNAMS( 8 ) = 'CLINK'
            IVARNAMS( 9 ) = 'CVTYPE'

        ELSE IF( CATEGORY .EQ. 'POINT' ) THEN
            NINVARR = 5

        END IF

C.........  Allocate memory for and read in required inventory characteristics
        CALL RDINVCHR( CATEGORY, ENAME, SDEV, NSRC, NINVARR, IVARNAMS )

C.........  Reset TPFLAG if average day emissions are being used since
C           we don't want to apply the monthly adjustment factors in this case.
        IF ( INVPIDX .EQ. 1 ) THEN
            DO S = 1, NSRC
                IF ( MOD( TPFLAG( S ), MTPRFAC ) .EQ. 0 ) THEN
                    TPFLAG( S ) = TPFLAG( S ) / MTPRFAC
                END IF
            END DO
        END IF

C.........  Build unique lists of SCCs per SIC from the inventory arrays
        CALL GENUSLST

C.........  Define the minimum and maximum time zones in the inventory
        TZMIN = MINVAL( TZONES )
        TZMAX = MAXVAL( TZONES )

C.........  Adjust TZMIN and TZMAX for possibility of daylight savings
        TZMIN = TZMIN - 1
        TZMAX = TZMAX + 1

C.........  Read special files...
C.........  Read region codes file
        IF( USEEXPGEO() ) THEN
            CALL RDGEOCODES( NINVIFIP, INVCFIP )
        ELSE
            CALL RDSTCY( CDEV, NINVIFIP, INVCFIP )
        END IF

C.........  Populate filter for sources that use daylight time
        CALL SETDAYLT

C.........  Output hourly emissions in local time
        IF( LTFLAG ) THEN
            FLTRDAYL = 0         ! reset daily light saving value to zero
            TZONE  = 0           ! reset output time zone to zero
            TZONES = 0           ! reset county-specific time zone to zero
            TZMIN  = 0
            TZMAX  = 0
        END IF

C.........  Read holidays file
        CALL RDHDAYS( HDEV, EARLST, LATEST )

C.........  Determine all of the variables to be output based on the
C           activities and input pollutants.
        CALL TMNAMUNT

C.........  Reset the number of all output variables as the number of pollutants
C           and emission types, instead of the number of pollutants and
C           activities
        NIPPA = NIPOL

C.........  Allocate memory for I/O pol names, activities, & emission types
C.........  Will be resetting EANAM to include the emission types instead
C           of the activities
        DEALLOCATE( EANAM )
        ALLOCATE( EANAM( NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EANAM', PROGNAME )
        ALLOCATE( ALLIN( NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ALLIN', PROGNAME )
        ALLOCATE( EAREAD2D( NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EAREAD2D', PROGNAME )
        EAREAD2D = ' '

C.........  Create 1-d arrays of I/O pol names
C.........  If pollutant is created by one or more activities, then give a
C           warning.
        N = 0
        DO I = 1, NIPOL

C.............  Look for pollutant in list of pollutants created by activities
C               First, make sure that EMTPOL has been allocated; it won't be
C               when not using emission factors
            IF( ALLOCATED( EMTPOL ) ) THEN
                J = INDEX1( EINAM( I ), NEPOL, EMTPOL )
            ELSE
                J = 0
            END IF

C.............  If pollutant created by activity, skip from this list, unless
C               pollutant is also part of the inventory pollutants
            IF( J .GT. 0 ) THEN
                L1 = LEN_TRIM( EINAM ( I ) )
                L2 = LEN_TRIM( EAREAD( I ) )
                MESG = 'WARNING: Pollutant "' // EINAM(I)( 1:L1 ) //
     &                 '" is explicitly in the inventory and' //
     &                 CRLF() // BLANK10 // 'it is also generated by '
     &                 // 'activity data.'
                CALL M3MSG2( MESG )
            END IF

            N = N + 1
            ALLIN( N ) = EAREAD( I )
            EANAM( N ) = EINAM ( I )

        END DO

C.........  Reset number of pollutants and emission types based on those used

        NIPPA = N

C.........  Read episode time period lists from PROCDATES.txt
C.........  or set NTPERIOD=1 and allocate RUNLEN(:)
C.........  Get episode settings from the Models-3 environment variables
C.........  when $GE_DAT/procdates.txt is not available for episode time periods

        CALL RDDATES( PFLAG, KDEV, NTPERIOD )

        IF ( .NOT.PFLAG ) THEN
            SDATE  = 0
            STIME  = 0
            NSTEPS = 1
            CALL GETM3EPI( TZONE, SDATE, STIME, TSTEP, NSTEPS )

            TSTEP = 10000           ! Only 1-hour time steps supported
            JYEAR = SDATE / 1000
            CALL DAYMON( SDATE, JMNTH, JDAYS )
            STDATE(1) = JDAYS + 100*( JMNTH + 100*JYEAR )
            STTIME(1) = STIME
            RUNLEN(1) = NSTEPS
            ITDATE(1) = SDATE
        END IF

        DO II = 1, NTPERIOD

            NDAYS = 0

C.............  Get episode settings from episode time periods file
            IF( PFLAG ) THEN

                WRITE( SCDATE,'(I8)' ) STDATE( II )
                JYEAR = STR2INT( SCDATE( 1:4 ) )
                JMNTH = STR2INT( SCDATE( 5:6 ) )
                JDAYS = STR2INT( SCDATE( 7:  ) )
                SDATE = JULIAN( JYEAR, JMNTH, JDAYS )
                ITDATE( II ) = JYEAR*1000 + SDATE

                SDATE = ITDATE( II )
                STIME = STTIME( II )
                TSTEP  = 10000  ! Only 1-hour time steps supported
                NSTEPS = RUNLEN ( II ) / TSTEP
             END IF

C............  Determine number of days in episode
             IF( II == 1 ) THEN
                 EARLST = SDATE
                 LATEST = SDATE
             END IF

C............  Earliest day is start time in maximum time zone
             EARLYDATE = SDATE
             EARLYTIME = STIME
             CALL NEXTIME( EARLYDATE, EARLYTIME,
     &                    -( TZMAX - TZONE )*10000 )

C............  Latest day is end time in minimum time zone
C............  Calculate the ending date and time
             EDATE = SDATE
             ETIME = STIME
             CALL NEXTIME( EDATE, ETIME, NSTEPS * 10000 )

             LATEDATE = EDATE
             LATETIME = ETIME
             CALL NEXTIME( LATEDATE, LATETIME,
     &                    -( TZMIN - TZONE )*10000 )

             NDAYS = SECSDIFF( EARLYDATE, 0, LATEDATE, 0 ) / ( 24*3600 )
             NDAYS = NDAYS + 1
             IF( NDAYS > TDMAX ) TDMAX = NDAYS

C..............  Check earliest and latest dates during episode
             IF( EARLYDATE < EARLST ) EARLST = EARLYDATE
             IF( LATEDATE  > LATEST ) LATEST = LATEDATE

        END DO ! end of entire time episode periodes loops


C.........  Read temporal-profile cross-reference file and profiles, and
C.........  put into tables
C.........  Only read entries for pollutants that are in the inventory.
C.........  Only read if not using uniform temporal profiles
C.........  Assign temporal profiles to sources.

        CALL PROCTPRO( NFLAG, METPROFLAG, PNAME )

        IF ( .NOT.NFLAG ) CALL ASGNTPRO()

C.........  Check requested episode against available emission factors

C.........  For day-specific data input...
        IF( DFLAG ) THEN

C.............  Get header description of day-specific input file
            IF( .NOT. DESC3( DNAME ) ) THEN
                CALL M3EXIT( PROGNAME, 0, 0,
     &                       'Could not get description of file "'
     &                       // DNAME( 1:LEN_TRIM( DNAME ) ) // '"', 2 )
            END IF

C.............  Allocate memory for pollutant pointer
            ALLOCATE( DYPNAM( NVARS3D ), STAT=IOS )
            CALL CHECKMEM( IOS, 'DYPNAM', PROGNAME )
            ALLOCATE( DYPDSC( NVARS3D ), STAT=IOS )
            CALL CHECKMEM( IOS, 'DYPDSC', PROGNAME )
            DYPNAM = ' '  ! array
            DYPDSC = ' '  ! array

CC.............  Set day-specific file dates, check dates, and report problems
            CALL PDSETUP( DNAME, EARLST, STIME, LATEST, ETIME, TZONE,
     &                    NIPPA, EANAM, NDYPOA, NDYSRC, EFLAG, DYPNAM,
     &                    DYPDSC )

        ENDIF

C.........  Allocate memory for reading day-specific emissions data
C.........  NDYSRC is initialized to zero in case DFLAG is false
        ALLOCATE( INDXD( NDYSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INDXD', PROGNAME )
        ALLOCATE( EMACD( NDYSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMACD', PROGNAME )

C.........  For hour-specific data input...
        IF( HFLAG ) THEN

C............. Get header description of hour-specific input file
            IF( .NOT. DESC3( HNAME ) ) THEN
                CALL M3EXIT( PROGNAME, 0, 0,
     &                       'Could not get description of file "'
     &                       // HNAME( 1:LEN_TRIM( HNAME ) ) // '"', 2 )
            ENDIF

C.............  Allocate memory for pollutant pointer
            ALLOCATE( HRPNAM( NVARS3D ), STAT=IOS )
            CALL CHECKMEM( IOS, 'HRPNAM', PROGNAME )
            ALLOCATE( HRPDSC( NVARS3D ), STAT=IOS )
            CALL CHECKMEM( IOS, 'HRPDSC', PROGNAME )
            HRPNAM = ' '  ! array
            HRPDSC = ' '  ! array

C.............  Set day-specific file dates, check dates, and report problems
            CALL PDSETUP( HNAME, EARLST, STIME, LATEST, ETIME, TZONE,
     &                    NIPPA, EANAM, NHRPOA, NHRSRC, EFLAG2, HRPNAM,
     &                    HRPDSC )

        ENDIF

C.........  Allocate memory for reading hour-specific emissions data
C.........  NHRSRC is initialized to 0 in case HFLAG is false
        ALLOCATE( INDXH( NHRSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INDXH', PROGNAME )
        ALLOCATE( EMACH( NHRSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMACH', PROGNAME )

        IF( EFLAG .OR. EFLAG2 ) THEN
            MESG = 'Problem with day- or hour-specific inputs'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ENDIF

C......  Get episode settings from episode time periods file

      DO II = 1, NTPERIOD

        LDATE = -1     ! initializing previous date

C.........  Initialize EAREAD2D every processing date
        DO I = 1, NIPPA
            EAREAD2D( I ) = EAREAD( I )
        END DO

C.........  It is important that all major arrays must be allocated by this
C           point because the next memory allocation step is going to pick a
C           data structure that will fit within the limits of the host.

C.........  Allocate memory, but allow flexibility in memory allocation
C           for second dimension.
C.........  The second dimension (the number of pollutants and emission types)
C            can be different depending on the memory available.
C.........  To determine the approproate size, first attempt to allocate memory
C           for all pollutants & emission types to start, and if this fails,
C           divide pollutants into even groups and try again.

        NGSZ = NIPPA            ! No. of pollutant & emis types in each group
        NGRP = 1                ! Number of groups

        DO

            ALLOCATE( TMAT ( NSRC, NGSZ, 24 ),
     &                EMAC ( NSRC, NGSZ )    ,
     &                EMACV( NSRC, NGSZ )    ,
     &                EMIST( NSRC, NGSZ )    ,
     &                EMFAC( NSRC, NGSZ )    , STAT=IOS9 )

            IF( IOS9 .GT. 0 ) THEN

                IF( NGSZ .EQ. 1 ) THEN
                    J = 8 * NSRC * 31    ! Assume 8-byte reals
                    WRITE( MESG,94010 )
     &                'Insufficient memory to run program.' //
     &                CRLF() // BLANK5 // 'Could not allocate ' //
     &                'pollutant-dependent block of', J, 'bytes.'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF

                NGRP = NGRP + 1
                NGSZ = ( NIPPA + NGRP - 1 ) / NGRP

                IF( ALLOCATED( TMAT  ) ) DEALLOCATE( TMAT )
                IF( ALLOCATED( EMAC  ) ) DEALLOCATE( EMAC )
                IF( ALLOCATED( EMACV ) ) DEALLOCATE( EMACV )
                IF( ALLOCATED( EMIST ) ) DEALLOCATE( EMIST )
                IF( ALLOCATED( EMFAC ) ) DEALLOCATE( EMFAC )

            ELSE
                EXIT

            END IF

        END DO

        WRITE( MESG, '( 2( A, I4, 2X ) )' )
     &     'Processing', NGRP, 'pollutant groups of size', NGSZ
        CALL M3MESG( MESG )

C.........  Allocate a few small arrays based on the size of the groups
C.........  NOTE that this has a small potential for a problem if these little
C           arrays exceed the total memory limit.
        IF( ALLOCATED( ALLIN2D ) ) DEALLOCATE( ALLIN2D )
        IF( ALLOCATED( EANAM2D ) ) DEALLOCATE( EANAM2D )
        IF( ALLOCATED( IPOL2D  ) ) DEALLOCATE( IPOL2D  )
        IF( ALLOCATED( LDSPOA  ) ) DEALLOCATE( LDSPOA  )
        IF( ALLOCATED( LHSPOA  ) ) DEALLOCATE( LHSPOA  )
        IF( ALLOCATED( LHPROF  ) ) DEALLOCATE( LHPROF  )

        ALLOCATE( ALLIN2D( NGSZ, NGRP ),
     &            EANAM2D( NGSZ, NGRP ),
     &             IPOL2D( NGSZ, NGRP ),
     &             LDSPOA( NGSZ ),
     &             LHSPOA( NGSZ ),
     &             LHPROF( NGSZ ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ALLIN2D...LHPROF', PROGNAME )

C.........  Create 2-d arrays of I/O pol names, activities, & emission types

        I = 0
        DO N = 1, NGRP
        DO J = 1, NGSZ
            I = I + 1
            IF ( I .LE. NIPPA ) THEN
                IPOL2D ( J,N ) = I
                EANAM2D( J,N ) = EANAM( I )
                ALLIN2D( J,N ) = ALLIN( I )
            ELSE
                IPOL2D ( J,N ) = 0
                EANAM2D( J,N ) = ' '
                ALLIN2D( J,N ) = ' '
            END IF
        END DO
        END DO

C......  Determine number of days in episode

        IF( PFLAG ) THEN
           SDATE  = ITDATE( II )
           STIME  = STTIME( II )
           TSTEP  = 10000  ! Only 1-hour time steps supported
           NSTEPS = RUNLEN ( II ) / TSTEP
        END IF

C.........  Earliest day is start time in maximum time zone
        EARLYDATE = SDATE
        EARLYTIME = STIME
        CALL NEXTIME( EARLYDATE, EARLYTIME,
     &               -( TZMAX - TZONE )*10000 )

C.........  Latest day is end time in minimum time zone
        EDATE = SDATE
        ETIME = STIME
        CALL NEXTIME( EDATE, ETIME, NSTEPS * 10000 )

        LATEDATE = EDATE
        LATETIME = ETIME
        CALL NEXTIME( LATEDATE, LATETIME,
     &               -( TZMIN - TZONE )*10000 )

        NDAYS = SECSDIFF( EARLYDATE, 0, LATEDATE, 0 ) / ( 24*3600 )
        NDAYS = NDAYS + 1

C.....  Allocate memory for emission factor arrays
        IF( ALLOCATED( TEMPEF ) ) DEALLOCATE ( TEMPEF )
        ALLOCATE( TEMPEF( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TEMPEF', PROGNAME )

        TEMPEF = 0.

C.........  Compare base year with episode and warn if not consistent
        IF( SDATE / 1000 .NE. BYEAR ) THEN

            WRITE( MESG,94010 ) 'WARNING: Inventory base year ', BYEAR,
     &             'is inconsistent with year ' // CRLF() // BLANK10 //
     &             'of episode start date', SDATE/1000
            CALL M3MSG2( MESG )

        ENDIF

C.........  Set up and open I/O API output file(s) ...
        CALL OPENTMP( II, ENAME, SDATE, STIME, TSTEP, NSTEPS, TZONE,
     &                NPELV, TNAME, PDEV, PFLAG )

C.........  Give a note if running for a projected year
        IF( PYEAR .GT. 0 ) THEN

            WRITE( MESG,94010 ) 'NOTE: Emissions based on projected '//
     &             'year', PYEAR
            CALL M3MSG2( MESG )

        END IF

C.........  Loop through pollutant/emission-type groups
        DO N = 1, NGRP

C.............  If this is the last group, NGSZ may be larger than actual
C               number of variables, so reset based on total number
            IF( N == NGRP ) THEN
                NGSZ = NIPPA - ( NGRP - 1 )*NGSZ
            END IF

C.............  Skip group if the first pollutant in group is blank (this
C               shouldn't happen, but it is happening, and it's easier to
C               make this fix).
            IF ( EANAM2D( 1,N ) .EQ. ' ' ) CYCLE

C.............  Write message stating the pols/emission-types being processed
            CALL POLMESG( NGSZ, EANAM2D( 1,N ) )

C.............  Set up logical arrays that indicate which pollutants/activities
C               are day-specific and which are hour-specific.
C.............  Also set flag for which hour-specific pollutants/activities
C               are actually diurnal profiles instead of emissions
            LDSPOA = .FALSE.   ! array
            DO I = 1, NDYPOA
                J = INDEX1( DYPNAM( I ), NGSZ,  EANAM2D( 1,N ) )
                IF ( J /= 0 ) LDSPOA( J ) = .TRUE.
            END DO

            LHSPOA = .FALSE.   ! array
            LHPROF = .FALSE.   ! array
            DO I = 1, NHRPOA
                J = INDEX1( HRPNAM( I ), NGSZ,  EANAM2D( 1,N ) )
                IF ( J /= 0 ) LHSPOA( J ) = .TRUE.

                CALL UPCASE( HRPDSC( I ) )
                K = INDEX( HRPDSC( I ), 'PROFILE' )
                IF( K .GT. 0 ) LHPROF( J ) = .TRUE.
            END DO

C.............  Initialize emissions, activities, and other arrays for this
C.............  pollutant/emission-type group

            TMAT  = 0.
            EMAC  = 0.
            EMACV = 0.
            EMIST = 0.

C.............  Read in pollutant emissions or activities from inventory for
C               current group

            DO I = 1, NGSZ

                EBUF = EANAM2D( I,N )
                CBUF = ALLIN2D( I,N )
                L1   = LEN_TRIM( CBUF )

C.................  Skip blanks that can occur when NGRP > 1
                IF ( CBUF .EQ. ' ' ) CYCLE

C.................  Read the emissions data in either map format
C                   or the old format.
                CALL RDMAPPOL( NSRC, 1, 1, CBUF, EMAC( 1,I ) )

C...............  If there are any missing values in the data, give an
C                 error to avoid problems in genhemis routine
                RTMP = MINVAL( EMAC( 1:NSRC,I ) )
                IF( RTMP .LT. 0 ) THEN
                    EFLAG = .TRUE.
                    MESG = 'ERROR: Missing or negative emission '//
     &                     'value(s) in inventory for "' //
     &                     CBUF( 1:L1 ) // '".'
                    CALL M3MSG2( MESG )
                END IF

C...............  If pollutant name is average day-based, remove the
C                 prefix from the input pollutant name
                K = INDEX1( CBUF, NIPPA, EAREAD )
                J = INDEX( CBUF, AVEDAYRT )
                IF( J .GT. 0 ) THEN
                    CBUF = CBUF( CPRTLEN3+1:L1 )
                    ALLIN2D( I,N ) = CBUF
                    EAREAD2D ( K ) = CBUF
                END IF

            END DO

C.............  Abort if error found
            IF( EFLAG ) THEN
                MESG = 'Problem with input data.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  For each time step and pollutant or emission type in current
C               group, generate hourly emissions, and write layer-1 emissions
C               file (or all data).
            JDATE = SDATE
            JTIME = STIME

            DO T = 1, NSTEPS        !  Loop through time steps for current pollutant group

                WRITE( MESG, '( A, I4, 2X, A, I9.7, A, I6.6 )' ) 
     &                'Processing pol-group', N, 'for', JDATE, ':', JTIME
                CALL M3MESG( MESG )

C.................  Adjust sources' time zones to account for daylight time...
C                   Subtract 1 if date is daylight time and TZONES is not already
C                   converted.  Add 1 if date is standard and TZONES has been converted.
C                   FLTRDAYL is a source-array of 0s and 1s to permit sources
C                   to not get daylight time conversion.

                IF( ISDSTIME( JDATE ) .AND. .NOT. DAYLIT ) THEN

                    DAYLIT = .TRUE.
                    TZONES = TZONES - 1 * FLTRDAYL   ! arrays

                ELSE IF( .NOT. ISDSTIME( JDATE ) .AND. DAYLIT ) THEN

                    DAYLIT = .FALSE.
                    TZONES = TZONES + 1 * FLTRDAYL   ! arrays

                END IF

C.................  Generate and write hourly emissions for current hour

                CALL GENHEMIS( N, NGRP, NGSZ, JDATE, JTIME, TZONE, DNAME, HNAME,
     &                  PNAME, ALLIN2D( 1,N ), EANAM2D( 1,N ), EAREAD2D,
     &                  LDATE )

                DO I = 1, NGSZ      !  Loop through pollutants/emission-types in this group
                                    !  Skip blanks that can occur when NGRP > 1
                    CBUF = EANAM2D( I,N )
                    IF ( CBUF .EQ. ' ' ) CYCLE

                    IF( .NOT. WRITESET( TNAME, CBUF, ALLFILES, JDATE,
     &                                  JTIME, EMIST( 1,I ) )    ) THEN
                        MESG = 'Could not write "' // TRIM( CBUF ) //
     &                         '" to file "' // TRIM( TNAME ) // '"'
                        CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )

                    END IF

                END DO  ! End loop on pollutants/emission-types I in this group

C.................  Advance the output date/time by one time step

                CALL NEXTIME( JDATE, JTIME, TSTEP )

            END DO      ! End loop on time steps T

         END DO         ! End loop on pollutant groups N

C.............  Write supplemental temporal profiles file

         CALL WRTSUP( PDEV, NSRC, NIPPA, EANAM )

         IF( .NOT. CLOSESET( TNAME ) ) THEN
             MESG = 'Could not close file "' // TRIM( TNAME ) // '".'
             CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
         END IF

         DEALLOCATE( TMAT, EMAC, EMACV, EMIST, EMFAC )

       END DO            ! End loop on time period II

       CALL M3EXIT( PROGNAME, 0, 0, ' ', 0 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats.............94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94050   FORMAT( A, 1X, I2.2, A, 1X, A, 1X, I6.6, 1X,
     &          A, 1X, I3.3, 1X, A, 1X, I3.3, 1X, A   )

        END PROGRAM TEMPORAL

