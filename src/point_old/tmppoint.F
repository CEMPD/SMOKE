
C Version "@(#)$Id$ $Source$ $Date$ 

        PROGRAM TMPPOINT

C***********************************************************************
C  program body starts at line 258
C
C  DESCRIPTION:
C       Produces source-level inventory-speciated temporally allocated
C       emissions, using the emissions-source vector produced by RAWPOINT
C       and the temporal profiles and cross-reference tables for point
C       sources (extracted as the points-only portions of the AIRS tables)
C
C  PRECONDITIONS REQUIRED:  
C       M3IO source-level point source emissions in canonical order
C       produced by RAWPOINT
C       Sorted, cut-down, only-active-point-sources XREF file.
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       Models-3 I/O; 
C       FIND2, FIND3, GETNUM, GETYN, DAYMON, TIME2SEC, TRIMLEN, WKDAY
C
C  REVISION  HISTORY:
C       Prototype  7/95 by CJC.
C
C       Modified  10/95 by CJC -- proper handling of time zones, 
C       daylight saving time, annualized vs. day-specific inventories;
C       uses TPFLAG, to determine source-dependent applicability of 
C       monthly, weekly profiles (applicable iff 
C       MOD( TPFLAG( source ), [M|W]PRFAC ) is zero
C       Optionally:  uses day-specific, hour-specific data from DAYPOINT,
C       uses elevated-sources file to generate time-stepped elevated
C       sources emissions file for input to plume rise.
C
C***********************************************************************

      IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'PTDIMS3.EXT'   !  point-source dimensioning parameters
        INCLUDE 'CHDIMS3.EXT'   !  emis chem parms (inventory + model)
        INCLUDE 'TMDIMS3.EXT'   !  emis temporal modeling parms
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.


C...........   EXTERNAL FUNCTIONS and their descriptions:

        INTEGER         FIND1, FIND2, FIND3, FIND4
        INTEGER         GETNUM
        LOGICAL         GETYN
        LOGICAL         ISDSTIME
        INTEGER         INDEX1
        INTEGER         LBLANK
        INTEGER         PROMPTFFILE
        CHARACTER*16    PROMPTMFILE
        INTEGER         RDTPROF
        INTEGER         STR2INT
        INTEGER         TIME2SEC
        INTEGER         TRIMLEN
        INTEGER         WKDAY           !  day of week (1...7)

      EXTERNAL      FIND1, FIND2, FIND3, FIND4, GETNUM, GETYN,
     &              INDEX1, ISDSTIME, LBLANK, PROMPTFFILE, PROMPTMFILE,
     &              RDTPROF, STR2INT, TIME2SEC, TRIMLEN, WKDAY
     
                                            
C...........   PARAMETERS and their descriptions

        INTEGER     MXREC
        PARAMETER ( MXREC = NPSRC * NIPOL )


C...........   LOCAL VARIABLES and their descriptions:
C...........   Point Sources input and output arrays

        INTEGER     IFIP  ( NPSRC )  !  source FIPS (county) ID
        INTEGER     ISCC  ( NPSRC )  !  source SCC
        INTEGER     IPLANT( NPSRC )  !  Plant ID
        INTEGER     ISTACK( NPSRC )  !  Stack ID
        INTEGER     TZONES( NPSRC )  !  time zones
        INTEGER     TPFLAG( NPSRC )  !  governs month, week profile use
        REAL        EMIS  ( NPSRC, NIPOL )  !  input inventory emissions.
        REAL        EMISV ( NPSRC, NIPOL )  !  day-corrected emissions.
        REAL        EMIST ( NPSRC, NIPOL )  !  timestepped output emssions
        INTEGER     MDEX  ( NPSRC, NIPOL )  ! monthly       profile subscript
        INTEGER     WDEX  ( NPSRC, NIPOL )  ! weeky         profile subscript
        INTEGER     DDEX  ( NPSRC, NIPOL )  ! wkday diurnal profile subscript
        INTEGER     EDEX  ( NPSRC, NIPOL )  ! wkend-diurnal profile subscript

        INTEGER     MREF  ( NIPOL )  ! monthly       xref
        INTEGER     WREF  ( NIPOL )  ! weeky         xref
        INTEGER     DREF  ( NIPOL )  ! wkday diurnal xref


C...........   Day-specific, hour-specific, elevated point source data:

        INTEGER       NDAYPT            !  current source-count
        INTEGER       INDXD( NPDAY )    !  subscript into netCDF inventory file
        REAL          EMISD( NPDAY, NIPOL )  !  emissions values (tons/yr)

        COMMON / DAYREC / NDAYPT, INDXD, EMISD

        INTEGER       NHRPT            !  current source-count
        INTEGER       INDXH( NPHR )    !  subscript into netCDF inventory file
        REAL          EMISH( NPHR, NIPOL )  !  emissions values (tons/yr)

        COMMON / HOURREC / NHRPT, INDXH, EMISH

        INTEGER       NELPT            !  current source-count
        INTEGER       INDXE( NPELV )    !  subscript into netCDF inventory file
        REAL          EMISE( NPELV, NIPOL )  !  emissions values (tons/yr)

        COMMON / ELEVREC / NELPT, INDXE, EMISE

C...........   Time-stepping Matrix.  
C.......   TMAT( TNDX(S,C), H ) is profile coeffs for source S, hour H
C.......   TMAT( C,H ) gives the coefficients for the effective 
C.......   temporal profile 7*(12*MONTH-code# + WKDAY-code#) + DIURNAL-code#

        REAL            TMAT( NPSRC, NIPOL, 24 )    ! coeff profiles


C.......   Actual-SCC  table

        INTEGER      NSCC
        INTEGER      ASCC( NPSRC )  !  actutally-occurring SCC table


C...........   Temporal Profile XREF tables:  give monthly, weekly, hourly 
C.......   profile-indexes for lookup into temporal profile tables.
C.......   Before sorting:

        INTEGER         NPRNA
        INTEGER         INDXA( MXREC )
        INTEGER         IFIPA( MXREC )
        INTEGER         ISCCA( MXREC )
        INTEGER         IPLTA( MXREC )
        INTEGER         ISTKA( MXREC )
        INTEGER         ICODA( MXREC )
        INTEGER         MPRNA( MXREC )
        INTEGER         WPRNA( MXREC )
        INTEGER         DPRNA( MXREC )

C...........   After sorting:
C.......   Default no-FIP, no-ACC profile

        INTEGER         MPR00( NIPOL )  !  month-in-year
        INTEGER         WPR00( NIPOL )  !  day-in-week
        INTEGER         DPR00( NIPOL )  !  hour-in-day (diurnal)

C.......   No-FIP, by-SCC profiles

        INTEGER         NPRN1
        INTEGER         ISCC1( NPSRC )
        INTEGER         MPRN1( NPSRC, NIPOL )
        INTEGER         WPRN1( NPSRC, NIPOL )
        INTEGER         DPRN1( NPSRC, NIPOL )

C.......   State-FIP, by-SCC profiles

        INTEGER         NPRN2
        INTEGER         IFIP2( NPSRC )
        INTEGER         ISCC2( NPSRC )
        INTEGER         MPRN2( NPSRC, NIPOL )
        INTEGER         WPRN2( NPSRC, NIPOL )
        INTEGER         DPRN2( NPSRC, NIPOL )

C.......   By-FIP,SCC profiles

        INTEGER         NPRN3
        INTEGER         IFIP3( NPSRC )
        INTEGER         ISCC3( NPSRC )
        INTEGER         MPRN3( NPSRC, NIPOL )
        INTEGER         WPRN3( NPSRC, NIPOL )
        INTEGER         DPRN3( NPSRC, NIPOL )

C.......   By-FIP,SCC,PLANT profiles

        INTEGER         NPRN4
        INTEGER         IFIP4( NPSRC )
        INTEGER         ISCC4( NPSRC )
        INTEGER         IPLT4( NPSRC )
        INTEGER         MPRN4( NPSRC, NIPOL )
        INTEGER         WPRN4( NPSRC, NIPOL )
        INTEGER         DPRN4( NPSRC, NIPOL )

C.......   By-FIP,SCC,PLANT,STACK profiles

        INTEGER         NPRN5
        INTEGER         IFIP5( NPSRC )
        INTEGER         ISCC5( NPSRC )
        INTEGER         IPLT5( NPSRC )
        INTEGER         ISTK5( NPSRC )
        INTEGER         MPRN5( NPSRC, NIPOL )
        INTEGER         WPRN5( NPSRC, NIPOL )
        INTEGER         DPRN5( NPSRC, NIPOL )


C...........   Temporal Profile arrays

        INTEGER         MONREF( NMCOD )     ! Monthly         XREF table
        INTEGER         WEKREF( NWCOD )     ! Weekly          XREF table
        INTEGER         DIUREF( NDCOD )     ! Weekday-diurnal XREF table
        INTEGER         ENDREF( NDCOD )     ! Weekend-diurnal XREF table

        REAL            MONFAC( 12, NMCOD )     ! Monthly
        REAL            WEKFAC(  7, NWCOD )     ! Weekly
        REAL            DIUFAC( 24, NDCOD )     ! Weekday-diurnal
        REAL            ENDFAC( 24, NDCOD )     ! Weekend-diurnal

C.......   Numbers of month, week, day, weekday, weekend profiles read

        INTEGER         NMON, NWEK, NDIU, NEND, M, W, D, H

C...........   Other local variables

        LOGICAL         EFLAG   !  error-flag
        LOGICAL         WFLAG   !  use weekly  profiles
        LOGICAL         MFLAG   !  use monthly profiles
        LOGICAL         DFLAG   !  day-specific  file available
        LOGICAL         HFLAG   !  hour-specific file available
        LOGICAL         UFLAG   !  generating upper-level output file
        LOGICAL         YFLAG   !  generating output xref file
        INTEGER         TZONE   !  output-file time zone
        INTEGER         TSTEP   !  output time step
        INTEGER         JDATE, JTIME
        INTEGER         NSTEPS
        INTEGER         CDEV    !  for actual-SCC file
        INTEGER         LDEV    !  unit number for log file
        INTEGER         XDEV    !  ...input cross-reference file
        INTEGER         PDEV    !  ...input profile file
        INTEGER         YDEV    !  ...optional    output cross-reference file
        INTEGER         UDEV    !  ...optional input elevated sources file
        INTEGER         IOS
        INTEGER         IREC
        INTEGER         C0, C1
        INTEGER         I, J, K, F, S, T, V
        INTEGER         FIP, SCC, PLT, STK, LID
        INTEGER         COD
        INTEGER         IMON, IWEK, IDIU
        INTEGER         LS, EC
        INTEGER         MONTH, DAY, LDATE

        CHARACTER*16    ENAME   !  logical name for point-source input file
        CHARACTER*16    DNAME   !  " day-specific  input file, or "NONE"
        CHARACTER*16    HNAME   !  " hour-specific input file, or "NONE"
        CHARACTER*16    TNAME   !  " timestepped (low-level)   output file
        CHARACTER*16    UNAME   !  " (upper-level) output file, or "NONE"
        CHARACTER*80    LINE    !  buffer for XREF input lines
        CHARACTER*256   MESG    !  buffer for M3EXIT() messages


C***********************************************************************
C   begin body of program TMPPOINT

        LDEV = INIT3()

        WRITE( *,92000 ) 
     &' ',
     &'Program TMPPOINT to take the point source emissions file from',
     &'RAWPOINT, the ACTUAL SCC file, the ASCII temporal profile and',
     &'cross-reference file, optional DAY-SPECIFIC and HOUR-SPECIFIC',
     &'point source file from PDPOINT, and produce temporally',
     &'allocated (time-stepped lower level) and optionally upper-level',
     &'emissions for the requested time period.',
     &' '
        WRITE( *,92000 ) 
     &'You will need to enter the logical names for the input and',
     &'output files (and to have set them prior to program launch,',
     &'using "setenv <logicalname> <pathname>").  You will be ',
     &'prompted for which species to allocate -- these should only ',
     &'be the chemical species, and not the FIP-related or ASC-related',
     &'variables.)', 
     &  ' '
        WRITE( *,92000 ) 
     &'You will also be prompted for the TIME PERIOD to be covered',
     &'by this run of the program, as specified by STARTING DATE',
     &'in format YYYYDDD (= 1000*year + Julian day n0 (1...365,6) )',
     &'STARTING TIME in format HHMMSS (= 10000*hr + 100*min + sec) )',
     &'and DURATION in format HHMMSS.',
     &  ' '
        WRITE( *,92000 ) 
     &'You may use END_OF-FILE (control-D) to quit the program',
     &'during logical-name entry.  Default responses are indicated',
     &'in brackets [LIKE THIS].',
     &' '

        IF ( .NOT. GETYN( 'Continue with program?', .TRUE. ) ) THEN
            CALL M3EXIT( 'TMPPOINT', 0, 0, 'Ending program', 2 )
        END IF

        DO  2  M = 1, NMCOD
        DO  1  I = 1, 12
            MONFAC( I,M ) = 1.0
1       CONTINUE
2       CONTINUE

        DO  4  W = 1, NWCOD
        DO  3  I = 1, 7
            WEKFAC( I,W ) = 1.0
3       CONTINUE
4       CONTINUE

        DO  6  D = 1, NDCOD
        DO  5  I = 1, 24
            DIUFAC( I,D ) = 1.0
            ENDFAC( I,D ) = 1.0
5       CONTINUE
6       CONTINUE

        DO  9  I = 1, NIPOL
            MPR00( I ) = IMISS3
            WPR00( I ) = IMISS3
            DPR00( I ) = IMISS3
            DO  7  S = 1, NPSRC
                MPRN1( S,I ) = IMISS3
                WPRN1( S,I ) = IMISS3
                DPRN1( S,I ) = IMISS3
                MPRN2( S,I ) = IMISS3
                WPRN2( S,I ) = IMISS3
                DPRN2( S,I ) = IMISS3
                MPRN3( S,I ) = IMISS3
                WPRN3( S,I ) = IMISS3
                DPRN3( S,I ) = IMISS3
                MPRN4( S,I ) = IMISS3
                WPRN4( S,I ) = IMISS3
                DPRN4( S,I ) = IMISS3
                MPRN5( S,I ) = IMISS3
                WPRN5( S,I ) = IMISS3
                DPRN5( S,I ) = IMISS3
7           CONTINUE
            DO  8  S = 1, NPELV
                EMISE( S,I ) = 0.0
8           CONTINUE
9       CONTINUE


C.......   Get file name; open input point sources, temporal cross-reference,
C.......   and temporal profiles files

        ENAME = PROMPTMFILE( 
     &          'Enter name for POINT SOURCE input file',
     &          FSREAD3, 'PNTS', 'TMPPOINT' )

        DNAME = PROMPTMFILE( 
     &          'Enter name for DAY-SPECIFIC file, or "NONE"',
     &          FSREAD3, 'PDAY', 'TMPPOINT' )
        DFLAG = ( DNAME( 1:5 ) .NE. 'NONE ' )

        HNAME = PROMPTMFILE( 
     &          'Enter name for HOUR-SPECIFIC file, or "NONE"',
     &          FSREAD3, 'NONE', 'TMPPOINT' )
        HFLAG = ( HNAME( 1:5 ) .NE. 'NONE ' )

        CDEV = PROMPTFFILE( 
     &           'Enter logical name for ACTUAL SCC file >> ',
     &           .TRUE., .TRUE., 'PSCC', 'TMPPOINT' )

        XDEV = PROMPTFFILE( 
     &           'Enter name for TEMPORAL XREF file',
     &           .TRUE., .TRUE., 'PTREF', 'TMPPOINT' )
        IF ( XDEV .LE. 0 ) THEN
            CALL M3EXIT( 'TMPPOINT', 0, 0,
     &                   'Could not open TEMPORAL XREF file', 2 )
        END IF

        PDEV = PROMPTFFILE( 
     &           'Enter name for TEMPORAL PROFILES file',
     &           .TRUE., .TRUE., 'PTPRO', 'TMPPOINT' )
        IF ( PDEV .LE. 0 ) THEN
            CALL M3EXIT( 'TMPPOINT', 0, 0,
     &                   'Could not open TEMPORAL PROFILES file', 2 )
        END IF

        UDEV = PROMPTFFILE( 
     &           'Enter name for ELEVATED SOURCES file, or "NONE"',
     &           .TRUE., .TRUE., 'PELV', 'TMPPOINT' )
        IF      ( UDEV .EQ. -2 ) THEN
            UFLAG = .FALSE.
        ELSE IF ( UDEV .GT. 0 ) THEN
            UFLAG = .TRUE.
        ELSE
            CALL M3EXIT( 'TMPPOINT', 0, 0, 
     &                   'Could not open ELEVATED SOURCES file', 2 )
        END IF


C.......   Get time period for output file:

        IF ( DFLAG ) THEN
            IF ( .NOT. DESC3( DNAME ) ) THEN
                CALL M3EXIT( 'TMPPOINT', 0, 0, 
     &                      'Could not get description of file "' 
     &                      // DNAME( 1:TRIMLEN( DNAME ) ) // '"', 2 )
            END IF
            JDATE  = SDATE3D
            NSTEPS = 24 * MXREC3D
        ELSE
            JDATE  = 1987200
            NSTEPS = 24
        END IF          !  if using day-specific info, or not

        IF ( ISDSTIME( JDATE ) ) THEN
            TZONE = 4
        ELSE 
            TZONE = 5
        END IF
        TZONE  = GETNUM( -12, 12, TZONE, 
     &           'Enter time zone (0 for GMT, 5 for EST, 4 for EDT)' )
        JDATE  = GETNUM( 0, 9999999, JDATE, 
     &                  'Enter starting date (YYYYDDD)' )
        JTIME  = GETNUM( 0, 999999, 0, 
     &                  'Enter starting time  (HHMMSS)' )
        TSTEP  = 10000
        NSTEPS = GETNUM( 1, 999999, NSTEPS,
     &                   'Enter output duration (hours)' )


C.......   Get list of output variables (shrink list in file description
C.......   data structure to just those variables desired):

        IF ( .NOT. DESC3( ENAME ) ) THEN
            CALL M3EXIT( 'TMPPOINT', 0, 0, 
     &                  'Could not get description of file "' 
     &                  // ENAME( 1:TRIMLEN( ENAME ) ) // '"', 2 )
        END IF


C.......   Determine whether input represents peak-day emissions or
C.......   annual total emissions:

        C0 =   INDEX( FDESC3D( 6 ), '"' )
        IF ( C0 .EQ. 0 ) THEN
            MESG = 'Bad field [' // FDESC3D( 6 )( 1:24 ) //
     &             '] in file description -- no leading quote'
            CALL M3EXIT( 'TMPPOINT', 0, 0, MESG, 2 )
        END IF
        C1 =  INDEX( FDESC3D( 6 )( C0+1 : 24 ) , '"' )
        IF ( C1 .EQ. 0 ) THEN
            MESG = 'Bad field [' // FDESC3D( 6 )( 1:24 ) //
     &             '] in file description -- no trailing quote'
            CALL M3EXIT( 'TMPPOINT', 0, 0, MESG, 2 )
        END IF

 
C...........   Compute offset of pollutant codes within input VNAME3D:

        DO  11  V = 1, NIPOL

            I = INDEX1( EINAM( V ), NVARS3D, VNAME3D )
            IF( I .LE. 0 ) THEN
                CALL M3EXIT( 'TMPPOINT', 0, 0, 
     &                       'Could not find variable "' // 
     &                       EINAM( V )( 1:TRIMLEN( EINAM( I ) ) )
     &                       // '" in "' // 
     &                       ENAME( 1:TRIMLEN( ENAME ) ) //  '"',
     &                       2 )
            END IF

            VNAME3D( V ) = EINAM( V )
            UNITS3D( V ) = 'tons/hour'
            VDESC3D( V ) = VDESC3D( I )
            VTYPE3D( V ) = M3REAL

11      CONTINUE

        NVARS3D = NIPOL
        SDATE3D = JDATE
        STIME3D = JTIME
        TSTEP3D = 10000
        WRITE( FDESC3D( 9 ),94010 ) 'Time zone ', TZONE


C...........   Read in emissions indexing data:  IFIP, ASC7, ASC3, 
C...........   IPLANT, ISTACK, TZONES

        IF ( .NOT. READ3( ENAME, 'IFIP', ALLAYS3,0,0, IFIP ) ) THEN
            CALL M3EXIT( 'TMPPOINT', 0, 0, 
     &                   'Error reading FIP from file' // ENAME, 2)
        END IF

        IF ( .NOT. READ3( ENAME, 'ISCC', ALLAYS3,0,0, ISCC ) ) THEN
            CALL M3EXIT( 'TMPPOINT', 0, 0, 
     &                   'Error reading ISCC from file' // ENAME, 2 )
        END IF

        IF ( .NOT. READ3( ENAME, 'IPLANT', ALLAYS3,0,0, IPLANT ) ) THEN
            CALL M3EXIT( 'TMPPOINT', 0, 0, 
     &                   'Error reading IPLANT from file' // ENAME, 2 )
        END IF

        IF ( .NOT. READ3( ENAME, 'ISTACK', ALLAYS3,0,0, ISTACK ) ) THEN
            CALL M3EXIT( 'TMPPOINT', 0, 0, 
     &                   'Error reading ISTACK from file' // ENAME, 2 )
        END IF

        IF ( .NOT. READ3( ENAME, 'TZONES', ALLAYS3,0,0, TZONES ) ) THEN
            CALL M3EXIT( 'TMPPOINT', 0, 0, 
     &                   'Error reading TZONES from file' // ENAME, 2 )
        END IF

        IF ( .NOT. READ3( ENAME, 'TPFLAG', ALLAYS3,0,0, TPFLAG ) ) THEN
            CALL M3EXIT( 'TMPPOINT', 0, 0, 
     &                   'Error reading TPFLAG from file' // ENAME, 2 )
        END IF


C.......   Get file names; open output files
C.......   Deal with optional elevated point source files

        YDEV = PROMPTFFILE( 
     &  'Enter name for XREF OUTPUT file, or "NONE"',
     &           .FALSE., .TRUE., 'NONE', 'TMPPOINT' )
        IF ( YDEV .EQ. -2 ) THEN
            YFLAG = .FALSE.
        ELSE IF ( YDEV .GT. 0 ) THEN
            YFLAG = .TRUE.
        ELSE
            CALL M3EXIT( 'TMPPOINT', 0, 0, 
     &                   'Could not open XREF OUTPUT file', 2 )
        END IF


        IF ( UFLAG ) THEN       ! generating upper-level point sources file
        
            MESG  = 'Enter name for output LOW-LEVEL EMIS file'
            TNAME = PROMPTMFILE( MESG, FSUNKN3, 'PTMP', 'TMPPOINT' )

            FTYPE3D = IDDATA3
            NROWS3D = NPELV

            FDESC3D( 1 ) = 'Upper-level Point Source emissions values.'
            DO  21  I = 2, MXDESC3
                FDESC3D( I ) = ' '
21          CONTINUE

            MESG  = 'Enter name for output ELEVATED EMIS file'
            UNAME = PROMPTMFILE( MESG, FSUNKN3, 'ETMP', 'TMPPOINT' )

            WRITE( *,92000 ) 
     &             ' ', 'Reading elevated sources file...', ' '

            EFLAG = .FALSE.     !  for finding invalid sources
            IREC  = 0
            S     = 0

22          CONTINUE

                IREC = IREC + 1
                READ( UDEV,*,END=23,IOSTAT=IOS ) FIP, SCC, PLT, STK

                IF ( IOS .NE. 0 ) THEN

                    EFLAG = .TRUE.
                    WRITE( MESG,94010 )
     &              'I/O error', IOS,
     &              'reading elevated sources file at line', IREC
                    CALL M3MESG( MESG )

                END IF

                I = FIND4( FIP, SCC, PLT, STK ,
     &                     NPSRC, IFIP, ISCC, IPLANT, ISTACK )

                IF ( I .LT. 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 )
     &                  'Invalid source ' //
     &                  'FIP',   FIP,
     &                  'SCC',   SCC,
     &                  'Plant', PLT,
     &                  'Stack', STK, 'in elevated sources file.'
                    CALL M3MESG( MESG )

                ELSEIF( S .GE. NPELV ) THEN 

                    EFLAG = .TRUE.
                    WRITE( MESG, 94010 )
     &                  'Max elevated sources size', NPELV,
     &                  'exceeded at line', IREC
                    CALL M3MESG( MESG )

                ELSE
                    S = S + 1
                    INDXE( S ) = I

                END IF          !  if s < npup

            GO TO  22       !  to head of elevated sources loop

23          CONTINUE    !  exit from loop reading elevated sources file

            IF ( EFLAG ) THEN
                MESG = 'Error reading elevated sources file'
                CALL M3EXIT( 'TMPPOINT', 0, 0, MESG, 2 )

            ENDIF

            NELPT = S

        ELSE

            MESG = 'Enter name for output EMIS file'
            TNAME = PROMPTMFILE( MESG, FSUNKN3, 'PTMP', 'TMPPOINT' )

            NELPT = 0

        END IF          !  if generating upper-level point sources file


C.......   Read the ACTUAL-SCC FILE

        WRITE( *,92000 ) ' ', 'Reading ACTUAL SCC file...', ' '

        IREC  =  0
        EFLAG = .FALSE.
        I     =  1
        LID   =  0
        ASCC( 1 ) = 0

33      CONTINUE        !  head of the CDEV-read loop

            READ( CDEV, *, END=34, IOSTAT=IOS ) SCC

            IREC = IREC + 1

            IF ( IOS .NE. 0 ) THEN

                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 
     &              'I/O error', IOS, 
     &              'reading ACTUAL SCC file at line', IREC
                CALL M3MESG( MESG )

            ELSE IF ( I .GE. NPSRC ) THEN 

                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 
     &              'Max SCC table size', NPSRC,
     &              'exceeded at line', IREC
                CALL M3MESG( MESG )

            ELSE IF ( LID .GE. SCC ) THEN 

                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 
     &              'SCC table out of order at line', IREC
                CALL M3MESG( MESG )

            ELSE

                I   = I + 1
                LID = SCC
                ASCC( I ) = SCC

            END IF

            GO TO  33   !  to head of loop

34      CONTINUE        !  end of the SDEV-read loop

        IF ( EFLAG ) THEN
            MESG = 'Error reading ACTUAL SCC file'
            CALL M3EXIT( 'TMPPOINT', 0, 0, MESG, 2 )

        END IF

        NSCC = I


C.......   Read temporal-profile cross-reference file.

        WRITE( *,92000 ) 
     &         ' ', 'Reading TEMPORAL XREF file...', ' '

        IREC  = 0
        K     = 0
        EFLAG = .FALSE.

44      CONTINUE

            IREC = IREC + 1
            READ( XDEV, 93000, END=45, IOSTAT=IOS ) LINE

            IF ( IOS .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 
     &              'I/O error', IOS, 
     &              'reading TEMPORAL XREF file at line', IREC
                CALL M3MESG( MESG )
                GO TO 44
            END IF

C.............  Ensure valid pollutant code
            COD = STR2INT( LINE( 32:36 ) )
            IF( COD .LE. 0 ) THEN
                I = 0
            ELSE
                I = FIND1( COD, NIPOL, PCODE )
                IF ( I .LT. 0 ) THEN
                    GO TO  44
                END IF
            END IF

 
            LS = LBLANK( LINE( 21:30 ) )     ! Number of leading spaces
            EC = TRIMLEN( LINE( 21:30 ) )  ! Size without ending blanks
 
C............. Convert character SCC field to integer SCC number
            IF( LS .EQ. 0 .AND. EC .EQ. 10 ) THEN  ! Area x-ref entry
                GO TO 44

            ELSEIF( LINE( 21:21) .GT. '9' ) THEN   ! Motor vehicle ASC code
                GO TO 44

            ELSE
                SCC = STR2INT( LINE( 21:30 ) )

            ENDIF

C.............  Allow for default SCC entries
            IF( SCC .GT. 0 ) THEN
                J = FIND1( SCC, NSCC, ASCC )  ! Specific SCC

            ELSE
                J = 1                         ! Default SCC

            ENDIF

            IMON = STR2INT( LINE( 38:40 ) )
            IWEK = STR2INT( LINE( 42:44 ) )
            IDIU = STR2INT( LINE( 45:48 ) )

C.............  Ensure SCC needed
            IF ( J .LE. 0 ) THEN
                GO TO 44

C.............  Ensure temporal profile numbers are valid
            ELSE IF ( IMON .LE. 0 ) THEN  ! Bad monthly profile code
 
                WRITE( MESG,94010 )
     &                 'Bad IMON=', IMON,
     &                 'at line', IREC,
     &                 'for FIP', FIP,
     &                 'and SCC', SCC
                CALL M3MESG( MESG )
                EFLAG = .TRUE.
                GO TO  44
 
            ELSE IF ( IWEK .LE. 0 ) THEN  ! Bad weekly profile code
 
                WRITE( MESG,94010 )
     &                 'Bad IWEK=', IWEK,
     &                 'at line', IREC,
     &                 'for FIP', FIP,
     &                 'and SCC', SCC
                CALL M3MESG( MESG )
                EFLAG = .TRUE.
                GO TO  44
 
            ELSE IF ( IDIU .LE. 0 ) THEN  ! Bad diurnal profile code
 
                WRITE( MESG,94010 )
     &                 'Bad IDIU=', IDIU,
     &                 'at line', IREC,
     &                 'for FIP', FIP,
     &                 'and SCC', SCC
                    CALL M3MESG( MESG )
                    EFLAG = .TRUE.
                    GO TO  44

            ELSEIF ( K .GE. MXREC ) THEN

                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 
     &                 'Max XREF table size', MXREC,
     &                 'exceeded at line', IREC
                CALL M3MESG( MESG )
                GO TO 44

            ELSE    !  else incremented K within bounds:

                K = K + 1
                INDXA( K ) = K
                ICODA( K ) = I
                IFIPA( K ) = STR2INT( LINE(  1: 5 ) )
                IPLTA( K ) = STR2INT( LINE(  7:11 ) )
                ISTKA( K ) = STR2INT( LINE( 13:15 ) )
                ISCCA( K ) = SCC
                MPRNA( K ) = IMON
                WPRNA( K ) = IWEK
                DPRNA( K ) = IDIU

            END IF      !  No oveflow, SCC foudn in list, no errors

        GO TO  44       !  to head of loop:  read lines from XREF file

45      CONTINUE        !  exit from loop:  read lines from XREF file

        IF ( EFLAG ) THEN
            MESG = 'Error reading TEMPORAL XREF file'
            CALL M3EXIT( 'TMPPOINT', 0, 0, MESG, 2 )
        END IF
        
        NPRNA = K

C.......   Sort and Process temporal cross-references, according to category

        WRITE( *,92000 ) 
     &      ' ', 'Processing TEMPORAL XREF file...', ' '
        
        CALL SORTI4( NPRNA, INDXA, IFIPA, ISCCA, IPLTA, ISTKA )

        EFLAG = .FALSE.         !  for detecting XREF errors
        NPRN1 = 0
        NPRN2 = 0
        NPRN3 = 0
        NPRN4 = 0
        NPRN5 = 0

        DO  66  J = 1, NPRNA
            
            K    = INDXA( J )
            FIP  = IFIPA( K )
            PLT  = IPLTA( K )
            STK  = ISTKA( K )
            SCC  = ISCCA( K )
            I    = ICODA( K )
            IMON = MPRNA( K )
            IWEK = WPRNA( K )
            IDIU = DPRNA( K )

            IF ( FIP .LE. 0 ) THEN  !  FIP-independent fallback profiles)

                IF ( SCC .LE. 0 ) THEN  !  ultimate fallback profiles )

                    IF ( I .EQ. 0 ) THEN
                        DO  46  I = 1, NIPOL
                            MPR00( I ) = IMON
                            WPR00( I ) = IWEK
                            DPR00( I ) = IDIU
46                      CONTINUE
                    ELSE
                        MPR00( I ) = IMON
                        WPR00( I ) = IWEK
                        DPR00( I ) = IDIU
                        DO  47  V = 1, NIPOL
                            IF ( MPR00( V ) .EQ. IMISS3 ) 
     &                           MPR00( V ) = IMON
                            IF ( WPR00( V ) .EQ. IMISS3 ) 
     &                           WPR00( V ) = IWEK
                            IF ( DPR00( V ) .EQ. IMISS3 ) 
     &                           DPR00( V ) = IDIU
47                      CONTINUE
                    END IF

                ELSE                    !  ASC-dependent-only profiles

                    IF ( NPRN1 .EQ. 0  .OR.
     &                   SCC .NE. ISCC1( NPRN1 ) ) THEN

                        NPRN1 = NPRN1 + 1
                        IF ( NPRN1 .LE. NPSRC ) THEN
                            ISCC1( NPRN1 ) = SCC
                        ENDIF

                    END IF      !  if new scc encountered

                    IF ( I .EQ. 0 ) THEN
                        DO  48  I = 1, NIPOL
                            MPRN1( NPRN1,I ) = IMON
                            WPRN1( NPRN1,I ) = IWEK
                            DPRN1( NPRN1,I ) = IDIU
48                      CONTINUE
                    ELSE
                        MPRN1( NPRN1,I ) = IMON
                        WPRN1( NPRN1,I ) = IWEK
                        DPRN1( NPRN1,I ) = IDIU
                        DO  49  V = 1, NIPOL
                            IF ( MPRN1( NPRN1,V ) .EQ. IMISS3 ) 
     &                           MPRN1( NPRN1,V ) = IMON
                            IF ( WPRN1( NPRN1,V ) .EQ. IMISS3 ) 
     &                           WPRN1( NPRN1,V ) = IWEK
                            IF ( DPRN1( NPRN1,V ) .EQ. IMISS3 ) 
     &                           DPRN1( NPRN1,V ) = IDIU
49                      CONTINUE
                    END IF

                END IF

            ELSE IF ( MOD( FIP,1000 ) .LE. 0 ) THEN     !  state-ASC dependent

                IF ( NPRN2 .EQ. 0                 .OR.
     &               FIP/1000 .NE. IFIP2( NPRN2 ) .OR.
     &               SCC      .NE. ISCC2( NPRN2 )      ) THEN

                    NPRN2 = NPRN2 + 1
                    IF ( NPRN2 .LE. NPSRC ) THEN
                        IFIP2( NPRN2 ) = FIP / 1000 !  state-code only
                        ISCC2( NPRN2 ) = SCC
                    END IF

                END IF

                IF ( I .EQ. 0 ) THEN
                    DO  50  I = 1, NIPOL
                        MPRN2( NPRN2,I ) = IMON
                        WPRN2( NPRN2,I ) = IWEK
                        DPRN2( NPRN2,I ) = IDIU
50                  CONTINUE
                ELSE
                    MPRN2( NPRN2,I ) = IMON
                    WPRN2( NPRN2,I ) = IWEK
                    DPRN2( NPRN2,I ) = IDIU
                    DO  51  V = 1, NIPOL
                        IF ( MPRN2( NPRN2,V ) .EQ. IMISS3 ) 
     &                          MPRN2( NPRN2,V ) = IMON
                        IF ( WPRN2( NPRN2,V ) .EQ. IMISS3 ) 
     &                          WPRN2( NPRN2,V ) = IWEK
                        IF ( DPRN2( NPRN2,V ) .EQ. IMISS3 ) 
     &                          DPRN2( NPRN2,V ) = IDIU
51                  CONTINUE
                END IF

            ELSE IF ( PLT .LE. 0 ) THEN !  by FIP,SCC only

                IF ( NPRN3 .EQ. 0              .OR.
     &               FIP   .NE. IFIP3( NPRN3 ) .OR.
     &               SCC   .NE. ISCC3( NPRN3 )      ) THEN

                    NPRN3 = NPRN3 + 1
                    IF ( NPRN3 .LE. NPSRC ) THEN
                        IFIP3( NPRN3 ) = FIP
                        ISCC3( NPRN3 ) = SCC
                    END IF

                END IF

                IF ( I .EQ. 0 ) THEN
                    DO  52  I = 1, NIPOL
                        MPRN3( NPRN3,I ) = IMON
                        WPRN3( NPRN3,I ) = IWEK
                        DPRN3( NPRN3,I ) = IDIU
52                  CONTINUE
                ELSE
                    MPRN3( NPRN3,I ) = IMON
                    WPRN3( NPRN3,I ) = IWEK
                    DPRN3( NPRN3,I ) = IDIU
                    DO  53  V = 1, NIPOL
                        IF ( MPRN3( NPRN3,V ) .EQ. IMISS3 ) 
     &                       MPRN3( NPRN3,V ) = IMON
                        IF ( WPRN3( NPRN3,V ) .EQ. IMISS3 ) 
     &                       WPRN3( NPRN3,V ) = IWEK
                        IF ( DPRN3( NPRN3,V ) .EQ. IMISS3 ) 
     &                       DPRN3( NPRN3,V ) = IDIU
53                  CONTINUE
                END IF

            ELSE IF ( STK .LE. 0 ) THEN !  profiles by FIP,SCC,PLANT

                IF ( NPRN4 .EQ. 0             .OR.
     &               FIP .NE. IFIP4( NPRN4 )  .OR.
     &               PLT .NE. IPLT4( NPRN4 )  .OR.
     &               SCC .NE. ISCC4( NPRN4 ) ) THEN

                    NPRN4 = NPRN4 + 1
                    IF ( NPRN4 .LE. NPSRC ) THEN
                        IFIP4( NPRN4 ) = FIP
                        IPLT4( NPRN4 ) = PLT
                        ISCC4( NPRN4 ) = SCC
                    END IF

                END IF

                IF ( I .EQ. 0 ) THEN
                    DO  54  I = 1, NIPOL
                        MPRN4( NPRN4,I ) = IMON
                        WPRN4( NPRN4,I ) = IWEK
                        DPRN4( NPRN4,I ) = IDIU
54                  CONTINUE
                ELSE
                    MPRN4( NPRN4,I ) = IMON
                    WPRN4( NPRN4,I ) = IWEK
                    DPRN4( NPRN4,I ) = IDIU
                    DO  55  V = 1, NIPOL
                        IF ( MPRN4( NPRN4,V ) .EQ. IMISS3 ) 
     &                       MPRN4( NPRN4,V ) = IMON
                        IF ( WPRN4( NPRN4,V ) .EQ. IMISS3 ) 
     &                       WPRN4( NPRN4,V ) = IWEK
                        IF ( DPRN4( NPRN4,V ) .EQ. IMISS3 ) 
     &                       DPRN4( NPRN4,V ) = IDIU
55                  CONTINUE
                END IF

            ELSE        !  select by FIP,SCC,PLANT,STACK:

                S = FIND4( FIP, SCC, PLT, STK,
     &                     NPSRC, IFIP, ISCC, IPLANT, ISTACK )

                IF ( S .LT. 0 ) THEN    !  irrelevant XREF entry.

                    GO TO  66

                ELSEIF ( NPRN5 .EQ. 0               .OR.
     &                   FIP   .NE. IFIP5( NPRN5 )  .OR.
     &                   PLT   .NE. IPLT5( NPRN5 )  .OR.
     &                   STK   .NE. ISTK5( NPRN5 )  .OR.
     &                   SCC   .NE. ISCC5( NPRN5 )       ) THEN

                    NPRN5 = NPRN5 + 1
                    IF ( NPRN5 .LE. NPSRC ) THEN
                        IFIP5( NPRN5 ) = FIP
                        IPLT5( NPRN5 ) = PLT
                        ISTK5( NPRN5 ) = STK
                        ISCC5( NPRN5 ) = SCC
                    END IF

                END IF

                IF ( I .EQ. 0 ) THEN
                    DO  56  I = 1, NIPOL
                        MPRN5( NPRN5,I ) = IMON
                        WPRN5( NPRN5,I ) = IWEK
                        DPRN5( NPRN5,I ) = IDIU
56                  CONTINUE
                ELSE
                    MPRN5( NPRN5,I ) = IMON
                    WPRN5( NPRN5,I ) = IWEK
                    DPRN5( NPRN5,I ) = IDIU
                    DO  57  V = 1, NIPOL
                        IF ( MPRN5( NPRN5,V ) .EQ. IMISS3 ) 
     &                       MPRN5( NPRN5,V ) = IMON
                        IF ( WPRN5( NPRN5,V ) .EQ. IMISS3 ) 
     &                       WPRN5( NPRN5,V ) = IWEK
                        IF ( DPRN5( NPRN5,V ) .EQ. IMISS3 ) 
     &                       DPRN5( NPRN5,V ) = IDIU
57                  CONTINUE
                END IF

            END IF      !  if not found, or fip=0, or ...

66      CONTINUE        !  end loop processing temporal XREFs

        IF( NPRN1 .GT. NPSRC ) THEN
            WRITE( MESG,94010 )
     &          'Max XREF count', NPSRC,
     &          'exceeded in XREF file: count=', NPRN1
            CALL M3EXIT( 'TMPPOINT', 0, 0, MESG, 2 )

        ELSE IF( NPRN2 .GT. NPSRC ) THEN
            WRITE( MESG,94010 )
     &          'Max XREF count', NPSRC,
     &          'exceeded in XREF file: count=', NPRN2
            CALL M3EXIT( 'TMPPOINT', 0, 0, MESG, 2 )

        ELSE IF( NPRN3 .GT. NPSRC ) THEN
            WRITE( MESG,94010 )
     &          'Max XREF count', NPSRC,
     &          'exceeded in XREF file: count=', NPRN3
            CALL M3EXIT( 'TMPPOINT', 0, 0, MESG, 2 )

        ELSE IF( NPRN4 .GT. NPSRC ) THEN
            WRITE( MESG,94010 )
     &          'Max XREF count', NPSRC,
     &          'exceeded in XREF file: count=', NPRN4
            CALL M3EXIT( 'TMPPOINT', 0, 0, MESG, 2 )

        ELSE IF( NPRN5 .GT. NPSRC ) THEN
            WRITE( MESG,94010 )
     &          'Max XREF count', NPSRC,
     &          'exceeded in XREF file: count=', NPRN5
            CALL M3EXIT( 'TMPPOINT', 0, 0, MESG, 2 )

        END IF

C...........   If requested, write out filled-in XREF file:

        IF ( YFLAG ) THEN

            WRITE( YDEV,93020, IOSTAT=IOS )
     &      ( 0, 0, -9, 0, 
     &        PCODE( V ), MPR00( V ), WPR00( V ), DPR00( V ), 
     &              V = 1, NIPOL ),
     &      ( ( 0, -999, -9, ISCC1( S ), 
     &          PCODE( V ), MPRN1( S,V ), WPRN1( S,V ), DPRN1( S,V ), 
     &              V = 1, NIPOL ), S = 1, NPRN1 ),
     &      ( ( 1000*IFIP2( S ), -999, -9, ISCC2( S ), 
     &          PCODE( V ), MPRN2( S,V ), WPRN2( S,V ), DPRN2( S,V ), 
     &              V = 1, NIPOL ), S = 1, NPRN2 ),
     &      ( ( IFIP3( S ), -999, -9, ISCC3( S ), 
     &          PCODE( V ), MPRN3( S,V ), WPRN3( S,V ), DPRN3( S,V ), 
     &              V = 1, NIPOL ), S = 1, NPRN3 ),
     &      ( ( IFIP4( S ), IPLT4( S ), -9, ISCC4( S ), 
     &          PCODE( V ), MPRN4( S,V ), WPRN4( S,V ), DPRN4( S,V ), 
     &              V = 1, NIPOL ), S = 1, NPRN4 ),
     &      ( ( IFIP5( S ), IPLT5( S ), ISTK5( S ), ISCC5( S ), 
     &          PCODE( V ), MPRN5( S,V ), WPRN5( S,V ), DPRN5( S,V ), 
     &              V = 1, NIPOL ), S = 1, NPRN5 )

            IF ( IOS .NE. 0 ) THEN
                WRITE( MESG,94010 ) 
     &              'I/O error', IOS, 
     &              'writing OUTPUT XREF file'
                CALL M3EXIT( 'TMPPOINT', 0, 0, MESG, 1 )
            END IF

            CLOSE( YDEV )

        END IF  !  if yflag (optional TEMPORAL XREF file being produced)


C...........   Read temporal-profiles file:  3 parts (monthly, weekly, diurnal)

        WRITE( *,92000 ) ' ', 'Reading TEMPORAL PROFILES file...', ' '

        NMON = RDTPROF( PDEV, 'MONTHLY', 12, NMCOD, MONREF, MONFAC )
        NWEK = RDTPROF( PDEV, 'WEEKLY',   7, NWCOD, WEKREF, WEKFAC )
        NDIU = RDTPROF( PDEV, 'WEEKDAY', 24, NDCOD, DIUREF, DIUFAC )
        NEND = RDTPROF( PDEV, 'WEEKEND', 24, NDCOD, ENDREF, ENDFAC )


C...........   Construct  -- mapping source number to effective
C...........   composite profile number.  Six-way procedure:
C...........       FIP-SCC-PLANT-STACK match; else 
C...........       FIP-SCC-PLANT match;       else 
C...........       FIP-SCC match;             else 
C...........       FIP match;                 else 
C...........       FIP-state-code match;      else 
C...........       default profile.

        EFLAG = .FALSE.         !  for detecting xref-mapping errors

        DO  81  S = 1, NPSRC

            MFLAG = ( MOD( TPFLAG( S ), MTPRFAC ) .EQ. 0 )
            WFLAG = ( MOD( TPFLAG( S ), WTPRFAC ) .EQ. 0 )

            IF ( .NOT. MFLAG ) THEN
                DO  71  V = 1, NIPOL
                    MDEX( S,V ) = 0
71              CONTINUE
            END IF

            IF ( .NOT. WFLAG ) THEN
                DO  72  V = 1, NIPOL
                    WDEX( S,V ) = 0
72              CONTINUE
            END IF

            FIP = IFIP  ( S )
            SCC = ISCC  ( S )
            PLT = IPLANT( S )
            STK = ISTACK( S )

C...........   Try to get FIP-SCC-PLANT-STACK match;

            F = FIND4( FIP, SCC, PLT, STK , 
     &                 NPRN5, IFIP5, ISCC5, IPLT5, ISTK5 )

            IF ( F .GT. 0 ) THEN

                DO  73  V = 1, NIPOL
                    MREF( V ) = MPRN5( F,V )
                    WREF( V ) = WPRN5( F,V )
                    DREF( V ) = DPRN5( F,V )
73              CONTINUE

                GO TO  79       !  to end of sources-loop

            END IF      !  if find4() succeeded

C...........   Try to get FIP-SCC-PLANT match;

            F = FIND3( FIP, SCC, PLT,
     &                 NPRN4, IFIP4, ISCC4, IPLT4 )

            IF ( F .GT. 0 ) THEN

                DO  74  V = 1, NIPOL
                    MREF( V ) = MPRN4( F,V )
                    WREF( V ) = WPRN4( F,V )
                    DREF( V ) = DPRN4( F,V )
74              CONTINUE

                GO TO  79       !  to end of sources-loop

            END IF      !  if find3() succeeded

C...........   Try to get FIP-SCC match;

            F = FIND2( FIP, SCC, NPRN3, IFIP3, ISCC3 )

            IF ( F .GT. 0 ) THEN

                DO  75  V = 1, NIPOL
                    MREF( V ) = MPRN3( F,V )
                    WREF( V ) = WPRN3( F,V )
                    DREF( V ) = DPRN3( F,V )
75              CONTINUE

                GO TO  79       !  to end of sources-loop

            END IF      !  if find2() succeeded: FIP,SCC case

C...........   Try to get state-SCC match;

            F = FIND2( FIP/1000, SCC, NPRN2, IFIP2, ISCC2 )

            IF ( F .GT. 0 ) THEN

                DO  76  V = 1, NIPOL
                    MREF( V ) = MPRN2( F,V )
                    WREF( V ) = WPRN2( F,V )
                    DREF( V ) = DPRN2( F,V )
76              CONTINUE

                GO TO  79       !  to end of sources-loop

            END IF      !  if find2() succeeded:  state,SCC case

C...........   Try to get SCC-only match;

            F = FIND1( SCC, NPRN1, ISCC1 )

            IF ( F .GT. 0 ) THEN

                DO  77  V = 1, NIPOL
                    MREF( V ) = MPRN1( F,V )
                    WREF( V ) = WPRN1( F,V )
                    DREF( V ) = DPRN1( F,V )
77              CONTINUE

            ELSE        ! else use default no-FIP, no-SCC profile

                DO  78  V = 1, NIPOL
                    MREF( V ) = MPR00( V )
                    WREF( V ) = WPR00( V )
                    DREF( V ) = DPR00( V )
78              CONTINUE

            END IF      !  if find3() FIP-and-ASC succeeded

C...............   Now look up these refs

79          CONTINUE

            DO  80  V = 1, NIPOL
                IF ( MFLAG ) THEN
                    MDEX( S,V ) = FIND1( MREF( V ), NMON, MONREF )
                    IF ( MDEX( S,V ) .LT. 0 ) THEN
                        WRITE( MESG,94010 ) 
     &                      'Profile MREF', MREF( V ),
     &                      'not found'
                        CALL M3MESG( MESG )
                        EFLAG = .TRUE.
                    END IF
                END IF

                IF ( WFLAG ) THEN
                    WDEX( S,V ) = FIND1( WREF( V ), NWEK, WEKREF )
                    IF ( WDEX( S,V ) .LT. 0 ) THEN
                        WRITE( MESG,94010 ) 
     &                      'Profile WREF', WREF( V ),
     &                      'not found'
                        CALL M3MESG( MESG )
                        EFLAG = .TRUE.
                    END IF
                END IF

                DDEX( S,V ) = FIND1( DREF( V ), NDIU, DIUREF )
                IF ( DDEX( S,V ) .LT. 0 ) THEN
                    WRITE( MESG,94010 ) 
     &                  'Profile DREF', DREF( V ), 'not found'
                         CALL M3MESG( MESG )
                     EFLAG = .TRUE.
                END IF

                IF ( NEND .EQ. 0 ) GO TO  80
                EDEX( S,V ) = FIND1( DREF( V ), NEND, ENDREF )
                IF ( EDEX( S,V ) .LT. 0 ) THEN
                    WRITE( MESG,94010 ) 
     &                  'Profile EREF', DREF( V ), 'not found'
                    CALL M3MESG( MESG )
                END IF
80          CONTINUE

81      CONTINUE        !  end loop on sources S, constructing TNDX

        IF ( EFLAG ) THEN
            MESG = 'Bad XREF: PROFILE matches encountered'
            CALL M3EXIT( 'TMPPOINT', 0, 0, MESG, 2 )
        END IF          ! if eflag:  error mapping sources to profiles


C...........   Construct TMAT -- array of effective composite 
C...........   profile coefficients

C...........   Read in emissions values.

        DO  88  V = 1, NIPOL

            IF ( .NOT. READ3( ENAME, EINAM( V ), ALLAYS3, 0, 0, 
     &                        EMIS( 1,V ) ) ) THEN
                CALL M3EXIT( 'TMPPOINT', 0, 0, 
     &                       'Error reading "' // 
     &                       EINAM( V )( 1:TRIMLEN( EINAM( V ) ) )//
     &                       '" from file "' // 
     &                       ENAME( 1:TRIMLEN( ENAME ) ) // '".', 2 )
            END IF

88      CONTINUE


C...........   If not using day-specific emissions, set NDAYPT to zero

        IF ( .NOT. DFLAG ) THEN

            NDAYPT = 0

            DO  99  V = 1, NIPOL
            DO  98  S = 1, NPSRC

                EMISV( S,V ) = EMIS( S,V )

98          CONTINUE
99          CONTINUE

        END IF          !  if not using day-specific emissions

C...........   For each time step, write out the point source emissions values:

        WRITE( *,92000 ) ' ', 'Writing out emissions values...', ' '

        LDATE = -1   !  initialization forces TMAT calculation on first iter.

        DO  199  T = 1, NSTEPS

C...............   If this is a new month, or transitions between
C...............   weekday and weekend, calculate a new TMAT


            IF ( LDATE .NE. JDATE ) THEN

                LDATE = JDATE

C...............   If DFLAG, prepare day-corrections:
                
                IF ( DFLAG ) THEN
                
                    IF ( .NOT. READ3( DNAME, ALLVAR3, ALLAYS3,
     &                                JDATE, JTIME, NDAYPT ) ) THEN

                        CALL M3EXIT( 'TMPPOINT', JDATE, JTIME,
     &                       'Could not read "' // 
     &                       DNAME( 1 : TRIMLEN( DNAME ) ) // '".', 2 )

                    END IF      !  if read3() failed on dname

                    DO  122  V = 1, NIPOL

                        DO  111  S = 1, NPSRC
                            EMISV( S,V ) = EMIS( S,V )
111                     CONTINUE

                        DO  112  S = 1, NDAYPT
                            EMISV( INDXD( S ),V ) = EMISD( S,V )
112                     CONTINUE

122                 CONTINUE

                END IF          ! if using day-specific emissions

C...............   Get day-of-week, month:
              
                CALL DAYMON( JDATE, MONTH, DAY )  !  get month
                DAY = WKDAY( JDATE )              !  get day-of-week

                IF ( NEND .EQ. 0 ) THEN           !  no weekend packet

                    WRITE( *,92011 ) 'Processing day', DAY
                    CALL MKTMAT( NPSRC, JDATE, MONTH, DAY,
     &                           TZONE, TZONES,
     &                           MONFAC, MDEX,
     &                           WEKFAC, WDEX,
     &                           DIUFAC, DDEX, 
     &                           NDAYPT, INDXD, TMAT )

                ELSE IF ( DAY .GE. 6 ) THEN     !  day is Saturday or Sunday
                                                !  in SAI's brain-damaged
                                                !  day-of-week conventions.
                    WRITE( *,92011 ) 'Processing weekend day', DAY
                    CALL MKTMAT( NPSRC, JDATE, MONTH, DAY,
     &                           TZONE, TZONES,
     &                           MONFAC, MDEX,
     &                           WEKFAC, WDEX,
     &                           ENDFAC, EDEX, 
     &                           NDAYPT, INDXD, TMAT )

                ELSE

                    WRITE( *,92011 ) 'Processing midweek day', DAY
                    CALL MKTMAT( NPSRC, JDATE,  MONTH, DAY, 
     &                           TZONE, TZONES,
     &                           MONFAC, MDEX,
     &                           WEKFAC, WDEX,
     &                           DIUFAC, DDEX, 
     &                           NDAYPT, INDXD, TMAT )

                END IF

            END IF      ! if jdate not ldate

C...............   If HFLAG, prepare hour-corrections:

            H = 1 + MOD( JTIME / 10000 , 24 )

            IF ( HFLAG ) THEN

                IF ( .NOT. READ3( HNAME, ALLVAR3, ALLAYS3,
     &                            JDATE, JTIME, NHRPT ) ) THEN

                    CALL M3EXIT( 'TMPPOINT', JDATE, JTIME,
     &                   'Could not read "' // 
     &                   HNAME( 1 : TRIMLEN( HNAME ) ) // '".', 2 )

                END IF  !  if read3() failed on dname

                DO  177  V = 1, NIPOL
                   
                    DO  174  S = 1, NPSRC
                        EMIST( S,V ) = EMISV( S,V ) * TMAT( S,V,H )
174                 CONTINUE
                   
                    DO  175  S = 1, NHRPT
                        EMIST( INDXH( S ),V ) = EMISH( S,V )
175                 CONTINUE

                    DO  176  S = 1, NELPT       !  elevated sources
                        I = INDXE( S )
                        EMISE( S,V ) = EMIST( I,V )
                        EMIST( I,V ) = 0.0
176                 CONTINUE

177             CONTINUE

             ELSE               ! else not hflag:

                DO  188  V = 1, NIPOL

                    DO  186  S = 1, NPSRC
                        EMIST( S,V ) = EMISV( S,V ) * TMAT( S,V,H )
186                 CONTINUE

                    DO  187  S = 1, NELPT       !  elevated sources
                        I = INDXE( S )
                        EMISE( S,V ) = EMIST( I,V )
                        EMIST( I,V ) = 0.0
187                 CONTINUE

188             CONTINUE

             END IF             ! if using hour-specific emissions or not

            IF ( .NOT.WRITE3( TNAME, 'ALL', JDATE, JTIME, EMIST ) ) THEN
                CALL M3EXIT( 'TMPPOINT', JDATE, JTIME,
     &                       'Could not write timestep to "' // 
     &                       TNAME( 1 : TRIMLEN( TNAME ) ) // '".', 2 )
            END IF

            IF ( UFLAG ) THEN   !  writing elevated sources file

                IF ( .NOT.WRITE3( UNAME, 'ALL', 
     &                            JDATE, JTIME, NELPT ) ) THEN
                    CALL M3EXIT( 'TMPPOINT', JDATE, JTIME,
     &                       'Could not write timestep to "' // 
     &                       UNAME( 1 : TRIMLEN( UNAME ) ) // '".', 2 )
                END IF

            END IF      !  if writing elevated sources file

            CALL NEXTIME( JDATE, JTIME, TSTEP )

199     CONTINUE

999     CONTINUE          !  exit program

        MESG = 'Normal completion of TMPPOINT'
        CALL M3EXIT( 'TMPPOINT', 0, 0, MESG, 0 )


C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )

92010   FORMAT ( 5X , A, I10 )

92011   FORMAT (/5X , A, I10 )

92020   FORMAT ( 5X , A, I9 )


C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

93010   FORMAT( A16 )

93020   FORMAT( I5.5, 1X, I4, 1X, I3, 1X, I10, 1X, I5, 3 I4 )


C...........   Internal buffering formats.............94xxx

94010   FORMAT( 10( A, :, I10, :, 2X ) )

        END

