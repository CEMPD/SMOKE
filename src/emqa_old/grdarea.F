
        PROGRAM GRDAREA

C***********************************************************************
C  program body starts at line  146
C
C  DESCRIPTION:
C       Use gridding-transform matrices from program GRDAMAT to grid
C       source-level area source emissions files.
C
C  PRECONDITIONS REQUIRED:  
C       M3IO source-level area source emissions in canonical order used
C       by GRDAMAT
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       Models-3 I/O; GETYN, TRIMLEN
C
C  REVISION  HISTORY:
C       Prototype  2/95 by CJC.
C
C       Revised 9/96 by CJC:  match with char-string pollutant names,
C	rather than numeric names, in loop 33
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE)
C                Modeling System
C File Version @(#)$Id$
C
C COPYRIGHT (C) 1996, MCNC--North Carolina Supercomputing Center
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
C Pathname:    $Source$
C Last updated: $Date$ 
C
C***************************************************************************/

      IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'ARDIMS3.EXT'   !  area-source dimensioning parameters
        INCLUDE 'CHDIMS3.EXT'   !  emis chem parms (inventory + model)
        INCLUDE 'GRDIMS3.EXT'   !  grid parameters
        INCLUDE 'TMDIMS3.EXT'   !  temporal parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.


C...........   EXTERNAL FUNCTIONS and their descriptions:

        LOGICAL         DSCGRID
        INTEGER         ENVINT
        LOGICAL         ENVYN
        LOGICAL         GETYN
        CHARACTER*10    HHMMSS
        INTEGER         INDEX1
        CHARACTER*14    MMDDYY
        INTEGER         PROMPTFFILE
        CHARACTER*16    PROMPTMFILE
        INTEGER         TRIMLEN
        INTEGER         WKDAY
        REAL            YR2DAY

        EXTERNAL        DSCGRID, ENVINT, ENVYN, GETYN, HHMMSS, 
     &                  INDEX1, MMDDYY, PROMPTFFILE, PROMPTMFILE, 
     &                  TRIMLEN, WKDAY, YR2DAY

C...........   PARAMETER

        CHARACTER*16    NONE 
        PARAMETER     ( NONE = 'NONE' )


C...........   LOCAL VARIABLES and their descriptions:
C...........   Area Sources input and output arrays

        REAL            EMISV( NASRC )  !  one variable of input emissions.
        REAL            EMISG( NGRID )  !  gridded output variable

C...........   Gridding Matrix

        INTEGER     NS( NGRID )
        INTEGER     IS( NMATX )
        REAL        CS( NMATX )

        COMMON  / GRIDMAT / NS, IS, CS

C...........   Inventory arrays
        INTEGER         IFIP  ( NASRC )   ! FIPs  code for each source
        INTEGER         INVYR ( NASRC )   ! Inventory year code for each source
        INTEGER         ISPTR ( NASRC )   ! Pointer for each state
        INTEGER         STATES( NASID )   ! Unique list of state IDs

C...........   Report arrays
        REAL            EMSUM( NASID, NIPOL ) ! Summed emissions by state

C...........   Other local variables
 
        INTEGER         BTIME   !  beginning time of day in run
        INTEGER         IOS
        INTEGER         I, K, L, S, T, V
        INTEGER         JDATE, JTIME, TSTEP
        INTEGER         LDATE, LTIME       ! previous date and time in loop
        INTEGER         LDEV
        INTEGER         LSTA               ! previous state code
        INTEGER         NSTA               ! actual number of states
        INTEGER         NSTEPS
        INTEGER         RDEV               ! output report unit number
        INTEGER         RTIME              ! time for reporting state totals
        INTEGER         STA                ! temporay state code
        INTEGER         VTYPE( NIPOL )

        LOGICAL         EFLAG   !  Error flag
        LOGICAL         OUTPUTFLAG
        LOGICAL         ST_REPORT  ! Output state-based emissions totals
 
        CHARACTER*4     OUTNAM  !  contains default  logical output name
        CHARACTER*7     REPNAM  !  contains default logical report output name
        CHARACTER*16    ANAME   !  logical name for area-source      input file
        CHARACTER*16    ENAME   !  logical name for emission source output file
        CHARACTER*16    GNAME   !  logical name for grid matrix      input file
        CHARACTER*16    INAME   !  logical name for inventory input file

        CHARACTER*16    SCRBUF  !  scratch name buffer
        CHARACTER*16    VNAME( NIPOL )
        CHARACTER*16    UNITS( NIPOL )
        CHARACTER*16    VDESC( NIPOL )
        CHARACTER*256   MESG


C***********************************************************************
C   begin body of program GRDAREA

        LDEV = INIT3()

        WRITE( *,92000 ) 
     &  ' ',
     &  'Program GRDAREA to take a sorted source level area source ',
     &  'emissions file and a gridding-transform matrix computed by',
     &  'program GRDAMAT, and produce gridded emissions.  You will be',
     &  'able to choose which variables are gridded; the output time',
     &  'step structure of the gridded file will be the same as the',
     &  'time step structure of the input area source emissions file.',
     &  ' '
        WRITE( *,92000 ) 
     &  'You will need to enter the logical names for the input and',
     &  'output files (and to have set them prior to program launch,',
     &  'using "setenv <logicalname> <pathname>").  Input files must',
     &  'have been sorted as indicated, prior to program execution.',
     &  'You will be prompted for which species to grid (and do not',
     &  'need to grid them all; in particular, you SHOULD NOT grid',
     &  'the FIP-related or ASC-related variables.)',
     &  ' '
        WRITE( *,92000 ) 
     &  'You may use END_OF-FILE (control-D) to quit the program',
     &  'during logical-name entry.  Default responses are indicated',
     &  'in brackets [LIKE THIS].',
     &  ' '

        IF ( .NOT. GETYN( 'Continue with program?', .TRUE. ) ) THEN
            CALL M3EXIT( 'GRDAREA', 0, 0 , 'Ending program.', 2 )
        END IF

C.......   Get environment variable settings
        ST_REPORT = ENVYN( 'ST_REPORT', 'Per-state report writing flag',
     &                      .TRUE., IOS )

        RTIME = ENVINT( 'ST_REPTIME', 
     &                  'Per-State summation reporting time (HHMMSS)',
     &                   230000, IOS )

        GNAME = PROMPTMFILE( 
     &          'Enter logical name for input GRID MATRIX file',
     &          FSREAD3, 'AGMAT', 'GRDAREA' )

        ANAME = PROMPTMFILE( 
     &          'Enter logical name for ANNUAL or TIME-STEPPED ' //
     &          'AREA SOURCE file',
     &          FSREAD3, 'ATMP', 'GRDAREA' )


C.......   Borrow much of file description from input area file.
C.......   Get list of output variables (shrink list in file description
C.......   data structure to just those variables desired):

        IF ( .NOT. DESC3( ANAME ) ) THEN
            CALL M3ERR( 'GRDAREA', 0, 0, 
     &                  'Could not get description of file ' // ANAME,
     &                  .TRUE. )

        ELSEIF( NROWS3D .NE. NASRC ) THEN
            WRITE( MESG, 94010 )
     &      'Dimension mismatch.  AREA SOURCES file:', NROWS3D,
     &      'program:', NASRC

        END IF

        JDATE = SDATE3D
        JTIME = STIME3D
        TSTEP = TSTEP3D

C.........  Based on input files and control status so far, set default
C........   output file names

C.........  Annual inventory file as input
        IF ( TSTEP .EQ. 0 ) THEN
            NSTEPS  = 1
            FDESC3D( 1 ) = 'Area gridded inventory pollutants'
            OUTNAM = 'AGRD'
            REPNAM = 'AGRDREP'
 
C.........  Hourly inventory file as input
        ELSE
            NSTEPS  = MXREC3D
            FDESC3D( 1 ) = 'Area gridded hourly inventory pollutants'
            OUTNAM = 'AGRT'
            REPNAM = 'AGRTREP'
 
        END IF

C.......   Read FIPS code for reporting purpose
        IF( ST_REPORT ) THEN
 
C............  If we've read in temporal file, must get inventory file name
            IF ( TSTEP .NE. 0 ) THEN
 
                INAME = PROMPTMFILE(
     &                  'Enter logical name for AREA INVENTORY file',
     &                  FSREAD3, 'AREA', 'CSGAREA' )
            ELSE
                INAME = ANAME
 
            ENDIF
 
            IF( .NOT. READ3( INAME, 'FIP', ALLAYS3, 0, 0, IFIP ) ) THEN
 
                CALL M3EXIT( 'SPCAMAT', 0, 0,
     &                       'Error reading FIP from file' // ANAME, 2 )
 
            END IF
 
            IF( .NOT.
     &          READ3( INAME, 'INVYR', ALLAYS3, 0, 0, INVYR ) ) THEN
 
                CALL M3EXIT( 'SPCAMAT', 0, 0,
     &               'Error reading INVYR from file' // ANAME, 2 )
 
            END IF
 
            RDEV  = PROMPTFFILE(
     &             'Enter logical name for TOTALS REPORT or "NONE"',
     &             .FALSE., .TRUE., REPNAM, 'CSGAREA' )
 
C............  Create state codes arrays
            LSTA = -9
            K    = 0
            DO 22 I = 1, NASRC
 
                STA = IFIP( I ) / 1000
 
                IF( STA .GT. LSTA ) THEN
 
                    K = K + 1
                    STATES( K ) = STA
                    LSTA        = STA
 
                ENDIF
 
                ISPTR ( I ) = K
 
22          CONTINUE
 
            NSTA = K
 
            IF( NSTA .GT. NASID ) THEN
                WRITE( MESG,94010 )
     &                'Number of states is ', NSTA,
     &                'but maximum (NASID) is', NASID
                CALL M3EXIT( 'CSGAREA', JDATE, JTIME, MESG, 2 )
            ENDIF
 
C............  Initialize emission summary array
            DO 28 I = 1, NASID
                DO 27 S = 1, NIPOL
                    EMSUM( I,S ) = 0.
27              CONTINUE
28          CONTINUE
 
        ENDIF

        IF ( .NOT.DSCGRID( GRDNM, SCRBUF, GDTYP3D, 
     &                     P_ALP3D, P_BET3D, P_GAM3D, XCENT3D, YCENT3D, 
     &                     XORIG3D, YORIG3D, XCELL3D, YCELL3D, 
     &                     NCOLS3D, NROWS3D, NTHIK3D ) ) THEN

            CALL M3EXIT( 'CSGAREA', 0, 0,
     &              'Could not get description of grid ' // GRDNM, 2 )

        END IF

C......... Store variable properties based on DESC3 call
        EFLAG = .FALSE.
        DO  31 V = 1, NIPOL
            L = INDEX1( EINAM( V ), NVARS3D, VNAME3D )
 
            IF( L .LE. 0 ) THEN
 
                EFLAG = .TRUE.
                MESG = 'ERROR: Inventory pollutant ' // EINAM( V ) //
     &                 ' not found in input file!'
                CALL M3MSG2( MESG )
 
            ELSE
                VNAME( V ) = VNAME3D( L )
                UNITS( V ) = UNITS3D( L )
                VTYPE( V ) = VTYPE3D( L )
                VDESC( V ) = VDESC3D( L )
            ENDIF
 
31      CONTINUE

        IF( EFLAG ) THEN
            CALL M3EXIT( 'GRDAREA', 0, 0, 'Bad input file', 2 )
 
        ENDIF

        NVARS3D = NIPOL
        FTYPE3D = GRDDED3   !  shares most of file-description with input file.
        NLAYS3D = 1
        GDNAM3D = GRDNM
        VGTYP3D = IMISS3
        VGTOP3D = BADVAL3

C.......... Set up output variable characteristics
        DO  33 V = 1, NIPOL
            VNAME3D( V ) = VNAME( V )
            UNITS3D( V ) = UNITS( V )
            VTYPE3D( V ) = VTYPE( V )
            VDESC3D( V ) = VDESC( V )
33      CONTINUE

C.......   Get file name; open output gridded area sources file

        ENAME = PROMPTMFILE( 
     &          'Enter logical name for OUTPUT EMIS file or "NONE"',
     &          FSUNKN3, OUTNAM, 'GRDAREA' )

        OUTPUTFLAG = ( ENAME .NE. NONE )


C.......   Read gridding-transform matrix:

        CALL M3MSG2( 'Reading GRIDDING matrix...' )

        IF ( .NOT. READ3( GNAME, 'ALL', 1, 0, 0, NS ) ) THEN

            CALL M3ERR( 'GRDAREA', 0, 0,
     &                  'Could not read gridding matrix from file "' //
     &                  GNAME( 1 : TRIMLEN( GNAME ) ) // '".',
     &                  .TRUE. )

        END IF

C.......   Write out the area source emissions values:

        CALL M3MSG2( 'Calculating gridded area source emissions...' )

        BTIME = JTIME
        LDATE = 0
        DO  199  T = 1, NSTEPS

C...............   If this is a new month, or new day, write message
            IF ( JDATE .GT. 0 .AND. LDATE .NE. JDATE  ) THEN
 
                MESG = 'Processing ' //
     &                 DAYS( WKDAY( JDATE ) ) // MMDDYY( JDATE )
                CALL M3MSG2( MESG( 1:TRIMLEN( MESG ) ) )
 
            END IF

C.............  Write to screen because WRITE3 only writes to LDEV
            IF( JDATE .GT. 0 ) WRITE( *, 93020 ) HHMMSS( JTIME )

            DO  133  V = 1, NIPOL

                IF( READ3( ANAME, EINAM( V ), ALLAYS3, 
     &                     JDATE, JTIME, EMISV ) ) THEN

                    CALL SMATVEC( NASRC, NGRID, NMATX, NS, IS, CS,
     &                            EMISV, EMISG )

                    IF( OUTPUTFLAG ) THEN
                        IF ( .NOT. WRITE3( ENAME, EINAM( V ), 
     &                                     JDATE, JTIME, EMISG ) ) THEN

                            CALL M3EXIT( 'GRDAREA', JDATE, JTIME,
     &                           'Could not write "' // 
     &                           EINAM( V )( 1: TRIMLEN( EINAM( V ))) //
     &                          '" to ' // ENAME, 2 )

                        END IF
                    END IF

                ELSE

                    CALL M3EXIT( 'GRDAREA', JDATE, JTIME, 
     &                           'Could not read "' // 
     &                           EINAM( V )( 1:TRIMLEN( EINAM( V ))) //
     &                           '" from ' // ANAME, 2 )

                END IF  !  if read3() succeeds or not

C.................  When reporting is switched on, calculate emissions
C.................  summaries by state and report to output QA file
C.................  EMSUM must be initialized
                IF( ST_REPORT ) THEN
 
C.....................  For annual values, convert to day-specific
                    IF( TSTEP .EQ. 0 ) THEN
 
                        DO 111 S = 1, NASRC
                            EMISV( S ) = YR2DAY( INVYR( S ) ) *
     &                                   EMISV( S )
111                     CONTINUE
 
                    ENDIF

                    CALL AGRDSUM( NASRC, NGRID, NMATX, NS, IS,
     &                            CS, ISPTR, EMISV, EMSUM( 1,V ) )
                ENDIF

133          CONTINUE   ! End loop on output variables

            LDATE = JDATE
            LTIME = JTIME

            CALL NEXTIME( JDATE, JTIME, TSTEP )

C.............  Write report if report flag is switched on and previous
C.............  time equals the report time
 
            IF( ST_REPORT .AND.
     &          ( LTIME .EQ. RTIME .OR. T .EQ. NSTEPS ) ) THEN
 
                IF( TSTEP .NE. 0 .AND. RDEV .GT. 0)
     &          WRITE( RDEV, 94010 ) DAYS( WKDAY( LDATE ) ) //
     &                               MMDDYY( LDATE ) // ' from',
     &                               BTIME, 'to', LTIME
 
                IF( TSTEP .EQ. 0 ) THEN
                    MESG = 'Area source emissions ' //
     &                     '(tons/ave inventory day) within grid ' //
     &                     GRDNM
                ELSE
                    MESG = 'Area source emissions ' //
     &                     '(tons/specific day) within grid ' // GRDNM
                ENDIF

                IF( RDEV .GT. 0) 
     &          CALL STATELST( RDEV  , MESG , NASID , NSTA , NIPOL,
     &                         'F7.1', EINAM, STATES, EMSUM         )

C.................  Reinitialize emissions summary
                DO 198 I = 1, NASID
                    DO 197 S = 1, NIPOL
                        EMSUM( I,S ) = 0.
197                 CONTINUE
198             CONTINUE
 
C.................  Set begining time for next day
                BTIME = JTIME
 
            ENDIF
 
199     CONTINUE

999   CONTINUE          !  exit program

      CALL M3EXIT( 'GRDAREA', 0, 0, 'Normal completion', 0 )


C******************  FORMAT  STATEMENTS   ******************************

C...........   Error and warning message formats..... 91xxx


C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )


C...........   Formatted file I/O formats............ 93xxx

93010   FORMAT( A16 )

93020   FORMAT( 8X, 'at time ', A8 )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT ( 10 ( A, :, I10, :, 2X ) )

C...........   Miscellaneous formats................. 95xxx

95000   FORMAT ( /5X , A , $ )          !  generic prompt format.


        END

