
        PROGRAM GRDELVPT

C***********************************************************************
C  program body starts at line 161
C
C  DESCRIPTION:
C       Use gridding-transform matrices from program GRDPMAT to grid
C       elevated point source emissions files.
C
C  PRECONDITIONS REQUIRED:  
C       M3IO source-level point source emissions in canonical order used
C       by RAWPOINT
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       Models-3 I/O; GETYN, PROMPTMFILE, TRIMLEN
C
C  REVISION  HISTORY:
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

      IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'PTDIMS3.EXT'   !  point-source dimensioning parameters
        INCLUDE 'CHDIMS3.EXT'   !  emis chem parms (inventory + model)
        INCLUDE 'GRDIMS3.EXT'   !  grid parameters
        INCLUDE 'TMDIMS3.EXT'   !  temporal parameters (for DAYS variable)
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

        EXTERNAL        DSCGRID, ENVINT, ENVYN, GETYN, HHMMSS, INDEX1, 
     &                  MMDDYY, PROMPTFFILE, PROMPTMFILE, TRIMLEN, 
     &                  WKDAY, YR2DAY

C...........   PARAMETER

        CHARACTER*16    NONE
        PARAMETER     ( NONE = 'NONE' )

C...........   LOCAL VARIABLES and their descriptions:
C...........   Point Sources input and output arrays (one variable at a time)

        REAL            EMISV( NPSRC )  !  input emissions.
        REAL            EMISG( NGRID )  !  gridded output emissions

C...........   Inventory arrays
        INTEGER         IFIP  ( NPSRC )   ! FIPs  code for each source
        INTEGER         ISPTR ( NPSRC )   ! Pointer for each state
        INTEGER         ICPTR ( NPSRC )   ! Pointer for each county
        INTEGER         CTYS  ( NPFIP )   ! Unique list of county IDs
        INTEGER         STATES( NPSID )   ! Unique list of state IDs

C...........   Report arrays

        REAL            EMSUM( NPSID, NIPOL ) ! Summed emissions by state
        REAL            CYSUM( NPFIP, NIPOL ) ! Summed emissions by county

C...........   LOCAL VARIABLES and their descriptions:
C...........   Point Sources input and output arrays (one variable at a time)

        INTEGER       NELVIN             !  current source-count
        INTEGER       INDXIN( NPELV )    !  subscript into netCDF inventory file
        REAL          EMISIN( NPELV, NIPOL )  !  emissions values (tons/yr)

        COMMON / ELEVIN / NELVIN, INDXIN, EMISIN

C...........   Gridding Matrix (all active coeffs are 1.0, since 
C...........   any point source "lives in" just one cell)

        INTEGER         NS( NGRID )
        INTEGER         IS( NMATX )

        COMMON  / GRIDMAT / NS, IS

C...........   Control Matrix

        REAL            CC( NPSRC, NIPOL )

C...........   Other local variables

        INTEGER         BTIME   !  beginning time of day in run
        INTEGER         IOS
        INTEGER         NSTA               ! actual number of states
        INTEGER         NCTY               ! actual number of counties
        INTEGER         CTY                ! temporary county code
        INTEGER         STA                ! temporay state code
   
        INTEGER         I, L, K, M, N, S, T, V
        INTEGER         LDATE, LTIME, JDATE, JTIME, TSTEP
        INTEGER         LDEV
        INTEGER         LSTA               ! previous state code
        INTEGER         LCTY               ! previous county code
        INTEGER         NSTEPS
        INTEGER         RDEV               ! output report unit number
        INTEGER         FDEV               ! state, county names unit number
 
        INTEGER         RTIME              ! time for reporting state totals
        INTEGER         VTYPE( NIPOL )

        LOGICAL         EFLAG, OUTPUTFLAG, ST_REPORT, CY_REPORT
        LOGICAL         CONTROLFLAG

        CHARACTER*5     OUTNAM  !  contains default  logical output name
        CHARACTER*16    PNAME   !  logical name for point-source      input file
        CHARACTER*16    GNAME   !  logical name for grid       matrix input file
        CHARACTER*16    CNAME   !  logical name for control    matrix input file
        CHARACTER*16    ENAME   !  logical name for emission output file
        CHARACTER*16    INAME   !  logical name for inventory input file

        CHARACTER*8     REPNAM  !  contains default logical report output name
        CHARACTER*16    SCRBUF  !  scratch name buffer
        CHARACTER*16    VNAME( NIPOL )
        CHARACTER*16    UNITS( NIPOL )
        CHARACTER*16    VDESC( NIPOL )
        CHARACTER*256   MESG


C***********************************************************************
C   begin body of program GRDELVPT

        LDEV = INIT3()

        CALL INITEM( LDEV )

        WRITE( *,92000 ) 
     &  ' ',
     &  'Program GRDELVPT to take a sorted source level point source ',
     &  'emissions file and a gridding-transform matrix computed by',
     &  'program GRDPMAT and produce gridded emissions.  The output ',
     &  'time step structure of the gridded file will be the same as ',
     &  'the time step structure of the input point source emissions',
     &  'file.',
     &  ' '
        WRITE( *,92000 ) 
     &  'You will need to enter the logical names for the input and',
     &  'output files (and to have set them prior to program launch,',
     &  'using "setenv <logicalname> <pathname>").',
     &  ' ',
     &  'You may use END_OF-FILE (control-D) to quit the program',
     &  'during logical-name entry.  Default responses are indicated',
     &  'in brackets [LIKE THIS].',
     &  ' '

        IF ( .NOT. GETYN( 'Continue with program?', .TRUE. ) ) THEN
            CALL M3EXIT( 'GRDELVPT', 0, 0, 'Ending Program.', 2 )
        END IF


C.......   Get environment variable settings
        ST_REPORT = ENVYN( 'ST_REPORT', 'Per-state report writing flag',
     &                      .TRUE., IOS )

        CY_REPORT = ENVYN( 'CY_REPORT', 'Per-state report writing flag',
     &                      .TRUE., IOS )

C.......  If county totals desired, get county names

        IF ( CY_REPORT )
     &          FDEV = PROMPTFFILE(
     &          'Enter name for STATE/COUNTY FIPS CODE file',
     &          .TRUE., .TRUE., 'STCYNAM', 'GRDELVPT' )

        RTIME = ENVINT( 'ST_REPTIME', 
     &                  'Per-State summation reporting time (HHMMSS)',
     &                  230000, IOS )

C.......   Get file names; open input gridding, speciation, 
C.......   optional control matrix files, point source file, 
C.......   optional output file

        GNAME = PROMPTMFILE( 
     &          'Enter logical name for GRID MATRIX file',
     &          FSREAD3, 'PGMAT', 'GRDELVPT' )

        CNAME = PROMPTMFILE(
     &          'Enter logical name for CONTROL MATRIX file ' //
     &          'or "NONE"',
     &          FSREAD3, 'PCMAT', 'GRDELVPT' )

        CONTROLFLAG = ( CNAME .NE. NONE )

        PNAME = PROMPTMFILE( 
     &          'Enter logical name for elevated TIME-STEPPED ' //
     &          'POINT SOURCE file',
     &          FSREAD3, 'ETMP', 'GRDELVPT' )

C.......   Build description of output file, and optionally open it:

        IF ( .NOT. DESC3( PNAME ) ) THEN
            CALL M3EXIT( 'GRDELVPT', 0, 0, 
     &          'Could not get description of file ' // PNAME, 2 )

        ELSEIF( NROWS3D .NE. NPELV ) THEN
            WRITE( MESG, 94010 )
     &      'Dimension mismatch.  Elv POINT SOURCES file:', NROWS3D,
     &      'program:', NPELV
            CALL M3EXIT( 'GRDELVPT', 0, 0, MESG, 2 )

        ENDIF

        JDATE = SDATE3D
        JTIME = STIME3D
        TSTEP = TSTEP3D

C.........  Based on input files and control status so far, set default
C........   output file names

        NSTEPS  = MXREC3D
        FDESC3D( 1 ) = 'Point gridded hourly inventory pollutants'
        OUTNAM = 'EPGRT'
        REPNAM = 'EPGRTREP'

C.......   Read FIPS code for reporting purpose
        IF( ST_REPORT .OR.  CY_REPORT ) THEN

C............  If we've read in temporal file, must get inventory file name
            IF ( TSTEP .NE. 0 ) THEN

                INAME = PROMPTMFILE(
     &                  'Enter logical name for POINT INVENTORY file',
     &                  FSREAD3, 'PNTS', 'GRDELVPT' )
            ELSE
                INAME = PNAME

            ENDIF

            IF( .NOT. READ3( INAME, 'IFIP', ALLAYS3, 0, 0, IFIP ) ) THEN

                CALL M3EXIT( 'GRDELVPT', 0, 0,
     &                       'Error reading FIP from file' // PNAME, 2 )

            END IF

            RDEV  = PROMPTFFILE(
     &             'Enter logical name for TOTALS REPORT or "NONE"',
     &             .FALSE., .TRUE., REPNAM, 'GRDELVPT' )


C............  Create state codes arrays
            LSTA = -9
            K    = 0
            DO 32 I = 1, NPSRC

                STA = IFIP( I ) / 1000

                IF( STA .GT. LSTA ) THEN

                    K = K + 1
                    STATES( K ) = STA
                    LSTA        = STA

                ENDIF

                ISPTR ( I ) = K

32          CONTINUE

            NSTA = K

            IF( NSTA .GT. NPSID ) THEN

                WRITE( MESG,94010 )
     &                'Number of states is ', NSTA,
     &                'but maximum (NPSID) is', NPSID
                CALL M3EXIT( 'GRDPOINT', JDATE, JTIME, MESG, 2 )
            ENDIF

C............  Initialize emission summary array
            DO 28 I = 1, NPSID
                DO 27 S = 1, NIPOL
                    EMSUM( I,S ) = 0.
27              CONTINUE
28          CONTINUE

            IF ( CY_REPORT ) THEN

C............  Create county codes arrays
              LCTY = -9
              K    = 0
              DO 42 I = 1, NPSRC

                CTY = IFIP( I )

                IF( CTY .GT. LCTY ) THEN

                    K = K + 1
                    CTYS( K ) = CTY
                    LCTY      = CTY
                ENDIF

                ICPTR ( I ) = K

42            CONTINUE

              NCTY = K

              IF( NCTY .GT. NPFIP ) THEN
                WRITE( MESG,94010 )
     &                'Number of counties is ', NCTY,
     &                'but maximum (NPFIP) is', NPFIP
                CALL M3EXIT( 'GRDELVPT', JDATE, JTIME, MESG, 2 )
              ENDIF

C............  Initialize cty emission summary array
              DO 38 I = 1, NPFIP
                DO 37 S = 1, NIPOL
                    CYSUM( I,S ) = 0.
37              CONTINUE
38            CONTINUE

            ENDIF

        ENDIF

        IF ( .NOT. DSCGRID( GRDNM, SCRBUF,
     &                      GDTYP3D, P_ALP3D, P_BET3D, P_GAM3D,
     &                      XCENT3D, YCENT3D, XORIG3D, YORIG3D,
     &                      XCELL3D, YCELL3D,
     &                      NCOLS3D, NROWS3D, NTHIK3D ) ) THEN

            CALL M3EXIT( 'GRDELVPT', 0, 0,
     &              'Could not get description of grid ' // GRDNM, 2 )

        END IF

C......... Store variable properties based on DESC3 call
        EFLAG = .FALSE.
        DO  22 V = 1, NIPOL
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

22      CONTINUE

        IF( EFLAG ) THEN
            CALL M3EXIT( 'GRDELVPT', 0, 0, 'Bad input file', 2 )

        ENDIF

        NVARS3D = NIPOL
        FTYPE3D = GRDDED3   !  shares most of file-description with input file.
        NLAYS3D = 1
        GDNAM3D = GRDNM
        VGTYP3D = IMISS3
        VGTOP3D = BADVAL3
        JDATE   = SDATE3D
        JTIME   = STIME3D
        TSTEP   = TSTEP3D

    
C.......... Set up output variable characteristics
        DO  33 V = 1, NIPOL
            VNAME3D( V ) = VNAME( V )
            UNITS3D( V ) = UNITS( V )
            VTYPE3D( V ) = VTYPE( V )
            VDESC3D( V ) = VDESC( V )
33      CONTINUE

        DO 44 I = 2, MXDESC3
            FDESC3D( I ) = ' '
44      CONTINUE

        ENAME = PROMPTMFILE( 
     &          'Enter logical name for OUTPUT EMIS file or "NONE"',
     &          FSUNKN3, OUTNAM, 'GRDELVPT' )

        OUTPUTFLAG = ( ENAME .NE. NONE )

C.......   Read gridding matrix:

        CALL M3MSG2( 'Reading GRIDDING matrix...' )

        IF ( .NOT. READ3( GNAME, 'ALL', 1, 0, 0, NS ) ) THEN

            CALL M3EXIT( 'GRDELVPT', 0, 0,
     &           'Could not read gridding matrix from file "' //
     &           GNAME( 1 : TRIMLEN( GNAME ) ) // '".', 2 )

        END IF      !  if read3() failed for gridding matrix


C.......   If CONTROLFLAG, read and apply control matrix:

        IF( CONTROLFLAG ) THEN

            CALL M3MSG2( 'Reading CONTROL matrix...' )

            IF ( .NOT. DESC3( CNAME ) ) THEN
                WRITE( MESG, '(2A)' )
     &          'Could not get description of file ', CNAME
                CALL M3EXIT( 'GRDELVPT', 0, 0, MESG, 2 )
            ELSE IF ( NROWS3D .NE. NPSRC ) THEN
                WRITE( MESG, 94010 )
     &          'Dimension mismatch.  C-MATRIX file:', NROWS3D,
     &          'program:', NPSRC
                CALL M3EXIT( 'GRDELVPT', 0, 0, MESG, 2 )
            END IF
            DO  66 S = 1, NIPOL

                SCRBUF = 'CONTROL_' // EINAM( S )( 1:8 )
                IF ( .NOT. READ3( CNAME, SCRBUF, 1, 0, 0, CC( 1,S ) )
     &             ) THEN

                    CALL M3EXIT( 'GRDELVPT', 0, 0,
     &                   'Could not read "' // SCRBUF( 1:12 ) //
     &                   '" from CONTROL FILE "' //
     &                   CNAME( 1 : TRIMLEN( CNAME ) ) // '".',
     &                   2 )

                END IF      !  if read3() failed for control matrix

66          CONTINUE

        END IF          !  if controlflag


C.......   Transform and write out point source emissions values:

        CALL M3MSG2( 'Calculating gridded elv pt source emissions...' )

        BTIME = JTIME
        LDATE = 0
        DO  199  T = 1, NSTEPS

C...............   If this is a new month, or new day, write message
            IF ( JDATE .GT. 0 .AND. LDATE. NE. JDATE ) THEN
 
                MESG = 'Processing ' //
     &                 DAYS( WKDAY( JDATE ) ) // MMDDYY( JDATE )
                CALL M3MSG2( MESG( 1:TRIMLEN( MESG ) ) )
 
            END IF
 
C.............  Write to screen because WRITE3 only writes to LDEV
            IF( JDATE .GT. 0 ) WRITE( *, 93020 ) HHMMSS( JTIME )

            IF(  READ3( PNAME, ALLVAR3 , ALLAYS3,
     &               JDATE, JTIME, NELVIN )  ) THEN

            ELSE                !  else read3() failed

                CALL M3EXIT( 'GRDELVPT', JDATE, JTIME,
     &                       'Could not read time step from ' // PNAME,
     &                       2 )

            END IF          !  if read3() succeeds or n

            DO  198  V = 1, NIPOL

C.......  Initialize emissions array

                    DO 200 N = 1, NPSRC

                      EMISV( N ) = 0.0

200                 CONTINUE

C.......  Find the elevated sources only

                    DO 210 N = 1, NELVIN

                       S = INDXIN( N )
                       EMISV( S ) = EMISIN( N, V )

210                 ENDDO

                    IF ( CONTROLFLAG) THEN

                     DO M = 1, NPSRC
                        EMISV( M ) = EMISV( M ) * CC( M, V )
                     ENDDO

                    ENDIF
 
                          
                    CALL PMATVEC( NPSRC, NGRID, NMATX, NS, IS,
     &                            EMISV, EMISG )

                    IF ( OUTPUTFLAG ) THEN

                      IF ( .NOT. WRITE3( ENAME, EINAM( V ),
     &                                   JDATE, JTIME, 
     &                                   EMISG ) ) THEN

       
                        CALL M3EXIT( 'GRDELVPT', JDATE, JTIME,
     &                      'Could not write "' // 
     &                      EINAM( V )( 1: TRIMLEN( EINAM( V ))) // 
     &                      '" to ' // ENAME, 2 )

                      END IF          !  if write3() failed
      
                    END IF


C.........................  When reporting is switched on, calculate emissions
C.........................  summaries by state and county and report to output 
C.........................   QA file
C.........................   CYSUM and EMSUM must be initialized
                    IF( ST_REPORT ) THEN

                            CALL PGRDSUM( NPSRC, NGRID, NMATX, NS, IS,
     &                                    ISPTR, EMISV,
     &                                    EMSUM( 1,V ) )
                    ENDIF

                    IF ( CY_REPORT ) THEN
                            CALL PGRDSUM( NPSRC, NGRID, NMATX, NS, IS,
     &                                    ICPTR, EMISV,
     &                                    CYSUM( 1,V ) )
                    ENDIF

198         CONTINUE            ! end loop on input variables

            LDATE = JDATE
            LTIME = JTIME

            CALL NEXTIME( JDATE, JTIME, TSTEP )

C.............  Write report if report flag is switched on and previous
C.............  time equals the report time

            IF ( ST_REPORT .OR. CY_REPORT ) THEN 

              IF ( LTIME .EQ. RTIME .OR. T .EQ. NSTEPS ) THEN


                IF( TSTEP .NE. 0 .AND. RDEV .GT. 0)
     &          WRITE( RDEV, 94010 ) DAYS( WKDAY( LDATE ) ) //
     &                               MMDDYY( LDATE ) // ' from',
     &                               BTIME, 'to', LTIME

                MESG = 'Elevated Point source emissions ' //
     &                     '(tons/specific day) within grid ' // GRDNM

                IF ( ST_REPORT ) THEN 

                  IF( RDEV .GT. 0)
     &            CALL STATELST( RDEV  , MESG , NPSID , NSTA , NIPOL,
     &                         'F7.1', EINAM, STATES, EMSUM         )

C.................  Reinitialize emissions summary
                  DO 298 I = 1, NPSID
                    DO 297 S = 1, NIPOL
                        EMSUM( I,S ) = 0.
297                 CONTINUE
298               CONTINUE

                ENDIF

                IF ( CY_REPORT ) THEN

                  IF( RDEV .GT. 0)
     &            CALL CYLST( RDEV, FDEV , MESG , NPFIP , NCTY , NIPOL,
     &                         'F7.2', EINAM, CTYS, CYSUM         )

C.................  Reinitialize emissions summary
                  DO 398 I = 1, NPFIP
                    DO 397 S = 1, NIPOL
                        CYSUM( I,S ) = 0.
397                 CONTINUE
398               CONTINUE

                ENDIF

C.................  Set begining time for next day
                BTIME = JTIME

           ENDIF

          ENDIF 
199     CONTINUE          !  end loop on time steps

999     CONTINUE          !  exit program

        CALL M3EXIT( 'GRDELVPT', 0, 0, 
     &               'Normal completion of PROGRAM GRDELVPT', 0 )


C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )

C...........   Formatted file I/O formats............ 93xxx

93020   FORMAT( 8X, 'at time ', A8 )

C...........   Internal buffering formats............ 94xxx
 
94010   FORMAT ( 10 ( A, :, I10, :, 2X ) )

        END

