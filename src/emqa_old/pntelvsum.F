
        PROGRAM PNTELVSUM

C***********************************************************************
C  program body starts at line 145
C
C  DESCRIPTION:
C       Report only on elevated point source emissions 
C       summed by state and county.  Report created using
C       ETMP created by CSGELEV which is used to support
C       generating UAM, UAM-V or CAMx emission input files. 
C       This program will typically be used only when 
C       creating UAM-V or CAMx ready emissions files.
C        
C  PRECONDITIONS REQUIRED:  
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       Models-3 I/O; 
C       ENVYN, FIND1, GETYN, MMDDYY, PROMPTMFILE, PROMPTFFILE, TRIMLEN, 
C       YR2DAY, WKDAY
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
        INCLUDE 'TMDIMS3.EXT'   !  Temporal params
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS and their descriptions:

        LOGICAL         ENVYN
        INTEGER         FIND1
        LOGICAL         GETYN
        CHARACTER*14    MMDDYY
        INTEGER         PROMPTFFILE
        CHARACTER*16    PROMPTMFILE
        INTEGER         TRIMLEN
        REAL            YR2DAY
        INTEGER         WKDAY

        EXTERNAL        ENVYN, FIND1, GETYN, MMDDYY, PROMPTFFILE, 
     &                  PROMPTMFILE, TRIMLEN, YR2DAY
                                            
C...........   PARAMETERS and their descriptions
        INTEGER MAXST, MAXCO ! Max number of states and counties
        PARAMETER ( MAXST = 60, MAXCO = 3500 )
        
        INTEGER NSTPOL, NCOPOL
        PARAMETER ( NSTPOL = NPSID * NIPOL, NCOPOL = NPFIP * NIPOL )

C...........   LOCAL VARIABLES and their descriptions:
C...........   Point Sources input and output arrays

        REAL         COSUM   ( NPFIP, NIPOL ) / NCOPOL * 0.0 /
        REAL         STSUM   ( NPSID, NIPOL ) / NSTPOL * 0.0 /

C...........   LOCAL VARIABLES and their descriptions:
C...........   Point Sources input and output arrays (one variable at a time)

        INTEGER       NELVIN             !  current source-count
        INTEGER       INDXIN( NPELV )    !  subscript into netCDF inventory file
        REAL          EMISIN( NPELV, NIPOL )  !  emissions values (tons/yr)

        COMMON / ELEVIN / NELVIN, INDXIN, EMISIN

        INTEGER      COIDX   ( NPSRC )  !  Index for summing src data to county
        INTEGER      ICOFIP  ( MAXCO )  !  All county FIPS IDs
        INTEGER      ICOLIST ( NPFIP )  !  List of counties in PTMP file
        INTEGER      IFIP    ( NPSRC )
        INTEGER      ISTFIP  ( MAXST )  !  All state FIPS IDs
        INTEGER      ISTLIST ( NPSID )  !  List of states in PTMP file
        INTEGER      STIDX   ( NPSRC )  !  Index for summing src data to state

        CHARACTER*20 STATES  ( MAXST )  !  All state names
        CHARACTER*20 COUNTIES( MAXCO )  !  All county names

C...........   Logical names and unit number

        INTEGER         FDEV    !  unit number of state and county codes/names
        INTEGER         LDEV    !  unit number for log file
        INTEGER         ODEV    ! unit number of output file(s)
 
        CHARACTER*16    ENAME   !  logical name for point-src inventory file
        CHARACTER*16    PNAME   !  logical name for point-src hourly file
        CHARACTER*16    INAME   !  logical name for file to read emissions

C...........   Other local variables

        REAL            EMIS

        INTEGER         ACT_ST, ACT_CO  ! Actual number of states and counties
                                        ! in the input point source file
        INTEGER         I, IDX, J, K, J1, J2, S, T
        INTEGER         IOS             !  I/O status
        INTEGER         FIPNUM          !  Buffer for reading the FIPS code
        INTEGER         JDATE, JTIME    ! current processing date and time
        INTEGER         LASTCO, LASTST  ! The last state and county FIPS numbers
        INTEGER         LDATE, LTIME     ! The previous iteration's date/time
                                        ! encountered when summing the emis.
        INTEGER         NSRCS
        INTEGER         NSTEPS          ! number of time steps for processing
        INTEGER         SDATE, STIME    ! start date and time
        INTEGER         TOT_ST, TOT_CO  ! Total number of states and counties 
                                        ! in the FIPS file.
        INTEGER         TSTEP           ! Time step (HHMMSS) 
 
        LOGICAL         EFLAG           ! Errors flag
        LOGICAL         HRLYTOT         ! true iff PSUMHRLY e.v. set to Y

        CHARACTER*16    NAMBUF
        CHARACTER*20    FIPNAME
        CHARACTER*256   MESG    !  buffer for M3EXIT() messages


C***********************************************************************
C   begin body of program PNTELVSUM

        LDEV = INIT3()

        CALL INITEM( LDEV )

        IF ( .NOT. GETYN( 'Continue with program?', .TRUE. ) ) THEN
            CALL M3EXIT( 'PNTELVSUM', 0, 0, 'Ending Program', 2 )
        END IF

        EFLAG = .FALSE.

C.......   Get file name; open input point source file

        ENAME = PROMPTMFILE( 
     &          'Enter name for POINT SOURCE INVENTORY file',
     &          FSREAD3, 'PNTS', 'PNTELVSUM' )

        IF ( .NOT. DESC3( ENAME ) ) THEN
            CALL M3EXIT( 'PNTELVSUM', 0, 0, 
     &                  'Could not get description of file "' 
     &                  // ENAME( 1:TRIMLEN( ENAME ) ) // '"', 2 )
        END IF

        NSRCS = NROWS3D
        IF ( NSRCS .GT. NPSRC ) THEN
            WRITE ( MESG, 94010 )
     &              'Number of sources (NROWS3D=', NSRCS,
     &              ') exceeds program limits (NPSRC=' , NPSRC, ')'
            CALL M3EXIT( 'PNTELVSUM', 0, 0, MESG, 2 )
        END IF

        PNAME = PROMPTMFILE( 
     &          'Enter name for elevated POINT SOURCE HOURLY file',
     &          FSREAD3, 'ETMP', 'PNTELVSUM' )

C..............  Reset default name of file from which to read the emissions

        INAME = PNAME

C..............  Prompt for hour- or day- totals
        HRLYTOT = ENVYN( 'PSUMHRLY', 'Per-hour report writing flag',
     &                       .FALSE., IOS )

        IF ( .NOT. DESC3( INAME ) ) THEN
            CALL M3EXIT( 'PNTELVSUM', 0, 0, 
     &                  'Could not get description of file "' 
     &                  // INAME( 1:TRIMLEN( INAME ) ) // '"', 2 )
        END IF

        NSTEPS = MAX( 1, MXREC3D )  ! min 1 in case input is time-independent
        SDATE  = SDATE3D
        STIME  = STIME3D
        TSTEP  = TSTEP3D

        FDEV = PROMPTFFILE( 
     &          'Enter name for STATE/COUNTY FIPS CODE file',
     &          .TRUE., .TRUE., 'STCYNAM', 'PNTELVSUM' )

C.........  Open  output file
        IF( HRLYTOT ) THEN          ! Hourly report from ETMP
            NAMBUF = 'EPSUMRPT_HR'
        ELSE
            NAMBUF = 'EPSUMRPT_DAY'  ! Daily report from ETMP
        ENDIF

        ODEV = PROMPTFFILE( 
     &              'Enter name for EMISSIONS SUM REPORT file',
     &              .FALSE., .TRUE., NAMBUF, 'PNTELVSUM' )

C.........  Read state and county names file and create tables
        TOT_ST = 0
        TOT_CO = 0
30      CONTINUE ! Begin "infinite" loop; until EOF.
            READ ( FDEV, 93010, END = 40 ) FIPNUM, FIPNAME
            IF ( FIPNUM .NE. 0 ) THEN  ! If FIPNUM=0, we just read a blank line
                IF ( MOD ( FIPNUM, 1000 ) .EQ. 0 ) THEN
                    TOT_ST = TOT_ST + 1
                    IF ( TOT_ST .LE. MAXST ) THEN
                        ISTFIP ( TOT_ST ) = FIPNUM / 1000
                        STATES ( TOT_ST ) = FIPNAME
                    END IF

                ELSE
                    TOT_CO = TOT_CO + 1
                    IF ( TOT_CO .LE. MAXCO ) THEN
                        ICOFIP ( TOT_CO ) = FIPNUM
                        COUNTIES ( TOT_CO ) = FIPNAME
                    END IF

                END IF
            END IF
        GO TO 30
40      CONTINUE

        IF( TOT_ST .GT. MAXST ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 
     &           'Number of states in STCYNAM (', TOT_ST,
     &           ') exceeds program limits (MAXST=', MAXST, ')'
            CALL M3MSG2( MESG )
        END IF

        IF( TOT_CO .GT. MAXCO ) THEN
            EFLAG = .TRUE.
            WRITE ( MESG, 94010 )
     &            'Number of counties in STCYNAM (', TOT_CO,
     &            ') exceeds program limits (MAXCO=', MAXCO, ')'
            CALL M3MSG2( MESG )
        END IF

        IF( EFLAG ) THEN
            MESG = 'Must recompile program with larger dimensions'
            CALL M3EXIT( 'PNTELVSUM', 0, 0, MESG, 2 )
        ENDIF
        
C.........  Read point source inventory characteristics

        IF ( .NOT. READ3( ENAME, 'IFIP', ALLAYS3,0,0, IFIP ) ) THEN
            CALL M3EXIT( 'REPTPRO', 0, 0, 
     &                   'Error reading FIP from file' // ENAME, 2)
        END IF

C.........  Create index for grouping source emissions by state and county

        LASTCO = 0
        LASTST = 0
        ACT_ST = 0
        ACT_CO = 0

        DO 60 S = 1, NSRCS


C.............  We've moved to a new FIP code (i.e., we've switched counties),
C.............    so bump the number of counties and record the FIP code.

            IF( IFIP( S ) .NE. LASTCO ) THEN
                LASTCO = IFIP ( S )
                ACT_CO = ACT_CO + 1

                IF ( ACT_CO .LE. NPFIP ) THEN
                    ICOLIST ( ACT_CO ) = LASTCO
                ENDIF

C................. If the current state code is not equal to the last state
C.................   code, then bump the number of states and record the code.
                IF ( ( IFIP( S ) / 1000 ) .NE. LASTST ) THEN
                    LASTST = IFIP ( S ) / 1000
                    ACT_ST = ACT_ST + 1

                    IF ( ACT_ST .LE. NPSID ) THEN
                        ISTLIST( ACT_ST ) = LASTST
                    ENDIF

                END IF
            END IF

            STIDX( S ) = ACT_ST
            COIDX( S ) = ACT_CO
            
60      CONTINUE

        IF ( ACT_CO .GT. NPFIP ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 
     &          'Number of counties in point source file (',
     &          ACT_CO, ') exceeds program limits (NPFIP=', NPFIP, ')'
            CALL M3MSG2( MESG )
        END IF

        IF ( ACT_ST .GT. NPSID ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 
     &          'Number of states in point source file (',
     &          ACT_ST, ') exceeds program limits (NPSID=', NPSID, ')'
            CALL M3MSG2( MESG )
        ENDIF

        IF( EFLAG ) THEN
            MESG = 'Must recompile program with larger dimensions'
            CALL M3EXIT( 'PNTELVSUM', 0, 0, MESG, 2 )
        ENDIF

C.........  Loop through time steps (could be 1) and write output
        JDATE = SDATE
        JTIME = STIME
        LDATE = -9
        DO 222 T = 1, NSTEPS

C.............  Read emissions from file.  Read indiv. so NSRCS need not = NPSRC
            DO 50 I = 1, NIPOL

               IF( .NOT. READ3( INAME, ALLVAR3, ALLAYS3, JDATE,
     &                          JTIME, NELVIN )  ) THEN

                   CALL M3EXIT( 'PNTELVSUM', 0, 0, 
     &                       'Error reading ' // EINAM ( I ) //
     &                       ' from file' // ENAME, 2 )

               ELSEIF ( NELVIN .NE. NPELV ) THEN

                   WRITE( MESG, 94010 )
     &             'Dimension mismatch. Elv POINT SOURCES file:', 
     &             NELVIN, 'program:', NPELV
                   CALL M3EXIT( 'PNTELVSUM', 0, 0, MESG, 2 )

               END IF
 
50          CONTINUE

C.............  Initialize summed emissions arrays
            IF( JDATE .NE. LDATE ) THEN
                DO 67 I = 1, NIPOL
                    DO 61 J = 1, NPSID
                        STSUM( J,I ) = 0.
61                  CONTINUE
                    DO 65 J = 1, NPFIP
                        COSUM( J,I ) = 0.
65                  CONTINUE
67              CONTINUE
            ENDIF

C.............  Sum emissions for this hour and all species by state and county
            DO 88 S = 1, NELVIN
                DO 70 I = 1, NIPOL

                    K = INDXIN( S ) 
 
                    J1 = STIDX( K )
                    J2 = COIDX( K )
                    EMIS = EMISIN( S,I ) 
                    STSUM( J1,I ) = STSUM( J1,I ) + EMIS
                    COSUM( J2,I ) = COSUM( J2,I ) + EMIS

70              CONTINUE
88          CONTINUE

C.............  Increment date and time here so can check LDATE
            LDATE = JDATE
            LTIME = JTIME
            CALL NEXTIME( JDATE, JTIME, TSTEP )

C.............  Write state summary emissions
            IF( ( .NOT. HRLYTOT .AND.
     &            ( JDATE .NE. LDATE .OR. T .EQ. NSTEPS ) ) .OR. ! daily ETMP tot
     &          HRLYTOT                              ) THEN ! hourly ETMP total

C.................  Output headers for report
                IF( HRLYTOT ) THEN        ! Hourly totals from ETMP
                    WRITE( ODEV,93020 ) DAYS( WKDAY( LDATE ) ) //
     &                  MMDDYY( LDATE ), LTIME,
     &                  ( EINAM ( I ), I=1,NIPOL )

                ELSE   ! Daily totals from ETMP
                    WRITE( ODEV,93023 ) DAYS( WKDAY( LDATE ) ) //
     &                  MMDDYY( LDATE ) , ( EINAM ( I ), I=1,NIPOL )

                ENDIF

                DO 121 J = 1, ACT_ST
 
                    IDX = FIND1( ISTLIST ( J ), TOT_ST, ISTFIP )
                    IF( IDX .GT. 0 ) THEN
            	        FIPNAME = STATES( IDX )
                    ELSE
                        FIPNAME = '*** no state name'
                    ENDIF

                    WRITE( ODEV, 93055 ) 
     &                  ISTLIST( J ), ( STSUM( J,I ), I = 1, NIPOL ),
     &                  FIPNAME

121             CONTINUE

C.................  Write county break header        
                WRITE( ODEV,93030 ) 

C.................  Write county summary emissions
                LASTST = 0
                DO 181 J = 1, ACT_CO

                    IF ( ( ICOLIST ( J ) / 1000 ) .NE. LASTST ) THEN
                        LASTST = ICOLIST ( J ) / 1000
                        IDX = FIND1 ( LASTST, TOT_ST, ISTFIP )
                        IF ( IDX .GT. 0 ) THEN
                            FIPNAME = STATES ( IDX )
                        ELSE
                            FIPNAME = '*** no state name'
                        ENDIF
 
                        WRITE( ODEV,93045 )
     &                      FIPNAME, ( EINAM ( I ), I=1,NIPOL )
                    END IF

                    IDX = FIND1 ( ICOLIST ( J ), TOT_CO, ICOFIP )
                    IF ( IDX .GT. 0 ) THEN
            	        FIPNAME = COUNTIES ( IDX )
                    ELSE
                        FIPNAME = '*** no county name'
                    ENDIF

                    WRITE( ODEV, 93075 ) 
     &                  ICOLIST ( J ), ( COSUM ( J, I ), I = 1, NIPOL ),
     &                  FIPNAME

181             CONTINUE


C.................  Write space between one hour and the next in output files
                WRITE( ODEV, 94010 ) ' '

            ENDIF  ! whether to output anything for this hour or not

222     CONTINUE
        
        MESG = 'Normal completion of PNTELVSUM'
        CALL M3EXIT( 'PNTELVSUM', 0, 0, MESG, 0 )


C******************  FORMAT  STATEMENTS   ******************************


C...........   Formatted file I/O formats............ 93xxx

93010   FORMAT ( 2X, I5, 1X, A )

93020   FORMAT( 'Elv. pt source emissions in inventory (tons/hour) for ' / 
     &          A ' at ', I6.6, / 
     &          23('-'), ' By State ', 23('-') /
     &          'STATE' / 'CODE     ', <NIPOL>(A10), 'STATE NAME' )

93023   FORMAT( 'Elv. pt source emissions in inventory (tons/day) for ' /
     &          A /
     &          23('-'), ' By State ', 23('-') /
     &          'STATE' / 'CODE     ', <NIPOL>(A10), 'STATE NAME' )

93025   FORMAT( 'Elv. pt source emissions in inventory (tons/day)' //
     &          23('-'), ' By State ', 23('-') /
     &          'STATE' / 'CODE     ', <NIPOL>(A10), 'STATE NAME' )

93030   FORMAT( / 23('-'), ' By County ', 23('-') )

93045   FORMAT( /, 19X, A, /
     &          'FIPs' / 'CODE      ', <NIPOL>(A10), 'COUNTY NAME' )

93055   FORMAT ( I2.2, 2X, <NIPOL>(F10.2), 2X, A )

93070   FORMAT ( I5.5, F10.2, 2X, A )

93075   FORMAT ( I5.5, <NIPOL>(F10.2), 2X, A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT ( 10( A, :, I10, :, 2X ) )

        END
