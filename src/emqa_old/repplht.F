
        PROGRAM REPPLHT

C***********************************************************************
C  program body starts at line 127
C
C  DESCRIPTION:
C       Report minimum and maximum plume heights for utilities with 
C       excessive NOX emissions.
C
C  PRECONDITIONS REQUIRED:  
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       Models-3 I/O; 
C       GETREAL, ENVREAL, PROMPTMFILE, PROMPTFFILE, TRIMLEN, LBLANK, 
C       ENVINT, GETNUM
C
C  REVISION  HISTORY:
C       Prototype  11/96 by Chad Barker
C
C
C***********************************************************************
C 
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C  
C COPYRIGHT (C) 1996, Chad Barker, Turner Engineering
C All Rights Reserved
C  
C Chad Barker, Senior Computer Scientist
C Timothy L. Turner, Owner
C Turner Engineering
C 1310 Holly Creek Lane
C Chapel Hill NC 27516
C (919) 929 6851
C  
C barker@mcnc.org
C  
C Pathname: $Source$
C Last updated: $Date$ 
C  
C**********************************************************************

      IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'PTDIMS3.EXT'   !  point-source dimensioning parameters
        INCLUDE 'CHDIMS3.EXT'   !  emis chem parms (inventory + model)
        INCLUDE 'TMDIMS3.EXT'   !  emis temporal modeling parms
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.
        INCLUDE 'VGRD3.EXT'     !  emis layer params


C...........   EXTERNAL FUNCTIONS and their descriptions:

        REAL            GETREAL, ENVREAL
        CHARACTER*16    PROMPTMFILE
        INTEGER         LBLANK, TRIMLEN, PROMPTFFILE, ENVINT, GETNUM

        EXTERNAL        GETREAL, ENVREAL, PROMPTMFILE, TRIMLEN
        EXTERNAL        LBLANK, PROMPTFFILE, ENVINT, GETNUM
     
                                            
C...........   PARAMETERS and their descriptions
        INTEGER MAXDAYS
        PARAMETER ( MAXDAYS = 7 )

        INTEGER MINVAL, MAXVAL
        PARAMETER ( MINVAL = -1, MAXVAL = 999 )

        INTEGER MAXLARGE, LRGPOL
        PARAMETER ( MAXLARGE = 10000, LRGPOL = NIPOL * MAXLARGE )

C...........   LOCAL VARIABLES and their descriptions:
C...........   Point Sources input and output arrays

        INTEGER     IFIP    ( NPSRC )  !  source FIPS (county) ID
        INTEGER     ISCC    ( NPSRC )  !  source SCC
        INTEGER     IPLANT  ( NPSRC )  !  Plant ID
        INTEGER     ISTACK  ( NPSRC )  !  Stack ID
        REAL        STKHT   ( NPSRC )
        REAL        STKDM   ( NPSRC )
        REAL        STKTK   ( NPSRC )
        REAL        STKVE   ( NPSRC )
        REAL        LFRAC   ( NPSRC, NLAYS ) 
        REAL        NOXSUM  ( NPSRC ) / NPSRC * 0.0 /
        INTEGER     TOOMUCH ( NPSRC ) / NPSRC * 0 /
        REAL        SPECBUF ( NPSRC, NIPOL )
        REAL        SPECSUM ( MAXLARGE, NIPOL ) / LRGPOL * 0.0 /
        INTEGER     MINLAYER( MAXLARGE, MAXDAYS )
        INTEGER     MAXLAYER( MAXLARGE, MAXDAYS )

C...........   Other local variables

        REAL            MAXNOX
        INTEGER         PTSTEP, HTSTEP   !  input time step
        INTEGER         PDATE, PTIME, HDATE, HTIME
        INTEGER         CURDATE, CURTIME, LASTDATE
        INTEGER         ENDDATE, ENDTIME, NDAYS, CURDAY
        INTEGER         NLAYERS, MAXLAY, PNSTEPS, HNSTEPS, NTIMES
        INTEGER         NPTS
        INTEGER         STATUS

        INTEGER         LDEV    !  unit number for log file
        INTEGER         OUNIT, MUNIT

        INTEGER         I, J, TS, LAYER, CURSRC, INDEX
        LOGICAL         HITMAX
        CHARACTER*3     ISMAX

        CHARACTER*16    ENAME   !  logical name for point-source input file
        CHARACTER*16    HNAME   !  " hour-specific input file
        CHARACTER*16    PNAME   !  " plume-rise input file
        CHARACTER*10    BUF1,BUF2
        CHARACTER*16    SPECNM ( NIPOL ) ! Species variable names

        CHARACTER*256   MESG    !  buffer for M3EXIT() messages


C***********************************************************************
C   begin body of program REPPLHT

        LDEV = INIT3()

C.......   Get file name; open input point sources

        ENAME = PROMPTMFILE( 
     &          'Enter name for POINT SOURCE input file',
     &          FSREAD3, 'PNTS', 'REPPLHT' )

        IF ( .NOT. DESC3( ENAME ) ) THEN
            CALL M3EXIT( 'REPPLHT', 0, 0, 
     &                  'Could not get description of file "' 
     &                  // ENAME( 1:TRIMLEN( ENAME ) ) // '"', 2 )
        END IF


        NPTS = NROWS3D
        IF ( NPTS .GT. NPSRC ) THEN
            WRITE ( BUF1, 600 ) NPTS
            WRITE ( BUF2, 600 ) NPSRC
            CALL M3EXIT( 'REPPLHT', 0, 0, 
     &                   'Number of sources (NROWS3D='  // 
     &                   BUF1 ( LBLANK(BUF1):TRIMLEN(BUF1) ) //
     &                   ') exceeds program limits (NPSRC=' //
     &                   BUF2 ( LBLANK(BUF2):TRIMLEN(BUF2) ) // 
     &                   ')', 2 )
        END IF



        HNAME = PROMPTMFILE( 
     &          'Enter name for POINT-SOURCES HOURLY EMISSIONS file',
     &          FSREAD3, 'PTMP', 'REPPLHT' )

        IF ( .NOT. DESC3( HNAME ) ) THEN
            CALL M3EXIT( 'REPTPRO', 0, 0, 
     &                   'Could not get description of file "' 
     &                   // HNAME( 1:TRIMLEN( HNAME ) ) // '"', 2 )
        END IF
        IF ( NROWS3D .NE. NPTS ) THEN
            WRITE ( BUF1, 600 ) NROWS3D
            WRITE ( BUF2, 600 ) NPTS
            CALL M3EXIT( 'REPPLHT', 0, 0, 
     &                   'Number of sources in file ' // HNAME //
     &                   ' (NROWS3D='  // 
     &                   BUF1 ( LBLANK(BUF1):TRIMLEN(BUF1) ) //
     &                   ') does not match number of sources in ' //
     &                   'point source file ' // ENAME // ' (NROWS3D='//
     &                   BUF2 ( LBLANK(BUF2):TRIMLEN(BUF2) ) // 
     &                   ')', 2 )
        END IF

        HDATE = SDATE3D
        HTIME = STIME3D
        HTSTEP = TSTEP3D
        HNSTEPS = MXREC3D

        DO I = 1, NIPOL
            SPECNM ( I ) = VNAME3D ( I )
C            WRITE (*,730) VNAME3D ( I )
C            IF ( VNAME3D ( I ) .EQ. 'NOX' ) THEN
C                NOXVAR = I
C            END IF
        END DO




        PNAME = PROMPTMFILE( 
     &         'Enter name for POINT-SOURCES PLUME RISE FRACTIONS file',
     &          FSREAD3, 'PLAY', 'REPPLHT' )


        IF ( .NOT. DESC3( PNAME ) ) THEN
            CALL M3EXIT( 'REPPLHT', 0, 0, 
     &                   'Could not get description of file "' 
     &                   // PNAME( 1:TRIMLEN( PNAME ) ) // '"', 2 )
        END IF
        IF ( NROWS3D .NE. NPTS ) THEN
            WRITE ( BUF1, 600 ) NROWS3D
            WRITE ( BUF2, 600 ) NPTS
            CALL M3EXIT( 'REPPLHT', 0, 0, 
     &                   'Number of sources in file ' // PNAME //
     &                   ' (NROWS3D='  // 
     &                   BUF1 ( LBLANK(BUF1):TRIMLEN(BUF1) ) //
     &                   ') does not match number of sources in ' //
     &                   'point source file ' // ENAME // ' (NROWS3D='//
     &                   BUF2 ( LBLANK(BUF2):TRIMLEN(BUF2) ) // 
     &                   ')', 2 )
        END IF

        NLAYERS = NLAYS3D

        PDATE = SDATE3D
        PTIME = STIME3D
        PTSTEP = TSTEP3D
        PNSTEPS = MXREC3D - 1

        CALL M3MSG2( 'Sources with NOx totals above the ' //
    &                'cutoff will be included in reports...' )

        MAXNOX = ENVREAL ( 'MAXNOX', 'NOX cutoff point', 1.0, STATUS )
        MAXNOX = GETREAL ( 0.0, 99999.0, MAXNOX, 
     &                     'Enter NOx emissions cutoff [tons/day]' )

        MAXLAY = ENVINT ( 'MAXLAY', 'Max plume rise layer', NLAYERS, 
     &                    STATUS )
        MAXLAY = GETNUM ( 1, NLAYERS, MAXLAY, 
     &                    'Enter maximum plume rise layer' )

C.......... Open the report output files.  PTPDPLH stands for Point-sources
C.......... RePort Daily PLume Height.
        OUNIT = PROMPTFFILE( 
     &          'Enter name for PLUME RISE REPORT file',
     &          .FALSE., .TRUE., 'PRPDPLH', 'REPPLHT' )

        MUNIT = PROMPTFFILE( 
     &          'Enter name for PLUME RISE MAX-HEIGHT REPORT file',
     &          .FALSE., .TRUE., 'PLRMXHT', 'REPPLHT' )



C.......... Get time period of input files and check the number of sources
C.......... Determine how many days are covered in the plume rise fractions
C.......... file.
        ENDDATE = PDATE
        ENDTIME = PTIME
        CALL NEXTIME ( ENDDATE, ENDTIME, ( PTSTEP * PNSTEPS ) )
        NDAYS = ENDDATE - PDATE + 1
        IF ( NDAYS .GT. MAXDAYS ) THEN
            WRITE ( BUF1, 600 ) NDAYS
            WRITE ( BUF2, 600 ) MAXDAYS
            CALL M3EXIT( 'REPPLHT', 0, 0, 
     &                   'Number of days (NDAYS='  // 
     &                   BUF1 ( LBLANK(BUF1):TRIMLEN(BUF1) ) //
     &                   ') exceeds program limits (MAXDAYS=' //
     &                   BUF2 ( LBLANK(BUF2):TRIMLEN(BUF2) ) // 
     &                   ')', 2 )
        END IF


C.......... Set the min and max layer values for each source for each day.
        DO CURDAY = 1, MAXDAYS
            DO J = 1, MAXLARGE
                MINLAYER ( J, CURDAY ) = MAXVAL
                MAXLAYER ( J, CURDAY ) = MINVAL
            END DO
        END DO


C...........   Read in emissions indexing data:  IFIP, ISCC, 
C...........   IPLANT, ISTACK

        IF ( .NOT. READ3( ENAME, 'IFIP', ALLAYS3,0,0, IFIP ) ) THEN
            CALL M3EXIT( 'REPPLHT', 0, 0, 
     &                   'Error reading FIP from file' // ENAME, 2)
        END IF

        IF ( .NOT. READ3( ENAME, 'ISCC', ALLAYS3,0,0, ISCC ) ) THEN
            CALL M3EXIT( 'REPPLHT', 0, 0, 
     &                   'Error reading ISCC from file' // ENAME, 2 )
        END IF

        IF ( .NOT. READ3( ENAME, 'IPLANT', ALLAYS3,0,0, IPLANT ) ) THEN
            CALL M3EXIT( 'REPPLHT', 0, 0, 
     &                   'Error reading IPLANT from file' // ENAME, 2 )
        END IF

        IF ( .NOT. READ3( ENAME, 'ISTACK', ALLAYS3,0,0, ISTACK ) ) THEN
            CALL M3EXIT( 'REPPLHT', 0, 0, 
     &                   'Error reading ISTACK from file' // ENAME, 2 )
        END IF

        IF ( .NOT. READ3( ENAME, 'STKHT', ALLAYS3,0,0, STKHT ) ) THEN
            CALL M3EXIT( 'REPPLHT', 0, 0, 
     &                   'Error reading STKHT from file' // ENAME, 2 )
        END IF

        IF ( .NOT. READ3( ENAME, 'STKDM', ALLAYS3,0,0, STKDM ) ) THEN
            CALL M3EXIT( 'REPPLHT', 0, 0, 
     &                   'Error reading STKDM from file' // ENAME, 2 )
        END IF

        IF ( .NOT. READ3( ENAME, 'STKTK', ALLAYS3,0,0, STKTK ) ) THEN
            CALL M3EXIT( 'REPPLHT', 0, 0, 
     &                   'Error reading STKTK from file' // ENAME, 2 )
        END IF

        IF ( .NOT. READ3( ENAME, 'STKVE', ALLAYS3,0,0, STKVE ) ) THEN
            CALL M3EXIT( 'REPPLHT', 0, 0, 
     &                   'Error reading STKVE from file' // ENAME, 2 )
        END IF


        CURDATE = HDATE
        CURTIME = HTIME
        LASTDATE = CURDATE

C.......... Find the sources whose NOx emissions exceed the threshold on 
C.......... any given day.
        DO TS = 1, HNSTEPS
C.............. Read NOx for the current date and time.
            IF ( .NOT. READ3 ( HNAME, 'NOX', ALLAYS3, CURDATE, 
     &                         CURTIME, SPECBUF ( 1, 1 ) ) ) THEN
                CALL M3EXIT ( 'REPTPRO', 0, 0, 
     &                        'Error reading NOX from file ' // 
     &                        HNAME, 2 )
            END IF

C.............. Sum the current hour's NOx into the total for all sources
            DO J = 1, NPTS
                NOXSUM ( J ) = NOXSUM ( J ) + SPECBUF ( J, 1 )
            END DO

C.............. Increment to the next time step.
            CALL NEXTIME ( CURDATE, CURTIME, HTSTEP )


C.............. If we've gone to a new day (or we're done), then check the total
C.............. NOx for the "previous" day to see if it exceeds the maximum.
            IF ( ( CURDATE .NE. LASTDATE ) .OR. 
     &           ( TS .EQ. HNSTEPS ) ) THEN

               DO J = 1, NPTS
C...................... TOOMUCH(J)>0 if this source has exceeded the maximum
C...................... for any other day.
                    IF ( ( NOXSUM ( J ) .GT. MAXNOX ) .AND.
     &                   ( TOOMUCH ( J ) .EQ. 0 ) ) THEN

C.......................... Check to see if we have too many large sources.
                        IF ( CURSRC .GE. MAXLARGE ) THEN
                            WRITE ( BUF1, 600 ) CURSRC
                            WRITE ( BUF2, 600 ) MAXLARGE
                            CALL M3EXIT( 'REPTPRO', 0, 0, 
     &                          'Number of large sources (' // 
     &                          BUF1 ( LBLANK(BUF1):TRIMLEN(BUF1) ) //
     &                          ') exceeds program limits (MAXLARGE=' //
     &                          BUF2 ( LBLANK(BUF2):TRIMLEN(BUF2) ) // 
     &                          ')', 2 )
                        END IF

C.......................... TOOMUCH(J) becomes an index into the SPECSUM and
C.......................... NOXBYTS arrays.
                        CURSRC = CURSRC + 1
                        TOOMUCH ( J ) = CURSRC
                    END IF

                    NOXSUM ( J ) = 0.0
                END DO
            END IF

        END DO

        CURDATE = HDATE
        CURTIME = HTIME

C.......... Now do the loop again, but this time, for each large source, 
C.......... compute the sum for all emission species and store the NOx 
C.......... emissions for each time step.
        DO TS = 1, HNSTEPS
            IF ( .NOT. READ3 ( HNAME, ALLVAR3, ALLAYS3, CURDATE, 
     &                         CURTIME, SPECBUF ) ) THEN
                CALL M3EXIT ( 'REPTPRO', 0, 0, 
     &                        'Error reading ALLVARS from file ' //
     &                        HNAME, 2 )
            END IF

            DO I = 1, NIPOL
                DO J = 1, NPTS
                    INDEX = TOOMUCH ( J )
                    IF ( INDEX .GT. 0 ) THEN
                        SPECSUM ( INDEX, I ) = 
     &                      SPECSUM ( INDEX, I ) + SPECBUF ( J, I )
                    END IF
                END DO
            END DO

            CALL NEXTIME ( CURDATE, CURTIME, HTSTEP )
        END DO


C.......... Start writing the report headers.
        WRITE (OUNIT, 720) 'Plume rise minimum and maximum height and '
        WRITE (OUNIT, 730) 'plume rise fractions for'
        WRITE (OUNIT, 700) 
     &          'point sources producing > ',MAXNOX,' tons/day of NOX'
        WRITE (OUNIT, 710) 
     &          'The maximum possible emissions layer is: ', NLAYERS
        WRITE (OUNIT, 720) 'The minimum and maximum layers are reported'
        WRITE (OUNIT, 730) ' for the following days: '
        WRITE (OUNIT, 730) '   DATE  BEG TIME  END TIME'
        WRITE (OUNIT, 730) '-------  --------  --------'

        WRITE (MUNIT, 730) 
     &          'Plume rise stack and emissions information for'
        WRITE (MUNIT, 700) 
     &          'point sources producing > ',MAXNOX,' tons/day of NOX'
        WRITE (MUNIT, 705)
     $          'with plume rise maximum layer >= ', MAXLAY
        WRITE (MUNIT, 710) 
     &          'The maximum possible emissions layer is: ', NLAYERS
        WRITE (MUNIT, 720) 'The minimum and maximum layers are reported'
        WRITE (MUNIT, 730) ' for the following days: '
        WRITE (MUNIT, 730) '   DATE  BEG TIME  END TIME'
        WRITE (MUNIT, 730) '-------  --------  --------'

C.......... Write the dates and time-ranges for the report period.
        CURDATE = PDATE
        CURTIME = PTIME
        WRITE (OUNIT, 740) CURDATE, CURTIME
        WRITE (MUNIT, 740) CURDATE, CURTIME
        DO WHILE ( CURDATE .LT. ENDDATE )
            NTIMES = ( ( 240000 - CURTIME ) / PTSTEP ) - 1
            CALL NEXTIME ( CURDATE, CURTIME, PTSTEP * NTIMES )
            WRITE (OUNIT, 745) CURTIME
            WRITE (MUNIT, 745) CURTIME
            CALL NEXTIME ( CURDATE, CURTIME, PTSTEP )
            WRITE (OUNIT, 740) CURDATE, CURTIME
            WRITE (MUNIT, 740) CURDATE, CURTIME
        END DO
        WRITE (OUNIT, 745) ENDTIME
        WRITE (MUNIT, 745) ENDTIME
        
C        PRINT *, 'Finished writing date and time ranges.'

C.......... Write the report headers.
        WRITE (OUNIT, 770) 'MIN', 'MAX'
        WRITE (OUNIT, 720) '    ID   FIP      SCC  PLANT  STACK'
        DO I = 1, NIPOL
            WRITE (OUNIT,760) SPECNM ( I )( 1:TRIMLEN ( SPECNM ( I ) ) )
        END DO
        WRITE (OUNIT, 730) ' LAYER LAYER'
        WRITE (OUNIT, 720) '------ ----- -------- ------ ------'
        DO I = 1, NIPOL
            WRITE ( OUNIT, 720 ) ' --------'
        END DO
        WRITE (OUNIT, 730) ' ----- -----'

        WRITE (MUNIT, 900) '[m]', '[m]', '[K]', '[m/s]'
        WRITE (MUNIT, 930) '    ID   FIP      SCC  PLANT  STACK',
     &          '   STKHT   STKDM   STKTK   STKVE'
        DO I = 1, NIPOL
            WRITE (MUNIT,910) SPECNM ( I )( 1:TRIMLEN ( SPECNM ( I ) ) )
        END DO
        WRITE (MUNIT, 950) '------ ----- -------- ------ ------',
     &          ' ------- ------- ------- -------'
        DO I = 1, NIPOL
            WRITE ( MUNIT, 940 ) ' ---------'
        END DO
        WRITE ( MUNIT, 920 ) ' '


C.......... First, determine the minimum and maximum plume heights
        CURDATE = PDATE
        CURTIME = PTIME
        LASTDATE = CURDATE
        CURDAY = 1
        
        DO TS = 1, PNSTEPS
C.................. Read LFRAC for the current date, time, and layer.
            IF ( .NOT. READ3 ( PNAME, 'LFRAC', ALLAYS3, CURDATE, 
     &                         CURTIME, LFRAC ) ) THEN
                CALL M3EXIT ( 'REPPLHT', 0, 0, 
     &                        'Error reading LFRAC from file ' // 
     &                        PNAME, 2 )
            END IF

C            PRINT '( A, I6.5 )', 'Just read layer for time ', CURTIME
            
            DO LAYER = 1, NLAYERS
C.................. Check to see if this is the min or max layer so far.
                DO J = 1, NPTS
                    INDEX = TOOMUCH ( J )
                    IF ( ( INDEX .GT. 0 ) .AND. 
     &                   ( LFRAC ( J, LAYER ) .GT. 0.0 ) ) THEN
                        MINLAYER ( INDEX, CURDAY ) = 
     &                          MIN( LAYER, 
     &                               MINLAYER ( INDEX, CURDAY ) )
                        MAXLAYER ( INDEX, CURDAY ) = 
     &                          MAX( LAYER, 
     &                               MAXLAYER ( INDEX, CURDAY ) )
                    END IF
                END DO
C                PRINT '( A, I3 )', 'Did MIN/MAX for layer ', LAYER
            END DO

C.............. Increment to the next time step.
            CALL NEXTIME ( CURDATE, CURTIME, PTSTEP )

            IF ( CURDATE .NE. LASTDATE ) THEN
                LASTDATE = CURDATE
                CURDAY = CURDAY + 1
            END IF

        END DO

C        PRINT *, 'Done reading and min/maxing layers'

C.......... For each large source, write the identifying info as well as
C.......... the min and max plume heights for each day.
        DO J = 1, NPTS
            INDEX = TOOMUCH ( J )
            IF ( INDEX .GT. 0 ) THEN
                WRITE(OUNIT,800) 
     &                    J, IFIP(J), ISCC(J), IPLANT(J), ISTACK(J),
     &                    ( SPECSUM ( INDEX, I ), I = 1, NIPOL )

                HITMAX = .FALSE.

                DO CURDAY = 1, NDAYS
                    IF ( MINLAYER ( INDEX, CURDAY ) .GT. NLAYERS ) THEN
                        MINLAYER ( INDEX, CURDAY ) = 0
                    END IF
                    IF ( MAXLAYER ( INDEX, CURDAY ) .LT. 1 ) THEN
                        MAXLAYER ( INDEX, CURDAY ) = 0
                    END IF 
                    IF ( MAXLAYER ( INDEX, CURDAY ) .GE. MAXLAY ) THEN
                        HITMAX = .TRUE.
                        ISMAX = '***'
                    ELSE
                	    ISMAX = '   '
                    END IF
                    IF ( CURDAY .EQ. 1 ) THEN
                        WRITE (OUNIT, 805) MINLAYER ( INDEX, CURDAY ), 
     &                                     MAXLAYER ( INDEX, CURDAY ), 
     &                                     ISMAX
                    ELSE
                        WRITE (OUNIT, 810) MINLAYER ( INDEX, CURDAY ), 
     &                                     MAXLAYER ( INDEX, CURDAY ), 
     &                                     ISMAX
                    END IF
            	END DO

            	IF ( HITMAX ) THEN
                    WRITE(MUNIT,820) 
     &                    J, IFIP(J), ISCC(J), IPLANT(J), ISTACK(J),
     &                    STKHT(J), STKDM(J), STKTK(J), STKVE(J),
     &                    ( SPECSUM ( INDEX, I ), I = 1, NIPOL )
                END IF
            END IF
        END DO


999     CONTINUE          !  exit program

        MESG = 'Normal completion of REPPLHT'
        CALL M3EXIT( 'REPPLHT', 0, 0, MESG, 0 )


C******************  FORMAT  STATEMENTS   ******************************

600     FORMAT ( I10 )
700     FORMAT ( 3X, A, F6.3, A )
705     FORMAT ( 3X, A, I3, A )
710     FORMAT ( /, A, I3, / )
720     FORMAT ( A, $ )
730     FORMAT ( A )
740     FORMAT ( I7, 4X, I6.5, $ )
745     FORMAT ( 4X, I6.5 )
760     FORMAT ( A9, $ )
770     FORMAT ( /, 65X, A, 3X, A )
800     FORMAT ( I6, 1X, I5.5, 1X, I8.8, 1X, I6, 1X, I6, <NIPOL>F9.3, $)
805     FORMAT ( 1X, I5, 1X, I5, 1X, A3 )
810     FORMAT ( 63X, I5, 1X, I5, 1X, A3 )
820     FORMAT ( I6, 1X, I5.5, 1X, I8.8, 1X, I6, 1X, I6, 1X,
     &           F7.2, 1X, F7.2, 1X, F7.2, 1X, F7.2, <NIPOL>F10.4 )
900     FORMAT ( /, 39X, A, 5X, A, 5X, A, 4X, A )
910     FORMAT ( A10, $ )
920     FORMAT ( A )
930     FORMAT ( A, A, $ )
940     FORMAT ( A, $ )
950     FORMAT ( /, A, A, $ )
        END
