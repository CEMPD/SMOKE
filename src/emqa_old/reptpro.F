
        PROGRAM REPTPRO

C***********************************************************************
C  program body starts at line 111
C
C  DESCRIPTION:
C       Report temporal profiles of utilities with excessive NOX emissions.
C
C  PRECONDITIONS REQUIRED:  
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       Models-3 I/O; 
C       GETREAL, ENVREAL, PROMPTMFILE, PROMPTFFILE, TRIMLEN, LBLANK
C
C  REVISION  HISTORY:
C       Prototype  11/96 by Chad Barker
C
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1998, Turner Engineering & MCNC - NCSC
C All Rights Reserved
C
C See file COPYRIGHT for conditions of use.
C
C Turner Engineering               Environmental Programs Group
C 105 Sudbury Lane                 MCNC--North Carolina Supercomputing Center
C Chapel Hill NC 27516             P.O. Box 12889
C (919) 929 6851                   Research Triangle Park, NC  27709-2889
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
        INCLUDE 'TMDIMS3.EXT'   !  emis temporal modeling parms
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.


C...........   EXTERNAL FUNCTIONS and their descriptions:

        REAL            GETREAL, ENVREAL
        CHARACTER*16    PROMPTMFILE
        INTEGER         LBLANK, TRIMLEN, PROMPTFFILE

        EXTERNAL        GETREAL, ENVREAL, PROMPTMFILE, TRIMLEN
        EXTERNAL        LBLANK, PROMPTFFILE
     
                                            
C...........   PARAMETERS and their descriptions
        INTEGER MAXLARGE, MAXSTEPS
        PARAMETER ( MAXLARGE = 10000, MAXSTEPS = 168 ) ! 168 hrs = 7 days

        INTEGER LRGPOL
        PARAMETER ( LRGPOL = NIPOL * MAXLARGE )
        
C...........   LOCAL VARIABLES and their descriptions:
C...........   Point Sources input and output arrays

        INTEGER     IFIP    ( NPSRC )  !  source FIPS (county) ID
        INTEGER     ISCC    ( NPSRC )  !  source SCC
        INTEGER     IPLANT  ( NPSRC )  !  Plant ID
        INTEGER     ISTACK  ( NPSRC )  !  Stack ID
        INTEGER     TOOMUCH ( NPSRC ) / NPSRC * 0 /
        REAL        SPECBUF ( NPSRC, NIPOL )
        REAL        NOXSUM  ( NPSRC ) / NPSRC * 0.0 /
        REAL        SPECSUM ( MAXLARGE, NIPOL ) / LRGPOL * 0.0 /
        REAL        NOXBYTS ( MAXLARGE, MAXSTEPS )

C...........   Other local variables

        REAL            MAXNOX
        INTEGER         TSTEP   !  input time step
        INTEGER         JDATE, JTIME
        INTEGER         CURDATE, CURTIME, LASTDATE, CURSRC
        INTEGER         NSTEPS, NTIMES, INDEX
        INTEGER         NPTS, NOXVAR
        INTEGER         STATUS, LBL

        INTEGER         LDEV    !  unit number for log file
        INTEGER         OUNIT

        INTEGER         I, J, TS

        CHARACTER*16    ENAME   !  logical name for point-source input file
        CHARACTER*16    HNAME   !  " hour-specific input file
        CHARACTER*10    BUF1,BUF2

        CHARACTER*256   MESG    !  buffer for M3EXIT() messages


C***********************************************************************
C   begin body of program REPTPRO

        LDEV = INIT3()

        CALL INITEM( LDEV )

C.......   Get file name; open input point sources

        ENAME = PROMPTMFILE( 
     &          'Enter name for POINT SOURCE input file',
     &          FSREAD3, 'PNTS', 'REPTPRO' )

        HNAME = PROMPTMFILE( 
     &          'Enter name for POINT-SOURCES HOURLY EMISSIONS file',
     &          FSREAD3, 'PTMP', 'REPTPRO' )

        MAXNOX = ENVREAL ( 'MAXNOX', 'NOX cutoff point', 1.0, STATUS )
        MAXNOX = GETREAL ( 0.0, 99999.0, MAXNOX, 
     &                     'Enter cutoff point for NOX emissions' )

C.......... Open the report output file.
        OUNIT = PROMPTFFILE( 
     &          'Enter name for TEMPORAL PROFILES REPORT file',
     &          .FALSE., .TRUE., 'PRPHNOX', 'REPTPRO' )



        IF ( .NOT. DESC3( ENAME ) ) THEN
            CALL M3EXIT( 'REPTPRO', 0, 0, 
     &                  'Could not get description of file "' 
     &                  // ENAME( 1:TRIMLEN( ENAME ) ) // '"', 2 )
        END IF


        NPTS = NROWS3D
        IF ( NPTS .GT. NPSRC ) THEN
            WRITE ( BUF1, 600 ) NPTS
            WRITE ( BUF2, 600 ) NPSRC
            CALL M3EXIT( 'REPTPRO', 0, 0, 
     &                   'Number of sources (NROWS3D='  // 
     &                   BUF1 ( LBLANK(BUF1):TRIMLEN(BUF1) ) //
     &                   ') exceeds program limits (NPSRC=' //
     &                   BUF2 ( LBLANK(BUF2):TRIMLEN(BUF2) ) // 
     &                   ')', 2 )
        END IF


C.......... Get time period of input file and check the number of time steps
        IF ( .NOT. DESC3( HNAME ) ) THEN
            CALL M3EXIT( 'REPTPRO', 0, 0, 
     &                   'Could not get description of file "' 
     &                   // HNAME( 1:TRIMLEN( HNAME ) ) // '"', 2 )
        END IF

        JDATE = SDATE3D
        JTIME = STIME3D
        TSTEP = TSTEP3D
        NSTEPS = MXREC3D
        IF ( NSTEPS .GT. MAXSTEPS ) THEN
            WRITE ( BUF1, 600 ) NSTEPS
            WRITE ( BUF2, 600 ) MAXSTEPS
            CALL M3EXIT( 'REPTPRO', 0, 0, 
     &                   'Number of time steps (MXREC3D='  // 
     &                   BUF1 ( LBLANK(BUF1):TRIMLEN(BUF1) ) //
     &                   ') exceeds program limits (MAXSTEPS=' //
     &                   BUF2 ( LBLANK(BUF2):TRIMLEN(BUF2) ) // 
     &                   ')', 2 )
        END IF

C.......... Find out which of the policy variables is NOx
        DO I = 1, NIPOL
C            WRITE (*,730) VNAME3D ( I )
            IF ( VNAME3D ( I ) .EQ. 'NOX' ) THEN
                NOXVAR = I
            END IF
        END DO

C...........   Read in emissions indexing data:  IFIP, ISCC, 
C...........   IPLANT, ISTACK

        IF ( .NOT. READ3( ENAME, 'IFIP', ALLAYS3,0,0, IFIP ) ) THEN
            CALL M3EXIT( 'REPTPRO', 0, 0, 
     &                   'Error reading FIP from file' // ENAME, 2)
        END IF

        IF ( .NOT. READ3( ENAME, 'ISCC', ALLAYS3,0,0, ISCC ) ) THEN
            CALL M3EXIT( 'REPTPRO', 0, 0, 
     &                   'Error reading ISCC from file' // ENAME, 2 )
        END IF

        IF ( .NOT. READ3( ENAME, 'IPLANT', ALLAYS3,0,0, IPLANT ) ) THEN
            CALL M3EXIT( 'REPTPRO', 0, 0, 
     &                   'Error reading IPLANT from file' // ENAME, 2 )
        END IF

        IF ( .NOT. READ3( ENAME, 'ISTACK', ALLAYS3,0,0, ISTACK ) ) THEN
            CALL M3EXIT( 'REPTPRO', 0, 0, 
     &                   'Error reading ISTACK from file' // ENAME, 2 )
        END IF


        CURDATE = JDATE
        CURTIME = JTIME
        LASTDATE = CURDATE
        CURSRC = 0
        NTIMES = 0

C.......... Start writing the report headers.
        WRITE (OUNIT, 730) 'Temporal profiles of NOX emissions for'
        WRITE (OUNIT, 700) 
     &          'point sources producing > ',MAXNOX,' tons/day of NOX'
        WRITE (OUNIT, 720) 
     &   'The emissions are reported for the following dates and times:'
        WRITE (OUNIT, 740) CURDATE
        LBL = 750
C.......... First, find the point sources that exceed the maximum NOX
C.......... emissions for any day.

        DO TS = 1, NSTEPS
            ! Write current time in the header
            IF ( LBL .EQ. 750 ) THEN ! for some reason, I couldn't get 
                                     ! WRITE(OUNIT,LBL) ... to work.
                WRITE (OUNIT,750) CURTIME 
            ELSE
                WRITE (OUNIT,760) CURTIME ! Write current time in the header
                LBL = 750                  ! Reset the format to the default
            ENDIF
            NTIMES = NTIMES + 1           ! Count the times we've written

C.............. Read NOx for the current date and time.
            IF ( .NOT. READ3 ( HNAME, 'NOX', ALLAYS3, CURDATE, 
     &                         CURTIME, SPECBUF (1, 1) ) ) THEN
                CALL M3EXIT ( 'REPTPRO', 0, 0, 
     &                        'Error reading NOX from file ' // 
     &                        HNAME, 2 )
            END IF

C.............. Sum the current hour's NOx into the total for all sources
            DO J = 1, NPTS
                NOXSUM ( J ) = NOXSUM ( J ) + SPECBUF ( J, 1 )
            END DO

C.............. Increment to the next time step.
            CALL NEXTIME ( CURDATE, CURTIME, TSTEP )


C.............. If we've gone to a new day (or we're done), then check the total
C.............. NOx for the "previous" day to see if it exceeds the maximum.
            IF ( ( CURDATE .NE. LASTDATE ) .OR. ( TS .EQ. NSTEPS ) )THEN
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
                END DO

C.................. Reset NOXSUM, for the next day.
                DO J = 1, NPTS
                    NOXSUM ( J ) = 0.0
                END DO

C.................. If we're not at the end, then write the new date
C.................. at the beginning of the next line.
                IF ( TS .LT. NSTEPS ) THEN
                    WRITE (OUNIT, 740) CURDATE
                    LASTDATE = CURDATE
                    NTIMES = 0
                END IF
            ELSE IF ( NTIMES .EQ. 8 ) THEN
C.................. If we've put eight times on the current line, then 
C.................. move to the beginning of the next line.
                LBL = 760
                NTIMES = 0
            END IF

        END DO


        CURDATE = JDATE
        CURTIME = JTIME

C.......... Now do the loop again, but this time, for each large source, 
C.......... compute the sum for all emission species and store the NOx 
C.......... emissions for each time step.
        DO TS = 1, NSTEPS
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

            DO J = 1, NPTS
                INDEX = TOOMUCH ( J )
                IF ( INDEX .GT. 0 ) THEN
                    NOXBYTS ( INDEX, TS ) = SPECBUF ( J, NOXVAR )
                END IF
            END DO

            CALL NEXTIME ( CURDATE, CURTIME, TSTEP )
        END DO


C.......... Write the report headers.
        WRITE (OUNIT, 710) 
     &          '    ID    FIP       SCC   PLANT   STACK',
     &          ( VNAME3D ( I )( 1:TRIMLEN ( VNAME3D ( I ) ) ), 
     &            I = 1, NIPOL )
        WRITE (OUNIT, 720) 
     &          '------  -----  --------  ------  ------'
        DO I = 1, NIPOL
            WRITE ( OUNIT, 720 ) '  ---------'
        END DO
        WRITE ( OUNIT, 730 ) ''


C.......... For each large source, write the identifying info as well as
C.......... the NOx emissions for each time step.
        DO J = 1, NPTS
            INDEX = TOOMUCH ( J )
            IF ( INDEX .GT. 0 ) THEN
                WRITE(OUNIT,800) 
     &                    J, IFIP(J), ISCC(J), IPLANT(J), ISTACK(J),
     &                    ( SPECSUM ( INDEX, I ), I = 1, NIPOL )

        		CURDATE = JDATE
        		CURTIME = JTIME
        		LASTDATE = 0
        		NTIMES = 8

                DO TS = 1, NSTEPS
                    IF ( ( NTIMES .EQ. 8 ) .OR. 
     &                   ( CURDATE .NE. LASTDATE ) ) THEN
                        NTIMES = 0
                        IF ( CURDATE .EQ. LASTDATE ) THEN
                            WRITE (OUNIT, 805) NOXBYTS ( INDEX, TS )
                        ELSE
                            LASTDATE = CURDATE
                            WRITE (OUNIT, 810) 
     &                            CURDATE, NOXBYTS ( INDEX, TS )
                        END IF
                    ELSE
                        WRITE (OUNIT, 820) NOXBYTS ( INDEX, TS )
                    END IF
                    NTIMES = NTIMES + 1

                    CALL NEXTIME ( CURDATE, CURTIME, TSTEP )
            	END DO

                WRITE (OUNIT, 830)
            ENDIF
        END DO


999     CONTINUE          !  exit program

        MESG = 'Normal completion of REPTPRO'
        CALL M3EXIT( 'REPTPRO', 0, 0, MESG, 0 )


C******************  FORMAT  STATEMENTS   ******************************

600     FORMAT ( I10 )
700     FORMAT ( 3X, A, F6.3, A, / )
710     FORMAT ( //, A, <NIPOL>A11 )
720     FORMAT ( A, $ )
730     FORMAT ( A )
740     FORMAT ( /, I7, $ )
750     FORMAT ( 1X, I6.6, $ )
760     FORMAT ( /, 8X, I6.6, $ )
800     FORMAT ( I6, 2X, I5.5, 2X, I8.8, 2X, I6, 2X, I6, <NIPOL>F11.4,$)
805     FORMAT ( /, 15X, F7.4, $ )
810     FORMAT ( /, 6X, I7, 2X, F7.4, $ )
820     FORMAT ( 1X, F7.4, $ )
830     FORMAT ( / )
        END

