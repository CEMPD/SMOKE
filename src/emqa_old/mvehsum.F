
        PROGRAM MVEHSUM

C***********************************************************************
C  program body starts at line 145
C
C  DESCRIPTION:
C       Report mobile-source emissions summed by vehicle type and emissions
C       type for selected counties.
C
C  PRECONDITIONS REQUIRED:  
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       Models-3 I/O; 
C       PROMPTMFILE, PROMPTFFILE, TRIMLEN, LBLANK, 
C       GETREAL, ENVREAL, FIND1
C
C  REVISION  HISTORY:
C       Prototype  7/97 by Marc Houyoux
C
C
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1998, MCNC--North Carolina Supercomputing Center
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

        INCLUDE 'MBDIMS3.EXT'   !  mobile-source dimensioning parameters
        INCLUDE 'CHDIMS3.EXT'   !  emis chem parms (inventory + model)
        INCLUDE 'TMDIMS3.EXT'   !  Temporal emis parms (needs DAYS)
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.


C...........   EXTERNAL FUNCTIONS and their descriptions:

        CHARACTER*2     CRLF
        REAL            ENVREAL
        INTEGER         FIND1
        INTEGER         GETNUM
        CHARACTER*10    HHMMSS
        INTEGER         LBLANK
        CHARACTER*14    MMDDYY
        INTEGER         PROMPTFFILE
        CHARACTER*16    PROMPTMFILE
        INTEGER         TRIMLEN
        INTEGER         WKDAY

        EXTERNAL        ENVREAL, FIND1, GETNUM, HHMMSS, LBLANK, 
     &                  MMDDYY, PROMPTFFILE, PROMPTMFILE, TRIMLEN,
     &                  WKDAY
     
                                            
C...........   PARAMETERS and their descriptions
        INTEGER MAXST, MAXCO ! Max number of states and counties
        PARAMETER ( MAXST = 60, MAXCO = 3500 )

        INTEGER MYMXSTA, MYMXFIP ! Local version of vars in MBDIMS3.EXT
        PARAMETER ( MYMXSTA = 60, MYMXFIP = 3500 )
        
        CHARACTER*5  BLANK5
        INTEGER      MXDAYS    ! Maximum no. of days
        INTEGER      MXOFIP    ! Maximum no. output FIPs codes
        PARAMETER  ( BLANK5 = ' ',
     &               MXDAYS = 10 ,
     &               MXOFIP = 50   )

C...........   LOCAL VARIABLES and their descriptions:
C...........   Point Sources input and output arrays

        INTEGER      ISTFIP  ( MAXST )  !  All state FIPS IDs
        CHARACTER*20 STATES  ( MAXST )  !  All state names
        INTEGER      ICOFIP  ( MAXCO )  !  All county FIPS IDs
        CHARACTER*20 COUNTIES( MAXCO )  !  All county names
        INTEGER      TOT_ST, TOT_CO     ! Total number of states and counties 
                                        ! in the FIPS file.

        INTEGER      IFIP    ( NMSRC )
        INTEGER      EIDX    ( NMSRC )
        
        REAL         TMAT ( NMSRC, NVTYPE, MXEMIS ) ! For reading one time step

        INTEGER      INDEXA( MXOFIP )
        INTEGER      OFIPA ( MXOFIP )
        INTEGER      OFIP  ( MXOFIP )
        INTEGER      BSRC  ( MXOFIP )
        INTEGER      ESRC  ( MXOFIP )

        REAL         EMSUM( MXOFIP, MXEMIS, NVTYPE, MXDAYS )

C...........   Other local variables

        INTEGER         NSRCS

        INTEGER         LDEV    !  unit number for log file
        INTEGER         ODEV
        INTEGER         FUNIT

        INTEGER         D, I, J, K, S, V
        INTEGER         FIP
        INTEGER         ICNT
        INTEGER         JDATE, JTIME
        INTEGER         LDATE
        INTEGER         NETYPE
        INTEGER         NOFIP   !  number of output FIPs codes
        INTEGER         NDAYS, NS
        INTEGER         PFIP    !  FIPs code from previous iteration
        INTEGER         SDATE, STIME, TSTEP, NSTEPS
        CHARACTER*7     EMTNAM( MXEMIS ) ! Scenario-specific emission type names
        CHARACTER*16    ENAME   !  logical name for mobile-source input file
        CHARACTER*16    HNAME   !   " for mobile-source emissions input file
        CHARACTER*10    BUF1,BUF2
        INTEGER         FIPNUM  !  Buffer for reading the FIPS code
        CHARACTER*20    FIPNAME !  Buffer for reading the state/county name

        CHARACTER*30    FMTBUF  !  format buffer
        CHARACTER*256   MESG    !  buffer for M3EXIT() messages


C***********************************************************************
C   begin body of program MVEHSUM

        LDEV = INIT3()

        CALL INITEM( LDEV )

C.......   Get file name; open input area source

        ENAME = PROMPTMFILE( 
     &          'Enter name for MOBILE SOURCE input file',
     &          FSREAD3, 'MOBL', 'MVEHSUM' )

        IF ( .NOT. DESC3( ENAME ) ) THEN
            CALL M3EXIT( 'MVEHSUM', 0, 0, 
     &                  'Could not get description of file "' 
     &                  // ENAME( 1:TRIMLEN( ENAME ) ) // '"', 2 )
        END IF

        NSRCS = NROWS3D
        IF ( NSRCS .GT. NMSRC ) THEN
            WRITE ( BUF1, 600 ) NSRCS
            WRITE ( BUF2, 600 ) NMSRC
            CALL M3EXIT( 'MVEHSUM', 0, 0, 
     &                   'Number of sources (NROWS3D='  // 
     &                   BUF1 ( LBLANK(BUF1):TRIMLEN(BUF1) ) //
     &                   ') exceeds program limits (NMSRC=' //
     &                   BUF2 ( LBLANK(BUF2):TRIMLEN(BUF2) ) // 
     &                   ')', 2 )
        END IF

        HNAME = PROMPTMFILE( 
     &          'Enter name for MOBILE SOURCE EMISSIONS input file',
     &          FSREAD3, 'MTMP', 'MVEHSUM' )

        IF ( .NOT. DESC3( HNAME ) ) THEN
            CALL M3EXIT( 'MVEHSUM', 0, 0, 
     &                  'Could not get description of file "' 
     &                  // HNAME( 1:TRIMLEN( HNAME ) ) // '"', 2 )
        END IF

        SDATE = SDATE3D
        STIME = STIME3D
        TSTEP = TSTEP3D
        NSTEPS = MXREC3D
        NETYPE = NVARS3D
        IF ( NSRCS .NE. NROWS3D ) THEN
            WRITE ( BUF1, 600 ) NSRCS
            WRITE ( BUF2, 600 ) NROWS3D
            CALL M3EXIT( 'MVEHSUM', 0, 0, 
     &                   'Number of sources in MOBL (NROWS3D='  // 
     &                   BUF1 ( LBLANK(BUF1):TRIMLEN(BUF1) ) //
     &                   ') does not match number in MTMP (NROWS3D=' //
     &                   BUF2 ( LBLANK(BUF2):TRIMLEN(BUF2) ) // 
     &                   ')', 2 )
        END IF

C.........  Store variable names

        DO 50 I = 1, NETYPE
            EMTNAM( I ) = VNAME3D( I )
50      CONTINUE

        FUNIT = PROMPTFFILE( 
     &          'Enter name for STATE/COUNTY FIPS CODE file',
     &          .TRUE., .TRUE., 'STCYNAM', 'MVEHSUM' )


C..........  Read FIPs code for generating output report

        I = 0
70      CONTINUE
            FIP   = 0
            FIP   = GETNUM( 0, 99999, FIP, 
     &                'Enter next FIPs code [XXXXX] or 0 to end' )

            IF( FIP .EQ. 0 ) GO TO 82  ! to end of loop

cmh NOTE: later add ability to do state codes
cmh NOTE: later add assurances that state/FIPs codes are all unique and
cmh       that the list is filled out when a state code is given

            I = I + 1

            IF( I .LE. MXOFIP ) THEN
                 OFIPA ( I ) = FIP
                 INDEXA( I ) = I
            ENDIF

        GO TO 70  ! To head of loop

82      CONTINUE

        IF( I .EQ. 0 ) THEN

            MESG = 'No FIPs codes entered for generating output'
            CALL M3EXIT( 'MVEHSUM', 0, 0, MESG, 2 )

        ELSEIF( I .GT. MXOFIP ) THEN

            WRITE( MESG, 94010 ) 
     &             'Maximum output FIPs code (MXOFIP)=', MXOFIP,
     &             CRLF() // BLANK5 //
     &             'Actual output FIPs code          =', I, 
     &             CRLF() // BLANK5 //
     &             'Reset MXOFIP and recompile.'
            CALL M3EXIT( 'MVEHSUM', 0, 0, MESG, 2 )

        ENDIF

        NOFIP = I

C..........  Prompt for output file name

        ODEV = PROMPTFFILE( 
     &             'Enter name for detailed report output file',
     &             .FALSE., .TRUE., 'MDTLREP', 'MVEHSUM' )

        ! Read complete list of state & County fips codes and names
        TOT_ST = 0
        TOT_CO = 0
90      CONTINUE ! Begin "infinite" loop; until EOF.
            READ ( FUNIT, 700, END = 91 ) FIPNUM, FIPNAME
            IF ( FIPNUM .NE. 0 ) THEN  ! If FIPNUM=0, we just read a blank line
                IF ( MOD ( FIPNUM, 1000 ) .EQ. 0 ) THEN
                    TOT_ST = TOT_ST + 1
                    IF ( TOT_ST .GT. MAXST ) THEN
                        WRITE ( BUF1, 600 ) TOT_ST
                        WRITE ( BUF2, 600 ) MAXST
                        CALL M3EXIT( 'MVEHSUM', 0, 0, 
     &                      'Number of states in STCYNAM (' // 
     &                      BUF1 ( LBLANK(BUF1):TRIMLEN(BUF1) ) //
     &                      ') exceeds program limits (MAXST=' //
     &                      BUF2 ( LBLANK(BUF2):TRIMLEN(BUF2) ) // 
     &                      ')', 2 )
                    END IF

                    ISTFIP ( TOT_ST ) = FIPNUM / 1000
                    STATES ( TOT_ST ) = FIPNAME

                ELSE
                    TOT_CO = TOT_CO + 1
                    IF ( TOT_CO .GT. MAXCO ) THEN
                        WRITE ( BUF1, 600 ) TOT_CO
                        WRITE ( BUF2, 600 ) MAXCO
                        CALL M3EXIT( 'MVEHSUM', 0, 0, 
     &                      'Number of counties in STCYNAM (' // 
     &                      BUF1 ( LBLANK(BUF1):TRIMLEN(BUF1) ) //
     &                      ') exceeds program limits (MAXCO=' //
     &                      BUF2 ( LBLANK(BUF2):TRIMLEN(BUF2) ) // 
     &                      ')', 2 )
                    END IF

                    ICOFIP ( TOT_CO ) = FIPNUM
                    COUNTIES ( TOT_CO ) = FIPNAME

                END IF
            END IF
        GO TO 90
91      CONTINUE

C.........  Read variables from MOBL

        IF ( .NOT. READ3 ( ENAME, 'IFIP', ALLAYS3, 0, 0, IFIP ) ) THEN
            CALL M3EXIT( 'MVEHSUM', 0, 0, 
     &                   'Error reading IFIP from file' // ENAME, 2)
        END IF

c       IF ( .NOT. READ3 ( ENAME, 'VMT', ALLAYS3, 0, 0, VMT ) ) THEN
c           CALL M3EXIT( 'MVEHSUM', 0, 0, 
c    &                   'Error reading VMT from file' // ENAME, 2)
c       END IF

C.........  Sort output FIPs code list
        CALL SORTI1( NOFIP, INDEXA, OFIPA )

        DO 300 I = 1, NOFIP
            OFIP( I ) = OFIPA( INDEXA( I ) )
300     CONTINUE

        ICNT = 0
        DO 320 S = 1, NMSRC
            FIP = IFIP( S )
            K = FIND1( FIP, NOFIP, OFIP )
            IF( K .GT. 0 ) THEN
                EIDX( S ) = K
            ELSE
                EIDX( S ) = 0
            ENDIF
320     CONTINUE

C.................  Initialize summation variables

        DO 440 I = 1, MXOFIP
            DO 430 J = 1, MXEMIS
                DO 420 V = 1, NVTYPE
                    DO 410 D = 1, MXDAYS

                        EMSUM( I, J, V, D ) = 0.

410                 CONTINUE
420             CONTINUE
430         CONTINUE
440     CONTINUE

C.........  Read MTMP summing all the data by output FIPs codes for each 
C.........  emissions type and vehicle class for each day.

        JDATE = SDATE
        JTIME = STIME
        LDATE = -1
        D = 0
        DO 510 K = 1, NSTEPS

            IF ( .NOT. READ3 ( HNAME, ALLVAR3, ALLAYS3, JDATE, 
     &                         JTIME, TMAT ) ) THEN
                CALL M3EXIT ( 'MVEHSUM', 0, 0, 
     &                        'Error reading ALLVARS from file ' // 
     &                        HNAME, 2 )
            END IF

C.............  Perform once-per-day steps
            IF( JDATE .NE. LDATE ) THEN 

                D = D + 1

                MESG = 'Processing ' //
     &                 DAYS( WKDAY( JDATE ) ) // MMDDYY( JDATE )
                CALL M3MSG2( MESG( 1:TRIMLEN( MESG ) ) )

            ENDIF

C.............  Don't want spaces here, so no M3MSG2
            WRITE( LDEV, 810 ) '        at time ' // HHMMSS( JTIME )

            PFIP = 0
            DO 460 S = 1, NMSRC

                I   = EIDX( S )
                FIP = IFIP( S )

                IF( I .GT. 0 ) THEN
                    DO 450 J = 1, MXEMIS
                        DO 445 V = 1, NVTYPE

                            EMSUM( I, J, V, D ) = EMSUM( I, J, V, D ) + 
     &                                             TMAT( S, V, J )

445                     CONTINUE
450                 CONTINUE

                    IF( OFIP( I ) .NE. PFIP ) THEN
                        PFIP      = OFIP( I )
                        BSRC( I ) = S
                        ESRC( I ) = S

                    ELSE
                        ESRC( I ) = ESRC( I ) +  1
                    ENDIF

                ENDIF

460         CONTINUE

            LDATE = JDATE
            CALL NEXTIME ( JDATE, JTIME, TSTEP ) ! Next time step.

510     CONTINUE

        NDAYS = D

C.........  Write report of emissions totals 

        WRITE( ODEV, 93020 ) ( EMTNAM( J ), J = 1, NETYPE )
        WRITE( FMTBUF, '( "(26X,",I2.2,"(1X,A9))" )' ) NETYPE
        WRITE( ODEV, FMTBUF ) ( '[ton/day]', J = 1, NETYPE )

        DO 550 I = 1, NOFIP

            NS = ESRC( I ) - BSRC( I ) + 1

            WRITE( ODEV, 93010 ) I, OFIP( I ), BSRC( I ), ESRC( I ), NS

            DO 540 D = 1, NDAYS

                DO 530 V = 1, NVTYPE
                    WRITE( ODEV, 93030 ) D, OFIP( I ), VTYPE3( V ),  
     &                   ( EMSUM( I, J, V, D ), J = 1,NETYPE )

530             CONTINUE
540         CONTINUE

550     CONTINUE


999     CONTINUE          !  exit program

        MESG = 'Normal completion of MVEHSUM'
        CALL M3EXIT( 'MVEHSUM', 0, 0, MESG, 0 )


C******************  FORMAT  STATEMENTS   ******************************

600     FORMAT ( I10 )
700     FORMAT ( 2X, I5, 1X, A )
800     FORMAT ( A, / )
810     FORMAT ( A )
820     FORMAT ( A10, $ )
830     FORMAT ( I2.2, 5X, 3F14.2, 4X, A )
840     FORMAT ( A, $ )
850     FORMAT ( //, A, / )
860     FORMAT ( /, 19X, A )
870     FORMAT ( I5.5, 2X, 3F14.2, 4X, A )
880     FORMAT ( F14.2, 2X, F14.2, 2X, F14.2 )
890     FORMAT ( A, F12.2, A )

93010   FORMAT ( I5, 7X, I6.5, ' FIPS', /, 6X, I6, ' Src1',/,  
     &           6X, I6, ' Src2', /, 6X, I6, ' Nsrc' )

93020   FORMAT ( 26X, 20( 1X, A9 ) )


93030   FORMAT ( 6X, 'Day=', I2, 1X, I5.5, 1X, A5, 2X, 20( 1X, E9.3 ) )

94010   FORMAT( 10 (A, :, 1X, I6, :, 1X ) )

        END
