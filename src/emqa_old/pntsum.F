
        PROGRAM PNTSUM

C***********************************************************************
C  program body starts at line 113
C
C  DESCRIPTION:
C       Report point source emissions summed by state and county.
C
C  PRECONDITIONS REQUIRED:  
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       Models-3 I/O; 
C       GETREAL, ENVREAL, PROMPTMFILE, PROMPTFFILE, TRIMLEN, LBLANK, 
C       ENVINT, FIND1
C
C  REVISION  HISTORY:
C       Prototype  2/97 by Chad Barker
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
C***********************************************************************

      IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'PTDIMS3.EXT'   !  point-source dimensioning parameters
        INCLUDE 'CHDIMS3.EXT'   !  emis chem parms (inventory + model)
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.


C...........   EXTERNAL FUNCTIONS and their descriptions:

        REAL            GETREAL, ENVREAL
        CHARACTER*16    PROMPTMFILE
        INTEGER         LBLANK, TRIMLEN, PROMPTFFILE, ENVINT, FIND1

        EXTERNAL        GETREAL, ENVREAL, PROMPTMFILE, TRIMLEN
        EXTERNAL        LBLANK, PROMPTFFILE, ENVINT, FIND1
     
                                            
C...........   PARAMETERS and their descriptions
        INTEGER MAXST, MAXCO ! Max number of states and counties
        PARAMETER ( MAXST = 60, MAXCO = 3500 )
        
        INTEGER NSTPOL, NCOPOL
        PARAMETER ( NSTPOL = NPSID * NIPOL, NCOPOL = NPFIP * NIPOL )

C...........   LOCAL VARIABLES and their descriptions:
C...........   Point Sources input and output arrays

        INTEGER      ISTFIP  ( MAXST )  !  All state FIPS IDs
        CHARACTER*20 STATES  ( MAXST )  !  All state names
        INTEGER      ICOFIP  ( MAXCO )  !  All county FIPS IDs
        CHARACTER*20 COUNTIES( MAXCO )  !  All county names
        INTEGER      ISTLIST ( NPSID )  !  List of states in PTMP file
        INTEGER      ICOLIST ( NPFIP )  !  List of counties in PTMP file
        REAL         STSUM   ( NPSID, NIPOL ) / NSTPOL * 0.0 /
        REAL         COSUM   ( NPFIP, NIPOL ) / NCOPOL * 0.0 /
        INTEGER      IFIP    ( NPSRC )
        INTEGER      INVYR   ( NPSRC )
        REAL         SPECBUF ( NPSRC, NIPOL )

C...........   Other local variables

        INTEGER         NSRCS

        INTEGER         LDEV    !  unit number for log file
        INTEGER         OUNIT, FUNIT

        REAL            YR2DAY
        INTEGER         I, J, K, INDEX
        INTEGER         LASTCO, LASTST  ! The last state and county FIPS numbers
                                        ! encountered when summing the emmis.
        INTEGER         TOT_ST, TOT_CO  ! Total number of states and counties 
                                        ! in the FIPS file.
        INTEGER         ACT_ST, ACT_CO  ! Actual number of states and counties
                                        ! in the PNTS source file

        CHARACTER*16    ENAME   !  logical name for point-source input file
        CHARACTER*10    BUF1,BUF2
        INTEGER         FIPNUM  !  Buffer for reading the FIPS code
        CHARACTER*20    FIPNAME !  Buffer for reading the state/county name

        CHARACTER*256   MESG    !  buffer for M3EXIT() messages


C***********************************************************************
C   begin body of program PNTSUM

        LDEV = INIT3()

        CALL INITEM( LDEV )

C.......   Get file name; open input PNTS source

        ENAME = PROMPTMFILE( 
     &          'Enter name for POINT SOURCE input file',
     &          FSREAD3, 'PNTS', 'REPTPRO' )

        IF ( .NOT. DESC3( ENAME ) ) THEN
            CALL M3EXIT( 'PNTSUM', 0, 0, 
     &                  'Could not get description of file "' 
     &                  // ENAME( 1:TRIMLEN( ENAME ) ) // '"', 2 )
        END IF

        NSRCS = NROWS3D
        IF ( NSRCS .GT. NPSRC ) THEN
            WRITE ( BUF1, 600 ) NSRCS
            WRITE ( BUF2, 600 ) NPSRC
            CALL M3EXIT( 'PNTSUM', 0, 0, 
     &                   'Number of sources (NROWS3D='  // 
     &                   BUF1 ( LBLANK(BUF1):TRIMLEN(BUF1) ) //
     &                   ') exceeds program limits (NPSRC=' //
     &                   BUF2 ( LBLANK(BUF2):TRIMLEN(BUF2) ) // 
     &                   ')', 2 )
        END IF

        FUNIT = PROMPTFFILE( 
     &          'Enter name for STATE/COUNTY FIPS CODE file',
     &          .TRUE., .TRUE., 'STCYNAM', 'PNTSUM' )

        OUNIT = PROMPTFFILE( 
     &          'Enter name for EMISSIONS SUM REPORT file',
     &          .FALSE., .TRUE., 'PSUMRPT', 'PNTSUM' )


        TOT_ST = 0
        TOT_CO = 0
30      CONTINUE ! Begin "infinite" loop; until EOF.
            READ ( FUNIT, 700, END = 40 ) FIPNUM, FIPNAME
            IF ( FIPNUM .NE. 0 ) THEN  ! If FIPNUM=0, we just read a blank line
                IF ( MOD ( FIPNUM, 1000 ) .EQ. 0 ) THEN
                    TOT_ST = TOT_ST + 1
                    IF ( TOT_ST .GT. MAXST ) THEN
                        WRITE ( BUF1, 600 ) TOT_ST
                        WRITE ( BUF2, 600 ) MAXST
                        CALL M3EXIT( 'PNTSUM', 0, 0, 
     &                      'Number of states in STCYNAM (' // 
     &                      BUF1 ( LBLANK(BUF1):TRIMLEN(BUF1) ) //
     &                      ') exceeds program limits (MAXST=' //
     &                      BUF2 ( LBLANK(BUF2):TRIMLEN(BUF2) ) // 
     &                      ')', 2 )
                    END IF

                    ISTFIP ( TOT_ST ) = FIPNUM / 1000
                    STATES ( TOT_ST ) = FIPNAME
C                    PRINT '( A, I2, A, A, A, I2, A )', 
C     &                'Got state #', TOT_ST, ': ',
C     &                STATES ( TOT_ST ), ' (FIP=', ISTFIP ( TOT_ST ), ')'
                ELSE
                    TOT_CO = TOT_CO + 1
                    IF ( TOT_CO .GT. MAXCO ) THEN
                        WRITE ( BUF1, 600 ) TOT_CO
                        WRITE ( BUF2, 600 ) MAXCO
                        CALL M3EXIT( 'PNTSUM', 0, 0, 
     &                      'Number of counties in STCYNAM (' // 
     &                      BUF1 ( LBLANK(BUF1):TRIMLEN(BUF1) ) //
     &                      ') exceeds program limits (MAXCO=' //
     &                      BUF2 ( LBLANK(BUF2):TRIMLEN(BUF2) ) // 
     &                      ')', 2 )
                    END IF

                    ICOFIP ( TOT_CO ) = FIPNUM
                    COUNTIES ( TOT_CO ) = FIPNAME
C                    PRINT '( A, I4, A, A, A, I5, A )', 
C     &                'Got county #', TOT_CO, ': ',
C     &                COUNTIES ( TOT_CO ), ' (FIP=', ICOFIP ( TOT_CO ),')'
                END IF
            END IF
        GO TO 30
40      CONTINUE

        IF ( .NOT. READ3( ENAME, 'IFIP', ALLAYS3,0,0, IFIP ) ) THEN
            CALL M3EXIT( 'REPTPRO', 0, 0, 
     &                   'Error reading FIP from file' // ENAME, 2)
        END IF

        IF ( .NOT. READ3 ( ENAME, 'INVYR', ALLAYS3,0,0, INVYR ) ) THEN
            CALL M3EXIT( 'REPTPRO', 0, 0, 
     &                   'Error reading INVYR from file' // ENAME, 2)
        END IF

        DO 50 I = 1, NIPOL
C           PRINT *, 'Reading ' // EINAM (I) // '...'
            IF ( .NOT. READ3 ( ENAME, EINAM ( I ), ALLAYS3, 0, 0, 
     &                         SPECBUF ( 1, I ) ) ) THEN
                CALL M3EXIT( 'PNTSUM', 0, 0, 
     &                       'Error reading ' // EINAM ( I ) //
     &                       ' from file' // ENAME, 2 )
            END IF
50      CONTINUE

        LASTCO = 0
        LASTST = 0
        ACT_ST = 0
        ACT_CO = 0

        DO 60 J = 1, NSRCS
            K = INVYR( J )
            IF ( MOD( K,4 ) .NE. 0 ) THEN           !  nonleap years
                YR2DAY = 1.0 / 365.0
            ELSE IF ( MOD( K,100 ) .NE. 0 ) THEN    !  noncentury leap years
                YR2DAY = 1.0 / 366.0
            ELSE IF ( MOD( K,400 ) .NE. 0 ) THEN    !  century nonleap years
                YR2DAY = 1.0 / 365.0
            ELSE                                    !  leap centuries
                YR2DAY = 1.0 / 366.0
            END IF

            IF ( IFIP ( J ) .NE. LASTCO ) THEN
                ! We've moved to a new FIP code (i.e., we've switched counties),
                ! so bump the number of counties and record the FIP code.
                LASTCO = IFIP ( J )
                ACT_CO = ACT_CO + 1
                IF ( ACT_CO .GT. NPFIP ) THEN
                    WRITE ( BUF1, 600 ) ACT_CO
                    WRITE ( BUF2, 600 ) NPFIP
                    CALL M3EXIT( 'PNTSUM', 0, 0, 
     &                  'Number of counties in PNTS (' // 
     &                  BUF1 ( LBLANK(BUF1):TRIMLEN(BUF1) ) //
     &                  ') exceeds program limits (NPFIP=' //
     &                  BUF2 ( LBLANK(BUF2):TRIMLEN(BUF2) ) // 
     &                  ')', 2 )
                END IF

                ICOLIST ( ACT_CO ) = LASTCO
                
                ! If the current state code is not equal to the last state
                ! code, then bump the number of states and record the code.
                IF ( ( IFIP ( J ) / 1000 ) .NE. LASTST ) THEN
                    LASTST = IFIP ( J ) / 1000
                    ACT_ST = ACT_ST + 1
                    IF ( ACT_ST .GT. NPSID ) THEN
                        WRITE ( BUF1, 600 ) ACT_ST
                        WRITE ( BUF2, 600 ) NPSID
                        CALL M3EXIT( 'PNTSUM', 0, 0, 
     &                      'Number of states in PNTS (' // 
     &                      BUF1 ( LBLANK(BUF1):TRIMLEN(BUF1) ) //
     &                      ') exceeds program limits (NPSID=' //
     &                      BUF2 ( LBLANK(BUF2):TRIMLEN(BUF2) ) // 
     &                      ')', 2 )
                    END IF

                    ISTLIST ( ACT_ST ) = LASTST
                END IF
            END IF
            
            DO 70 I = 1, NIPOL
                STSUM ( ACT_ST, I ) = STSUM ( ACT_ST, I ) + 
     &                                ( SPECBUF ( J, I ) * YR2DAY )
                COSUM ( ACT_CO, I ) = COSUM ( ACT_CO, I ) + 
     &                                ( SPECBUF ( J, I ) * YR2DAY )
70          CONTINUE
60      CONTINUE


        WRITE ( OUNIT, 800 ) 'Point source emissions (tons/day)'
        WRITE ( OUNIT, 800 ) 
     &        '----------------------- By State -----------------------'
        WRITE ( OUNIT, 810 ) 'STATE'
        WRITE ( OUNIT, 840 ) 'CODE'
        DO 80 I = 1, NIPOL
            WRITE ( OUNIT, 820 ) EINAM ( I )( 1:TRIMLEN ( EINAM ( I ) ))
80      CONTINUE
        WRITE ( OUNIT, 810 ) '  STATE NAME'

        DO 90 J = 1, ACT_ST
            INDEX = FIND1 ( ISTLIST ( J ), TOT_ST, ISTFIP )
            IF ( INDEX .GT. 0 ) THEN
            	FIPNAME = STATES ( INDEX )
            ELSE
                FIPNAME = '*** no state name'
            ENDIF

            WRITE ( OUNIT, 830 ) 
     &           ISTLIST ( J ), ( STSUM ( J, I ), I = 1, NIPOL ),
     &           FIPNAME
90      CONTINUE
        

        WRITE ( OUNIT, 850 ) '----------------------- By County ',
     &                       '-----------------------'

        LASTST = 0
        DO 100 J = 1, ACT_CO
            IF ( ( ICOLIST ( J ) / 1000 ) .NE. LASTST ) THEN
                LASTST = ICOLIST ( J ) / 1000
                INDEX = FIND1 ( LASTST, TOT_ST, ISTFIP )
                IF ( INDEX .GT. 0 ) THEN
                	FIPNAME = STATES ( INDEX )
                ELSE
                    FIPNAME = '*** no state name'
                ENDIF
        
                WRITE ( OUNIT, 860 ) FIPNAME
                WRITE ( OUNIT, 810 ) 'FIPs'
                WRITE ( OUNIT, 840 ) 'CODE '
                DO 110 I = 1, NIPOL
                    WRITE ( OUNIT, 820 ) 
     &                    EINAM ( I )( 1:TRIMLEN ( EINAM ( I ) ) )
110             CONTINUE
                WRITE ( OUNIT, 810 ) '  COUNTY NAME'
            END IF

            INDEX = FIND1 ( ICOLIST ( J ), TOT_CO, ICOFIP )
            IF ( INDEX .GT. 0 ) THEN
            	FIPNAME = COUNTIES ( INDEX )
            ELSE
                FIPNAME = '*** no county name'
            ENDIF

            WRITE ( OUNIT, 870 ) 
     &           ICOLIST ( J ), ( COSUM ( J, I ), I = 1, NIPOL ),
     &           FIPNAME
100     CONTINUE
        
999     CONTINUE          !  exit program

        MESG = 'Normal completion of PNTSUM'
        CALL M3EXIT( 'PNTSUM', 0, 0, MESG, 0 )


C******************  FORMAT  STATEMENTS   ******************************

600     FORMAT ( I10 )
700     FORMAT ( 2X, I5, 1X, A )
800     FORMAT ( A, / )
810     FORMAT ( A )
820     FORMAT ( A10, $ )
830     FORMAT ( I2.2, 2X, <NIPOL>F10.2, 2X, A )
840     FORMAT ( A, $ )
850     FORMAT ( //, A, A )
860     FORMAT ( /, 19X, A )
870     FORMAT ( I5.5, <NIPOL>F10.2, 2X, A )
        END
