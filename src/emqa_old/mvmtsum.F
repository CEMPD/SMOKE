
        PROGRAM MVMTSUM

C***********************************************************************
C  program body starts at line 123
C
C  DESCRIPTION:
C       Report mobile VMT summed by state and county.
C
C  PRECONDITIONS REQUIRED:  
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       Models-3 I/O; 
C       PROMPTMFILE, PROMPTFFILE, TRIMLEN, LBLANK, 
C       GETREAL, ENVREAL, FIND1
C
C  REVISION  HISTORY:
C       Prototype  1/97 by Chad Barker
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

        INCLUDE 'MBDIMS3.EXT'   !  mobile-source dimensioning parameters
        INCLUDE 'CHDIMS3.EXT'   !  emis chem parms (inventory + model)
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.


C...........   EXTERNAL FUNCTIONS and their descriptions:

        CHARACTER*16    PROMPTMFILE
        REAL            GETREAL, ENVREAL
        INTEGER         LBLANK, TRIMLEN, PROMPTFFILE, FIND1

        EXTERNAL        PROMPTMFILE, TRIMLEN, GETREAL
        EXTERNAL        LBLANK, PROMPTFFILE, ENVREAL, FIND1
     
                                            
C...........   PARAMETERS and their descriptions
        INTEGER MAXST, MAXCO ! Max number of states and counties
        PARAMETER ( MAXST = 60, MAXCO = 3500 )

        INTEGER MYMXSTA, MYMXFIP ! Local version of vars in MBDIMS3.EXT
        PARAMETER ( MYMXSTA = 60, MYMXFIP = 3500 )

C...........   LOCAL VARIABLES and their descriptions:
C...........   Point Sources input and output arrays

        INTEGER      ISTFIP  ( MAXST )  !  All state FIPS IDs
        CHARACTER*20 STATES  ( MAXST )  !  All state names
        INTEGER      ICOFIP  ( MAXCO )  !  All county FIPS IDs
        CHARACTER*20 COUNTIES( MAXCO )  !  All county names
        INTEGER      TOT_ST, TOT_CO     ! Total number of states and counties 
                                        ! in the FIPS file.

        ! Raw variables read from the MOBL file
        INTEGER      IFIP    ( NMSRC )
        INTEGER      ILINK   ( NMSRC )
        REAL         VMT     ( NMSRC )

        INTEGER      ISTLIST ( MYMXSTA )  ! List of states in MOBL file
        INTEGER      ICOLIST ( MYMXFIP )  ! List of counties in MOBL file
        INTEGER      ACT_ST, ACT_CO     ! Actual number of states and counties
                                        ! in the MOBL file

        ! Variables/arrays to hold various sums.
        REAL         INVLINK / 0 /, INVNLINK / 0 /
        REAL         STLINK  ( MYMXSTA )  / MYMXSTA * 0.0 /
        REAL         STNLINK ( MYMXSTA )  / MYMXSTA * 0.0 /
        REAL         COLINK  ( MYMXFIP )  / MYMXFIP * 0.0 /
        REAL         CONLINK ( MYMXFIP )  / MYMXFIP * 0.0 /

C...........   Other local variables

        INTEGER         NSRCS

        INTEGER         LDEV    !  unit number for log file
        INTEGER         OUNIT, FUNIT

        INTEGER         J, INDEX, STATUS
        REAL            MINVMT, TOTVMT
        INTEGER         LASTCO, LASTST  ! The last state and county FIPS numbers
                                        ! encountered when summing the VMT.

        CHARACTER*16    ENAME   !  logical name for mobile-source input file
        CHARACTER*10    BUF1,BUF2
        INTEGER         FIPNUM  !  Buffer for reading the FIPS code
        CHARACTER*20    FIPNAME !  Buffer for reading the state/county name

        CHARACTER*256   MESG    !  buffer for M3EXIT() messages


C***********************************************************************
C   begin body of program MVMTSUM

        LDEV = INIT3()

        CALL INITEM( LDEV )

C.......   Get file name; open input area source

        ENAME = PROMPTMFILE( 
     &          'Enter name for MOBILE SOURCE input file',
     &          FSREAD3, 'MOBL', 'MVMTSUM' )

        IF ( .NOT. DESC3( ENAME ) ) THEN
            CALL M3EXIT( 'MVMTSUM', 0, 0, 
     &                  'Could not get description of file "' 
     &                  // ENAME( 1:TRIMLEN( ENAME ) ) // '"', 2 )
        END IF

        NSRCS = NROWS3D
        IF ( NSRCS .GT. NMSRC ) THEN
            WRITE ( BUF1, 600 ) NSRCS
            WRITE ( BUF2, 600 ) NMSRC
            CALL M3EXIT( 'MVMTSUM', 0, 0, 
     &                   'Number of sources (NROWS3D='  // 
     &                   BUF1 ( LBLANK(BUF1):TRIMLEN(BUF1) ) //
     &                   ') exceeds program limits (NMSRC=' //
     &                   BUF2 ( LBLANK(BUF2):TRIMLEN(BUF2) ) // 
     &                   ')', 2 )
        END IF

        FUNIT = PROMPTFFILE( 
     &          'Enter name for STATE/COUNTY FIPS CODE file',
     &          .TRUE., .TRUE., 'STCYNAM', 'MVMTSUM' )

        MINVMT = ENVREAL ( 'MINVMT', 'Minimum VMT', 0, STATUS )
        MINVMT = GETREAL ( 0, 999999999.99, MINVMT, 
     &                     'Enter minimum VMT for counties in report' )

        OUNIT = PROMPTFFILE( 
     &          'Enter name for VMT SUM REPORT file',
     &          .FALSE., .TRUE., 'MVMTRPT', 'MVMTSUM' )


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
                        CALL M3EXIT( 'MVMTSUM', 0, 0, 
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
                        CALL M3EXIT( 'MVMTSUM', 0, 0, 
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

        IF ( .NOT. READ3 ( ENAME, 'IFIP', ALLAYS3, 0, 0, IFIP ) ) THEN
            CALL M3EXIT( 'MVMTSUM', 0, 0, 
     &                   'Error reading IFIP from file' // ENAME, 2)
        END IF

        IF ( .NOT. READ3 ( ENAME, 'ILINK', ALLAYS3, 0, 0, ILINK ) ) THEN
            CALL M3EXIT( 'MVMTSUM', 0, 0, 
     &                   'Error reading ILINK from file' // ENAME, 2)
        END IF

        IF ( .NOT. READ3 ( ENAME, 'VMT', ALLAYS3, 0, 0, VMT ) ) THEN
            CALL M3EXIT( 'MVMTSUM', 0, 0, 
     &                   'Error reading VMT from file' // ENAME, 2)
        END IF


        LASTCO = 0
        LASTST = 0
        ACT_ST = 0
        ACT_CO = 0

        DO 60 J = 1, NSRCS
            IF ( IFIP ( J ) .NE. LASTCO ) THEN
                ! We've moved to a new FIP code (i.e., we've switched counties),
                ! so bump the number of counties and record the FIP code.
                LASTCO = IFIP ( J )
                ACT_CO = ACT_CO + 1
                IF ( ACT_CO .GT. MYMXFIP ) THEN
                    WRITE ( BUF1, 600 ) ACT_CO
                    WRITE ( BUF2, 600 ) MYMXFIP
                    CALL M3EXIT( 'MVMTSUM', 0, 0, 
     &                  'Number of counties in MOBL (' // 
     &                  BUF1 ( LBLANK(BUF1):TRIMLEN(BUF1) ) //
     &                  ') exceeds program limits (MYMXFIP=' //
     &                  BUF2 ( LBLANK(BUF2):TRIMLEN(BUF2) ) // 
     &                  ')', 2 )
                END IF

                ICOLIST ( ACT_CO ) = LASTCO
                
                ! If the current state code is not equal to the last state
                ! code, then bump the number of states and record the code.
                IF ( ( IFIP ( J ) / 1000 ) .NE. LASTST ) THEN
                    LASTST = IFIP ( J ) / 1000
                    ACT_ST = ACT_ST + 1
                    IF ( ACT_ST .GT. MYMXSTA ) THEN
                        WRITE ( BUF1, 600 ) ACT_ST
                        WRITE ( BUF2, 600 ) MYMXSTA
                        CALL M3EXIT( 'MVMTSUM', 0, 0, 
     &                      'Number of states in AREA (' // 
     &                      BUF1 ( LBLANK(BUF1):TRIMLEN(BUF1) ) //
     &                      ') exceeds program limits (MYMXSTA=' //
     &                      BUF2 ( LBLANK(BUF2):TRIMLEN(BUF2) ) // 
     &                      ')', 2 )
                    END IF

                    ISTLIST ( ACT_ST ) = LASTST
                END IF
            END IF

            IF ( ILINK ( J ) .EQ. 0 ) THEN
                STNLINK ( ACT_ST ) = STNLINK ( ACT_ST ) + VMT ( J )
                CONLINK ( ACT_CO ) = CONLINK ( ACT_CO ) + VMT ( J )
                INVNLINK = INVNLINK + VMT ( J )
            ELSE
                STLINK ( ACT_ST ) = STLINK ( ACT_ST ) + VMT ( J )
                COLINK ( ACT_CO ) = COLINK ( ACT_CO ) + VMT ( J )
                INVLINK = INVLINK + VMT ( J )
            END IF
60      CONTINUE


        WRITE ( OUNIT, 800 ) 
     &      'Mobile source vehicle-miles-traveled (miles/day)'
        WRITE ( OUNIT, 800 ) 
     &      '---------------------- Inventory -----------------------'
        WRITE ( OUNIT, 810 ) 
     &      '      NON-LINK            LINK           TOTAL'
        WRITE ( OUNIT, 880 ) INVNLINK, INVLINK, INVNLINK + INVLINK

        WRITE ( OUNIT, 850 ) 
     &      '----------------------- By State -----------------------'
        WRITE ( OUNIT, 810 ) 'STATE'
        WRITE ( OUNIT, 810 ) 
     &      'CODE         NON-LINK          LINK         TOTAL' // 
     &      '    STATE NAME'

        DO 90 J = 1, ACT_ST
            INDEX = FIND1 ( ISTLIST ( J ), TOT_ST, ISTFIP )
            IF ( INDEX .GT. 0 ) THEN
            	FIPNAME = STATES ( INDEX )
            ELSE
                FIPNAME = '*** no state name'
            ENDIF

            WRITE ( OUNIT, 830 ) 
     &           ISTLIST ( J ), STNLINK ( J ), STLINK ( J ),
     &           ( STNLINK ( J ) + STLINK ( J ) ), FIPNAME
90      CONTINUE
        

        WRITE ( OUNIT, 850 ) 
     &      '----------------------- By County -----------------------'
        WRITE ( OUNIT, 890 ) 
     &      'Totals for counties with VMT > ', MINVMT, ' miles/day'

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
                WRITE ( OUNIT, 810 ) 
     &           'CODE         NON-LINK          LINK         TOTAL' //
     &           '    COUNTY NAME'
            END IF

            TOTVMT = CONLINK ( J ) + COLINK ( J )
            IF ( TOTVMT .GE. MINVMT ) THEN
                INDEX = FIND1 ( ICOLIST ( J ), TOT_CO, ICOFIP )
                IF ( INDEX .GT. 0 ) THEN
                	FIPNAME = COUNTIES ( INDEX )
                ELSE
                    FIPNAME = '*** no county name'
                ENDIF

                WRITE ( OUNIT, 870 ) 
     &               ICOLIST ( J ), CONLINK ( J ), COLINK ( J ),
     &               TOTVMT, FIPNAME
            END IF
100     CONTINUE
        
999     CONTINUE          !  exit program

        MESG = 'Normal completion of MVMTSUM'
        CALL M3EXIT( 'MVMTSUM', 0, 0, MESG, 0 )


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
        END
