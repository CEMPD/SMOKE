        PROGRAM MOBILSUM

C***********************************************************************
C  program body starts at line 152
C
C  DESCRIPTION:
C       Report mobile-source emissions summed by state and county.
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
C
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1997, Chad Barker, Turner Engineering
C All Rights Reserved
C
C See file COPYRIGHT for conditions of use.
C
C Chad Barker, Senior Computer Scientist
C Timothy L. Turner, Owner
C Turner Engineering
C 1310 Holly Creek Lane
C Chapel Hill NC 27516
C (919) 929 6851
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
        
        INTEGER EMISSTA, EMISFIP
        PARAMETER ( EMISSTA = MXEMIS * MYMXSTA, 
     &              EMISFIP = MXEMIS * MYMXFIP )

        REAL GM2TON
        PARAMETER ( GM2TON = 1.102311E-06 )

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
        
        ! Raw variables read from MTMP file
        REAL         TMAT ( NMSRC, NVTYPE, MXEMIS ) ! For reading one time step

        ! Arrays to map raw variables (above) to State and County lists (below)
        INTEGER      STMAP   ( NMSRC )
        INTEGER      COMAP   ( NMSRC )

        INTEGER      ISTLIST ( MYMXSTA )  ! List of states in MOBL file
        INTEGER      ICOLIST ( MYMXFIP )  ! List of counties in MOBL file
        INTEGER      ACT_ST, ACT_CO     ! Actual number of states and counties
                                        ! in the MOBL file

        ! Variables/arrays to hold various sums.
        REAL         TOTCOVMT ( MYMXFIP )  / MYMXFIP * 0.0 /
        REAL         INVLINK  ( MXEMIS ) / MXEMIS * 0 /
        REAL         INVNLINK ( MXEMIS ) / MXEMIS * 0 /
        REAL         STLINK   ( MXEMIS, MYMXSTA ) / EMISSTA * 0.0 /
        REAL         STNLINK  ( MXEMIS, MYMXSTA ) / EMISSTA * 0.0 /
        REAL         COLINK   ( MXEMIS, MYMXFIP ) / EMISFIP * 0.0 /
        REAL         CONLINK  ( MXEMIS, MYMXFIP ) / EMISFIP * 0.0 /


        CHARACTER*16 INVNAMES ( MXEMIS )
        INTEGER      NUMINV / 0 /
        
        ! Map between TMAT vars and INVNAMES
        INTEGER      INVMAP ( MXEMIS ) / MXEMIS * 0 / 

C...........   Other local variables

        INTEGER         NSRCS

        INTEGER         LDEV    !  unit number for log file
        INTEGER         OUNIT ( MXEMIS )
        INTEGER         FUNIT

        INTEGER         JDATE, JTIME, TSTEP, NSTEPS
        INTEGER         NUMVARS, NUMLAYS
        INTEGER         CURDATE, CURTIME
        INTEGER         I, J, K, L, LOC, STATUS, POSN
        REAL            MINVMT
        INTEGER         LASTCO, LASTST  ! The last state and county FIPS numbers
                                        ! encountered when summing the VMT.

        CHARACTER*16    ENAME   !  logical name for mobile-source input file
        CHARACTER*16    HNAME   !   " for mobile-source emissions input file
        CHARACTER*10    BUF1,BUF2
        INTEGER         FIPNUM  !  Buffer for reading the FIPS code
        CHARACTER*20    FIPNAME !  Buffer for reading the state/county name
        CHARACTER*16    VARBUF

        CHARACTER*256   MESG    !  buffer for M3EXIT() messages


C***********************************************************************
C   begin body of program MOBILSUM

        LDEV = INIT3()

C.......   Get file name; open input area source

        ENAME = PROMPTMFILE( 
     &          'Enter name for MOBILE SOURCE input file',
     &          FSREAD3, 'MOBL', 'MOBILSUM' )

        IF ( .NOT. DESC3( ENAME ) ) THEN
            CALL M3EXIT( 'MOBILSUM', 0, 0, 
     &                  'Could not get description of file "' 
     &                  // ENAME( 1:TRIMLEN( ENAME ) ) // '"', 2 )
        END IF

        NSRCS = NROWS3D
        IF ( NSRCS .GT. NMSRC ) THEN
            WRITE ( BUF1, 600 ) NSRCS
            WRITE ( BUF2, 600 ) NMSRC
            CALL M3EXIT( 'MOBILSUM', 0, 0, 
     &                   'Number of sources (NROWS3D='  // 
     &                   BUF1 ( LBLANK(BUF1):TRIMLEN(BUF1) ) //
     &                   ') exceeds program limits (NMSRC=' //
     &                   BUF2 ( LBLANK(BUF2):TRIMLEN(BUF2) ) // 
     &                   ')', 2 )
        END IF

        HNAME = PROMPTMFILE( 
     &          'Enter name for MOBILE SOURCE EMISSIONS input file',
     &          FSREAD3, 'MTMP', 'MOBILSUM' )

        IF ( .NOT. DESC3( HNAME ) ) THEN
            CALL M3EXIT( 'MOBILSUM', 0, 0, 
     &                  'Could not get description of file "' 
     &                  // HNAME( 1:TRIMLEN( HNAME ) ) // '"', 2 )
        END IF

        JDATE = SDATE3D
        JTIME = STIME3D
        TSTEP = TSTEP3D
        NSTEPS = MXREC3D
        NUMVARS = NVARS3D
        NUMLAYS = NLAYS3D
        IF ( NSRCS .NE. NROWS3D ) THEN
            WRITE ( BUF1, 600 ) NSRCS
            WRITE ( BUF2, 600 ) NROWS3D
            CALL M3EXIT( 'MOBILSUM', 0, 0, 
     &                   'Number of sources in MOBL (NROWS3D='  // 
     &                   BUF1 ( LBLANK(BUF1):TRIMLEN(BUF1) ) //
     &                   ') does not match number in MTMP (NROWS3D=' //
     &                   BUF2 ( LBLANK(BUF2):TRIMLEN(BUF2) ) // 
     &                   ')', 2 )
        END IF

        ! Find the unique inventory pollutant names
        DO 500 I = 1, NUMVARS
            ! Look for the underscore in the name
            POSN = INDEX ( VNAME3D ( I ), '_' )
            IF ( POSN .LE. 0 ) THEN
                INVMAP ( I ) = 0
            ELSE
                ! Take the part of the name beyond the underscore
                VARBUF = VNAME3D ( I ) 
     &                      ( ( POSN + 1 ):LEN ( VNAME3D ( I ) ) )

                ! If this is a new inventory name, then record it.
                IF ( ( NUMINV .EQ. 0 ) .OR. 
     &               ( VARBUF .NE. INVNAMES ( NUMINV ) ) ) THEN

                    NUMINV = NUMINV + 1
                    INVNAMES ( NUMINV ) = VARBUF
                END IF

                ! This map will be used to translate variables in the TMAT
                ! array to the appropriate columns in ST(CO)LINK(NLINK) arrays.
                INVMAP ( I ) = NUMINV
            END IF

500     CONTINUE


        FUNIT = PROMPTFFILE( 
     &          'Enter name for STATE/COUNTY FIPS CODE file',
     &          .TRUE., .TRUE., 'STCYNAM', 'MOBILSUM' )

        MINVMT = ENVREAL ( 'MINVMT', 'Minimum VMT', 0, STATUS )
        MINVMT = GETREAL ( 0, 999999999.99, MINVMT, 
     &                     'Enter minimum VMT for counties in report' )


        ! Prompt for each of the output report file names
        DO 490 I = 1, NUMINV
            OUNIT ( I ) = PROMPTFFILE( 
     &                  'Enter name for MOBILE ' // 
     &                  INVNAMES( I )( 1:TRIMLEN ( INVNAMES ( I ) ) ) //
     &                  ' SUM REPORT file', .FALSE., .TRUE., 'M' // 
     &                  INVNAMES( I )( 1:TRIMLEN ( INVNAMES ( I ) ) ) //
     &                  'RPT', 'MOBILSUM' )
490     CONTINUE

        ! Read complete list of state & County fips codes and names
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
                        CALL M3EXIT( 'MOBILSUM', 0, 0, 
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
                        CALL M3EXIT( 'MOBILSUM', 0, 0, 
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
        GO TO 30
40      CONTINUE


        ! Read variables from MOBL

        IF ( .NOT. READ3 ( ENAME, 'IFIP', ALLAYS3, 0, 0, IFIP ) ) THEN
            CALL M3EXIT( 'MOBILSUM', 0, 0, 
     &                   'Error reading IFIP from file' // ENAME, 2)
        END IF

        IF ( .NOT. READ3 ( ENAME, 'ILINK', ALLAYS3, 0, 0, ILINK ) ) THEN
            CALL M3EXIT( 'MOBILSUM', 0, 0, 
     &                   'Error reading ILINK from file' // ENAME, 2)
        END IF

        IF ( .NOT. READ3 ( ENAME, 'VMT', ALLAYS3, 0, 0, VMT ) ) THEN
            CALL M3EXIT( 'MOBILSUM', 0, 0, 
     &                   'Error reading VMT from file' // ENAME, 2)
        END IF



        ! Figure out the unique list of states and counties in the MOBL file.
        ! Also compute the number of VMT for each county.
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
                    CALL M3EXIT( 'MOBILSUM', 0, 0, 
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
                        CALL M3EXIT( 'MOBILSUM', 0, 0, 
     &                      'Number of states in AREA (' // 
     &                      BUF1 ( LBLANK(BUF1):TRIMLEN(BUF1) ) //
     &                      ') exceeds program limits (MYMXSTA=' //
     &                      BUF2 ( LBLANK(BUF2):TRIMLEN(BUF2) ) // 
     &                      ')', 2 )
                    END IF

                    ISTLIST ( ACT_ST ) = LASTST
                END IF
            END IF

            TOTCOVMT ( ACT_CO ) = TOTCOVMT ( ACT_CO ) + VMT ( J )

            STMAP ( J ) = ACT_ST
            COMAP ( J ) = ACT_CO
60      CONTINUE


        ! Read MTMP summing all the data by state and county, making sure
        ! we sum all "layers" (vehicle type) and time steps.
        CURDATE = JDATE
        CURTIME = JTIME
        DO 510 K = 1, NSTEPS
            IF ( .NOT. READ3 ( HNAME, ALLVAR3, ALLAYS3, CURDATE, 
     &                         CURTIME, TMAT ) ) THEN
                CALL M3EXIT ( 'MOBILSUM', 0, 0, 
     &                        'Error reading ALLVARS from file ' // 
     &                        HNAME, 2 )
            END IF
 
            WRITE( MESG,94010 ) 'Processing time step', K
            CALL M3MESG( MESG )
 
            DO 520 I = 1, NUMVARS

                DO 530 L = 1, NUMLAYS
                    DO 540 J = 1, NSRCS
                        IF ( ILINK ( J ) .EQ. 0 ) THEN
                            INVNLINK ( INVMAP ( I ) ) = 
     &                         INVNLINK ( INVMAP ( I ) ) + 
     &                         TMAT ( J, L, I )
                            STNLINK ( INVMAP ( I ), STMAP ( J ) ) = 
     &                         STNLINK ( INVMAP ( I ), STMAP ( J ) ) + 
     &                         TMAT ( J, L, I )
                            CONLINK ( INVMAP ( I ), COMAP ( J ) ) = 
     &                         CONLINK ( INVMAP ( I ), COMAP ( J ) ) + 
     &                         TMAT ( J, L, I )
                        ELSE
                            INVLINK ( INVMAP ( I ) ) = 
     &                         INVLINK ( INVMAP ( I ) ) + 
     &                         TMAT ( J, L, I )
                            STLINK ( INVMAP ( I ), STMAP ( J ) ) = 
     &                         STLINK ( INVMAP ( I ), STMAP ( J ) ) + 
     &                         TMAT ( J, L, I )
                            COLINK ( INVMAP ( I ), COMAP ( J ) ) = 
     &                         COLINK ( INVMAP ( I ), COMAP ( J ) ) + 
     &                         TMAT ( J, L, I )
                        END IF
540                 CONTINUE
530             CONTINUE
520         CONTINUE

            ! Increment to the next time step.
            CALL NEXTIME ( CURDATE, CURTIME, TSTEP )
510     CONTINUE


        ! Now generate the reports, one for each inventory pollutant.
        DO 110 I = 1, NUMINV
            WRITE ( OUNIT ( I ), 800 ) 
     &          'Mobile source ' // 
     &          INVNAMES ( I )( 1:TRIMLEN ( INVNAMES ( I ) ) ) // 
     &          ' emissions totals (tons/day)'
            WRITE ( OUNIT ( I ), 800 ) 
     &          '-------------------- ' // 
     &          INVNAMES ( I )( 1:TRIMLEN ( INVNAMES ( I ) ) ) // 
     &          ' Inventory ---------------------'
            WRITE ( OUNIT ( I ), 810 ) 
     &          '      NON-LINK            LINK           TOTAL'
            WRITE ( OUNIT ( I ), 880 ) 
     &          INVNLINK ( I ) * GM2TON, INVLINK ( I ) * GM2TON, 
     &          ( INVNLINK ( I ) + INVLINK ( I ) ) * GM2TON

            WRITE ( OUNIT ( I ), 850 ) 
     &          '--------------------- ' // 
     &          INVNAMES ( I )( 1:TRIMLEN ( INVNAMES ( I ) ) ) // 
     &          ' By State ---------------------'
            WRITE ( OUNIT ( I ), 810 ) 'STATE'
            WRITE ( OUNIT ( I ), 810 ) 
     &          'CODE         NON-LINK          LINK         TOTAL' // 
     &          '    STATE NAME'

            DO 90 J = 1, ACT_ST
                LOC = FIND1 ( ISTLIST ( J ), TOT_ST, ISTFIP )
                IF ( LOC .GT. 0 ) THEN
                	FIPNAME = STATES ( LOC )
                ELSE
                    FIPNAME = '*** no state name'
                ENDIF

                WRITE ( OUNIT ( I ), 830 ) 
     &              ISTLIST ( J ), STNLINK ( I, J ) * GM2TON, 
     &              STLINK ( I, J ) * GM2TON,
     &              ( STNLINK ( I, J ) + STLINK ( I, J ) ) * GM2TON, 
     &              FIPNAME
90          CONTINUE ! End of state loop
        

            WRITE ( OUNIT ( I ), 850 ) 
     &          '--------------------- ' // 
     &          INVNAMES ( I )( 1:TRIMLEN ( INVNAMES ( I ) ) ) // 
     &          ' By County ---------------------'
            WRITE ( OUNIT ( I ), 890 ) 
     &          'Totals for counties with VMT > ', MINVMT, ' miles/day'

            LASTST = 0
            DO 100 J = 1, ACT_CO
                IF ( ( ICOLIST ( J ) / 1000 ) .NE. LASTST ) THEN
                    LASTST = ICOLIST ( J ) / 1000
                    LOC = FIND1 ( LASTST, TOT_ST, ISTFIP )
                    IF ( LOC .GT. 0 ) THEN
                    	FIPNAME = STATES ( LOC )
                    ELSE
                        FIPNAME = '*** no state name'
                    ENDIF
        
                    WRITE ( OUNIT ( I ), 860 ) FIPNAME
                    WRITE ( OUNIT ( I ), 810 ) 'FIPs'
                    WRITE ( OUNIT ( I ), 810 ) 
     &                  'CODE         NON-LINK          LINK' //
     &                  '         TOTAL    COUNTY NAME'
                END IF

                IF ( TOTCOVMT ( J ) .GE. MINVMT ) THEN
                    LOC = FIND1 ( ICOLIST ( J ), TOT_CO, ICOFIP )
                    IF ( LOC .GT. 0 ) THEN
                    	FIPNAME = COUNTIES ( LOC )
                    ELSE
                        FIPNAME = '*** no county name'
                    ENDIF
    
                    WRITE ( OUNIT ( I ), 870 ) 
     &                  ICOLIST ( J ), CONLINK ( I, J ) * GM2TON, 
     &                  COLINK ( I, J ) * GM2TON, 
     &                  ( CONLINK ( I, J ) + COLINK ( I, J ) ) * GM2TON, 
     &                  FIPNAME
                END IF
100         CONTINUE ! End of county loop

110     CONTINUE ! End of emissions loop
        
999     CONTINUE          !  exit program

        MESG = 'Normal completion of MOBILSUM'
        CALL M3EXIT( 'MOBILSUM', 0, 0, MESG, 0 )


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

94010   FORMAT( 10 (A, :, 1X, I5, :, 1X ) )

        END
