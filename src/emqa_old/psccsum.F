
        PROGRAM PSCCSUM

C***********************************************************************
C  program body starts at line 122
C
C  DESCRIPTION:
C       Report point source emissions summed by SCC code.
C
C  PRECONDITIONS REQUIRED:  
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       Models-3 I/O; 
C       GETREAL, ENVREAL, PROMPTMFILE, PROMPTFFILE, TRIMLEN, LBLANK, 
C       ENVINT, FIND2, SORT2
C
C  REVISION  HISTORY:
C       Prototype  3/97 by Chad Barker
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
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.


C...........   EXTERNAL FUNCTIONS and their descriptions:

        REAL            GETREAL, ENVREAL
        CHARACTER*16    PROMPTMFILE
        INTEGER         LBLANK, TRIMLEN, PROMPTFFILE, ENVINT, FIND1

        EXTERNAL     GETREAL, ENVREAL, PROMPTMFILE, TRIMLEN
        EXTERNAL     LBLANK, PROMPTFFILE, ENVINT, FIND1, SORTI1, SORTI2
     
                                            
C...........   PARAMETERS and their descriptions
        INTEGER MAXSCC, NSCCPOL
        PARAMETER ( MAXSCC = 4000, NSCCPOL = NPSCC * NIPOL )

        INTEGER MAXST, MAXDESC
        PARAMETER ( MAXST = 60, MAXDESC = 40 )

C...........   LOCAL VARIABLES and their descriptions:
C...........   Point Sources input and output arrays

        INTEGER       ISTFIP   ( MAXST )  !  All state FIPS IDs
        CHARACTER*20  STATES   ( MAXST )  !  All state names
        INTEGER       IALLSCC  ( MAXSCC )  !  All SCC codes
        CHARACTER*150 SCCDESC  ( MAXSCC )  !  All SCC names
        LOGICAL       HAS_SUM  ( MAXSCC ) / MAXSCC * .FALSE. /
        INTEGER       ISCCLIST ( NPSCC )
        REAL          SPECSUM  ( NPSCC, NIPOL ) / NSCCPOL * 0.0 /
        INTEGER       SCCINDEX ( NPSRC )
        INTEGER       IFIPST   ( NPSRC )
        INTEGER       ISCC2    ( NPSRC )
        INTEGER       IFIP     ( NPSRC )
        INTEGER       ISCC     ( NPSRC )
        INTEGER       INVYR    ( NPSRC )
        REAL          SPECBUF  ( NPSRC, NIPOL )
        REAL          YR2DAY   ( NPSRC )

C...........   Other local variables

        INTEGER         NSRCS

        INTEGER         LDEV    !  unit number for log file
        INTEGER         OUNIT, SUNIT, FUNIT, RUNIT

        REAL            YR2DAY
        INTEGER         I, J, K, POSN, LOC, L2
        INTEGER         LASTSCC   ! The last SCC number
                                  ! encountered when summing the emmis.
        INTEGER         TOT_SCC   ! Total number of SCC codes/names 
                                  ! in the SCCDESC file.
        INTEGER         ACT_SCC   ! Actual number of SCC codes
                                  ! in the area source file
        INTEGER         TOT_ST, LASTST
        INTEGER         FIPNUM  !  Buffer for reading the FIPS code
        CHARACTER*20    FIPNAME !  Buffer for reading the state/county name

        CHARACTER*16    ENAME   !  logical name for point-source input file
        CHARACTER*10    BUF1,BUF2
        INTEGER         INSCC
        CHARACTER*150   SCCNAME !  Buffer for reading the SCC name
        CHARACTER*45    DESC

        CHARACTER*256   MESG    !  buffer for M3EXIT() messages

C***********************************************************************
C   begin body of program PSCCSUM

        LDEV = INIT3()

        CALL INITEM( LDEV ) 

C.......   Get file name; open input area source

        ENAME = PROMPTMFILE( 
     &          'Enter name for POINT SOURCE input file',
     &          FSREAD3, 'PNTS', 'PSCCSUM' )

        IF ( .NOT. DESC3( ENAME ) ) THEN
            CALL M3EXIT( 'PSCCSUM', 0, 0, 
     &                  'Could not get description of file "' 
     &                  // ENAME( 1:TRIMLEN( ENAME ) ) // '"', 2 )
        END IF

        NSRCS = NROWS3D
        IF ( NSRCS .GT. NPSRC ) THEN
            WRITE ( BUF1, 600 ) NSRCS
            WRITE ( BUF2, 600 ) NPSRC
            CALL M3EXIT( 'PSCCSUM', 0, 0, 
     &                   'Number of sources (NROWS3D='  // 
     &                   BUF1 ( LBLANK(BUF1):TRIMLEN(BUF1) ) //
     &                   ') exceeds program limits (NPSRC=' //
     &                   BUF2 ( LBLANK(BUF2):TRIMLEN(BUF2) ) // 
     &                   ')', 2 )
        END IF

        SUNIT = PROMPTFFILE( 
     &          'Enter name for SCC DESCRIPTION file',
     &          .TRUE., .TRUE., 'SCCDESC', 'PSCCSUM' )

        FUNIT = PROMPTFFILE( 
     &          'Enter name for STATE/COUNTY FIPS CODE file',
     &          .TRUE., .TRUE., 'STCYNAM', 'PSCCSUM' )

        OUNIT = PROMPTFFILE( 
     &          'Enter name for EMISSIONS SUM REPORT file',
     &          .FALSE., .TRUE., 'PSCCRPT', 'PSCCSUM' )

        RUNIT = PROMPTFFILE( 
     &          'Enter name for STATE SUM REPORT file',
     &          .FALSE., .TRUE., 'PSTRPT', 'PSCCSUM' )



        TOT_SCC = 0
30      CONTINUE ! Begin "infinite" loop; until EOF.
            READ ( SUNIT, '(I8, 2X, A)', END = 40 ) INSCC, SCCNAME
C            PRINT '( I8.8, A, A )', INSCC, ': ', SCCNAME
            TOT_SCC = TOT_SCC + 1
            IF ( TOT_SCC .GT. MAXSCC ) THEN
                WRITE ( BUF1, 600 ) TOT_SCC
                WRITE ( BUF2, 600 ) MAXSCC
                CALL M3EXIT( 'PSCCSUM', 0, 0, 
     &                      'Number of SCC codes in SCCDESC (' // 
     &                      BUF1 ( LBLANK(BUF1):TRIMLEN(BUF1) ) //
     &                      ') exceeds program limits (MAXSCC=' //
     &                      BUF2 ( LBLANK(BUF2):TRIMLEN(BUF2) ) // 
     &                      ')', 2 )
            END IF

            IALLSCC ( TOT_SCC ) = INSCC
            SCCDESC ( TOT_SCC ) = SCCNAME
C            PRINT '( A, I4, A, A, A, I8.8, A )', 
C     &            'Got SCC #', TOT_SCC, ': ',
C     &            SCCDESC ( TOT_SCC ), ' (SCC=', 
C     &            IALLSCC ( TOT_SCC ), ')'
        GO TO 30
40      CONTINUE


        TOT_ST = 0
160     CONTINUE ! Begin "infinite" loop; until EOF.
            READ( FUNIT, '( 2X, I5, 1X, A )', END=170 ) FIPNUM, FIPNAME
            IF ( FIPNUM .NE. 0 ) THEN  ! If FIPNUM=0, we just read a blank line
                IF ( MOD ( FIPNUM, 1000 ) .EQ. 0 ) THEN
                    TOT_ST = TOT_ST + 1
                    IF ( TOT_ST .GT. MAXST ) THEN
                        WRITE ( BUF1, 600 ) TOT_ST
                        WRITE ( BUF2, 600 ) MAXST
                        CALL M3EXIT( 'PSCCSUM', 0, 0, 
     &                      'Number of states in STCYNAM (' // 
     &                      BUF1 ( LBLANK(BUF1):TRIMLEN(BUF1) ) //
     &                      ') exceeds program limits (MAXST=' //
     &                      BUF2 ( LBLANK(BUF2):TRIMLEN(BUF2) ) // 
     &                      ')', 2 )
                    END IF

                    ISTFIP ( TOT_ST ) = FIPNUM / 1000
                    STATES ( TOT_ST ) = FIPNAME
C                    PRINT '( A, I2, A, A, A, I2, A )', 
C     &                  'Got state #', TOT_ST, ': ',
C     &                  STATES ( TOT_ST ), ' (FIP=', 
C     &                  ISTFIP ( TOT_ST ), ')'
                END IF
            END IF
        GO TO 160
170      CONTINUE


        IF ( .NOT. READ3 ( ENAME, 'IFIP', ALLAYS3, 0, 0, IFIP ) ) THEN
            CALL M3EXIT( 'PSCCSUM', 0, 0, 
     &                   'Error reading IFIP from file' // ENAME, 2)
        END IF

        IF ( .NOT. READ3 ( ENAME, 'ISCC', ALLAYS3, 0, 0, ISCC ) ) THEN
            CALL M3EXIT( 'PSCCSUM', 0, 0, 
     &                   'Error reading ISCC from file' // ENAME, 2)
        END IF

        IF ( .NOT. READ3 ( ENAME, 'INVYR', ALLAYS3, 0, 0, INVYR ) ) THEN
            CALL M3EXIT( 'PSCCSUM', 0, 0, 
     &                   'Error reading INVYR from file' // ENAME, 2)
        END IF

C        PRINT *, 'Reading emissions...'

        DO 50 I = 1, NIPOL
C           PRINT *, 'Reading ' // EINAM (I) // '...'
            IF ( .NOT. READ3 ( ENAME, EINAM ( I ), ALLAYS3, 0, 0, 
     &                         SPECBUF ( 1, I ) ) ) THEN
                CALL M3EXIT( 'PSCCSUM', 0, 0, 
     &                       'Error reading ' // EINAM ( I ) //
     &                       ' from file' // ENAME, 2 )
            END IF
50      CONTINUE

C        PRINT *, 'Sorting ISCC...'

        DO 45 J = 1, NSRCS
            SCCINDEX ( J ) = J
45      CONTINUE
        CALL SORTI1 ( NSRCS, SCCINDEX, ISCC )


C        PRINT *, 'Summing emissions...'

        DO 150 J = 1, NSRCS
            K = INVYR ( J )
            IF ( MOD( K,4 ) .NE. 0 ) THEN           !  nonleap years
                YR2DAY ( J ) = 1.0 / 365.0
            ELSE IF ( MOD( K,100 ) .NE. 0 ) THEN    !  noncentury leap years
                YR2DAY ( J ) = 1.0 / 366.0
            ELSE IF ( MOD( K,400 ) .NE. 0 ) THEN    !  century nonleap years
                YR2DAY ( J ) = 1.0 / 365.0
            ELSE                                    !  leap centuries
                YR2DAY ( J ) = 1.0 / 366.0
            END IF
150     CONTINUE

        LASTSCC = 0
        ACT_SCC = 0

        DO 60 J = 1, NSRCS
            POSN = SCCINDEX ( J )

            IF ( ISCC ( POSN ) .NE. LASTSCC ) THEN
                ! We've moved to a new SCC code...
                LASTSCC = ISCC ( POSN )

                ACT_SCC = ACT_SCC + 1
                IF ( ACT_SCC .GT. NPSCC ) THEN
                    WRITE ( BUF1, 600 ) ACT_SCC
                    WRITE ( BUF2, 600 ) NPSCC
                    CALL M3EXIT( 'PSCCSUM', 0, 0, 
     &                      'Number of SCC codes in PNTS (' // 
     &                      BUF1 ( LBLANK(BUF1):TRIMLEN(BUF1) ) //
     &                      ') exceeds program limits (NPSCC=' //
     &                      BUF2 ( LBLANK(BUF2):TRIMLEN(BUF2) ) // 
     &                      ')', 2 )
                END IF
                ISCCLIST ( ACT_SCC ) = ISCC ( POSN )
            END IF
            
            DO 70 I = 1, NIPOL
                SPECSUM ( ACT_SCC, I ) = SPECSUM ( ACT_SCC, I ) + 
     &                        ( SPECBUF ( POSN, I ) * YR2DAY ( POSN ) )
70          CONTINUE
60      CONTINUE


        ! Write the report of emissions summed by SCC code

        WRITE ( OUNIT, 800 ) 'Point source emissions (tons/day)'
        WRITE ( OUNIT, 850 ) '------------ ',
     &           'By SCC Code ------------'
        WRITE ( OUNIT, 810 ) 'SCC'
        WRITE ( OUNIT, 840 ) 'CODE      '
        DO 80 I = 1, NIPOL
            WRITE ( OUNIT, 820 ) EINAM ( I )( 1:TRIMLEN ( EINAM ( I ) ))
80      CONTINUE
        WRITE ( OUNIT, 810 ) ' '

        DO 90 J = 1, ACT_SCC
            WRITE ( OUNIT, 830 ) ISCCLIST ( J ), 
     &            ( SPECSUM ( J, I ), I = 1, NIPOL )
90      CONTINUE
        


        ! Write the report of emissions summed by Tier II SCC by state
        
        DO 140 J = 1, NSRCS
            ! We'll only need the state part of the FIP code.  So get rid 
            ! of the rest now.
            IFIPST ( J ) = IFIP ( J ) / 1000

            ! Also we'll only need the Tier II part of the SCC, so get rid
            ! of the lower four digits.
            ISCC2 ( J ) = ISCC ( J ) / 10000

            SCCINDEX ( J ) = J
140     CONTINUE
        CALL SORTI2 ( NSRCS, SCCINDEX, IFIPST, ISCC2 )

        DO 95 I = 1, NIPOL
            SPECSUM ( 1, I ) = 0.0
95      CONTINUE

        ! Write report header.
        WRITE ( RUNIT, 800 ) 'Point source emissions (tons/day)'
        WRITE ( RUNIT, 880 ) '------------------- ',
     &       'Tier II SCC by State -------------------'

        LASTST = 0
        LASTSCC = 0

        DO 100 J = 1, NSRCS
        	POSN = SCCINDEX ( J )
            IF ( IFIPST ( POSN ) .NE. LASTST ) THEN
                ! Lookup state name
                LASTST = IFIPST ( POSN )
                LOC = FIND1 ( LASTST, TOT_ST, ISTFIP )
                IF ( LOC .GT. 0 ) THEN
                	FIPNAME = STATES ( LOC )
                ELSE
                    FIPNAME = '*** no state name'
                ENDIF

                ! Write section header
                WRITE ( RUNIT, 860 ) LASTST, FIPNAME
                WRITE ( RUNIT, 810 ) 'SCC'
                WRITE ( RUNIT, 840 ) 'CODE  '
                DO 130 I = 1, NIPOL
                    WRITE ( RUNIT, 820 ) 
     &                    EINAM ( I )( 1:TRIMLEN ( EINAM ( I ) ) )
130             CONTINUE
                WRITE ( RUNIT, 810 ) '  TIER II SCC DESCRIPTION'
            END IF
            
            IF ( ISCC2 ( POSN ) .NE. LASTSCC ) THEN
                ! Look up new SCC Description & get teir II part.
                LASTSCC = ISCC2 ( POSN )

                ! Find this code in the master SCC list
                LOC = FIND1 ( ISCC ( POSN ), TOT_SCC, IALLSCC )

                IF ( LOC .LE. 0 ) THEN
                    SCCNAME = '*** no SCC description available'
                    L2 = 32 

                ELSE
                    L2 = INDEX ( SCCDESC ( LOC ), ',' ) - 1

                    SCCNAME = SCCDESC ( LOC )

                    IF ( L2 .LE. 0 ) THEN
                    	SCCNAME = '*** no SCC description available'
                        L2 = 32 
                    END IF 
                END IF

                DESC = SCCNAME ( 1 : MIN( L2, 45 ) )
            END IF 
            
            ! Sum the emissions 
            DO 110 I = 1, NIPOL
                SPECSUM ( 1, I ) = SPECSUM ( 1, I ) + 
     &                         ( SPECBUF ( POSN, I ) * YR2DAY ( POSN ) )
110         CONTINUE

            IF ( ( J .EQ. NSRCS ) .OR.
     &           ( IFIPST ( SCCINDEX ( J + 1 ) ) .NE. LASTST ) .OR.
     &           ( ISCC2 ( SCCINDEX ( J + 1 ) ) .NE. LASTSCC ) ) THEN
                ! Write the sum to the report
                WRITE ( RUNIT, 870 ) LASTSCC,
     &               ( SPECSUM ( 1, I ), I = 1, NIPOL ), DESC

                DO 120 I = 1, NIPOL
                    SPECSUM ( 1, I ) = 0.0
120             CONTINUE
            END IF
100     CONTINUE


999     CONTINUE          !  exit program

        MESG = 'Normal completion of PSCCSUM'
        CALL M3EXIT( 'PSCCSUM', 0, 0, MESG, 0 )


C******************  FORMAT  STATEMENTS   ******************************

600     FORMAT ( I10 )
700     FORMAT ( A10, 1X, A40 )
800     FORMAT ( A, / )
810     FORMAT ( A )
820     FORMAT ( A9, $ )
830     FORMAT ( I8.8, 2X, <NIPOL>F9.2 )
840     FORMAT ( A, $ )
850     FORMAT ( A, A, / )
860     FORMAT ( //, 15X, 'State Code ', I2, ':  ', A )
870     FORMAT ( I4, 2X, <NIPOL>F9.2, 2X, A )
880     FORMAT ( A, A )
        END
