        PROGRAM PLANTSUM

C***********************************************************************
C  program body starts at line 123
C
C  DESCRIPTION:
C       Report point source emissions summed by plant within state.
C       Report the Top N NOx producing plants, and report the plants
C       with NOx or VOC > x.
C
C  PRECONDITIONS REQUIRED:  
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       Models-3 I/O; 
C       GETREAL, ENVREAL, PROMPTMFILE, PROMPTFFILE, TRIMLEN, LBLANK, 
C       ENVINT, FIND2, SORTR1, SORTI2
C
C  REVISION  HISTORY:
C       Prototype  3/97 by Chad Barker
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

        INCLUDE 'PTDIMS3.EXT'   !  point-source dimensioning parameters
        INCLUDE 'CHDIMS3.EXT'   !  emis chem parms (inventory + model)
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.


C...........   EXTERNAL FUNCTIONS and their descriptions:

        REAL            GETREAL, ENVREAL
        CHARACTER*16    PROMPTMFILE
        INTEGER         LBLANK, TRIMLEN, PROMPTFFILE, ENVINT
        INTEGER         GETNUM, FIND1

        EXTERNAL     GETREAL, ENVREAL, PROMPTMFILE, TRIMLEN, GETNUM
        EXTERNAL     LBLANK, PROMPTFFILE, ENVINT, FIND1, SORTI2, SORTR1
     
                                            
C...........   PARAMETERS and their descriptions
        INTEGER NSRCPOL
        PARAMETER ( NSRCPOL = NPSRC * NIPOL )

        INTEGER MAXST
        PARAMETER ( MAXST = 60 )

C...........   LOCAL VARIABLES and their descriptions:
C...........   Point Sources input and output arrays

        INTEGER       ISTFIP   ( MAXST )  !  All state FIPS IDs
        CHARACTER*20  STATES   ( MAXST )  !  All state names
        INTEGER       SORTINDX ( NPSRC )
        INTEGER       IFIPST   ( NPSRC )  ! The state part of IFIP
        INTEGER       IFIP     ( NPSRC )
        INTEGER       INVYR    ( NPSRC )
        INTEGER       IPLANT   ( NPSRC )  !  Plant ID
        REAL          SPECBUF  ( NPSRC, NIPOL )
        REAL          YR2DAY   ( NPSRC )
        
        ! Summarized lists of states, plants and sum of emissions
        INTEGER       STLIST   ( NPSRC )
        INTEGER       PLANTLST ( NPSRC )
        REAL          SPECSUM  ( NPSRC, NIPOL ) / NSRCPOL * 0.0 /
        INTEGER       ACTPLANT ! Number of entries in the above lists
        
        ! These are used to sort the above lists
        REAL          TMPSUM   ( NPSRC )


C...........   Other local variables

        INTEGER         NSRCS
        INTEGER         STATUS, LISTSIZE
        REAL            MAXEMIS
        INTEGER         NOXVAR, VOCVAR

        INTEGER         LDEV    !  unit number for log file
        INTEGER         OUNIT, FUNIT, RUNIT

        INTEGER         I, J, K, L, POSN, LOC
        INTEGER         TOT_ST, LASTST
        INTEGER         LASTPLANT
        INTEGER         FIPNUM  !  Buffer for reading the FIPS code
        CHARACTER*20    FIPNAME !  Buffer for reading the state/county name

        CHARACTER*16    ENAME   !  logical name for point-source input file
        CHARACTER*10    BUF1,BUF2

        CHARACTER*256   MESG    !  buffer for M3EXIT() messages


C***********************************************************************
C   begin body of program PLANTSUM

        LDEV = INIT3()

C.......   Get file name; open input area source

        ENAME = PROMPTMFILE( 
     &          'Enter name for POINT SOURCE input file',
     &          FSREAD3, 'PNTS', 'PLANTSUM' )

        IF ( .NOT. DESC3( ENAME ) ) THEN
            CALL M3EXIT( 'PLANTSUM', 0, 0, 
     &                  'Could not get description of file "' 
     &                  // ENAME( 1:TRIMLEN( ENAME ) ) // '"', 2 )
        END IF

        NSRCS = NROWS3D
        IF ( NSRCS .GT. NPSRC ) THEN
            WRITE ( BUF1, 600 ) NSRCS
            WRITE ( BUF2, 600 ) NPSRC
            CALL M3EXIT( 'PLANTSUM', 0, 0, 
     &                   'Number of sources (NROWS3D='  // 
     &                   BUF1 ( LBLANK(BUF1):TRIMLEN(BUF1) ) //
     &                   ') exceeds program limits (NPSRC=' //
     &                   BUF2 ( LBLANK(BUF2):TRIMLEN(BUF2) ) // 
     &                   ')', 2 )
        END IF

        DO I = 1, NIPOL
C            WRITE (*,730) VNAME3D ( I )
            IF ( EINAM ( I ) .EQ. 'NOX' ) THEN
                NOXVAR = I
            ELSE IF ( ( EINAM ( I ) .EQ. 'VOC' ) .OR.
     &                ( EINAM ( I ) .EQ. 'ROG' ) ) THEN
                VOCVAR = I
            END IF
        END DO

        FUNIT = PROMPTFFILE( 
     &          'Enter name for STATE/COUNTY FIPS CODE file',
     &          .TRUE., .TRUE., 'STCYNAM', 'AREASUM' )

        LISTSIZE = ENVINT ( 'STLSTSZ', 'Max plants per state', NPSRC, 
     &                      STATUS )
        LISTSIZE = GETNUM ( 1, NPSRC, LISTSIZE, 
     &     'Enter maximum number of plants per state group in report' )

        MAXEMIS = ENVREAL ( 'MAXEMIS', 'NOX cutoff point', 1.0, STATUS )
        MAXEMIS = GETREAL ( 0.0, 99999.0, MAXEMIS, 
     &             'Enter cutoff point for NOX or ' // 
     &             EINAM ( VOCVAR )( 1:TRIMLEN ( EINAM ( VOCVAR ) ) ) //
     &             ' emissions' )

        OUNIT = PROMPTFFILE( 
     &          'Enter name for TOP NOX PLANTS REPORT file',
     &          .FALSE., .TRUE., 'PTOPRPT', 'PLANTSUM' )

        RUNIT = PROMPTFFILE( 
     &          'Enter name for LARGE PLANT EMISSIONS REPORT file',
     &          .FALSE., .TRUE., 'PLRGRPT', 'PLANTSUM' )



        TOT_ST = 0
160     CONTINUE ! Begin "infinite" loop; until EOF.
            READ( FUNIT, '( 2X, I5, 1X, A )', END=170 ) FIPNUM, FIPNAME
            IF ( FIPNUM .NE. 0 ) THEN  ! If FIPNUM=0, we just read a blank line
                IF ( MOD ( FIPNUM, 1000 ) .EQ. 0 ) THEN
                    TOT_ST = TOT_ST + 1
                    IF ( TOT_ST .GT. MAXST ) THEN
                        WRITE ( BUF1, 600 ) TOT_ST
                        WRITE ( BUF2, 600 ) MAXST
                        CALL M3EXIT( 'AREASUM', 0, 0, 
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
170     CONTINUE


        IF ( .NOT. READ3 ( ENAME, 'IFIP', ALLAYS3, 0, 0, IFIP ) ) THEN
            CALL M3EXIT( 'PLANTSUM', 0, 0, 
     &                   'Error reading IFIP from file' // ENAME, 2)
        END IF

        IF ( .NOT. READ3 ( ENAME, 'IPLANT', ALLAYS3,0,0, IPLANT ) ) THEN
            CALL M3EXIT( 'PLANTSUM', 0, 0, 
     &                   'Error reading IPLANT from file' // ENAME, 2)
        END IF

        IF ( .NOT. READ3 ( ENAME, 'INVYR', ALLAYS3, 0, 0, INVYR ) ) THEN
            CALL M3EXIT( 'PLANTSUM', 0, 0, 
     &                   'Error reading INVYR from file' // ENAME, 2)
        END IF

C        PRINT *, 'Reading emissions...'

        DO 50 I = 1, NIPOL
C           PRINT *, 'Reading ' // EINAM (I) // '...'
            IF ( .NOT. READ3 ( ENAME, EINAM ( I ), ALLAYS3, 0, 0, 
     &                         SPECBUF ( 1, I ) ) ) THEN
                CALL M3EXIT( 'PLANTSUM', 0, 0, 
     &                       'Error reading ' // EINAM ( I ) //
     &                       ' from file' // ENAME, 2 )
            END IF
50      CONTINUE

C        PRINT *, 'Sorting IFIPST, IPLANT...'

        DO 45 J = 1, NSRCS
            ! We'll only need the state part of the FIP code.  So get rid 
            ! of the rest now.
            IFIPST ( J ) = IFIP ( J ) / 1000

            SORTINDX ( J ) = J
45      CONTINUE
        CALL SORTI2 ( NSRCS, SORTINDX, IFIPST, IPLANT )


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

        LASTST = 0
        LASTPLANT = 0
        ACTPLANT = 0

        DO 60 J = 1, NSRCS
            POSN = SORTINDX ( J )

            ! We're assuming that a plant number will not be the same
            ! when we move from one state to the next.
            IF ( ( IPLANT ( POSN ) .NE. LASTPLANT ) .OR. 
     &           ( IFIPST ( POSN ) .NE. LASTST ) ) THEN
                ! We've moved to a new PLANT...
                LASTPLANT = IPLANT ( POSN )
                LASTST = IFIPST ( POSN )

                ACTPLANT = ACTPLANT + 1
                PLANTLST ( ACTPLANT ) = LASTPLANT
                STLIST ( ACTPLANT ) = LASTST
            END IF
            
            DO 70 I = 1, NIPOL
                SPECSUM ( ACTPLANT, I ) = SPECSUM ( ACTPLANT, I ) + 
     &                        ( SPECBUF ( POSN, I ) * YR2DAY ( POSN ) )
70          CONTINUE
60      CONTINUE


        ! Write the report of Top N NOx producing plants per state

        WRITE ( OUNIT, 810 ) 'Point source emissions (tons/day)'
        WRITE ( OUNIT, 890 ) '    Top', LISTSIZE, 'NOx producing plants'
        WRITE ( OUNIT, 850 ) '----------------- ',
     &           'By State -----------------'


        ! Loop through every state in STLIST
		J = 1
200     CONTINUE
            LASTST = STLIST ( J )
            LOC = FIND1 ( LASTST, TOT_ST, ISTFIP )
            IF ( LOC .GT. 0 ) THEN
                FIPNAME = STATES ( LOC )
            ELSE
                FIPNAME = '*** no state name'
            ENDIF

            ! Write section header
            WRITE ( OUNIT, 860 ) LASTST, FIPNAME
            WRITE ( OUNIT, 840 ) ' PLANT  '
            DO 80 I = 1, NIPOL
                WRITE ( OUNIT, 820 ) 
     &                    EINAM ( I )( 1:TRIMLEN ( EINAM ( I ) ) )
80          CONTINUE
            WRITE ( OUNIT, 810 ) ' '

            ! Create a list of NOX sums for the current state 
            K = 0
220         CONTINUE
                K = K + 1
                TMPSUM ( K ) = SPECSUM ( J, NOXVAR )
                SORTINDX ( K ) = K

                J = J + 1
                IF ( ( J .GT. ACTPLANT ) .OR. 
     &               ( STLIST ( J ) .NE. LASTST ) ) THEN
                    GO TO 230
                END IF
            GO TO 220
230         CONTINUE

            ! Sort the list just created
            CALL SORTR1 ( K, SORTINDX, TMPSUM )
            
            ! Write the Top LISTSIZE entries by reading the list in reverse
            DO 240 L = K, MAX( K - LISTSIZE + 1 , 1 ), -1
                POSN = SORTINDX ( L ) + J - K - 1
                WRITE ( OUNIT, 830 ) PLANTLST ( POSN ), 
     &                ( SPECSUM ( POSN, I ), I = 1, NIPOL )
240         CONTINUE

            IF ( J .GT. ACTPLANT ) THEN
                GO TO 210
            END IF
        GO TO 200
210     CONTINUE



        ! Write the report of Large emissions sources
        
        WRITE ( RUNIT, 810 ) 'Point source emissions (tons/day)'
        WRITE ( RUNIT, 900 ) '   All plants producing >', MAXEMIS, 
     &      'tons/day of NOX or ', 
     &      EINAM ( VOCVAR )( 1:TRIMLEN ( EINAM ( VOCVAR ) ) )
        WRITE ( RUNIT, 880 ) '----------------- ',
     &       'By State -----------------'

        LASTST = 0

        DO 100 J = 1, ACTPLANT
            IF ( STLIST ( J ) .NE. LASTST ) THEN
                ! Lookup state name
                LASTST = STLIST ( J )
                LOC = FIND1 ( LASTST, TOT_ST, ISTFIP )
                IF ( LOC .GT. 0 ) THEN
                	FIPNAME = STATES ( LOC )
                ELSE
                    FIPNAME = '*** no state name'
                ENDIF

                ! Write section header
                WRITE ( RUNIT, 860 ) LASTST, FIPNAME
                WRITE ( RUNIT, 840 ) ' PLANT  '
                DO 130 I = 1, NIPOL
                    WRITE ( RUNIT, 820 ) 
     &                    EINAM ( I )( 1:TRIMLEN ( EINAM ( I ) ) )
130             CONTINUE
                WRITE ( RUNIT, 810 ) ' '
            END IF

            IF ( ( SPECSUM ( J, NOXVAR ) .GT. MAXEMIS ) .OR.
     &           ( SPECSUM ( J, VOCVAR ) .GT. MAXEMIS ) ) THEN
                ! Write the sum to the report
                WRITE ( RUNIT, 830 ) PLANTLST ( J ),
     &               ( SPECSUM ( J, I ), I = 1, NIPOL )
            END IF
100     CONTINUE


999     CONTINUE          !  exit program

        MESG = 'Normal completion of PLANTSUM'
        CALL M3EXIT( 'PLANTSUM', 0, 0, MESG, 0 )


C******************  FORMAT  STATEMENTS   ******************************

600     FORMAT ( I10 )
700     FORMAT ( A10, 1X, A40 )
800     FORMAT ( A, / )
810     FORMAT ( A )
820     FORMAT ( A12, $ )
830     FORMAT ( I6, 2X, <NIPOL>F12.5 )
840     FORMAT ( A, $ )
850     FORMAT ( A, A, / )
860     FORMAT ( //, 8X, 'State Code ', I2, ':  ', A, / )
880     FORMAT ( A, A )
890     FORMAT ( A, 1X, I3, 1X, A, / )
900     FORMAT ( A, F5.1, 1X, A, A, / )
        END
