
        PROGRAM ASCTSUM

C***********************************************************************
C  program body starts at line 113
C
C  DESCRIPTION:
C       Report area emissions summed by state and county.
C
C  PRECONDITIONS REQUIRED:  
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       Models-3 I/O; 
C       GETREAL, ENVREAL, PROMPTMFILE, PROMPTFFILE, TRIMLEN, LBLANK, 
C       ENVINT, FIND2
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

        INCLUDE 'ARDIMS3.EXT'   !  point-source dimensioning parameters
        INCLUDE 'CHDIMS3.EXT'   !  emis chem parms (inventory + model)
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.


C...........   EXTERNAL FUNCTIONS and their descriptions:

        INTEGER         ENVINT
        REAL            ENVREAL
        INTEGER         FIND2
        REAL            GETREAL
        INTEGER         PROMPTFFILE
        CHARACTER*16    PROMPTMFILE
        INTEGER         LBLANK
        INTEGER         TRIMLEN

        EXTERNAL        ENVINT, ENVREAL, FIND2, GETREAL, PROMPTFFILE, 
     &                  PROMPTMFILE, LBLANK, TRIMLEN
                                            
C...........   PARAMETERS and their descriptions
        INTEGER MAXASCT, NASCPOL
        PARAMETER ( MAXASCT = 3000, NASCPOL = NASCC * NIPOL )

C...........   LOCAL VARIABLES and their descriptions:
C...........   Point Sources input and output arrays

        INTEGER      IALLASC7 ( MAXASCT )  !  All ASCT codes
        INTEGER      IALLASC3 ( MAXASCT )
        CHARACTER*40 ASCTNAMES( MAXASCT )  !  All ASCT names
        LOGICAL      HAS_SUM  ( MAXASCT ) / MAXASCT * .FALSE. /
        INTEGER      IASC7LIST( NASCC )
        INTEGER      IASC3LIST( NASCC )
        REAL         SPECSUM  ( NASCC, NIPOL ) / NASCPOL * 0.0 /
        INTEGER      IASC7    ( NASRC )
        INTEGER      IASC3    ( NASRC )
        INTEGER      INVYR    ( NASRC )
        INTEGER      ASCINDEX ( NASRC )
        REAL         SPECBUF  ( NASRC, NIPOL )

C...........   Other local variables

        INTEGER         NSRCS

        INTEGER         LDEV    !  unit number for log file
        INTEGER         OUNIT, FUNIT

        REAL            YR2DAY
        INTEGER         I, J, K, INDEX
        INTEGER         LASTASC7, LASTASC3 ! The last  ASCT numbers
                                           ! encountered when summing the emmis.
        INTEGER         TOT_ASCT  ! Total number of ASCT codes/names 
                                  ! in the ASCTNAM file.
        INTEGER         ACT_ASCT  ! Actual number of ASCT codes
                                  ! in the area source file

        CHARACTER*16    ENAME   !  logical name for point-source input file
        CHARACTER*10    BUF1,BUF2
        INTEGER         INASC7, INASC3  
        CHARACTER*10    INASC    !  Buffer for reading the ASCT code
        CHARACTER*40    ASCTNAME !  Buffer for reading the ASCT name

        CHARACTER*256   MESG    !  buffer for M3EXIT() messages

C***********************************************************************
C   begin body of program ASCTSUM

        LDEV = INIT3()

        CALL INITEM( LDEV )

C.......   Get file name; open input area source

        ENAME = PROMPTMFILE( 
     &          'Enter name for AREA SOURCE input file',
     &          FSREAD3, 'AREA', 'ASCTSUM' )

        IF ( .NOT. DESC3( ENAME ) ) THEN
            CALL M3EXIT( 'ASCTSUM', 0, 0, 
     &                  'Could not get description of file "' 
     &                  // ENAME( 1:TRIMLEN( ENAME ) ) // '"', 2 )
        END IF

        NSRCS = NROWS3D
        IF ( NSRCS .GT. NASRC ) THEN
            WRITE ( BUF1, 600 ) NSRCS
            WRITE ( BUF2, 600 ) NASRC
            CALL M3EXIT( 'ASCTSUM', 0, 0, 
     &                   'Number of sources (NROWS3D='  // 
     &                   BUF1 ( LBLANK(BUF1):TRIMLEN(BUF1) ) //
     &                   ') exceeds program limits (NASRC=' //
     &                   BUF2 ( LBLANK(BUF2):TRIMLEN(BUF2) ) // 
     &                   ')', 2 )
        END IF

        FUNIT = PROMPTFFILE( 
     &          'Enter name for ASCT DESCRIPTION file',
     &          .TRUE., .TRUE., 'ASCTNAM', 'ASCTSUM' )

        OUNIT = PROMPTFFILE( 
     &          'Enter name for EMISSIONS SUM REPORT file',
     &          .FALSE., .TRUE., 'ASCTRPT', 'ASCTSUM' )



        TOT_ASCT = 0
30      CONTINUE ! Begin "infinite" loop; until EOF.
            READ ( FUNIT, 700, END = 40 ) INASC, ASCTNAME
C            PRINT *, INASC // ': ' // ASCTNAME
            TOT_ASCT = TOT_ASCT + 1
            IF ( TOT_ASCT .GT. MAXASCT ) THEN
                WRITE ( BUF1, 600 ) TOT_ASCT
                WRITE ( BUF2, 600 ) MAXASCT
                CALL M3EXIT( 'ASCTSUM', 0, 0, 
     &                      'Number of ASCT codes in ASCTNAM (' // 
     &                      BUF1 ( LBLANK(BUF1):TRIMLEN(BUF1) ) //
     &                      ') exceeds program limits (MAXASCT=' //
     &                      BUF2 ( LBLANK(BUF2):TRIMLEN(BUF2) ) // 
     &                      ')', 2 )
            END IF

            READ ( INASC ( 1:7 ), '(I7)' ) INASC7
            READ ( INASC ( 8:10 ), '(I3)' ) INASC3
            IALLASC7 ( TOT_ASCT ) = INASC7
            IALLASC3 ( TOT_ASCT ) = INASC3
            ASCTNAMES ( TOT_ASCT ) = ASCTNAME
C            PRINT '( A, I2, A, A, A, I7.7, I3.3, A )', 
C     &            'Got ASCT #', TOT_ASCT, ': ',
C     &            ASCTNAMES ( TOT_ASCT ), ' (ASC=', 
C     &            IALLASC7 ( TOT_ASCT ), IALLASC3 ( TOT_ASCT ), ')'
        GO TO 30
40      CONTINUE

        IF ( .NOT. READ3 ( ENAME, 'ASC7', ALLAYS3,0,0, IASC7 ) ) THEN
            CALL M3EXIT( 'ASCTSUM', 0, 0, 
     &                   'Error reading IASC7 from file' // ENAME, 2)
        END IF

        IF ( .NOT. READ3 ( ENAME, 'ASC3', ALLAYS3,0,0, IASC3 ) ) THEN
            CALL M3EXIT( 'ASCTSUM', 0, 0, 
     &                   'Error reading IASC3 from file' // ENAME, 2)
        END IF

        IF ( .NOT. READ3 ( ENAME, 'INVYR', ALLAYS3,0,0, INVYR ) ) THEN
            CALL M3EXIT( 'AREASUM', 0, 0, 
     &                   'Error reading INVYR from file' // ENAME, 2)
        END IF


        DO 50 I = 1, NIPOL
C           PRINT *, 'Reading ' // EINAM (I) // '...'
            IF ( .NOT. READ3 ( ENAME, EINAM ( I ), ALLAYS3, 0, 0, 
     &                         SPECBUF ( 1, I ) ) ) THEN
                CALL M3EXIT( 'ASCTSUM', 0, 0, 
     &                       'Error reading ' // EINAM ( I ) //
     &                       ' from file' // ENAME, 2 )
            END IF
50      CONTINUE

        DO 45 J = 1, NSRCS
            ASCINDEX ( J ) = J
45      CONTINUE
        CALL SORTI2 ( NSRCS, ASCINDEX, IASC7, IASC3 )


        LASTASC3 = 0
        LASTASC7 = 0
        ACT_ASCT = 0

        DO 60 J = 1, NSRCS
            INDEX = ASCINDEX ( J )

            K = INVYR( INDEX )
            IF ( MOD( K,4 ) .NE. 0 ) THEN           !  nonleap years
                YR2DAY = 1.0 / 365.0
            ELSE IF ( MOD( K,100 ) .NE. 0 ) THEN    !  noncentury leap years
                YR2DAY = 1.0 / 366.0
            ELSE IF ( MOD( K,400 ) .NE. 0 ) THEN    !  century nonleap years
                YR2DAY = 1.0 / 365.0
            ELSE                                    !  leap centuries
                YR2DAY = 1.0 / 366.0
            END IF

            IF ( ( IASC3 ( INDEX ) .NE. LASTASC3 ) .OR. 
     &           ( IASC7 ( INDEX ) .NE. LASTASC7 ) ) THEN
                ! We've moved to a new ASCT code...
                LASTASC7 = IASC7 ( INDEX )
                LASTASC3 = IASC3 ( INDEX )

                ACT_ASCT = ACT_ASCT + 1
                IF ( ACT_ASCT .GT. NASCC ) THEN
                    WRITE ( BUF1, 600 ) ACT_ASCT
                    WRITE ( BUF2, 600 ) NASCC
                    CALL M3EXIT( 'ASCTSUM', 0, 0, 
     &                      'Number of ASCT codes in AREA (' // 
     &                      BUF1 ( LBLANK(BUF1):TRIMLEN(BUF1) ) //
     &                      ') exceeds program limits (NASCC=' //
     &                      BUF2 ( LBLANK(BUF2):TRIMLEN(BUF2) ) // 
     &                      ')', 2 )
                END IF
                IASC7LIST ( ACT_ASCT ) = IASC7 ( INDEX )
                IASC3LIST ( ACT_ASCT ) = IASC3 ( INDEX )
            END IF
            
            DO 70 I = 1, NIPOL
                SPECSUM ( ACT_ASCT, I ) = SPECSUM ( ACT_ASCT, I ) + 
     &                                 ( SPECBUF ( INDEX, I ) * YR2DAY )
70          CONTINUE
60      CONTINUE


        WRITE ( OUNIT, 800 ) 'Area source emissions (tons/day)'
        WRITE ( OUNIT, 850 ) '------------------------------ ',
     &            'By ASCT Code ------------------------------'
        WRITE ( OUNIT, 810 ) 'ASCT'
        WRITE ( OUNIT, 840 ) 'CODE       '
        DO 80 I = 1, NIPOL
            WRITE ( OUNIT, 820 ) EINAM ( I )( 1:TRIMLEN ( EINAM ( I ) ))
80      CONTINUE
        WRITE ( OUNIT, 810 ) '  ASCT DESCRIPTION'

        DO 90 J = 1, ACT_ASCT
            ! Find this code in the master ASCT list
            INDEX = FIND2 ( IASC7LIST ( J ), IASC3LIST ( J ), TOT_ASCT, 
     &                      IALLASC7, IALLASC3 )

            IF ( INDEX .LE. 0 ) THEN
                ASCTNAME = "*** no ASCT name available"
            ELSE
                ASCTNAME = ASCTNAMES ( INDEX )
                IF ( ASCTNAME .EQ. ' ' ) THEN
                	ASCTNAME = "*** no ASCT name available"
                END IF 
            END IF

            WRITE ( OUNIT, 830 ) IASC7LIST ( J ), IASC3LIST ( J ), 
     &            ( SPECSUM ( J, I ), I = 1, NIPOL ), ASCTNAME
90      CONTINUE
        

        
999     CONTINUE          !  exit program

        MESG = 'Normal completion of ASCTSUM'
        CALL M3EXIT( 'ASCTSUM', 0, 0, MESG, 0 )


C******************  FORMAT  STATEMENTS   ******************************

600     FORMAT ( I10 )
700     FORMAT ( A10, 1X, A40 )
800     FORMAT ( A, / )
810     FORMAT ( A )
820     FORMAT ( A9, $ )
830     FORMAT ( I7.7, I3.3, 1X, <NIPOL>F9.2, 2X, A )
840     FORMAT ( A, $ )
850     FORMAT ( A, A, / )

        END
