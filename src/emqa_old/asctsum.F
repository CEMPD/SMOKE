
        PROGRAM ASCTSUM

C***********************************************************************
C  program body starts at line 113
C
C  DESCRIPTION:
C       Report area emissions summed by ASCT and state.
C
C  PRECONDITIONS REQUIRED:  
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       Models-3 I/O; 
C       GETREAL, PROMPTMFILE, PROMPTFFILE, TRIMLEN, LBLANK, 
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

        CHARACTER*2     CRLF
        LOGICAL         ENVYN
        INTEGER         FIND2
        REAL            GETREAL
        INTEGER         PROMPTFFILE
        CHARACTER*16    PROMPTMFILE
        INTEGER         LBLANK
        INTEGER         TRIMLEN
        REAL            YR2DAY

        EXTERNAL        CRLF, ENVYN, FIND2, GETREAL, PROMPTFFILE, 
     &                  PROMPTMFILE, LBLANK, TRIMLEN, YR2DAY
                                            
C...........   LOCAL PARAMETERS 
        INTEGER     MAXASCT, NASCPOL
        CHARACTER*5 BLANK5

        PARAMETER ( MAXASCT = 3000
     &            , NASCPOL = NASID * NASCC * NIPOL
     &            , BLANK5  = ' '
     &            )

C...........   LOCAL VARIABLES and their descriptions:

C...........   Area Sources inventory characteristics
        INTEGER      SINDX    ( NASRC )  ! Sorting index
        INTEGER      IFIPS    ( NASRC )
        INTEGER      STIDS    ( NASRC )
        INTEGER      IASC7    ( NASRC )
        INTEGER      IASC3    ( NASRC )
        INTEGER      INVYR    ( NASRC )

C...........   Area Sources inventory emissions
        REAL         SPECBUF  ( NASRC, NIPOL )

C...........   ASCTs and FIPS found in the inventory for summing
        INTEGER      STATLIST( NASID )  ! list of state codes
        INTEGER      NASCT   ( NASID )  ! no. ASCTs in inventory
        INTEGER      IASC7LST( NASID, NASCC )
        INTEGER      IASC3LST( NASID, NASCC )

C...........   ASCT names table
        INTEGER      IALLASC7 ( MAXASCT )  !  All ASCT codes
        INTEGER      IALLASC3 ( MAXASCT )
        CHARACTER*40 ASCTNAMES( MAXASCT )  !  All ASCT names

C...........   Output emissions
        REAL         SPECSUM  ( NASID, NASCC, NIPOL ) / NASCPOL * 0.0 /

C...........   File names and unit numbers
        INTEGER         LDEV    !  unit number for log file
        INTEGER         OUNIT, FUNIT
        CHARACTER*16    ENAME   !  logical name for point-source input file

C...........   Other local variables

        INTEGER         F, I, J, K, S, N   ! counters and indices

        INTEGER         ID3, ID7, STA ! tmp source characteristics
        INTEGER         INASC7, INASC3 
        INTEGER         IOS                ! i/o status
        INTEGER         LASTASC7, LASTASC3 ! Previous ASCTs in loop
        INTEGER         LSTA               ! Previous state code in loop
        INTEGER         MAXN               ! max number of ASCTs by state
        INTEGER         NSRCS              ! No. sources
        INTEGER         NSTATE             ! No. states in inventory
        INTEGER         TOT_ASCT           ! No. of ASCT codes in  ASCTNAM file
        INTEGER         YR                 ! tmp inventory year

        LOGICAL         SFLAG              ! true: output by state codes
        LOGICAL         EFLAG              ! true: error found
        DATA            EFLAG / .FALSE. / 

        CHARACTER*10    BUF1,BUF2
        CHARACTER*10    INASC    !  Buffer for reading the ASCT code
        CHARACTER*40    ASCTNAME !  Buffer for reading the ASCT name

        CHARACTER*256   MESG    !  buffer for M3EXIT() messages
        CHARACTER*256   FMT

C***********************************************************************
C   begin body of program ASCTSUM

        LDEV = INIT3()

        CALL INITEM( LDEV )

C.......   Get environment variable settings
        MESG = 'Sum area sources by ASCT and state codes'
        SFLAG = ENVYN ( 'ASCTSUM_BYSTATE', MESG, .TRUE., IOS )

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


C.........  Read in list of ASCs and the descriptions
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

C.........  Read in inventory characteristics

        IF ( .NOT. READ3 ( ENAME, 'FIP', ALLAYS3,0,0, IFIPS ) ) THEN
            CALL M3EXIT( 'ASCTSUM', 0, 0, 
     &                   'Error reading FIP from file' // ENAME, 2)
        END IF

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

C.........  Read in inventory emissions
        DO 50 I = 1, NIPOL
            IF ( .NOT. READ3 ( ENAME, EINAM ( I ), ALLAYS3, 0, 0, 
     &                         SPECBUF ( 1, I ) ) ) THEN
                CALL M3EXIT( 'ASCTSUM', 0, 0, 
     &                       'Error reading ' // EINAM ( I ) //
     &                       ' from file' // ENAME, 2 )
            END IF
50      CONTINUE

C.........  Initialize index for sorting inventory and state codes
        DO 45 J = 1, NSRCS
            SINDX ( J ) = J
            STIDS ( J ) = IFIPS( J ) / 1000
45      CONTINUE

        IF( SFLAG ) THEN
            CALL SORTI3 ( NSRCS, SINDX, STIDS, IASC7, IASC3 )

C.........  If output is not by state, resort the inventory by ASCT
        ELSE 
            CALL SORTI2 ( NSRCS, SINDX, IASC7, IASC3 )
        END IF

        F = 0
        LSTA = -9
        MAXN = 0
        DO 60 J = 1, NSRCS

            S = SINDX( J )

            YR  = INVYR( S )
            STA = STIDS( S ) 
            ID7 = IASC7( S )
            ID3 = IASC3( S )

            IF( J .EQ. 1 .OR.
     &        ( SFLAG .AND. STA .NE. LSTA ) ) THEN
                F = F + 1
                N = 0
                LASTASC3 = 0
                LASTASC7 = 0
                LSTA = STA
            END IF

C.............  If we've moved to a new ASCT code...
            IF ( ( ID7 .NE. LASTASC7 ) .OR. 
     &           ( ID3 .NE. LASTASC3 ) ) THEN

                LASTASC7 = ID7
                LASTASC3 = ID3

                N = N + 1
                IF ( N .LE. NASCC .AND. F .LE. NASID ) THEN
                    IASC7LST ( F,N ) = ID7
                    IASC3LST ( F,N ) = ID3
                    STATLIST  ( F )  = STA 
                ELSE
                    IF( N .GT. MAXN ) MAXN = N
                    GO TO 60
                END IF

            END IF

            DO 70 I = 1, NIPOL
                SPECSUM ( F,N,I ) = SPECSUM ( F,N,I ) + 
     &                              SPECBUF ( S,I ) * YR2DAY( YR )
70          CONTINUE

            NASCT( F ) = N

60      CONTINUE

        IF( MAXN .GT. NASCC ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 
     &             'ERROR: Number of ASCT codes in AREA file: ', MAXN,
     &             CRLF() // BLANK5 //
     &             'exceeds program limits, NASCC=', NASCC
            CALL M3MSG2( MESG )

        END IF

        IF( F .GT. NASID ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 
     &             'ERROR: Number of state codes in AREA file: ', F,
     &             CRLF() // BLANK5 //
     &             'exceeds program limits, NASID=', NASID
            CALL M3MSG2( MESG )

        ELSE
            NSTATE = F
        END IF

        IF( EFLAG ) THEN
            MESG = 'Problem with program dimensions'
            CALL M3EXIT( 'ASCTSUM', 0, 0, MESG, 2 )
        END IF

C.........  Output emissions by State code and/or ASCT code  

        WRITE ( OUNIT, 800 ) 'Area source emissions (tons/day)'
        WRITE ( OUNIT, 850 ) '------------------------------ ',
     &            'By ASCT Code ------------------------------'
        WRITE ( OUNIT, 810 ) 'ASCT'
        WRITE ( OUNIT, 840 ) 'CODE       '
        DO 80 I = 1, NIPOL
            WRITE ( OUNIT, 820 ) EINAM ( I )( 1:TRIMLEN ( EINAM ( I ) ))
80      CONTINUE
        WRITE ( OUNIT, 810 ) '  ASCT DESCRIPTION'

        DO 101 F = 1, NSTATE 

            IF( SFLAG ) THEN
                WRITE( OUNIT, 93050 ) 'For state code: ', STATLIST( F )
                WRITE ( OUNIT, 850 ) '-------------------------'
            ENDIF

            DO 90 J = 1, NASCT( F )

C................. Find this code in the master ASCT list to get the name
                K = FIND2 ( IASC7LST ( F,J ), IASC3LST ( F,J ), 
     &                      TOT_ASCT, IALLASC7, IALLASC3 )

                IF ( K .LE. 0 ) THEN
                    ASCTNAME = "*** no ASCT name available"
                ELSE
                    ASCTNAME = ASCTNAMES ( K )
                    IF ( ASCTNAME .EQ. ' ' ) THEN
                	ASCTNAME = "*** no ASCT name available"
                    END IF 
                END IF

                WRITE ( FMT, 831 ) NIPOL
                WRITE ( OUNIT, FMT ) IASC7LST ( F,J ), IASC3LST ( F,J ),
     &                ( SPECSUM ( F,J,I ), I = 1, NIPOL ), ASCTNAME

90          CONTINUE
        
101     CONTINUE
        
999     CONTINUE          !  exit program

C......... Normal completion
        CALL M3EXIT( 'ASCTSUM', 0, 0, ' ', 0 )


C******************  FORMAT  STATEMENTS   ******************************

600     FORMAT ( I10 )
700     FORMAT ( A10, 1X, A40 )
800     FORMAT ( A, / )
810     FORMAT ( A )
820     FORMAT ( A9, $ )
C830    FORMAT ( I7.7, I3.3, 1X, <NIPOL>F9.2, 2X, A )
831     FORMAT ( '( I7.7, I3.3, 1X, ', I8, 'F9.2, 2X, A )' )
840     FORMAT ( A, $ )
850     FORMAT ( A, A, / )

93050   FORMAT ( A, I2.2 )

94010   FORMAT ( 10( A, :, I8, :, 1X ) )

        END
