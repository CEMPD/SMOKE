        PROGRAM AOICSUM

C***********************************************************************
C  program body starts at line 150
C
C  DESCRIPTION:
C       Report emissions summed by state and county.
C
C  PRECONDITIONS REQUIRED:  
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       Models-3 I/O; 
C       GETREAL, ENVREAL, PROMPTMFILE, PROMPTFFILE, TRIMLEN, LBLANK, 
C       ENVINT, FIND1, FIND2, SORTR1, GETNUM
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

        INCLUDE 'ARDIMS3.EXT'   !  point-source dimensioning parameters
        INCLUDE 'CHDIMS3.EXT'   !  emis chem parms (inventory + model)
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.


C...........   EXTERNAL FUNCTIONS and their descriptions:

        REAL            GETREAL, ENVREAL
        CHARACTER*16    PROMPTMFILE
        INTEGER         LBLANK, TRIMLEN, PROMPTFFILE, ENVINT
        INTEGER         FIND1, FIND2, GETNUM

        EXTERNAL        GETREAL, ENVREAL, PROMPTMFILE, TRIMLEN
        EXTERNAL        LBLANK, PROMPTFFILE, ENVINT, FIND1, FIND2
        EXTERNAL        SORTR1, GETNUM
     
                                            
C...........   PARAMETERS and their descriptions
        INTEGER MAXST, MAXCO ! Max number of states and counties
        PARAMETER ( MAXST = 60, MAXCO = 3500 )
        
        INTEGER MAXOIC, MAXASCT
        PARAMETER ( MAXOIC = 100, MAXASCT = 3000 )

        INTEGER NSTPOL, NCOPOL, NSTOIC, NCOOIC
        PARAMETER ( NSTPOL = MAXOIC * NASID * NIPOL, 
     &              NCOPOL = MAXOIC * NAFIP * NIPOL,
     &              NSTOIC = MAXOIC * NASID, NCOOIC = MAXOIC * NAFIP )

C...........   LOCAL VARIABLES and their descriptions:
C...........   Area Sources input and output arrays

        INTEGER      ISTFIP   ( MAXST )  !  All state FIPS IDs
        CHARACTER*20 STATES   ( MAXST )  !  All state names
        INTEGER      ICOFIP   ( MAXCO )  !  All county FIPS IDs
        CHARACTER*20 COUNTIES ( MAXCO )  !  All county names
        INTEGER      IALLASC7 ( MAXASCT ) !  All ASCT codes
        INTEGER      IALLASC3 ( MAXASCT )
        CHARACTER*40 ASCTNAMES( MAXASCT ) !  All ASCT names

        INTEGER      ISTLIST ( NASID )  !  List of states in AREA file
        INTEGER      ICOLIST ( NAFIP )  !  List of counties in AREA file
        INTEGER      IOICLIST( MAXOIC ) !  List of OIC codes in ASCC file

        LOGICAL      IN_RPT  ( MAXOIC ) / MAXOIC * .FALSE. / 
                            !  Indicates which OIC codes are reportable
        LOGICAL      HASSTSUM( MAXOIC, NASID ) / NSTOIC * .FALSE. /
        LOGICAL      HASCOSUM( MAXOIC, NAFIP ) / NCOOIC * .FALSE. /
        REAL         STSUM   ( MAXOIC, NASID, NIPOL ) / NSTPOL * 0.0 /
        REAL         COSUM   ( MAXOIC, NAFIP, NIPOL ) / NCOPOL * 0.0 /
        
        INTEGER      SORTINDEX ( NAFIP )
        REAL         TMPEMIS    ( NAFIP )

        INTEGER      IFIP    ( NASRC )
        INTEGER      IASC7   ( NASRC )
        INTEGER      IASC3   ( NASRC )
        INTEGER      INVYR   ( NASRC )
        REAL         SPECBUF ( NASRC, NIPOL )

C...........   Other local variables

        INTEGER         NSRCS

        INTEGER         LDEV    !  unit number for log file
        INTEGER         STUNIT, COUNIT, FUNIT, AUNIT, SUNIT, LUNIT

        REAL            YR2DAY
        INTEGER         I, J, K, INDEX, STATUS, POSN
        INTEGER         NOX_VAR / 0 /
        INTEGER         VOC_VAR / 0 /
        INTEGER         EMIS_VAR
        CHARACTER*1     RPT_NUM
        INTEGER         LASTCO, LASTST  ! The last state and county FIPS numbers
                                        ! encountered when summing the emmis.
        INTEGER         LASTOIC, CUROIC
        INTEGER         TOT_ST, TOT_CO  ! Total number of states and counties 
                                        ! in the FIPS file.
        INTEGER         TOT_OIC, TOT_ASCT
        INTEGER         COLISTSZ, COWRITTEN

        INTEGER         ACT_ST, ACT_CO  ! Actual number of states and counties
                                        ! in the area source file

        CHARACTER*16    ENAME   !  logical name for point-source input file
        CHARACTER*10    BUF1,BUF2
        INTEGER         FIPNUM  !  Buffer for reading the FIPS code
        CHARACTER*20    FIPNAME !  Buffer for reading the state/county name
        CHARACTER*20    FIPNAME2
        INTEGER         INASC7, INASC3, INOIC  
        CHARACTER*10    INASC    !  Buffer for reading the ASCT code
        CHARACTER*40    ASCTNAME !  Buffer for reading the ASCT name

        CHARACTER*256   MESG    !  buffer for M3EXIT() messages


C***********************************************************************
C   begin body of program REPPLHT

        LDEV = INIT3()

C.......   Get file name; open input area source

        ENAME = PROMPTMFILE( 
     &          'Enter name for AREA SOURCE input file',
     &          FSREAD3, 'AREA', 'AOICSUM' )

        IF ( .NOT. DESC3( ENAME ) ) THEN
            CALL M3EXIT( 'AOICSUM', 0, 0, 
     &                  'Could not get description of file "' 
     &                  // ENAME( 1:TRIMLEN( ENAME ) ) // '"', 2 )
        END IF

        NSRCS = NROWS3D
        IF ( NSRCS .GT. NASRC ) THEN
            WRITE ( BUF1, 600 ) NSRCS
            WRITE ( BUF2, 600 ) NASRC
            CALL M3EXIT( 'AOICSUM', 0, 0, 
     &                   'Number of sources (NROWS3D='  // 
     &                   BUF1 ( LBLANK(BUF1):TRIMLEN(BUF1) ) //
     &                   ') exceeds program limits (NASRC=' //
     &                   BUF2 ( LBLANK(BUF2):TRIMLEN(BUF2) ) // 
     &                   ')', 2 )
        END IF

        FUNIT = PROMPTFFILE( 
     &          'Enter name for STATE/COUNTY FIPS CODE file',
     &          .TRUE., .TRUE., 'STCYNAM', 'AOICSUM' )

        AUNIT = PROMPTFFILE( 
     &          'Enter name for OIC DESCRIPTION file',
     &          .TRUE., .TRUE., 'ASCTNAM', 'AOICSUM' )

        LUNIT = PROMPTFFILE( 
     &  'Enter name for LIST OF OIC CODES to limit reports (or "NONE")',
     &          .TRUE., .TRUE., 'OICLIMIT', 'AOICSUM' )

        SUNIT = PROMPTFFILE( 
     &          'Enter name for AREA ASCT CODES file',
     &          .TRUE., .TRUE., 'ASCC', 'AOICSUM' )

        COLISTSZ = ENVINT ( 'COLSTSZ', 'Max counties per OIC', MAXCO, 
     &                      STATUS )
        COLISTSZ = GETNUM ( 1, MAXCO, COLISTSZ, 
     &     'Enter maximum number of counties per OIC group in report' )

        STUNIT = PROMPTFFILE( 
     &          'Enter name for STATE EMISSIONS RANKING REPORT file',
     &          .FALSE., .TRUE., 'AOSTRPT', 'AOICSUM' )

        COUNIT = PROMPTFFILE( 
     &          'Enter name for COUNTY EMISSIONS RANKING REPORT file',
     &          .FALSE., .TRUE., 'AOCORPT', 'AOICSUM' )


        ! Read in all state and county fips codes and names.
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
                        CALL M3EXIT( 'AOICSUM', 0, 0, 
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
                        CALL M3EXIT( 'AOICSUM', 0, 0, 
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


        ! Read in all ASCT codes and names
        TOT_ASCT = 0
35      CONTINUE ! Begin "infinite" loop; until EOF.
            READ ( AUNIT, 710, END = 45 ) INASC, ASCTNAME
C            PRINT *, INASC // ': ' // ASCTNAME
            TOT_ASCT = TOT_ASCT + 1
            IF ( TOT_ASCT .GT. MAXASCT ) THEN
                WRITE ( BUF1, 600 ) TOT_ASCT
                WRITE ( BUF2, 600 ) MAXASCT
                CALL M3EXIT( 'AOICSUM', 0, 0, 
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
        GO TO 35
45      CONTINUE

        LASTOIC = 0
        TOT_OIC = 0

        ! Find the unique complement of OIC codes within the ASCC file
55      CONTINUE
            READ ( SUNIT, '(A)', END = 65 ) INASC

            READ ( INASC ( 1:4 ), '(I4)' ) INOIC
            IF ( INOIC .NE. LASTOIC ) THEN
                TOT_OIC = TOT_OIC + 1
                IF ( TOT_OIC .GT. MAXOIC ) THEN
                    WRITE ( BUF1, 600 ) TOT_OIC
                    WRITE ( BUF2, 600 ) MAXOIC
                    CALL M3EXIT( 'AOICSUM', 0, 0, 
     &                      'Number of OIC codes in ASCC (' // 
     &                      BUF1 ( LBLANK(BUF1):TRIMLEN(BUF1) ) //
     &                      ') exceeds program limits (MAXOIC=' //
     &                      BUF2 ( LBLANK(BUF2):TRIMLEN(BUF2) ) // 
     &                      ')', 2 )
                END IF

                LASTOIC = INOIC
                IOICLIST ( TOT_OIC ) = INOIC
            END IF
        GO TO 55
65      CONTINUE


        ! Set the flags that determine which OIC Codes we'll report on.
        IF ( LUNIT .LT. 0 ) THEN
            DO 200 I = 1, MAXOIC
                IN_RPT ( I ) = .TRUE.
200         END DO
        ELSE
210         CONTINUE
                READ ( LUNIT, '(I)', END = 220 ) INOIC
                INDEX = FIND1 ( INOIC, TOT_OIC, IOICLIST )
                IF ( INDEX .GT. 0 ) THEN
                    IN_RPT ( INDEX ) = .TRUE.
                ELSE
                    WRITE ( BUF1, 600 ) INOIC
                    CALL M3WARN ( 'AOICSUM', 0, 0, 
     &                   'Requested OIC code (' // BUF1 // 
     &                   ') not in file "ASCC"' )
                END IF
            GO TO 210
220         CONTINUE
        END IF



        IF ( .NOT. READ3 ( ENAME, 'FIP', ALLAYS3,0,0, IFIP ) ) THEN
            CALL M3EXIT( 'AOICSUM', 0, 0, 
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
            CALL M3EXIT( 'AOICSUM', 0, 0, 
     &                   'Error reading INVYR from file' // ENAME, 2)
        END IF

        DO 50 I = 1, NIPOL
C           PRINT *, 'Reading ' // EINAM (I) // '...'
            IF ( .NOT. READ3 ( ENAME, EINAM ( I ), ALLAYS3, 0, 0, 
     &                         SPECBUF ( 1, I ) ) ) THEN
                CALL M3EXIT( 'AOICSUM', 0, 0, 
     &                       'Error reading ' // EINAM ( I ) //
     &                       ' from file' // ENAME, 2 )
            END IF
            
            IF ( ( VOC_VAR .EQ. 0 ) .AND. 
     &           ( ( EINAM ( I ) .EQ. 'VOC' ) 
     &             .OR. EINAM ( I ) .EQ. 'ROG' ) ) THEN
                VOC_VAR = I
            END IF
            IF ( ( NOX_VAR .EQ. 0 ) .AND. 
     &           ( EINAM ( I ) .EQ. 'NOX' ) ) THEN
                NOX_VAR = I
            END IF
50      CONTINUE
C        PRINT '(A, I2, A, A )', 'VOC Variable is ', VOC_VAR, ' or ', 
C     &        EINAM ( VOC_VAR )

        LASTCO = 0
        LASTST = 0
        LASTOIC = 0
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
                IF ( ACT_CO .GT. NAFIP ) THEN
                    WRITE ( BUF1, 600 ) ACT_CO
                    WRITE ( BUF2, 600 ) NAFIP
                    CALL M3EXIT( 'AOICSUM', 0, 0, 
     &                  'Number of counties in AREA (' // 
     &                  BUF1 ( LBLANK(BUF1):TRIMLEN(BUF1) ) //
     &                  ') exceeds program limits (NAFIP=' //
     &                  BUF2 ( LBLANK(BUF2):TRIMLEN(BUF2) ) // 
     &                  ')', 2 )
                END IF

                ICOLIST ( ACT_CO ) = LASTCO
                
                ! If the current state code is not equal to the last state
                ! code, then bump the number of states and record the code.
                IF ( ( IFIP ( J ) / 1000 ) .NE. LASTST ) THEN
                    LASTST = IFIP ( J ) / 1000
                    ACT_ST = ACT_ST + 1
                    IF ( ACT_ST .GT. NASID ) THEN
                        WRITE ( BUF1, 600 ) ACT_ST
                        WRITE ( BUF2, 600 ) NASID
                        CALL M3EXIT( 'AOICSUM', 0, 0, 
     &                      'Number of states in AREA (' // 
     &                      BUF1 ( LBLANK(BUF1):TRIMLEN(BUF1) ) //
     &                      ') exceeds program limits (NASID=' //
     &                      BUF2 ( LBLANK(BUF2):TRIMLEN(BUF2) ) // 
     &                      ')', 2 )
                    END IF

                    ISTLIST ( ACT_ST ) = LASTST
                END IF
            END IF

            ! Check to see if we've moved to another OIC.  If so, find out 
            ! where it is in the list.
            CUROIC = IASC7 ( J ) / 1000
            IF ( CUROIC .NE. LASTOIC ) THEN
                LASTOIC = CUROIC
                INDEX = FIND1 ( CUROIC, TOT_OIC, IOICLIST )
            END IF

            HASSTSUM ( INDEX, ACT_ST ) = .TRUE.
            HASCOSUM ( INDEX, ACT_CO ) = .TRUE.

            DO 70 I = 1, NIPOL
                STSUM ( INDEX, ACT_ST, I ) = 
     &                  STSUM ( INDEX, ACT_ST, I ) + 
     &                  ( SPECBUF ( J, I ) * YR2DAY )
                COSUM ( INDEX, ACT_CO, I ) = 
     &                  COSUM ( INDEX, ACT_CO, I ) + 
     &                  ( SPECBUF ( J, I ) * YR2DAY )
70          CONTINUE
60      CONTINUE


        WRITE ( STUNIT, 810 ) 
     &     'Area source emissions (tons/day) by OIC group and State:'
        WRITE ( STUNIT, 810 ) 'Report includes:'
        WRITE ( STUNIT, 810 ) '  (1) States ranked by ' // 
     &     EINAM ( NOX_VAR )( 1:TRIMLEN ( EINAM ( NOX_VAR ) ) ) // 
     &     ' emissions totals'
        WRITE ( STUNIT, 810 ) '  (2) States ranked by ' // 
     &     EINAM ( VOC_VAR )( 1:TRIMLEN ( EINAM ( VOC_VAR ) ) ) // 
     &     ' emissions totals'

        WRITE ( COUNIT, 810 )
     &     'Area source emissions (tons/day) by OIC group and County:'
        WRITE ( COUNIT, 810 ) 'Report includes:'
        WRITE ( COUNIT, 890 ) '  (1) Top ', COLISTSZ, 
     &     ' counties ranked by ' // 
     &     EINAM ( NOX_VAR )( 1:TRIMLEN ( EINAM ( NOX_VAR ) ) ) // 
     &     ' emissions totals'
        WRITE ( COUNIT, 890 ) '  (2) Top ', COLISTSZ, 
     &     ' counties ranked by ' // 
     &     EINAM ( VOC_VAR )( 1:TRIMLEN ( EINAM ( VOC_VAR ) ) ) // 
     &     ' emissions totals'



        ! We're going to "loop" through each variable of interest.  We'll
        ! start with NOX_VAR, then at the end of the "DO 80..." loop, if
        ! the EMIS_VAR is NOX_VAR, we'll change it to VOC_VAR and loop again.
        EMIS_VAR = NOX_VAR
        RPT_NUM = "1"
        
75      CONTINUE

        WRITE ( STUNIT, 800 ) 'Report (' // RPT_NUM // 
     $     '): States ranked by ' // 
     &     EINAM ( EMIS_VAR )( 1:TRIMLEN ( EINAM ( EMIS_VAR ) ) ) // 
     &     ' emissions totals'
        WRITE ( COUNIT, 800 ) 'Report (' // RPT_NUM // 
     $     '): Counties ranked by ' // 
     &     EINAM ( EMIS_VAR )( 1:TRIMLEN ( EINAM ( EMIS_VAR ) ) ) // 
     &     ' emissions totals'


        ! Loop thru each OIC code.
        DO 80 K = 1, TOT_OIC
            ! Do we want a report on this OIC code?
            IF ( .NOT. IN_RPT ( K ) ) THEN
                GO TO 80
            END IF

            ! Find the description for this OIC code.
            INASC7 = IOICLIST ( K ) * 1000
            INASC3 = 0
            INDEX = FIND2 ( INASC7, INASC3, TOT_ASCT, 
     &                      IALLASC7, IALLASC3 )
            IF ( INDEX .LE. 0 ) THEN
                ASCTNAME = "*** no OIC name available"
            ELSE
                ASCTNAME = ASCTNAMES ( INDEX )
                IF ( ASCTNAME .EQ. ' ' ) THEN
                	ASCTNAME = "*** no OIC name available"
                END IF 
            END IF

            ! Write the STATE header
            WRITE ( STUNIT, 880 ) 
     &        '---- OIC = ',
     &        IOICLIST ( K ), ': ', ASCTNAME ( 1:TRIMLEN ( ASCTNAME ) ),
     &        ' -------------'
            WRITE ( STUNIT, 810 ) 'STATE'
            WRITE ( STUNIT, 840 ) 'CODE'
            DO 90 I = 1, NIPOL
                WRITE ( STUNIT, 820 ) 
     &                 EINAM ( I )( 1:TRIMLEN ( EINAM ( I ) ))
90          CONTINUE
            WRITE ( STUNIT, 810 ) '  STATE NAME'


            ! Write the COUNTY header
            WRITE ( COUNIT, 880 ) 
     &        '---- OIC = ',
     &        IOICLIST ( K ), ': ', ASCTNAME ( 1:TRIMLEN ( ASCTNAME ) ),
     &        ' -------------'
            WRITE ( COUNIT, 810 ) 'FIPs'
            WRITE ( COUNIT, 840 ) 'CODE '
            DO 100 I = 1, NIPOL
                WRITE ( COUNIT, 820 ) 
     &                 EINAM ( I )( 1:TRIMLEN ( EINAM ( I ) ))
100         CONTINUE
            WRITE ( COUNIT, 810 ) '  STATE NAME            COUNTY NAME'


            ! Order the states by increasing VOC emissions totals
            DO 110 J = 1, ACT_ST
                SORTINDEX ( J ) = J
                TMPEMIS ( J ) = STSUM ( K, J, EMIS_VAR )
110         CONTINUE
            CALL SORTR1 ( ACT_ST, SORTINDEX, TMPEMIS )

            ! Loop thru the states in reverse order.
            DO 120 J = ACT_ST, 1, -1
                POSN = SORTINDEX ( J )

                ! If this state never had a value for this OIC, don't
                ! bother to write it out.
                IF ( HASSTSUM ( K, POSN ) ) THEN
                    ! Lookup the state's name
                    INDEX = FIND1 ( ISTLIST ( POSN ), TOT_ST, ISTFIP )
                    IF ( INDEX .GT. 0 ) THEN
                    	FIPNAME = STATES ( INDEX )
                    ELSE
                        FIPNAME = '*** no state name'
                    ENDIF

                    ! Finally, write the data to the report.
                    WRITE ( STUNIT, 830 ) 
     &                  ISTLIST ( POSN ), 
     &                  ( STSUM ( K, POSN, I ), I = 1, NIPOL ),
     &                  FIPNAME
                END IF
120         CONTINUE


            ! Order the COUNTIES by increasing VOC emissions totals
            DO 130 J = 1, ACT_CO
                SORTINDEX ( J ) = J
                TMPEMIS ( J ) = COSUM ( K, J, EMIS_VAR )
130         CONTINUE
            CALL SORTR1 ( ACT_CO, SORTINDEX, TMPEMIS )

            ! Loop thru the COUNTIES in reverse order.
            COWRITTEN = 0
            DO 140 J = ACT_CO, 1, -1
                POSN = SORTINDEX ( J )

                ! If this state never had a value for this OIC, don't
                ! bother to write it out.
                IF ( HASCOSUM ( K, POSN ) ) THEN
                    ! Lookup the state's name
                    INDEX = FIND1 ( ICOLIST ( POSN ) / 1000, 
     &                              TOT_ST, ISTFIP )
                    IF ( INDEX .GT. 0 ) THEN
                    	FIPNAME = STATES ( INDEX )
                    ELSE
                        FIPNAME = '*** no state name'
                    ENDIF

                    ! Lookup the county's name
                    INDEX = FIND1 ( ICOLIST ( POSN ), TOT_CO, ICOFIP )
                    IF ( INDEX .GT. 0 ) THEN
                    	FIPNAME2 = COUNTIES ( INDEX )
                    ELSE
                        FIPNAME2 = '*** no county name'
                    ENDIF

                    ! Finally, write the data to the report.
                    WRITE ( COUNIT, 835 ) 
     &                  ICOLIST ( POSN ), 
     &                  ( COSUM ( K, POSN, I ), I = 1, NIPOL ),
     &                  FIPNAME, FIPNAME2

                    ! Bump the count of counties written, and check to see
                    ! if we've exceeded the maximum list size set by the user.
                    COWRITTEN = COWRITTEN + 1
                    IF ( COWRITTEN .GE. COLISTSZ ) THEN
                        GO TO 150
                    END IF
                END IF
140         CONTINUE
150         CONTINUE      

80      CONTINUE  ! end of DO for each OIC

        ! End of "loop" for each emissions variable of interest.
        IF ( EMIS_VAR .EQ. NOX_VAR ) THEN
            ! Set up for the next report.
            EMIS_VAR = VOC_VAR
            RPT_NUM = "2"
            GO TO 75 ! go back to the beginning of the "loop"
        END IF

999     CONTINUE          !  exit program

        MESG = 'Normal completion of AOICSUM'
        CALL M3EXIT( 'AOICSUM', 0, 0, MESG, 0 )


C******************  FORMAT  STATEMENTS   ******************************

600     FORMAT ( I10 )
700     FORMAT ( 2X, I5, 1X, A )
710     FORMAT ( A10, 1X, A40 )
800     FORMAT ( //, A, / )
810     FORMAT ( A )
820     FORMAT ( A10, $ )
830     FORMAT ( I2.2, 2X, <NIPOL>F10.2, 2X, A )
835     FORMAT ( I5.5, <NIPOL>F10.2, 2X, A, 2X, A )
840     FORMAT ( A, $ )
850     FORMAT ( //, A, A )
860     FORMAT ( /, 19X, A )
870     FORMAT ( I5.5, <NIPOL>F10.2, 2X, A )
880     FORMAT ( /, A, I4.4, A, A, A, / )
890     FORMAT ( A, I4, A )
        END
