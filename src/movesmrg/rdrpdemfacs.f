
        SUBROUTINE RDRPDEMFACS( REFIDX, MONTH )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C       Reads the emission factor for the given county and month
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:  none
C
C  REVISION  HISTORY:
C     04/10: Created by C. Seppanen
C     04/11: Modified by B.H. Baek
C
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2004, Environmental Modeling for Policy Development
C All Rights Reserved
C 
C Carolina Environmental Program
C University of North Carolina at Chapel Hill
C 137 E. Franklin St., CB# 6116
C Chapel Hill, NC 27599-6116
C 
C smoke@unc.edu
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C***********************************************************************

C.........  MODULES for public variables
C.........  This module is used for reference county information
        USE MODMBSET, ONLY: MCREFIDX

C.........  This module contains data structures and flags specific to Movesmrg
        USE MODMVSMRG, ONLY: MRCLIST, MVFILDIR, EMPOLIDX,
     &                       NEMTEMPS, EMTEMPS, RPDEMFACS,
     &                       NMVSPOLS, MVSPOLNAMS

C.........  This module contains the major data structure and control flags
        USE MODMERGE, ONLY: NSMATV, TSVDESC, NMSPC, NIPPA, EMNAM, EANAM

C.........  This module contains the lists of unique source characteristics
        USE MODLISTS, ONLY: NINVSCC, INVSCC

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'MVSCNST3.EXT'  !  MOVES contants

C...........   EXTERNAL FUNCTIONS and their descriptions:
        LOGICAL       BLKORCMT
        LOGICAL       CHKINT
        LOGICAL       CHKREAL
        INTEGER       FINDC
        INTEGER       INDEX1
        INTEGER       STR2INT
        REAL          STR2REAL
        CHARACTER(2)  CRLF

        EXTERNAL BLKORCMT, CHKINT, CHKREAL, FINDC,
     &           INDEX1, STR2INT, STR2REAL, CRLF

C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT(IN) :: REFIDX       ! ref. county index
        INTEGER, INTENT(IN) :: MONTH        ! current processing month

C...........   Local parameters
        INTEGER, PARAMETER :: NSEG = 200    ! number of segments
        INTEGER, PARAMETER :: NNONPOL = 8   ! number of non-pollutant fields in file

C...........   Local allocatable arrays
        CHARACTER(50),  ALLOCATABLE :: SEGMENT( : )    ! parsed input line
        LOGICAL,        ALLOCATABLE :: LINVSCC( : )    ! check inv SCC availability in lookup table

C...........   Other local variables
        INTEGER     I, J, L, LJ, L1, N, P, V  ! counters and indexes
        INTEGER     IOS         ! error status
        INTEGER  :: IREC = 0    ! record counter
        INTEGER  :: TDEV = 0    ! tmp. file unit
        INTEGER     SPDBIN      ! speed bin
        INTEGER     SCCIDX
        INTEGER     TMPIDX
        INTEGER     NSCC        ! no of processing SCCs
        INTEGER  :: TOGIDX = 0  ! index of TOG pollutant
        INTEGER  :: NHTOGIDX = 0  ! index of NONHAPTOG pollutant

        REAL        TMPVAL      ! temperature value
        REAL        PTMP        ! previous temperature value
        REAL        EMVAL       ! emission factor value

        LOGICAL     FOUND       ! true: header record was found
        LOGICAL     SKIPSCC     ! true: current SCC is not in inventory
        LOGICAL  :: EFLAG = .FALSE.
 
        CHARACTER(PLSLEN3)  SVBUF     ! tmp speciation name buffer
        CHARACTER(IOVLEN3)  CPOL      ! tmp pollutant buffer
        CHARACTER(IOVLEN3)  CSPC      ! tmp species buffer
        CHARACTER(SCCLEN3)  TSCC      ! current SCC
        CHARACTER(SCCLEN3)  PSCC      ! previous SCC
        CHARACTER(FIPLEN3)  CFIP      ! current ref county 
        CHARACTER(10000)    LINE          ! line buffer
        CHARACTER(100)      FILENAME      ! tmp. filename
        CHARACTER(200)      FULLFILE      ! tmp. filename with path
        CHARACTER(300)      MESG          ! message buffer

        CHARACTER(16) :: PROGNAME = 'RDRPDEMFACS'    ! program name

C***********************************************************************
C   begin body of subroutine RDRPDEMFACS

C.........  Open emission factors file based on MRCLIST file
        FILENAME = TRIM( MRCLIST( REFIDX, MONTH ) )
        
        IF( FILENAME .EQ. ' ' ) THEN
            WRITE( MESG, 94010 ) 'ERROR: No emission factors file ' //
     &        'for reference county ' //  MCREFIDX( REFIDX,1 ) // 
     &        ' and fuel month', MONTH
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        FULLFILE = TRIM( MVFILDIR ) // FILENAME
        OPEN( TDEV, FILE=FULLFILE, STATUS='OLD', IOSTAT=IOS )
        IF( IOS .NE. 0 ) THEN
            MESG = 'ERROR: Could not open emission factors file ' //
     &        FILENAME
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ELSE
            MESG = 'Reading emission factors file ' //
     &        FILENAME
            CALL M3MESG( MESG )
        END IF

C.........  Allocate memory to parse lines
        ALLOCATE( SEGMENT( NSEG ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SEGMENT', PROGNAME )

C.........  Read header line to get list of pollutants in file
        FOUND = .FALSE.
        IREC = 0
        NEMTEMPS = 0
        DO
        
            READ( TDEV, 93000, END=100, IOSTAT=IOS ) LINE
            
            IREC = IREC + 1
            
            IF( IOS .NE. 0 ) THEN
                WRITE( MESG, 94010 ) 'I/O error', IOS,
     &            'reading emission factors file ' //
     &            FILENAME // ' at line', IREC
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF
            
C.............  Check for header line
            IF( LINE( 1:12 ) .EQ. 'NUM_TEMP_BIN' ) THEN
                LJ = LEN_TRIM( LINE )
                NEMTEMPS = STR2INT( LINE( 13:LJ ) )

            ELSE IF( LINE( 1:15 ) .EQ. 'MOVESScenarioID' ) THEN
                FOUND = .TRUE.

                SEGMENT = ' '  ! array
                CALL PARSLINE( LINE, NSEG, SEGMENT )

C.................  Count number of pollutants
                NMVSPOLS = 0     ! Assume that both TOG/NONHAPTOG are existed
                DO J = NNONPOL + 1, NSEG 
                
                    IF( SEGMENT( J ) .NE. ' ' ) THEN
                        NMVSPOLS = NMVSPOLS + 1
                        IF( SEGMENT( J ) == 'NONHAPTOG' ) THEN
                            NHTOGIDX = INDEX1( 'NONHAPTOG', NIPPA, EANAM )
                        ELSE IF( SEGMENT( J ) == 'TOG' ) THEN
                            TOGIDX = INDEX1( 'TOG', NIPPA, EANAM )
                        END IF
                    ELSE
                        EXIT
                    END IF
                
                END DO
 
                IF( ALLOCATED( MVSPOLNAMS ) ) DEALLOCATE( MVSPOLNAMS )
                ALLOCATE( MVSPOLNAMS( NMVSPOLS ), STAT=IOS )
                CALL CHECKMEM( IOS, 'MVSPOLNAMS', PROGNAME )
                MVSPOLNAMS = ''

C.................  Store pollutant names                
                DO J = 1, NMVSPOLS
                    MVSPOLNAMS( J ) = SEGMENT( NNONPOL + J )
                END DO

                EXIT

            END IF

        END DO

100     CONTINUE

        REWIND( TDEV )

        IF( .NOT. FOUND ) THEN
            MESG = 'ERROR: Missing header line in emission factors file'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        IF( NMVSPOLS == 0 ) THEN
            MESG = 'ERROR: Emission factors file does not contain ' //
     &             'any pollutants'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        IF( TOGIDX == 0 .AND. NHTOGIDX == 0 ) THEN
            MESG = 'ERROR: Emission factors file does not contain ' //
     &             'data for both TOG/NONHAPTOG pollutants'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF
        
C.........  Build pollutant mapping table
C.........  Find emission pollutant in list of pollutants
        DO V = 1, NIPPA 

            CPOL = TRIM( EANAM( V ) )

            J = INDEX1( CPOL, NMVSPOLS, MVSPOLNAMS )
            IF( J .LE. 0 ) THEN
                MESG = 'WARNING: Emission factors file does not ' //
     &            'contain requested inventory pollutant '//TRIM( CPOL )
                CALL M3MESG( MESG )
            ELSE
                EMPOLIDX( V ) = J

            END IF

        END DO

C.............  Find model species in list of pollutants
        DO V = 1, NMSPC

            CSPC = EMNAM( V )
            J = INDEX1( CSPC, NMVSPOLS, MVSPOLNAMS )
            IF( J .LE. 0 ) THEN
                MESG = 'WARNING: Emission factors file does not ' //
     &            'contain requested model species ' // TRIM( CSPC )
                CALL M3MESG( MESG )
            ELSE
                EMPOLIDX( NIPPA + V ) = J
            END IF

        END DO

C.........  Error message
        IF( EFLAG ) THEN
            MESG = 'Error occurred while reading emission factors file'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Allocate memory to parse lines
        DEALLOCATE( SEGMENT )
        ALLOCATE( SEGMENT( NNONPOL + NMVSPOLS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SEGMENT', PROGNAME )

C.........  Read through file to determine maximum number of temperatures

C.........  Assumptions:
C             File will contain data for all SCCs in the inventory and
C               all 16 speed bins.
C             Lines are sorted by:
C                 temperature
C                 SCC (matching sorting of INVSCC)
C                 speed bin
C             Each SCC will have data for same set of emission temperatures, and
C               pollutants.

C.........  Limitations:
C             Program doesn't know if emission factors file is missing values.

C.........  Expected columns:
C MOVESScenarioID,yearID,monthID,FIPS,SCCsmoke,smokeProcID,avgSpeedBinID,temperature,relHumidity,THC,CO ...

        IF( NEMTEMPS .EQ. 0 ) THEN
        
          IREC = 0
          NEMTEMPS = 0
          PTMP = -999
          DO
        
            READ( TDEV, 93000, END=200, IOSTAT=IOS ) LINE
            
            IREC = IREC + 1
            
            IF( IOS .NE. 0 ) THEN
                WRITE( MESG, 94010 ) 'I/O error', IOS,
     &            'reading emission factors file ' //
     &            FILENAME // ' at line', IREC
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  Skip blank or comment lines
            IF( BLKORCMT( LINE ) ) CYCLE

C.............  Skip header line
            IF( LINE( 1:12 ) .EQ. 'NUM_TEMP_BIN' ) CYCLE
            IF( LINE( 1:15 ) .EQ. 'MOVESScenarioID' ) CYCLE

C.............  Parse line into segments
            CALL PARSLINE( LINE, NNONPOL + NMVSPOLS, SEGMENT )

C.............  Check that county matches requested county
            IF( .NOT. CHKINT( SEGMENT( 4 ) ) ) THEN
                WRITE( MESG, 94010 ) 'ERROR: Bad reference county ' //
     &            'FIPS code at line', IREC, 'of emission factors ' //
     &            'file.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF
 
            CFIP = SEGMENT( 4 )
            CALL PADZERO( CFIP )
            
            IF( CFIP .NE. MCREFIDX( REFIDX,1 ) ) THEN
                WRITE( MESG, 94010 ) 'ERROR: Reference county ' //
     &            'at line', IREC, 'of emission factors file ' //
     &            'does not match county listed in MRCLIST file.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  Check that fuel month matches requested month
            IF( .NOT. CHKINT( SEGMENT( 3 ) ) ) THEN
                WRITE( MESG, 94010 ) 'ERROR: Bad fuel month ' //
     &            'at line', IREC, 'of emission factors file.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            IF( STR2INT( SEGMENT( 3 ) ) .NE. MONTH ) THEN
                WRITE( MESG, 94010 ) 'ERROR: Fuel month at line',
     &            IREC, 'of emission factors file does not match ' //
     &            'fuel month listed in MRCLIST file.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF
            
C.............  Check temperature value
            IF( .NOT. CHKREAL( SEGMENT( 7 ) ) ) THEN
                WRITE( MESG, 94010 ) 'ERROR: Bad temperature value ' //
     &            'at line', IREC, 'of emission factors file.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF
            
            TMPVAL = STR2REAL( SEGMENT( 7 ) )
            IF( TMPVAL .NE. PTMP ) THEN

C.................  Check that temperatures are sorted
                IF( TMPVAL .LT. PTMP ) THEN
                    WRITE( MESG, 94010 ) 'ERROR: Temperature value ' //
     &                'at line', IREC, 'of emission factors file is ' //
     &                'smaller than previous temperature.'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF
            
                NEMTEMPS = NEMTEMPS + 1
                PTMP = TMPVAL
            END IF

          END DO

200       CONTINUE
        
          REWIND( TDEV )
        END IF
 
C.........  Allocate memory to store emission factors
        IF( ALLOCATED( RPDEMFACS ) ) THEN
            DEALLOCATE( RPDEMFACS )
        END IF
        ALLOCATE( RPDEMFACS( NINVSCC, 16, NEMTEMPS, NMVSPOLS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'RPDEMFACS', PROGNAME )
        RPDEMFACS = 0.  ! array

C.........  Allocate memory to store temperature values
        IF( ALLOCATED( EMTEMPS ) ) THEN
            DEALLOCATE( EMTEMPS )
        END IF
        ALLOCATE( EMTEMPS( NEMTEMPS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMTEMPS', PROGNAME )
        ALLOCATE( LINVSCC( NINVSCC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LINVSCC', PROGNAME )
        EMTEMPS = 0.  ! array
        LINVSCC = .FALSE.

C.........  Read and store emission factors
        IREC = 0
        PSCC = ' '
        SCCIDX = 0
        PTMP = -999
        TMPIDX = 0
        NSCC = 0
        DO
            READ( TDEV, 93000, END=300, IOSTAT=IOS ) LINE
            
            IREC = IREC + 1
            
            IF( IOS .NE. 0 ) THEN
                WRITE( MESG, 94010 ) 'I/O error', IOS,
     &            'reading emission factors file ' //
     &            FILENAME // ' at line', IREC
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  Skip blank or comment lines
            IF( BLKORCMT( LINE ) ) CYCLE

C.............  Skip header line
            IF( LINE( 1:12 ) .EQ. 'NUM_TEMP_BIN' ) CYCLE
            IF( LINE( 1:15 ) .EQ. 'MOVESScenarioID' ) CYCLE

C.............  Parse line into segments
            CALL PARSLINE( LINE, NNONPOL + NMVSPOLS, SEGMENT )

C.............  Set SCC index for current line
            TSCC = TRIM( SEGMENT( 5 ) )
            CALL PADZERO( TSCC )
            IF( TSCC .NE. PSCC ) THEN
                SKIPSCC = .FALSE.
            
                SCCIDX = SCCIDX + 1
                IF( SCCIDX .GT. NINVSCC ) THEN
                    SCCIDX = 1
                END IF
                
C.................  Make sure the SCC at this line is what we expect
                IF( TSCC .EQ. INVSCC( SCCIDX ) ) THEN
                    LINVSCC( SCCIDX ) = .TRUE.

                ELSE
C.....................  If the SCC is in the inventory, it means that the 
C                       emissions factors are out of order and that's a
C                       problem; otherwise, for SCCs in the emission 
C                       factors file that aren't in the inventory, we
C                       can skip this line
                    J = FINDC( TSCC, NINVSCC, INVSCC )
                    IF( J .LE. 0 ) THEN

C.........................  Emission factor SCC is not in the inventory
                        SKIPSCC = .TRUE.
                        SCCIDX = SCCIDX - 1
                        IF( SCCIDX .LT. 1 ) THEN
                            SCCIDX = NINVSCC
                        END IF
                    ELSE
                        WRITE( MESG, 94010 ) 'Expected SCC ' // 
     &                    TRIM( INVSCC( SCCIDX ) ) // ' in emission ' //
     &                    'factors file at line', IREC
                        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                    END IF

                END IF
                
                PSCC = TSCC
            END IF
            
            IF( SKIPSCC ) CYCLE
            NSCC = NSCC + 1

C.............  Set speed bin for current line            
            SPDBIN = STR2INT( SEGMENT( 6 ) )

C.............  Set temperature index for current line            
            TMPVAL = STR2REAL( SEGMENT( 7 ) )
            IF( TMPVAL .NE. PTMP ) THEN
                TMPIDX = TMPIDX + 1
                IF( TMPIDX .GT. NEMTEMPS ) THEN
                    WRITE( MESG, 94010 ) 'Unexpected temperature value ' //
     &                'in emission factors file at line', IREC
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF
                
                PTMP = TMPVAL
            END IF
            
            EMTEMPS( TMPIDX ) = TMPVAL

C.............  Store emission factors for each pollutant            
            DO P = 1, NMVSPOLS

                EMVAL = STR2REAL( SEGMENT( NNONPOL + P ) )

                IF( EMVAL < 0 ) THEN    ! reset negative EF to zero
                    EMVAL = 0.0
                    WRITE( MESG, 94010 ) 'WARNING: Resetting negative ' //
     &                ' emission factor to zero at line', IREC
                    CALL M3MESG( MESG )
                END IF                
                      
                RPDEMFACS( SCCIDX, SPDBIN, TMPIDX, P ) = EMVAL
            END DO

        END DO

300     CONTINUE        

C.........  Error message when inventory SCC is missing in the lookup table
        DO I = 1, NINVSCC
            IF( .NOT. LINVSCC( I ) ) THEN
                MESG = 'WARNING: Following SCC "' // INVSCC( I ) //
     &              '" is missing in this emission factors file'
                CALL M3MESG( MESG )
            END IF
        END DO

        IF( NSCC < 1 ) THEN
            EFLAG = .TRUE.
            MESG = 'ERROR: No inventory SCC is found in emission factors file'
            CALL M3MSG2( MESG )
        END IF

        IF( EFLAG )  CALL M3EXIT( PROGNAME, 0, 0, '', 2 )

        CLOSE( TDEV )
        
        DEALLOCATE( SEGMENT, LINVSCC )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )  
      
C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )
        
        END SUBROUTINE RDRPDEMFACS
