
        CHARACTER(LEN=280) FUNCTION WRSPDVMT( FREESPD, ARTSPD, SPDDIR, 
     &                                        SPDFLAG, ADJINV, ADJHOUR )

C***********************************************************************
C  function body starts at line 77
C
C  DESCRIPTION:
C       Writes the SPEED VMT output file
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:  none
C
C  REVISION  HISTORY:
C     10/01: Created by C. Seppanen
C
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2002, MCNC Environmental Modeling Center
C All Rights Reserved
C
C See file COPYRIGHT for conditions of use.
C
C Environmental Modeling Center
C MCNC
C P.O. Box 12889
C Research Triangle Park, NC  27709-2889
C
C smoke@emc.mcnc.org
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C***********************************************************************

C.........  MODULES for public variables
C...........  This module is for mobile-specific data
        USE MODMOBIL, ONLY: NSPDPROF, SPDNUMS, SPDPROFS

        IMPLICIT NONE
        
C...........   INCLUDES
        INCLUDE 'M6CNST3.EXT'  !  MOBILE6 constants
        
C...........   EXTERNAL FUNCTIONS and their descriptions:
        INTEGER        JUNIT
        CHARACTER*2    CRLF
        INTEGER        FIND1
        
        EXTERNAL  JUNIT, CRLF, FIND1

C...........   FUNCTION ARGUMENTS
        REAL,             INTENT (IN) :: FREESPD ! freeway speed
        REAL,             INTENT (IN) :: ARTSPD  ! arterial speed
        CHARACTER(LEN=*), INTENT (IN) :: SPDDIR  ! directory for writing speed files
        LOGICAL,          INTENT (IN) :: SPDFLAG ! true: speed profiles are available
        LOGICAL,          INTENT (IN) :: ADJINV  ! true: adjust inventory speeds for ramps
        LOGICAL,          INTENT (IN) :: ADJHOUR ! true: adjust hourly speeds for ramps

C...........   Local arrays
        REAL M6SPDPROF( 14 )                  ! M6 speed profile for a single hour

C...........   Other local variables
        INTEGER I,K                           ! counters and indices
        INTEGER IOS                           ! I/O status
        INTEGER M6HOUR                        ! hour adjusted for M6 (1 = 6 AM)
        INTEGER SDEV                          ! unit no. for speed file
        INTEGER PROFNUM                       ! speed profile number

        REAL    ADJFREESPD                    ! adjusted freeway speed

        LOGICAL :: FILEEXIST = .FALSE.        ! true: file already exists

        CHARACTER(LEN=6)   FREESPDSTR         ! freeway speed as string
        CHARACTER(LEN=6)   ARTSPDSTR          ! arterial speed as string
        CHARACTER(LEN=80)  FILENAME           ! filename of speed file
        CHARACTER(LEN=300) MESG               ! message buffer

        CHARACTER(LEN=16) :: PROGNAME = 'WRSPDVMT'   ! program name
        
C***********************************************************************
C   begin body of function WRSPDVMT

C.........  Make sure speed profiles are allowed
        IF( FREESPD < 0. .OR. ARTSPD < 0. ) THEN
            IF( .NOT. SPDFLAG ) THEN
                MESG = 'ERROR: Speed profiles were not read in, ' //
     &                 'but the SPDSUM file uses profiles'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF
        END IF

C.........  Create name of speed vmt file - if using values, put speed in name, 
C           otherwise, use 'pX' where X is the profile number
        IF( FREESPD >= 0. ) THEN
            WRITE( FREESPDSTR, '(F6.2)' ) FREESPD
        ELSE
            WRITE( FREESPDSTR, '(I6)' ) -INT( FREESPD )
            FREESPDSTR = 'p' // ADJUSTL( FREESPDSTR )
        END IF
            
        FREESPDSTR = ADJUSTL( FREESPDSTR )
        
        IF( ARTSPD >= 0. ) THEN
            WRITE( ARTSPDSTR, '(F6.2)' ) ARTSPD
        ELSE
            WRITE( ARTSPDSTR, '(I6)' ) -INT( ARTSPD )
            ARTSPDSTR = 'p' // ADJUSTL( ARTSPDSTR )
        END IF
            
        ARTSPDSTR  = ADJUSTL( ARTSPDSTR )

        IF( ( ADJINV .AND. FREESPD >= 0. ) .OR.
     &      ( ADJHOUR .AND. FREESPD < 0. ) ) THEN
            FILENAME = 'spd_adj'
        ELSE
            FILENAME = 'spd'
        END IF
        
        FILENAME = TRIM( FILENAME ) // TRIM( FREESPDSTR ) 
     &             // '_' // TRIM( ARTSPDSTR ) 
     &             // '.sv'

        WRSPDVMT = SPDDIR( 1:LEN_TRIM( SPDDIR ) ) // '/' 
     &             // FILENAME( 1:LEN_TRIM( FILENAME ) )

C.........  Check if the file already exists
        INQUIRE( FILE=WRSPDVMT, EXIST=FILEEXIST )
        
        IF( FILEEXIST ) THEN
            RETURN
        END IF
        
C.........  Otherwise, we need to create the file
        SDEV = JUNIT()

        OPEN( SDEV, FILE=WRSPDVMT, STATUS='NEW', IOSTAT=IOS )
            
        IF( IOS /= 0 ) THEN
            WRITE( MESG, 93000 ) 'ERROR: Could not create '
     &              // 'speed vmt file for freeway speed '
     &              // FREESPDSTR // ' and arterial speed '
     &              // ARTSPDSTR // '.'
            CALL M3MESG( MESG )
            
            WRSPDVMT = ''
            RETURN
        END IF

C.........  Write header line to file
        WRITE( SDEV, 93000 ) 'SPEED VMT'

C.........  If freeway speed is actual value, create single speed profile
        IF( FREESPD >= 0. ) THEN
        
C.............  If needed, adjust for ramps
            IF( ADJINV ) THEN
                ADJFREESPD = ADJUST_FREEWAY_SPD( FREESPD )
            ELSE
                ADJFREESPD = FREESPD
            END IF
                            
            CALL CALC_SPEED_PROF( ADJFREESPD, M6SPDPROF )

C.............  Write 24 copies of profile to file            
            DO I = 1, 24
                WRITE( SDEV, 93010 ) M6FREEWAY, I, M6SPDPROF
            END DO

C.........  Otherwise, match profile number and create separate M6 profiles for each hour
        ELSE
            PROFNUM = -FREESPD
            
            K = FIND1( PROFNUM, NSPDPROF, SPDNUMS )
            IF( K < 1 ) THEN
                WRITE( MESG,94010) 'ERROR: Could not find speed ' //
     &                 'profile', PROFNUM, 'in speed profiles file.'
                CALL M3MESG( MESG )
                WRSPDVMT = ''
                RETURN
            END IF
            
            DO I = 1,24
                M6HOUR = I+6    ! adjust for 1 = 6 AM
                IF( M6HOUR > 24 ) M6HOUR = M6HOUR - 24

C.................  If needed, adjust for ramps                
                IF( ADJHOUR ) THEN
                    ADJFREESPD = 
     &                  ADJUST_FREEWAY_SPD( SPDPROFS( K, M6HOUR ) )
                ELSE
                    ADJFREESPD = SPDPROFS( K, M6HOUR )
                END IF
                
                CALL CALC_SPEED_PROF( ADJFREESPD, M6SPDPROF )
                WRITE( SDEV, 93010 ) M6FREEWAY, I, M6SPDPROF
            END DO
        END IF
            
C.........  Repeat process for arterial speed
        IF( ARTSPD >= 0. ) THEN            
            CALL CALC_SPEED_PROF( ARTSPD,  M6SPDPROF )
         
            DO I = 1, 24
                WRITE( SDEV, 93010 ) M6ARTERIAL, I, M6SPDPROF
            END DO
        
        ELSE
            PROFNUM = -ARTSPD
            
            K = FIND1( PROFNUM, NSPDPROF, SPDNUMS )
            IF( K < 1 ) THEN
                WRITE( MESG,94010) 'ERROR: Could not find speed ' //
     &                 'profile', PROFNUM, 'in speed profiles file.'
                CALL M3MESG( MESG )
                WRSPDVMT = ''
                RETURN
            END IF
            
            DO I = 1,24
                M6HOUR = I+6    ! adjust for 1 = 6 AM 
                IF( M6HOUR > 24 ) M6HOUR = M6HOUR - 24
                CALL CALC_SPEED_PROF( SPDPROFS( K,M6HOUR ), M6SPDPROF )
                WRITE( SDEV, 93010 ) M6ARTERIAL, I, M6SPDPROF
            END DO
        END IF

C.........  Close the file        
        CLOSE( SDEV )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )
93010   FORMAT( I1, 1X, I2, 1X, 14( F6.4, 1X ) )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )
94020   FORMAT( 3( A, :, F6.2 ) )

C******************  INTERNAL SUBPROGRAMS  *****************************

        CONTAINS

C.............  This internal subroutine calculates a MOBILE6 speed profile
C               (uses 14 bins) based on a single average speed.

            SUBROUTINE CALC_SPEED_PROF( SPEED, SPDPROF )
        
C.................  Subprogram arguments
            REAL, INTENT (IN)    :: SPEED          ! speed to use
            REAL, INTENT (INOUT) :: SPDPROF( 14 )  ! speed profile

C.............  Local subprogram variables
            INTEGER I                             ! counter
            INTEGER UPPERBIN                      ! upper speed bin
            INTEGER LOWERBIN                      ! lower speed bin        
            
            REAL    UPPERSPD                      ! speed corresponding to UPPERBIN
            REAL    LOWERSPD                      ! speed corresponding to LOWERBIN            
            REAL    FRACTION                      ! fraction in LOWERBIN        

C.............................................................................

C.............  Initialize profile array
            SPDPROF = 0.

C.............  Set upper and lower speeds bracketing the actual speed
            IF( SPEED >= 65 ) THEN
                FRACTION = 0.
                UPPERBIN = 14
            ELSE IF( SPEED < 5 ) THEN
C.. NEED TO CHECK HOW MOBILE6 HANDLES THIS CASE
                FRACTION = 1.
                UPPERBIN = 2
            ELSE
                LOWERSPD = FLOOR( SPEED / 5 ) * 5
        
                UPPERSPD = LOWERSPD + 5 
        
C.................  Calculate fraction in lower bin
                FRACTION = ( ( 1 / SPEED    ) - ( 1 / UPPERSPD ) ) /
     &                     ( ( 1 / LOWERSPD ) - ( 1 / UPPERSPD ) )

C.................  Determine speed bin number
                UPPERBIN = ( UPPERSPD / 5 ) + 1
            END IF
        
            LOWERBIN = UPPERBIN - 1
        
C.............  Store fractions in appropriate speed bin        
            DO I = 1,14
                IF( I == LOWERBIN ) THEN
                    SPDPROF( I ) = FRACTION
                ELSE IF( I == UPPERBIN ) THEN
                    SPDPROF( I ) = 1 - FRACTION
                END IF
            END DO

            END SUBROUTINE CALC_SPEED_PROF       
        
C...............................................................................

C.............  This internal function adjusts a freeway speed for the effects
C               of travel on freeway ramps.

            REAL FUNCTION ADJUST_FREEWAY_SPD( SPEED )
            
C.................  Function arguments
            REAL, INTENT (IN) :: SPEED  ! freeway speed
            
C.............  Local function variables
            REAL DENOM           ! denominator of calculation
            
            DENOM = ( 1 / SPEED ) - ( RAMPVMT / RAMPSPEED )
            
            ADJUST_FREEWAY_SPD = ( 1 - RAMPVMT ) / DENOM
            
            END FUNCTION ADJUST_FREEWAY_SPD
            
        END FUNCTION WRSPDVMT
