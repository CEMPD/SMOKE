
        CHARACTER(LEN=80) FUNCTION WRSPDVMT( FREESPD, ARTSPD, SPDDIR )

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

C.........  This module is used for MOBILE6 setup information 
        USE MODMBSET
      
        IMPLICIT NONE
        
C...........   INCLUDES
        INCLUDE 'M6CNST3.EXT'  !  MOBILE6 constants
        
C...........   EXTERNAL FUNCTIONS and their descriptions:
        INTEGER        JUNIT
        CHARACTER*2    CRLF    
        
        EXTERNAL  JUNIT, CRLF

C...........   FUNCTION ARGUMENTS
        REAL,             INTENT (IN) :: FREESPD ! freeway speed
        REAL,             INTENT (IN) :: ARTSPD  ! arterial speed
        CHARACTER(LEN=*), INTENT (IN) :: SPDDIR  ! directory for writing speed files

C...........   Local arrays
        REAL FREEPROF( 14 )                ! speed profile for freeway
        REAL ARTPROF ( 14 )                !   "       "    "  arterial

C...........   Other local variables
        INTEGER I                             ! counters and indices
        INTEGER IOS                           ! I/O status
        INTEGER SDEV                          ! unit no. for speed file

        LOGICAL :: FILEEXIST = .FALSE.        ! true: file already exists

        CHARACTER(LEN=6)   FREESPDSTR         ! freeway speed as string
        CHARACTER(LEN=6)   ARTSPDSTR          ! arterial speed as string
        CHARACTER(LEN=20)  FILENAME           ! filename of speed file
        CHARACTER(LEN=300) MESG               ! message buffer

        CHARACTER(LEN=16) :: PROGNAME = 'WRSPDVMT'   ! program name
        
C***********************************************************************
C   begin body of function WRSPDVMT

C.........  Create speed profiles based on average values
        FREEPROF = 0.
        ARTPROF  = 0.
        
        CALL CALC_SPEED_PROF( FREESPD, FREEPROF )
        CALL CALC_SPEED_PROF( ARTSPD,  ARTPROF )

C.........  Create name of speed vmt file
        WRITE( FREESPDSTR, '(F6.2)' ) FREESPD
        FREESPDSTR = ADJUSTL( FREESPDSTR )
        
        WRITE( ARTSPDSTR,  '(F6.2)' ) ARTSPD
        ARTSPDSTR  = ADJUSTL( ARTSPDSTR )

        FILENAME = 'spd' // FREESPDSTR( 1:LEN_TRIM( FREESPDSTR ) ) 
     &             // '_' // ARTSPDSTR( 1:LEN_TRIM( ARTSPDSTR ) ) 
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
        
C.........  Write arrays to file
        DO I = 1, 24
            WRITE( SDEV, 93010 ) M6FREEWAY, I, FREEPROF
        END DO
        	
        DO I = 1, 24
            WRITE( SDEV, 93010 ) M6ARTERIAL, I, ARTPROF
        END DO

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
        
        END FUNCTION WRSPDVMT
