
        SUBROUTINE RDHOURTEMP( TNAME, NCOUNTY, JDATE, JTIME )

C***********************************************************************
C  subroutine body starts at line 74
C
C  DESCRIPTION:
C       Reads hourly temperatures, relative humidity values, and barometric
C       pressures from county file. Converts temperatures from K to F;
C       calculates average daily barometric pressure and converts from Pa
C       to inHG
C
C  PRECONDITIONS REQUIRED:
C       TKHOUR, RDHOUR, and BPHOUR must be allocated
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

C...........   This module is the derived meteorology data for emission factors
        USE MODMET, ONLY: TKHOUR, QVHOUR, BPHOUR, BPDAY, RHHOUR
        
        IMPLICIT NONE

C...........   INCLUDES:
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'CONST3.EXT'    !  physical and mathematical constants
                
C...........   EXTERNAL FUNCTIONS and their descriptions:
        REAL          CALCRELHUM
        CHARACTER(2)  CRLF    
        
        EXTERNAL  CALCRELHUM, CRLF

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(16), INTENT (IN) :: TNAME    ! logical name for meteorology file
        INTEGER,       INTENT (IN) :: NCOUNTY  ! no. counties in met file
        INTEGER,       INTENT (IN) :: JDATE    ! starting date (YYYYDDD)
        INTEGER,       INTENT (IN) :: JTIME    ! starting time (HHMMSS)

C...........   Local allocatable arrays

C...........   Other local variables
        INTEGER I, J, K                   ! counters and indices                     
        
        INTEGER IOS                       ! I/O status

        LOGICAL :: EFLAG      = .FALSE.   ! true: error found
        
        CHARACTER(300)     MESG     !  message buffer

        CHARACTER(16) :: PROGNAME = 'RDHOURTEMP'   ! program name
        
C***********************************************************************
C   begin body of subroutine RDHOURTEMP
        
C.........  Initialize data arrays
        TKHOUR = 0.
        QVHOUR = 0.
        BPHOUR = 0.
        BPDAY  = 0.
        RHHOUR = 0.

C.........  Loop through the time steps        
        DO I = 1, 24

C.............  Read temperature values        
            IF( .NOT. READ3( TNAME, 'TKCOUNTY', 1, JDATE, JTIME, 
     &                       TKHOUR( :,I ) ) ) THEN
                MESG = 'Could not read TKCOUNTY' //
     &                 ' from ' // TNAME 
                CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
            END IF
                   
C.............  Read relative humidity values
            IF( .NOT. READ3( TNAME, 'QVCOUNTY', 1, JDATE, JTIME,
     &                       QVHOUR( :,I ) ) ) THEN
                MESG = 'Could not read QVCOUNTY' //
     &                 ' from ' // TNAME
                CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
            END IF
            
C.............  Read barometric pressure values
            IF( .NOT. READ3( TNAME, 'BPCOUNTY', 1, JDATE, JTIME,
     &                       BPHOUR( :,I ) ) ) THEN
                MESG = 'Could not read BPCOUNTY' //
     &                 ' from ' // TNAME
                CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
            END IF
                   
            CALL NEXTIME( JDATE, JTIME, 10000 )
        
        END DO

C.........  Calculate average daily barometric pressure values
        DO I = 1, 24        
            BPDAY = BPDAY + BPHOUR( :,I )
        END DO

        BPDAY = BPDAY/24.
        
C.........  Calculate relative humidity values
        DO I = 1, 24
            DO J = 1, SIZE( BPDAY )
                RHHOUR( J,I ) = CALCRELHUM( TKHOUR( J,I ), BPDAY( J ), 
     &                                      QVHOUR( J,I ) )
            END DO
        END DO
        
C.........  Convert temps from Kelvin to Fahrenheit
        TKHOUR = CTOF * ( TKHOUR - CTOK ) + 32.

C.........  Convert pressure from Pa to inHG
        BPDAY = BPDAY * PA2INHG

        END SUBROUTINE RDHOURTEMP
        
