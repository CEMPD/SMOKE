
        SUBROUTINE AVERTEMP( NSRC, NCNTY, CNTYCODE, SRCARRAY,  
     &                       TSTEP, CNTYTEMP, CNTYQV, CNTYBP, NDAYSRC ) 

C***********************************************************************
C  subroutine body starts at line 78
C
C  DESCRIPTION:
C       Averages hourly meteorology data based on number of sources
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

C...........   MODULES for public variables

C...........   This module is the derived meteorology data for emission factors
        USE MODMET, ONLY: TKHOUR, QVHOUR, BPHOUR
        
        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters

C...........   EXTERNAL FUNCTIONS 
        CHARACTER*2  CRLF
        INTEGER      FIND1FIRST

        EXTERNAL     CRLF, FIND1FIRST

C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT   (IN) :: NSRC                   ! no. sources
        INTEGER, INTENT   (IN) :: NCNTY                  ! no. counties
        INTEGER, INTENT   (IN) :: CNTYCODE( NCNTY )      ! FIPS codes for counties
        INTEGER, INTENT   (IN) :: SRCARRAY( NSRC )       ! county codes for each source
        INTEGER, INTENT   (IN) :: TSTEP                  ! current time step
        REAL,    INTENT  (OUT) :: CNTYTEMP( NCNTY )      ! averaged temps by county
        REAL,    INTENT  (OUT) :: CNTYQV  ( NCNTY )      ! averaged mix. ratio by county
        REAL,    INTENT  (OUT) :: CNTYBP  ( NCNTY )      ! averaged baro. pressure by county
        INTEGER, INTENT(INOUT) :: NDAYSRC( NSRC,24 )     ! no. days to average over

C...........   Other local variables
        INTEGER I, J, K, L                ! counters and indices                     
        
        INTEGER IOS                       ! I/O status
        INTEGER CURRCNTY                  ! current county code
        INTEGER NUMSRC                    ! no. sources to be averaged

        REAL    TEMPSUM                   ! sum of temperatures
        REAL    QVSUM                     ! sum of mixing ratios
        REAL    BPSUM                     ! sum of barometric pressures

        LOGICAL :: INITIAL = .TRUE.       ! true: first time through routine

        CHARACTER(LEN=300) MESG           ! message buffer

        CHARACTER*16 :: PROGNAME = 'AVERTEMP' ! program name

C***********************************************************************
C   begin body of subroutine AVERTEMP
        
C.........  Loop through all counties
        DO I = 1, NCNTY
            NUMSRC = 0
            
            TEMPSUM = 0
            QVSUM   = 0
            BPSUM   = 0
        
            CURRCNTY = CNTYCODE( I )
            
            DO J = 1, NSRC
                IF( SRCARRAY( J ) /= CURRCNTY ) CYCLE

C.................  Skip sources with no days; this can happen when the
C                   gridding surrogates do not contain data for all counties
                IF( NDAYSRC( J,TSTEP ) == 0 ) CYCLE
                
                NUMSRC = NUMSRC + 1
                
                TEMPSUM = TEMPSUM + 
     &                    ( TKHOUR( J,TSTEP ) / NDAYSRC( J,TSTEP ) )
                QVSUM = QVSUM +
     &                    ( QVHOUR( J,TSTEP ) / NDAYSRC( J,TSTEP ) )
                BPSUM = BPSUM +
     &                    ( BPHOUR( J,TSTEP ) / NDAYSRC( J,TSTEP ) )
     
                TKHOUR( J,TSTEP ) = 0
                QVHOUR( J,TSTEP ) = 0
                BPHOUR( J,TSTEP ) = 0
                NDAYSRC( J,TSTEP ) = 0

            END DO
            
            IF( NUMSRC == 0 ) THEN
                WRITE( MESG,94010 ) 
     &                 'No valid meteorology data for reference county',
     &                 CURRCNTY, '. This is probably due to ',
     &                 'insufficient surrogate data.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            ELSE
                CNTYTEMP( I ) = TEMPSUM / NUMSRC
                CNTYQV  ( I ) = QVSUM   / NUMSRC
                CNTYBP  ( I ) = BPSUM   / NUMSRC
            END IF
        END DO

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I6, :, 1X ) )

94020   FORMAT( A, 4( 1X, F8.2, 1X, A ) )
 
        END SUBROUTINE AVERTEMP
