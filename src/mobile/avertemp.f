
        SUBROUTINE AVERTEMP( NSRC, NCNTY, CNTYCODE, SRCARRAY,  
     &                       TSTEP, HOURTEMP, CNTYTEMP, NDAYSRC ) 

C***********************************************************************
C  subroutine body starts at line 78
C
C  DESCRIPTION:
C       Averages hourly temperatures based on number of sources
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
        REAL,    INTENT(INOUT) :: HOURTEMP( NSRC, 24 )   ! hourly temps by source
        REAL,    INTENT  (OUT) :: CNTYTEMP( NCNTY )      ! averaged temps by county
        INTEGER, INTENT(INOUT) :: NDAYSRC( NSRC,24 )     ! no. days to average over

C...........   Other local variables
        INTEGER I, J, K, L                ! counters and indices                     
        
        INTEGER IOS                       ! I/O status
        INTEGER CURRCNTY                  ! current county code
        INTEGER NUMSRC                    ! no. sources to be averaged

        REAL    TEMPSUM                   ! sum of temperatures

        LOGICAL :: INITIAL = .TRUE.       ! true: first time through routine

        CHARACTER*16 :: PROGNAME = 'AVERTEMP' ! program name

C***********************************************************************
C   begin body of subroutine AVERTEMP
        
C.........  Loop through all counties
        DO I = 1, NCNTY
            NUMSRC = 0
            TEMPSUM = 0
        
            CURRCNTY = CNTYCODE( I )
            
            DO J = 1, NSRC
                IF( SRCARRAY( J ) /= CURRCNTY ) CYCLE
                
                NUMSRC = NUMSRC + 1
                TEMPSUM = TEMPSUM + 
     &                    ( HOURTEMP( J,TSTEP ) / NDAYSRC( J,TSTEP ) )
                HOURTEMP( J,TSTEP ) = 0
                NDAYSRC( J,TSTEP ) = 0

            END DO
            	
            CNTYTEMP( I ) = TEMPSUM / NUMSRC
        END DO

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I9, :, 1X ) )

94020   FORMAT( A, 4( 1X, F8.2, 1X, A ) )
 
        END SUBROUTINE AVERTEMP
