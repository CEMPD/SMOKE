
        SUBROUTINE AVERTEMP( NSRC, NSTEPS, NCOUNTY, SRCARRAY, 
     &                       HOURTEMP, COUNTYTEMP ) 

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters

C...........   EXTERNAL FUNCTIONS 
        CHARACTER*2  CRLF

        EXTERNAL     CRLF

C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT    (IN) :: NSRC                ! no. sources
        INTEGER, INTENT    (IN) :: NSTEPS              ! no. time steps
        INTEGER, INTENT    (IN) :: NCOUNTY             ! no. counties
        INTEGER, INTENT    (IN) :: SRCARRAY( NSRC )    ! county codes for each source
        REAL,    INTENT    (IN) :: HOURTEMP( NSRC, 0:NSTEPS-1 ) ! hourly temps by source
        REAL,    INTENT   (OUT) :: COUNTYTEMP( NCOUNTY,0:NSTEPS-1 ) ! hourly temps by county

C...........   Local arrays
        INTEGER IDX( NSRC )               ! index into sources array
        INTEGER COUNTIES( NCOUNTY )       ! array of valid county codes

C...........   Other local variables
        INTEGER I, J, K                   ! counters and indices                     
        
        INTEGER PREVCTY                   ! previous county code
        INTEGER CURRCTY                   ! current county code
        INTEGER COUNTYNUM                 ! current county number
        INTEGER NUMSRC                    ! no. sources to be averaged

        REAL    TEMPSUM                   ! sum of temperatures

        CHARACTER*16 :: PROGNAME = 'AVERTEMP' ! program name

C***********************************************************************
C   begin body of subroutine AVERTEMP

        COUNTYTEMP = 0.

C.........  Sort source array by county code
        DO I = 1, NSRC
            IDX( I ) = I
        END DO
        	
        CALL SORTI1( NSRC, IDX, SRCARRAY )

C.........  Create array of county codes
        PREVCTY = 0
        COUNTYNUM = 0
        
        DO I = 1, NSRC
            CURRCTY = SRCARRAY( IDX( I ) )
            
            IF( CURRCTY /= PREVCTY ) THEN
            	COUNTYNUM = COUNTYNUM + 1
                COUNTIES( COUNTYNUM ) = CURRCTY
            END IF
            
            PREVCTY = CURRCTY
        END DO
        
C.........  Loop through all counties
        DO I = 1, NCOUNTY
            CURRCTY = COUNTIES( I )
            
C.............  For current county, loop through all time steps            
            DO J = 0, NSTEPS - 1
                TEMPSUM = 0.
                NUMSRC  = 0
                
C.................  For current time step, loop through all sources
                DO K = 1, NSRC
                
C.....................  If current source is not in current county, skip it                
                    IF( SRCARRAY( K ) /= CURRCTY ) CYCLE
                    
                    NUMSRC = NUMSRC + 1 
                    TEMPSUM = TEMPSUM + HOURTEMP( K,J )
                END DO
                
                COUNTYTEMP( I,J ) = TEMPSUM / NUMSRC
                	
            END DO
            
        END DO

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I9, :, 1X ) )

94020   FORMAT( A, 4( 1X, F8.2, 1X, A ) )
 
        END SUBROUTINE AVERTEMP
