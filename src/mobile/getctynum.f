
        SUBROUTINE GETCTYNUM( NSRC, SRCARRAY, NCOUNTY )
                
        IMPLICIT NONE

C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT (IN)  :: NSRC           ! number of sources
        INTEGER, INTENT (IN)  :: SRCARRAY( NSRC ) ! array holding county codes for each source
        INTEGER, INTENT (OUT) :: NCOUNTY        ! number of unique counties

C...........   Local arrays
        INTEGER IDX( NSRC )               ! index into sources array

C...........   Other local variables
        INTEGER I, J                      ! counters and indices                     
        
        INTEGER PREVCTY                   ! previous county code
        INTEGER CURRCTY                   ! current county code

        CHARACTER*16 :: PROGNAME = 'GETCTYNUM'   ! program name

C***********************************************************************
C   begin body of subroutine GETCTYNUM        
        
C.........  Sort array by county code
        DO I = 1, NSRC
            IDX( I ) = I
        END DO
        	
        CALL SORTI1( NSRC, IDX, SRCARRAY )

        PREVCTY = 0
        NCOUNTY = 0
     
C.........  Loop through sources and count number of unique counties
        DO I = 1, NSRC
        
            CURRCTY = SRCARRAY( IDX( I ) )
        
C.........  Skip unused sources
            IF( CURRCTY == 0 ) CYCLE
        
            IF( CURRCTY /= PREVCTY ) THEN
                NCOUNTY = NCOUNTY + 1
            END IF

            PREVCTY = CURRCTY
        
        END DO        
        
        END SUBROUTINE GETCTYNUM
        