
        SUBROUTINE WRTIMEGR( TDEV, PERIOD )

C.........  MODULES for public variables
C.........  This module contains the inventory arrays
        USE MODSOURC

C.........  This module contains the information about the source category
        USE MODINFO
        
        USE MODMBSET
        
        IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        
C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT (IN) :: TDEV     ! time grouping file no.
        INTEGER, INTENT (IN) :: PERIOD   ! period to write file for (daily, weekly, etc.)

C...........   Other local variables         
        INTEGER I, J                      ! counters and indices

        CHARACTER(LEN=FIPLEN3) REFCOUNTY                 ! ref. county FIPS code
        CHARACTER(LEN=FIPLEN3) INVCOUNTY                 ! inv. county FIPS code
                
        CHARACTER*16 :: PROGNAME = 'WRTIMEGR'   ! program name

C***********************************************************************
C   begin body of subroutine WRTIMEGR

C.........  Loop through ref. counties in MVREFSORT
        DO I = 1, NREFC
            WRITE( REFCOUNTY, '(I6)' ) MVREFSORT( I,1 )
            CALL PADZERO( REFCOUNTY )

C.............  Check if time period for current ref. county matches the one
C               we are currently processing            
            IF( MVREFSORT( I,3 ) == PERIOD ) THEN
            	
C.................  Check if this ref. county is not spatially averaged            	
                IF( MVREFSORT( I,2 ) == 1 ) THEN                 

C.....................  Loop through inv. counties using this ref. county
                    J = MCREFIDX( I,2 )

                    DO
                        WRITE( INVCOUNTY, '(I6)' ) MCREFSORT( J,1 )
                        CALL PADZERO( INVCOUNTY )
                        
                        WRITE( TDEV,93010 ) INVCOUNTY, REFCOUNTY
                        
                        J = J + 1
                        
                        IF( J > NINVC ) EXIT
                        IF( I /= NREFC ) THEN
                            IF( J == MCREFIDX( I + 1,2 ) ) EXIT
                        END IF

                    END DO

C.................  Otherwise, write the ref. county number
                ELSE
                    WRITE( TDEV,93010 ) REFCOUNTY, REFCOUNTY
                END IF
                    
            END IF ! end time period check
        END DO ! end ref. county loop

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )  
93010   FORMAT( A6, 1X, A6 )  
        
        END SUBROUTINE WRTIMEGR
        