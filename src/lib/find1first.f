
        INTEGER FUNCTION FIND1FIRST( KEY, N, LIST )
    
        IMPLICIT NONE
        
C...........   EXTERNAL FUNCTIONS 
        INTEGER   FIND1
        
        EXTERNAL  FIND1   

C.........  Function arguments
        INTEGER, INTENT (IN) :: KEY        ! key to search for
        INTEGER, INTENT (IN) :: N          ! number of entries in LIST
        INTEGER, INTENT (IN) :: LIST( N )  ! table to be searched
        
C.........  Local function variables            
        INTEGER INDEX

C.........  Use FIND1 to get location of key        
        INDEX = FIND1( KEY, N, LIST )
        
C.........  If the key is found, search backward until the first entry is reached            
        IF( INDEX > 0 ) THEN
            DO
                IF( INDEX < 1 .OR. LIST( INDEX ) /= KEY ) EXIT
                INDEX = INDEX - 1
            END DO
            	
            INDEX = INDEX + 1	
        END IF

        FIND1FIRST = INDEX
        
        RETURN 
    
        END FUNCTION FIND1FIRST

C..........................................................
        
        INTEGER FUNCTION FINDR1FIRST( KEY, N, LIST )
    
        IMPLICIT NONE
        
C...........   EXTERNAL FUNCTIONS 
        INTEGER   FINDR1
        
        EXTERNAL  FINDR1   

C.........  Function arguments
        REAL,    INTENT (IN) :: KEY        ! key to search for
        INTEGER, INTENT (IN) :: N          ! number of entries in LIST
        REAL,    INTENT (IN) :: LIST( N )  ! table to be searched
        
C.........  Local function variables            
        INTEGER INDEX

C.........  Use FINDR1 to get location of key        
        INDEX = FINDR1( KEY, N, LIST )
        
C.........  If the key is found, search backward until the first entry is reached            
        IF( INDEX > 0 ) THEN
            DO
                IF( INDEX < 1 .OR. LIST( INDEX ) /= KEY ) EXIT
                INDEX = INDEX - 1
            END DO
            	
            INDEX = INDEX + 1	
        END IF

        FINDR1FIRST = INDEX
        
        RETURN 
    
        END FUNCTION FINDR1FIRST