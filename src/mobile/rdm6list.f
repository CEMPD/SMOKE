
        SUBROUTINE RDM6LIST( LDEV )

C.........  MODULES for public variables
                
        USE MODMBSET
                
        IMPLICIT NONE

C...........   EXTERNAL FUNCTIONS and their descriptions:
        INTEGER        GETFLINE
        
        EXTERNAL  GETFLINE

C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT (IN) :: LDEV     ! M6LIST file unit no.

C...........   Other local variables               

        INTEGER IOS                       ! I/O status        
        INTEGER NLINES                    ! number of lines

        CHARACTER*16 :: PROGNAME = 'RDM6LIST'   ! program name
        
C***********************************************************************
C   begin body of subroutine RDM6LIST

C.........  Get number of lines in file
        NLINES = GETFLINE( LDEV, 'M6LIST file' )
        
C.........  Allocate memory for storing the file
        ALLOCATE( M6LIST( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'M6LIST', PROGNAME )
        
C.........  Store lines of M6LIST file
        CALL RDLINES( LDEV, 'M6LIST file', NLINES, M6LIST )

        RETURN
        
        END SUBROUTINE RDM6LIST
        