
        SUBROUTINE WRM6HEADER( MDEV )
        
        IMPLICIT NONE
        
C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT (IN) :: MDEV    ! M6 input file unit no.        

C***********************************************************************
C   begin body of subroutine WRM6HEADER

        WRITE( MDEV,93000 ) 'MOBILE6 INPUT FILE :'
        WRITE( MDEV,93000 ) 'NO DESC OUTPUT     :'
        WRITE( MDEV,93000 ) 'DATABASE OUTPUT    :'
        WRITE( MDEV,93000 ) 'POLLUTANTS         : CO NOX HC'
        WRITE( MDEV,93000 ) 'PARTICULATES       :'
        WRITE( MDEV,93000 ) 'AIR TOXICS         :'
        WRITE( MDEV,93000 ) 'RUN DATA           :'
        WRITE( MDEV,93000 ) ' '

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )
        
        END SUBROUTINE WRM6HEADER