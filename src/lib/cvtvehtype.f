
       INTEGER FUNCTION CVTVEHTYPE( SMKVEH )
       
       IMPLICIT NONE

C.........  Function arguments
       INTEGER, INTENT (IN) :: SMKVEH   ! vehicle type in SMOKE code

C...........   LOCAL VARIABLES and their descriptions:
       LOGICAL :: EFLAG      = .FALSE.   ! true: error found

       CHARACTER*300          MESG      !  message buffer 
        
       CHARACTER*16  :: PROGNAME = 'CVTVEHTYPE' ! program name
        
C***********************************************************************
C   begin body of function CVTVEHTYPE

       SELECT CASE( SMKVEH )
       CASE( 100 )              ! light duty gasoline vehicles
           CVTVEHTYPE = 1
       CASE( 102 )              ! light duty gasoline trucks 1
           CVTVEHTYPE = 2
       CASE( 104 )              ! light duty gasoline trucks 2
           CVTVEHTYPE = 3
       CASE( 107 )              ! heavy duty gasoline vehicles
           CVTVEHTYPE = 4
       CASE( 3000 )             ! light duty diesel vehicles
           CVTVEHTYPE = 5
       CASE( 3006 )             ! light duty diesel trucks
           CVTVEHTYPE = 6
       CASE( 3007 )             ! heavy duty diesel vehicles
           CVTVEHTYPE = 7
       CASE( 108 )              ! motorcycles
           CVTVEHTYPE = 8
       CASE DEFAULT                   
           EFLAG = .TRUE.
            	
           WRITE( MESG, 94010 ) 'ERROR: Vehicle type ', SMKVEH,
     &                          'is not recognized.'
           CALL M3MESG( MESG )
       END SELECT

       IF( EFLAG ) THEN       
           MESG = 'Problem converting vehicle type'
           CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )       
       END IF

       RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )  
      
C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )
       
       END FUNCTION CVTVEHTYPE