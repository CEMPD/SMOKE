
       INTEGER FUNCTION CVTRDTYPE( SMKROAD, LASAFLAG )

C.........  MODULES for public variables       
       USE MODMBSET
       
       IMPLICIT NONE

C.........  Function arguments
       INTEGER, INTENT (IN) :: SMKROAD   ! road type in SMOKE code ( 1 - 19 )
       LOGICAL, INTENT (IN) :: LASAFLAG  ! true: treat local roads as arterial

C...........   LOCAL VARIABLES and their descriptions:
       LOGICAL :: EFLAG      = .FALSE.   ! true: error found

       CHARACTER*300          MESG      !  message buffer 
        
       CHARACTER*16  :: PROGNAME = 'CVTRDTYPE' ! program name
        
C***********************************************************************
C   begin body of function CVTRDTYPE

       SELECT CASE( SMKROAD )
       CASE( RURALINTERSTATE ) 
           CVTRDTYPE = FREEWAY
       CASE( RURALPRINCART:RURALMINORCOLL )    ! rural arterials through collectors
           CVTRDTYPE = ARTERIAL
       CASE( RURALLOCAL )
           IF( LASAFLAG ) THEN
               CVTRDTYPE = ARTERIAL
           ELSE
               CVTRDTYPE = LOCAL
           END IF
       CASE( URBANINTERSTATE:URBANFREEWAY )    ! urban interstate and freeway
           CVTRDTYPE = FREEWAY
       CASE( URBANPRINCART:URBANCOLL )         ! urban arterials through collectors
           CVTRDTYPE = ARTERIAL
       CASE( URBANLOCAL )
           IF( LASAFLAG ) THEN
               CVTRDTYPE = ARTERIAL
           ELSE
               CVTRDTYPE = LOCAL
           END IF
       CASE DEFAULT                   
           EFLAG = .TRUE.
            	
           WRITE( MESG, 94010 ) 'ERROR: Road type ', SMKROAD,
     &                          'is not recognized.'
           CALL M3MESG( MESG )
       END SELECT

       IF( EFLAG ) THEN       
           MESG = 'Problem converting road type'
           CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )       
       END IF

       RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )  
      
C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )
       
       END FUNCTION CVTRDTYPE