
        INTEGER FUNCTION EFPOSITION( SCEN, POL, VTYPE, ETYPE, 
     &                               FTYPE, HOUR )
    
        IMPLICIT NONE 

C.........  Function arguments
        INTEGER, INTENT (IN) :: SCEN        ! scenario
        INTEGER, INTENT (IN) :: POL         ! pollutant
        INTEGER, INTENT (IN) :: VTYPE       ! vehicle type
        INTEGER, INTENT (IN) :: ETYPE       ! emission process
        INTEGER, INTENT (IN) :: FTYPE       ! facility type
        INTEGER, INTENT (IN) :: HOUR        ! hour of day
        
C.........  Local function arrays
        INTEGER PER_FTYPE     ( 5 )      ! no. profiles per facility
        INTEGER PER_ETYPE     ( 8 )      ! no. profiles per emission process
        INTEGER PER_VTYPE     ( 8 )      ! no. profiles per vehicle type
        INTEGER PER_VTYPE_POL1( 8 )      ! no. profiles per vehicle type for pollutant 1
        INTEGER PER_POL       ( 3 )      ! no. profiles per pollutant
        
C.........  Local function variables
        INTEGER I, J
        INTEGER PER_SCEN                    ! no. profiles per scenario

        LOGICAL :: EFLAG = .FALSE.          ! true: error found
         
        CHARACTER(LEN=300)     MESG     !  message buffer

        CHARACTER*16 :: PROGNAME = 'EFPOSITION'   ! program name

C***********************************************************************
C   begin body of function EFPOSITION

C.........  Set up arrays with profile information
        PER_FTYPE      = (/1, 1, 1, 1, 1/)
        PER_ETYPE      = (/4, 1, 1, 1, 1, 4, 1, 1/)
        PER_VTYPE      = (/5, 5, 5, 5, 5, 5, 5, 5/)
        PER_VTYPE_POL1 = (/14, 14, 14, 14, 5, 5, 5, 10/)
        PER_POL        = (/81, 40, 40/)
        PER_SCEN       = 161

C.........  Check that input values are correct range
        IF( POL > 3 .OR. VTYPE > 8 .OR. ETYPE > 8 .OR. FTYPE > 5 .OR.
     &      HOUR > 24 ) THEN
     	    EFLAG = .TRUE.
     	    MESG = 'INTERNAL ERROR: Arguments out of range.'
     	    CALL M3MESG( MESG )
     	END IF

C.........  Check that only exhaust emission processes are used for
C           non-HC pollutants
     	IF( POL /= 1 .AND. ETYPE > 2 ) THEN
     	    EFLAG = .TRUE.
     	    MESG = 'INTERNAL ERROR: Cannot use non-exhaust emission ' //
     &             'processes with CO or NOx.'
            CALL M3MESG( MESG )
     	END IF

C.........  Check that no facility is used for non-road emission processes
        IF( ETYPE /= 1 .AND. ETYPE /= 6 ) THEN
            IF( FTYPE /= 5 ) THEN
                EFLAG = .TRUE.
                MESG = 'INTERNAL ERROR: Cannot use specified ' //
     &                 'facility type with non-road emission ' //
     &                 'process.'
                CALL M3MESG( MESG )
            END IF
        END IF

C.........  If error, print message and quit
        IF( EFLAG ) THEN
            MESG = 'Problem with input values.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Get to correct scenario start
        EFPOSITION = (SCEN - 1) * PER_SCEN

C.........  Account for previous pollutants        
        DO I = 1, POL - 1
            EFPOSITION = EFPOSITION + PER_POL( I )
        END DO

C.........  Account for previous vehicle types        	
        DO I = 1, VTYPE - 1
            IF( POL == 1 ) THEN
                EFPOSITION = EFPOSITION + PER_VTYPE_POL1( I )
            ELSE
                EFPOSITION = EFPOSITION + PER_VTYPE( I )
            END IF
        END DO

C.........  Account for previous emission processes        	
        DO I = 1, ETYPE - 1
            EFPOSITION = EFPOSITION + PER_ETYPE( I )
        END DO

C.........  Account for previous facilities only if facility is not 5        
        IF( FTYPE < 5 ) THEN	
            DO I = 1, FTYPE - 1
                EFPOSITION = EFPOSITION + PER_FTYPE( I )
            END DO
        END IF

C.........  Multiply by 24 hours        
        EFPOSITION = EFPOSITION * 24 + HOUR
        
        RETURN 
    
        END FUNCTION EFPOSITION

