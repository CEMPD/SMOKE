
        INTEGER FUNCTION EFPOSITION( SCEN, POL, VTYPE, ETYPE, 
     &                               FTYPE, HOUR, M6CNVT )
    
        IMPLICIT NONE 

C...........   INCLUDES:
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'M6CNST3.EXT'   !  Mobile6 constants

C.........  Function arguments
        INTEGER, INTENT (IN) :: SCEN   ! scenario
        INTEGER, INTENT (IN) :: POL    ! pollutant
        INTEGER, INTENT (IN) :: VTYPE  ! vehicle type
        INTEGER, INTENT (IN) :: ETYPE  ! emission process
        INTEGER, INTENT (IN) :: FTYPE  ! facility type
        INTEGER, INTENT (IN) :: HOUR   ! hour of day
        LOGICAL, INTENT (IN) :: M6CNVT ! true: convert M6 pollutant to SMOKE
        
C.........  Local function variables
        INTEGER I, J
        INTEGER POLIDX                 ! pollutant index
        INTEGER NUMVEH                 ! no. vehicle types for pollutant / process combo

        LOGICAL :: EFLAG = .FALSE.     ! true: error found
         
        CHARACTER(LEN=300)     MESG    !  message buffer

        CHARACTER*16 :: PROGNAME = 'EFPOSITION'   ! program name

C***********************************************************************
C   begin body of function EFPOSITION

        POLIDX = POL
         
C.........  Convert M6 pollutant number to pollutant index if needed
        IF( M6CNVT ) THEN
            IF( POL >= 7  ) POLIDX = POL - 3    ! SO4 - G25
            IF( POL >= 12 ) POLIDX = POL - 4    ! SO2 - ACRO
            IF( POL >= 28 ) POLIDX = POL - 10   ! O10 - G10
            IF( POL >= 34 ) POLIDX = POL - 13   ! B10 - T10
        END IF

C.........  Check that input values are correct range
        IF( POLIDX > MXM6POLS .OR. VTYPE > MXM6VTYP .OR. 
     &      ETYPE > MXM6EPR .OR. FTYPE > MXM6FACS .OR. HOUR > 24 ) THEN
     	    EFLAG = .TRUE.
     	    MESG = 'INTERNAL ERROR: Arguments out of range.'
     	    CALL M3MESG( MESG )
     	END IF

C.........  Check that no facility is used for non-road emission processes
        IF( ETYPE /= 1 .AND. ETYPE /= 6 .AND. 
     &      ETYPE /= 9 .AND. ETYPE /= 10 ) THEN
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
        EFPOSITION = (SCEN - 1) * NM6PROFS

C.........  Account for previous emission processes	
        DO I = 1, ETYPE - 1
            DO J = 1, MXM6POLS
                EFPOSITION = EFPOSITION + M6PER_EPPOL( I,J )
            END DO
        END DO

C.........  Account for previous pollutants        
        DO I = 1, POLIDX - 1
            EFPOSITION = EFPOSITION + M6PER_EPPOL( ETYPE,I )
        END DO

C.........  Account for previous vehicle types
        NUMVEH = M6PER_EPPOL( ETYPE, POLIDX ) / M6RDS_EP( ETYPE )

        DO I = 1, VTYPE - 1
            IF( I <= 4 ) THEN
                EFPOSITION = EFPOSITION + M6RDS_EP( ETYPE )
            ELSE
                IF( NUMVEH == 8 ) THEN
                    EFPOSITION = EFPOSITION + M6RDS_EP( ETYPE )
                END IF
            END IF
        END DO

C.........  Account for previous facilities only if facility is not 5        
        IF( FTYPE < 5 ) THEN	
            DO I = 1, FTYPE - 1
                EFPOSITION = EFPOSITION + 1
            END DO
        END IF

C.........  Multiply by 24 hours        
        EFPOSITION = EFPOSITION * 24 + HOUR
        
        RETURN 
    
        END FUNCTION EFPOSITION

