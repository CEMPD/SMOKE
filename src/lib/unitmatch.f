
        SUBROUTINE UNITMATCH( UNITBUF )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      This subroutine removes or adds trailing "s" to units to make them
C      consistent with the syntax all calling programs use for that unit.  
C      If a fraction of units is given, both the numerator and 
C      denominator are modified.
C
C  PRECONDITIONS REQUIRED:

C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created by M. Houyoux 11/2000
C
C**************************************************************************
C
C Project Title: EDSS Tools Library
C File: @(#)$Id$
C
C COPYRIGHT (C) 2004, Environmental Modeling for Policy Development
C All Rights Reserved
C
C Carolina Environmental Program
C University of North Carolina at Chapel Hill
C 137 E. Franklin St., CB# 6116
C Chapel Hill, NC 27599-6116
C
C smoke@unc.edu
C
C Pathname: $Source$
C Last updated: $Date$
C
C***************************************************************************

        IMPLICIT NONE

        INCLUDE 'IOSTRG3.EXT'   !  emissions constant parameters

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN OUT) :: UNITBUF

C...........   Local variables
        INTEGER                I, L2

        LOGICAL                FRACFLAG  ! true: units are a fraction

        CHARACTER(IOULEN3) CNUM    ! tmp numerator
        CHARACTER(IOULEN3) CDEN    ! tmp denominator (if any)
        CHARACTER(IOULEN3) PREFIX1 ! tmp numerator prefix (if any)
        CHARACTER(IOULEN3) PREFIX2 ! tmp denominator prefix (if any)

        CHARACTER(16) :: PROGNAME = 'UNITMATCH' ! program name

C***********************************************************************
C   begin body of subroutine UNITMATCH

C.........  Initialize prefixes
        PREFIX1 = ' '
        PREFIX2 = ' '

C.........  Separate out the numerator and denominator, if any
        I  = INDEX( UNITBUF, '/' )
        L2 = LEN_TRIM( UNITBUF )
        IF( I .GT. 0 ) THEN
            FRACFLAG = .TRUE.
            CNUM = UNITBUF( 1  :I-1 )
            CDEN = UNITBUF( I+1:L2  )
        ELSE
            FRACFLAG = .FALSE.
            CNUM = UNITBUF
            CDEN = ' '
        END IF

C.........  Separate out any leading adjustments (e.g., 10E6) from numerator
C.........  Make sure the numerator is left-justified
        I = INDEX( CNUM, '10E' )

        IF( I .GT. 0 ) THEN
            I = INDEX( CNUM, ' ' ) 

            IF( I .GT. 0 ) THEN
                PREFIX1 = CNUM( 1:I )
                CNUM    = UNITBUF( I+1:LEN_TRIM( CNUM ) )
                CNUM    = ADJUSTL( CNUM )
            END IF

        ELSE
            CNUM = ADJUSTL( CNUM )

        END IF

C.........  Separate out any leading adjustments (e.g., 10E6) from denominator
C.........  Make sure the numerator is left-justified
        IF( FRACFLAG ) THEN

            I = INDEX( CDEN, '10E' )

            IF( I .GT. 0 ) THEN
                I = INDEX( CDEN, ' ' ) 

                IF( I .GT. 0 ) THEN
                    PREFIX1 = CDEN( 1:I )
                    CDEN    = UNITBUF( I+1:LEN_TRIM( CDEN ) )
                    CDEN    = ADJUSTL( CDEN )
                END IF

            ELSE
                CDEN = ADJUSTL( CDEN )

            END IF

        END IF

C.........  Reset name of units for numerator
        CALL RENAME_UNITS( CNUM )

C.........  Reset name of units for denominator
        IF( FRACFLAG ) CALL RENAME_UNITS( CDEN )

C.........  Piece together output units
        IF( PREFIX1 .NE. ' ' ) THEN
            L2 = LEN_TRIM( PREFIX1 )
            CNUM = PREFIX1( 1:L2 ) // ' ' // CNUM
        END IF

        IF( PREFIX2 .NE. ' ' ) THEN
            L2 = LEN_TRIM( PREFIX2 )
            CNUM = PREFIX2( 1:L2 ) // ' ' // CDEN
        END IF

        IF( FRACFLAG ) THEN
            L2 = LEN_TRIM( CNUM )
            UNITBUF = CNUM( 1:L2 ) // '/' // CDEN

        ELSE
            UNITBUF = CNUM

        END IF

        RETURN

        CONTAINS

C.............  This internal subprogram changes the units to the names
C               to be consistent.
            SUBROUTINE RENAME_UNITS( BUFFER )

            CHARACTER(*), INTENT( IN OUT ) :: BUFFER

C.............  Convert units to SMOKE notation

            SELECT CASE( BUFFER )
            CASE( 'gms' )
                BUFFER = 'g'
            CASE( 'gm' )
                BUFFER = 'g'
            CASE( 'gram' )
                BUFFER = 'g'
            CASE( 'grams' )     
                BUFFER = 'g'
            CASE( 'kilogram' )  
                BUFFER = 'kg'
            CASE( 'kilograms' ) 
                BUFFER = 'kg'
            CASE( 'kgs' )       
                BUFFER = 'kg'
            CASE( 'ton' )       
                BUFFER = 'tons'
            CASE( 'meter' )     
                BUFFER = 'm'
            CASE( 'meters' )    
                BUFFER = 'm'
            CASE( 'feet' )      
                BUFFER = 'ft'
            CASE( 'mile' )      
                BUFFER = 'miles'
            CASE( 'mi' )        
                BUFFER = 'miles'
            CASE( 'mole' )      
                BUFFER = 'moles'
            CASE( 'gm mole' )   
                BUFFER = 'moles'
            CASE( 'gm moles' )  
                BUFFER = 'moles'
            CASE( 'gm-mole' )   
                BUFFER = 'moles'
            CASE( 'gm-moles' )   
                BUFFER = 'moles'
            CASE( 'g mole' )   
                BUFFER = 'moles'
            CASE( 'g moles' )  
                BUFFER = 'moles'
            CASE( 'g-mole' )   
                BUFFER = 'moles'
            CASE( 'g-moles' )   
                BUFFER = 'moles'
            CASE( 'years' )     
                BUFFER = 'yr'
            CASE( 'year' )      
                BUFFER = 'yr'
            CASE( 'yrs' )       
                BUFFER = 'yr'
            CASE( 'h' )         
                BUFFER = 'hr'
            CASE( 'hours' )     
                BUFFER = 'hr'
            CASE( 'hour' )      
                BUFFER = 'hr'
            CASE( 'hrs' )       
                BUFFER = 'hr'
            CASE( 'dy' )        
                BUFFER = 'day'
            CASE( 'dys' )       
                BUFFER = 'day'
            CASE( 'days' )      
                BUFFER = 'day'
            CASE( 'mins' )      
                BUFFER = 'min'
            CASE( 'mns' )       
                BUFFER = 'min'
            CASE( 'mn' )        
                BUFFER = 'min'
            CASE( 'minute' )    
                BUFFER = 'min'
            CASE( 'minutes' )   
                BUFFER = 'min'
            CASE( 'secs' )      
                BUFFER = 's'
            CASE( 'seconds' )   
                BUFFER = 's'
            CASE( 'second' )    
                BUFFER = 's'
            CASE( 'sec' )       
                BUFFER = 's'
            END SELECT

            END SUBROUTINE RENAME_UNITS

        END SUBROUTINE UNITMATCH

