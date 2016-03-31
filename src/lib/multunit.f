
        CHARACTER(*) FUNCTION MULTUNIT( UNIT1, UNIT2 )

C***********************************************************************
C  function body starts at line
C
C  DESCRIPTION:
C      This function combines the two units by comparing the numerators
C      and denominators and cancelling any of the same units. It works
C      only for one string in the numerator and in the denominator, but
C      could be enhanced to be more sophisticated. 
C
C  PRECONDITIONS REQUIRED:

C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created by M. Houyoux 10/99
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
C Last updated: %G 
C
C***************************************************************************

        IMPLICIT NONE

C.........  INCLUDES:
        INCLUDE 'IOSTRG3.EXT'   !  i/o API strings

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN) :: UNIT1    ! first unit
        CHARACTER(*), INTENT (IN) :: UNIT2    ! second unit

C...........   Other local variables
        INTEGER         K1, K2, L1, L2, LD1, LD2, LN1, LN2

        CHARACTER(IOULEN3) :: DEN1 = ' '
        CHARACTER(IOULEN3) :: DEN2 = ' '
        CHARACTER(IOULEN3) :: NUM1 = ' '
        CHARACTER(IOULEN3) :: NUM2 = ' '

        CHARACTER(16) :: PROGNAME = 'MULTUNIT' ! program name

C***********************************************************************
C   begin body of function MULTUNIT

C.........  Retrieve length of unit strings
        L1 = LEN_TRIM( UNIT1 )
        L2 = LEN_TRIM( UNIT2 )

C.........  Cases where one or more of the provided units are blank
        IF( L1 .LE. 0 .AND. L2 .LE. 0 ) THEN

            MULTUNIT = 'unitless'
            RETURN

        ELSE IF( L1 .LE. 0 ) THEN

            MULTUNIT = UNIT2
            RETURN

        ELSE IF( L2 .LE. 0 ) THEN

            MULTUNIT = UNIT1
            RETURN

        END IF

C.........  Determine positions of the divide-by symbol in the units
        K1 = INDEX( UNIT1, '/' )
        K2 = INDEX( UNIT2, '/' )

        IF( K1 .LE. 0 ) K1 = L1 + 1
        IF( K2 .LE. 0 ) K2 = L2 + 1

C.........  Define numerators, denominators, and their lengths
        NUM1 = UNIT1( 1 : K1-1 )
        IF( K1 .LT. L1 ) DEN1 = UNIT1( K1+1 : L1 )

        LN1  = LEN_TRIM( NUM1 )
        LD1  = LEN_TRIM( DEN1 )

        NUM2 = UNIT2( 1 : K2-1 )
        IF( K2 .LT. L2 ) DEN2 = UNIT2( K2+1 : L2 )
   
        LN2  = LEN_TRIM( NUM2 )
        LD2  = LEN_TRIM( DEN2 )

C.........  Set units by comparing numerators and denominators
        IF( ( NUM1 .EQ. DEN2 .AND. NUM2 .EQ. DEN1 ) .OR.
     &      ( NUM1 .EQ. DEN1 .AND. NUM2 .EQ. DEN2 )      ) THEN

            MULTUNIT = 'unitless'

        ELSE IF( NUM1 .EQ. DEN2 ) THEN

            MULTUNIT = NUM2( 1:LN2 )// '/'// DEN1

        ELSE IF( NUM2 .EQ. DEN1 ) THEN

            MULTUNIT = NUM1( 1:LN1 )// '/'// DEN2

        ELSE IF( NUM1 .EQ. DEN1 ) THEN

            MULTUNIT = UNIT2

        ELSE IF( NUM2 .EQ. DEN2 ) THEN

            MULTUNIT = UNIT1

        ELSE IF( NUM1 .EQ. '1' .AND. NUM2 .EQ. '1' ) THEN

            MULTUNIT = '1/'// DEN1( 1:LD1 )// '*'// DEN2

        ELSE IF( NUM1 .EQ. '1' ) THEN

            MULTUNIT = NUM2( 1:LN2 )// '/'// DEN1( 1:LD1 )// '*'// DEN2

        ELSE IF( NUM2 .EQ. '1' ) THEN

            MULTUNIT = NUM1( 1:LN1 )// '/'// DEN1( 1:LD1 )// '*'// DEN2

        ELSE IF( DEN1 .EQ. ' ' .AND. DEN2 .EQ. ' ' ) THEN

            MULTUNIT = NUM1( 1:LN1 )// '*'// NUM2( 1:LN2 )

        ELSE IF( DEN1 .EQ. ' ' ) THEN

            MULTUNIT = NUM1( 1:LN1 )// '*'// NUM2( 1:LN2 )// '/'// DEN2

        ELSE IF( DEN2 .EQ. ' ' ) THEN

            MULTUNIT = NUM1( 1:LN1 )// '*'// NUM2( 1:LN2 )// '/'// DEN1

        ELSE

            MULTUNIT = NUM1( 1:LN1 )// '*'// NUM2( 1:LN2 )// '/'// 
     &                 DEN1( 1:LD1 )// '*'// DEN2

        END IF

        RETURN

        END FUNCTION MULTUNIT

