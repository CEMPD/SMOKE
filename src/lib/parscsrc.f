
        SUBROUTINE PARSCSRC( INSTRING, MAXN, STARTS, ENDS, OUTCOL, 
     &                       NARR, STRARRAY )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      This subroutine splits a CSOURC entry into a character
C      array and returns only non-blank strings and the array length
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C
C**************************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
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

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN) :: INSTRING       ! input string
        INTEGER     , INTENT (IN) :: MAXN           ! max number of fields
        INTEGER     , INTENT (IN) :: STARTS( MAXN ) ! starting field pos'ns
        INTEGER     , INTENT (IN) :: ENDS  ( MAXN ) ! ending field pos'ns
        LOGICAL     , INTENT (IN) :: OUTCOL( MAXN ) ! true if column is valid 
        INTEGER     , INTENT(OUT) :: NARR           ! no. of non-blank chars
        CHARACTER(*), INTENT(OUT) :: STRARRAY( * )  ! output array of strings 

C...........   Other local variables
        INTEGER         I, J, L1, L2

        CHARACTER(300)  BUFFER 

        CHARACTER(16) :: PROGNAME = 'PARSCSRC' ! program name

C***********************************************************************
C   begin body of subroutine PARSCSRC

        J = 0
        DO I = 1, MAXN
            IF( OUTCOL( I ) ) THEN
                J = J + 1

C.................  Retrieve and left-justify contents of this part of string
                L1 = STARTS( I )
                L2 = ENDS  ( I )
                BUFFER = ADJUSTL( INSTRING( L1:L2 ) )

C.................  Convert missing entries to blanks 
                IF( BUFFER .EQ. EMCMISS3 ) BUFFER = ' '

                STRARRAY( J ) = BUFFER

            ENDIF
        ENDDO

        NARR = J

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE PARSCSRC

