
        SUBROUTINE BLDCSRC( CFIP , PLANT, CHAR1, CHAR2, 
     &                      CHAR3, CHAR4, CHAR5, CCOD , CSRC )

C***********************************************************************
C  subroutine body starts at line 97
C
C  DESCRIPTION:
C      This subroutine combines the source characteristics into a single
C      string.  It right justifies the fields and maintains the correct
C      segment lengths.  Note that if the results from this function are 
C      used in the calling program with concatenation, the programmer should
C      first use LEN_TRIM on the function result before doing the concatenation
C
C  PRECONDITIONS REQUIRED:
C      Values for at least CFIP, PLANT, CHAR1, and CCOD, with other CHAR* values
C      optional.
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C     Subroutines: Models-3 subroutines
C     Functions: Models-3 functions
C
C  REVISION  HISTORY:
C     Created by M. Houyoux 10/98
C
C**************************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1998, MCNC--North Carolina Supercomputing Center
C All Rights Reserved
C
C See file COPYRIGHT for conditions of use.
C
C Environmental Programs Group
C MCNC--North Carolina Supercomputing Center
C P.O. Box 12889
C Research Triangle Park, NC  27709-2889
C
C env_progs@mcnc.org
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C***************************************************************************

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(LEN=*), INTENT (IN) :: CFIP  ! string FIPS code
        CHARACTER(LEN=*), INTENT (IN) :: PLANT ! string plant code
        CHARACTER(LEN=*), INTENT (IN) :: CHAR1 ! source char 1
        CHARACTER(LEN=*), INTENT (IN) :: CHAR2 ! source char 2
        CHARACTER(LEN=*), INTENT (IN) :: CHAR3 ! source char 3
        CHARACTER(LEN=*), INTENT (IN) :: CHAR4 ! source char 4
        CHARACTER(LEN=*), INTENT (IN) :: CHAR5 ! source char 5
        CHARACTER(LEN=*), INTENT (IN) :: CCOD  ! string of int postn of pollutant
        CHARACTER(LEN=*), INTENT(OUT) :: CSRC  ! concatenated result

C...........   Other local variables
        INTEGER         L
        INTEGER         OUTLEN

        CHARACTER(LEN=FIPLEN3) B1
        CHARACTER(LEN=PLTLEN3) B2
        CHARACTER(LEN=CHRLEN3) B3
        CHARACTER(LEN=CHRLEN3) B4
        CHARACTER(LEN=CHRLEN3) B5
        CHARACTER(LEN=CHRLEN3) B6
        CHARACTER(LEN=CHRLEN3) B7
        CHARACTER(LEN=POLLEN3) B8

        CHARACTER*300   MESG 
        CHARACTER*300   STRING

        LOGICAL, SAVE :: FIRSTIME = .TRUE.

        CHARACTER*16 :: PROGNAME = 'BLDCSRC' ! program name

C***********************************************************************
C   begin body of subroutine BLDCSRC

C.........  Determine allocated length of string used for output
        OUTLEN = LEN( CSRC )

C.........  First copy the strings from the arrays of whatever length to
C           arrays that are the correct length.  Presumably, the calling
C           program has ensured that truncation won't be an issue.
        B1 = CFIP
        B2 = PLANT
        B3 = CHAR1
        B4 = CHAR2
        B5 = CHAR3
        B6 = CHAR4
        B7 = CHAR5
        B8 = CCOD

C.........  Now, right-justify using the correct lengths
        B1 = ADJUSTR( B1 )
        B2 = ADJUSTR( B2 )
        B3 = ADJUSTR( B3 )
        B4 = ADJUSTR( B4 )
        B5 = ADJUSTR( B5 )
        B6 = ADJUSTR( B6 )
        B7 = ADJUSTR( B7 )
        B8 = ADJUSTR( B8 )

        STRING = ' '
        STRING = B1 // B2 // B3 // B4 // B5 // B6 // B7 // B8
        L = LEN_TRIM( STRING )

        IF( L .GT. OUTLEN ) THEN 
            MESG= 'INTERNAL ERROR: Function BLDCSRC has been ' //
     &            'modified improperly'
            CALL M3MSG2( MESG )
            CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

        ELSE
            CSRC = STRING( 1:OUTLEN )

        ENDIF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

        END SUBROUTINE BLDCSRC

