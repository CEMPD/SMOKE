
        SUBROUTINE BLDCSRC( CFIP , STR2, STR3, STR4, 
     &                      STR5, STR6, STR7, CCOD , CSRC )

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
C      optional.  The width of the input strings *must* be consistent with
C      the right-justified width desired for output, which is defined by
C      the start and end positions of ARBEGL3/ARENDL3 for area sources,
C      MBBEGL3/MBENDL3 for mobile sources, and PTBEGL3/PTENDL3 for point 
C      sources.  These arrays are defined in EMSTRG3.EXT
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
        CHARACTER(*), INTENT (IN) :: CFIP  ! string FIPS code
        CHARACTER(*), INTENT (IN) :: STR2  ! ar:SCC, mb:SCC, pt: plant
        CHARACTER(*), INTENT (IN) :: STR3  ! ar:blk, mb:VTYPE, pt: plt char 1
        CHARACTER(*), INTENT (IN) :: STR4  ! ar:blk, mb:LINK, pt: plt char 2
        CHARACTER(*), INTENT (IN) :: STR5  ! ar:blk, mb:blk, pt: plt char 3
        CHARACTER(*), INTENT (IN) :: STR6  ! ar:blk, mb:blk, pt: plt char 4
        CHARACTER(*), INTENT (IN) :: STR7  ! ar:blk, mb:blk, pt: plt char 5
        CHARACTER(*), INTENT (IN) :: CCOD  ! string of int postn of pollutant
        CHARACTER(*), INTENT(OUT) :: CSRC  ! concatenated result

C...........   Other local variables
        INTEGER         L
        INTEGER         OUTLEN

        CHARACTER*300   MESG 
        CHARACTER*300   STRING

        LOGICAL, SAVE :: FIRSTIME = .TRUE.

        CHARACTER*16 :: PROGNAME = 'BLDCSRC' ! program name

C***********************************************************************
C   begin body of subroutine BLDCSRC

C.........  Determine allocated length of string used for output
        OUTLEN = LEN( CSRC )

C.........  Store right-justified entries in the output string

        STRING = ' '
        STRING = ADJUSTR( CFIP ) //
     &           ADJUSTR( STR2 ) //
     &           ADJUSTR( STR3 ) //
     &           ADJUSTR( STR4 ) //
     &           ADJUSTR( STR5 ) //
     &           ADJUSTR( STR6 ) //
     &           ADJUSTR( STR7 )
        STRING( POLPOS3:ALLLEN3 ) = ADJUSTR( CCOD )

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

