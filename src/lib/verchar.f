
        CHARACTER*16 FUNCTION VERCHAR( BUFFER )

C***********************************************************************
C  function body starts at line
C
C  DESCRIPTION:
C      This function returns the version number from a string filled in
C      by SCCS operations
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

C...........   EXTERNAL FUNCTIONS:
        INTEGER    TRIMLEN

        EXTERNAL   TRIMLEN

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(LEN=*)    BUFFER         ! string FIPS code

C...........   Other local variables
        INTEGER         J, L2

        CHARACTER*16 :: PROGNAME = 'VERCHAR' ! program name

C***********************************************************************
C   begin body of function VERCHAR

        L2 = TRIMLEN( BUFFER )
        VERCHAR = ADJUSTL( BUFFER( 1:L2 ) )

        IF( VERCHAR( 1:1 ) .EQ. '%' ) THEN
            VERCHAR = 'Dev'
        ELSE
            J  = INDEX  ( VERCHAR, ' ' )
            L2 = TRIMLEN( VERCHAR )
            VERCHAR = ADJUSTL( VERCHAR( J+1:L2 ) )
        ENDIF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

        END

