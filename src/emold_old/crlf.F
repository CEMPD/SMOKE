
        CHARACTER*2 FUNCTION CRLF()

C***********************************************************************
C  subroutine body starts at line  53
C
C  FUNCTION:  Return character string of character return and line feed
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:  none
C
C  REVISION  HISTORY:
C       prototype 1/97 by M Houyoux
C
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1996, MCNC--North Carolina Supercomputing Center
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
C****************************************************************************

        IMPLICIT NONE

        CHARACTER*2 VALUE
        SAVE        VALUE

        LOGICAL     FIRSTIME
        DATA        FIRSTIME / .TRUE. /
        SAVE        FIRSTIME

C***********************************************************************
C   begin body of subroutine CRLF

        IF( FIRSTIME ) THEN
            FIRSTIME = .FALSE.
            VALUE( 1:1 ) = CHAR( 13 ) 
            VALUE( 2:2 ) = CHAR( 10 ) 
        ENDIF

        CRLF = VALUE

        RETURN

        END
