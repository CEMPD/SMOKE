
        SUBROUTINE SETVPTR( NDREF, NDIN, VNAMREF, VNAMIN, 
     &                      IVPTR, STATUS )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C
C  PRECONDITIONS REQUIRED:  
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       CRLF, INDEX1, TRIMLEN
C
C  REVISION  HISTORY:
C       Copied from setvptr.F by M. Houyoux 1/99
CC
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1999, MCNC--North Carolina Supercomputing Center
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
C*************************************************************************

C.........  Subroutine paramaters and their descriptions
        INTEGER       NDREF            ! actual number of reference names (in)
        INTEGER       NDIN             ! actual number of test names (in)
        CHARACTER*(*) VNAMREF( NDREF ) ! reference names (in)
        CHARACTER*(*) VNAMIN ( NDIN )  ! test names (in)
        INTEGER       IVPTR  ( NDIN )  ! pointer from test to reference (out)
        INTEGER       STATUS           ! exit status

C.........  External functions
        INTEGER       INDEX1
        INTEGER       TRIMLEN

        EXTERNAL      INDEX1, TRIMLEN

C.........  Local variables
        INTEGER       L
        INTEGER       J, V

        CHARACTER*16  VBUF       

C***********************************************************************
C   begin body of subroutine SETVPTR

C.........  Initializations
        IVPTR  = 0 ! array
        STATUS = 0

        DO V = 1, NDIN

            VBUF = VNAMIN( V )
            L = TRIMLEN( VBUF )
            J = INDEX1 ( VBUF( 1:L ), NDREF, VNAMREF )

C.............  If test name found in reference list, then store position
            IF( J .GT. 0 ) THEN
                IVPTR( V ) = J
            ELSE
                STATUS = 1
            ENDIF

        ENDDO

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats.............94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END
