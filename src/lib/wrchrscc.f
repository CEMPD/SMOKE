
        SUBROUTINE WRCHRSCC( FDEV, NSRC, CSCC )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      This subroutine writes character string SCC codes
C
C  PRECONDITIONS REQUIRED:
C      Output temporal x-ref file opened on unit FDEV
C      Number of sources NPSRC defined correctly
C      Index sorted in order of increasing ISCC value and ISCC populated
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C      Subroutines: I/O API subroutines
C
C  REVISION  HISTORY:
C      Created 10/98 by M. Houyoux
C
C************************************************************************
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

C.........  SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN) :: FDEV          !  ASCII file unit
        INTEGER     , INTENT (IN) :: NSRC          !  actual source count
        CHARACTER(*), INTENT (IN) :: CSCC( NSRC )  !  unsorted SCCs

C...........   Sorting index
        INTEGER       INDX( NSRC )

C...........   Other local variables
        INTEGER       J, S

        CHARACTER*300 MESG                !  message buffer
        CHARACTER(LEN=SCCLEN3) TSCC, LSCC !  current and previous 10-digit SCC

        CHARACTER*16 :: PROGNAME = 'WRCHRSCC' !  program name

C***********************************************************************
C   begin body of subroutine WRCHRSCC

C.........  Initialize SCC sorting index     
        DO S = 1, NSRC
            INDX( S ) = S
        ENDDO

C.........  Sort all SCCs in the point sources inventory in increasing order
        CALL SORTIC( NSRC, INDX, CSCC )

        LSCC = '-9'
        DO S = 1, NSRC

            J = INDX( S )

            TSCC = CSCC( J )

            IF( TSCC .NE. LSCC ) THEN
                WRITE( FDEV, '(A)', ERR=6001 ) TSCC
                LSCC = TSCC
            ENDIF 

        ENDDO 

        RETURN

6001    MESG = 'ERROR writing SCC file'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93030   FORMAT( I8.8 )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE WRCHRSCC
