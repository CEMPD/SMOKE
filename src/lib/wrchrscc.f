
        SUBROUTINE WRPTSCC( FDEV, NPSRC, INDEX, ISCC )

C***********************************************************************
C  subroutine body starts at line 67
C
C  DESCRIPTION:
C      This subroutine writes a temporal x-ref file in SMOKE format
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
C****************************************************************************/
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
        INCLUDE 'PARMS3.EXT'    !  i/o api parameters

C.........  SUBROUTINE ARGUMENTS
        INTEGER      FDEV           !  ASCII file unit
        INTEGER      NPSRC          !  actual source count
        INTEGER      INDEX( NPSRC ) !  source FIPS (county) ID
        INTEGER      ISCC ( NPSRC ) !  source SCC

C...........   Other local variables
        INTEGER       J, S, SCC, LSCC

        CHARACTER*300 MESG             !  message buffer

        CHARACTER*16 :: PROGNAME = 'WRPTSCC' !  program name

C***********************************************************************
C   begin body of subroutine WRPTSCC

        LSCC = IMISS3
        DO S = 1, NPSRC

            J   = INDEX( S )
            SCC = ISCC ( J )

            IF( SCC .NE. LSCC ) THEN
                WRITE( FDEV, 93030, ERR=6001 ) SCC
                LSCC = SCC
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

        END
