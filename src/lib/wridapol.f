
        SUBROUTINE WRIDAPOL( CATEGORY, NSRC, BUFFER, VCNT, 
     &                       POLALL, STATUS )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      Write inventory pollutant-specific data for variables listed in VNAMES
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
        INTEGER         TRIMLEN

        EXTERNAL        TRIMLEN

C...........   SUBROUTINE ARGUMENTS
        CHARACTER*(*)   CATEGORY         ! Source category
        INTEGER         NSRC             ! Number of sources
        CHARACTER*(*)   BUFFER( NSRC )   ! Pollutant output buffer
        INTEGER         VCNT             ! Number of variables per pollutant
        REAL            POLDAT( NSRC,VCNT ) ! Pollutant-specific data
        INTEGER         STATUS           ! Exit status

C...........   Other local variables

        INTEGER         I, L1, L2, S     ! counters and indices

        CHARACTER*300   LINE

        CHARACTER*16 :: PROGNAME = 'WRIDAPOL' ! program name

C***********************************************************************
C   begin body of subroutine WRIDAPOL

        STATUS = 0

C.........  Get current length of buffer (NOTE: all sources are the same)
        L1 = TRIMLEN( BUFFER( 1 ) )
        L2 = LEN( BUFFER( 1 ) )

C.........  Append current pollutant to output buffer, w/ correct format
        IF( CATEGORY .EQ. 'AREA' ) THEN

            DO S = 1, NSRC

                WRITE( LINE, 94200, ERR=999 ) BUFFER( S )( 1:L1 ), 
     &               ( POLALL( S,I ), I = 1, VCNT )

                BUFFER( S ) = LINE( 1:L2 )

            ENDDO ! End loop over sources
 
        ELSEIF( CATEGORY .EQ. 'POINT' ) THEN

            DO S = 1, NSRC

                WRITE( LINE, 94210, ERR=999 ) BUFFER( S )( 1:L1 ), 
     &               ( POLALL( S,I ), I = 1, VCNT )

                BUFFER( S ) = LINE( 1:L2 )

            ENDDO ! End loop over sources

        ENDIF

        RETURN

999     MESG = 'ERROR: Problem writing to internal buffer "LINE" ' //
     &         'in subroutine ' // PROGNAME
        CALL M3MSG2( MESG )

        STATUS = 1

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94200   FORMAT( A, F10.4, F10.4, F11.4, F7.2, F3.0, F3.0 ) ! area

94210   FORMAT( A, F13.4, F13.4, F7.2, F3.0, F10.4, I3, I3 ) ! point

        END

