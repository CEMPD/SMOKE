
        SUBROUTINE WRSMET( NSRC, JDATE, JTIME, FNAME, VMIN, VMAX,
     &                     MIDX )
   
C***********************************************************************
C  subroutine WRSMET body starts at line < >
C
C  DESCRIPTION:
C      Write per-source meteorology data
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION HISTORY:
C
C***************************************************************************
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
C****************************************************************************

        IMPLICIT NONE

C...........   INCLUDES:
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations

C...........   SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN) :: NSRC           ! no. sources
        INTEGER     , INTENT (IN) :: JDATE          ! julian date
        INTEGER     , INTENT (IN) :: JTIME          ! time HHMMSS
        CHARACTER(*), INTENT (IN) :: FNAME          ! logical file name
        REAL        , INTENT (IN) :: VMIN( NSRC,4 )   ! daily min value
        REAL        , INTENT (IN) :: VMAX( NSRC,4 )   ! daily max value
        INTEGER     , INTENT (IN) :: MIDX( NSRC,4 ) ! min/max indices

C...........   Local variables
        CHARACTER*300   MESG    ! message buffer

        CHARACTER*16 :: PROGNAME = 'WRSMET' ! program name

C***********************************************************************
C   begin body of subroutine WRSMET

C.................  Write out min/max information for the current day

        IF( .NOT. WRITE3( FNAME, 'TKMIN', JDATE, JTIME, VMIN ) ) THEN 

            MESG = 'Could not write minimum data to "' //
     &              FNAME( 1:LEN_TRIM( FNAME ) ) //  '".'

            CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )

        END IF

        IF( .NOT. WRITE3( FNAME, 'TKMAX', JDATE, JTIME, VMAX ) ) THEN 

            MESG = 'Could not write maximum data to "' //
     &              FNAME( 1:LEN_TRIM( FNAME ) ) //  '".'

            CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )

        END IF

        IF( .NOT. WRITE3( FNAME, 'TMMI', JDATE, JTIME, MIDX ) ) THEN 

            MESG = 'Could not write min/max index to "' //
     &              FNAME( 1:LEN_TRIM( FNAME ) ) //  '".'

            CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )

        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I9, :, 1X ) )

        END SUBROUTINE WRSMET
