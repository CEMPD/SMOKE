
        SUBROUTINE OPENCTMP( PKTTYP, IDEV )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      This subroutine opens the temporary files that will contain the
C      indices to the control packet data tables
C
C      
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C
C
C****************************************************************************/
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
C***************************************************************************

        IMPLICIT NONE
        
C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS:
        INTEGER       GETEFILE

        EXTERNAL      GETEFILE

C...........   SUBROUTINE ARGUMENTS:

        CHARACTER(*), INTENT (IN) :: PKTTYP   ! packet type
        INTEGER     , INTENT(OUT) :: IDEV     ! logical file number


C...........   Other local variables

        INTEGER         IOS                   ! i/o status

        LOGICAL, SAVE :: FIRSTIME = .TRUE.  ! true: first time subroutine called
	
        CHARACTER*16    FILENM                ! file name
        CHARACTER*300   MESG                  ! message buffer
        CHARACTER, SAVE ::  PATHNM            ! path name for tmp file
        CHARACTER*16 :: PROGNAME = 'OPENCTMP' ! program name

C***********************************************************************
C   Begin body of subroutine OPENCTMP

        IF( FIRSTIME ) THEN

            FIRSTIME = .FALSE.

            MESG = 'Path where temporary control files will be written'
            CALL ENVSTR( 'TMP_CTL_PATH', MESG, '.', PATHNM, IOS )

        END IF

        SELECT CASE( PKTTYP )

           CASE ( 'CTG' )

              FILENM = PATHNM // '/.ctgtmp'
              IDEV = GETEFILE( FILENM, .FALSE., .TRUE., PROGNAME )

           CASE ( 'CONTROL' )

              FILENM = PATHNM // '/.ctltmp'
              IDEV = GETEFILE( FILENM, .FALSE., .TRUE., PROGNAME )

           CASE ( 'ALLOWABLE' )

              FILENM = PATHNM // '/.alwtmp'
              IDEV = GETEFILE( FILENM, .FALSE., .TRUE., PROGNAME )

           CASE ( 'ADD' )

              FILENM = PATHNM // '/.addtmp'
              IDEV = GETEFILE( FILENM, .FALSE., .TRUE., PROGNAME )

        END SELECT

        RETURN
        END SUBROUTINE OPENCTMP
