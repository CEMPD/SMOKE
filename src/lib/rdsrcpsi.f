
        SUBROUTINE RDSRCPSI( FNAME, ACT )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      This subroutine reads the PSIs per source file for 24 hours, for a
C      given activity.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
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

C.........  MODULES for public variables
C.........  This module contains emission factor tables and related
        USE MODEMFAC

C...........   This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file desc. data structures

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER*2     CRLF
        EXTERNAL        CRLF

C.........  SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN) :: FNAME     ! PSIs per source
        CHARACTER(*), INTENT (IN) :: ACT       ! activity

C.........  Other local variables
        INTEGER      I, L1, L2, L3, V       ! counters and indices
        INTEGER      IOS                    ! i/o status

        LOGICAL, SAVE :: FIRSTIME = .TRUE.  ! true: first time routine is called

        CHARACTER*300           MESG    ! message buffer
        CHARACTER(LEN=IOVLEN3)  VNAME   ! tmp inventory pollutant name

        CHARACTER*16 :: PROGNAME = 'RDSRCPSI' ! program name

C***********************************************************************
C   begin body of subroutine RDSRCPSI

C.........  For the first time the routine is called...
C.........  Allocate memory for the source-based PSIs
        IF( FIRSTIME ) THEN

            ALLOCATE( SRCPSI( NSRC,24 ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SRCPSI', PROGNAME )

            FIRSTIME = .FALSE.

        END IF

C.........  Retrieve file header
        IF ( .NOT. DESC3( FNAME ) ) THEN
            MESG = 'Could not get description of file ' // FNAME
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Find variable name for activity of interest
        DO V = 1, NVARS3D
            L1 = LEN_TRIM( ACT )
            I = INDEX( VDESC3D( V ), ACT( 1:L1 ) )
            IF( I .GT. 0 ) EXIT            ! Leave loop if activity name found
        END DO

C.........  Abort if activity not found in file
        IF( V .GT. NVARS3D ) THEN

            MESG = 'ERROR: Activity ' // ACT // ' not found in ' //
     &             'source PSIs file'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C.........  Otherwise, build variable name to read
        ELSE
            WRITE( VNAME, '("SRCPSI",I2.2)' ) V

        END IF

C.........  Read variable of interest
        IF( .NOT. READ3( FNAME, VNAME, ALLAYS3, 0, 0, SRCPSI ) ) THEN

            L1 = LEN_TRIM( VNAME )
            L2 = LEN_TRIM( ACT )
            L3 = LEN_TRIM( FNAME )

            WRITE( MESG,94010 ) 
     &             'Could not read '// VNAME(1:L1)// ' for activity '//
     &             ACT(1:L2) // 'from file "' // FNAME(1:L3) // '".'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
             
        END IF

C.........  Create unique PSI list
        DO H = 1, 24

C.............  Copy PSIs to tmp array
            TMPPSI( 1:NSRC ) = SRCPSI( 1:NSRC, H )

C.............  Create sorting index
            DO S = 1, NSRC
                INDX( S ) = S 
            END DO

C.............  Sort TMPPSI
            CALL SORTI1( NSRC, INDX, TMPPSI )

C.............  If a 

C.............  Loop through sources and flag those that

        END DO

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I10, :, 1X ) )

        END SUBROUTINE RDSRCPSI
