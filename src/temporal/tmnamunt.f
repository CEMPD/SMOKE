
        SUBROUTINE TMNAMUNT

C***********************************************************************
C  subroutine body starts at line 81
C
C  DESCRIPTION:
C       This program creates the temporal emissions output file variable names
C       and associated activities.  It also sets the units and conversion
C       factors for creating the output emission values.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created 10/99 by M. Houyoux
C
C*************************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2004, Environmental Modeling for Policy Development
C All Rights Reserved
C 
C Carolina Environmental Program
C University of North Carolina at Chapel Hill
C 137 E. Franklin St., CB# 6116
C Chapel Hill, NC 27599-6116
C 
C smoke@unc.edu
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C***************************************************************************

C.........  MODULES for public variables
C.........  This module contains the information about the source category
        USE MODINFO, ONLY: NIPPA, NIPOL, EANAM, EACNV, EAUNIT, EINAM

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'M6CNST3.EXT'   !  Mobile6 constants

C...........   EXTERNAL FUNCTIONS and their descriptions:
        INTEGER   INDEX1
        CHARACTER(IOULEN3) MULTUNIT
        REAL               UNITFAC 

        EXTERNAL     INDEX1, MULTUNIT, UNITFAC

C...........   Other local variables
        INTEGER         I, J, K, L, L2, M     !  counters and indices

        INTEGER         IOS               !  i/o status

        REAL            FAC1, FAC2        ! tmp conversion factors

        LOGICAL      :: FIXDESC = .FALSE. ! true: append info to description
        LOGICAL      :: EFLAG = .FALSE.   ! true: error found

        CHARACTER(16)      CURRUNIT   !  current unit
        CHARACTER(16)      CURRVNAME  !  current variable name
        CHARACTER(300)     MESG       !  message buffer
        CHARACTER(IOVLEN3) CBUF       !  tmp variable name

        CHARACTER(16) :: PROGNAME = 'TMNAMUNT' ! program name

C***********************************************************************
C   begin body of subroutine TMNAMUNT

C.........  Allocate memory for units conversions for inventory pollutants and
C           activities (stored in MODINFO)
        ALLOCATE( EACNV( NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EACNV', PROGNAME )

        EACNV  = 1.   ! array

C.........  Now loop through pollutants and create units and conversion factors
        DO I = 1, NIPOL

            M = INDEX1( EINAM( I ), NIPPA, EANAM )
            
            CBUF = EAUNIT ( M )
            FAC1 = UNITFAC( CBUF, 'tons', .TRUE. )
            FAC2 = UNITFAC( EAUNIT( M ), '1/yr', .FALSE. )

            IF ( FAC1 .LT. 0. ) FAC1 = 1.
            IF ( FAC2 .LT. 0. ) FAC2 = 1.

C.............  keep the orig miles/hr for Movesmrg to process MOVES lookup tables.
            IF( INDEX( CBUF,'miles' ) > 0 ) THEN
                EAUNIT( M ) = 'miles/hr'
            ELSE
                EAUNIT( M ) = 'tons/hr'
            END IF

            EACNV ( M ) = FAC1 / FAC2

        END DO

C.........  Abort if error was found
        IF ( EFLAG ) THEN
            MESG = 'Problem with emission types or emission factors'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE TMNAMUNT
