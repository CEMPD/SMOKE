
        SUBROUTINE ADJSHOUR( NSRCIN, NSTEPS, MINTEMP, MAXTEMP, DESC, 
     &                       HOURBYSRC ) 

C***********************************************************************
C  subroutine ADJSHOUR body starts at line 121
C
C  DESCRIPTION:
C      Check hourly values against minimum and maximum
C
C  PRECONDITIONS REQUIRED:
C      The min/max criteria and the values by source should be in the same 
C      units.
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
C COPYRIGHT (C) 2000, MCNC--North Carolina Supercomputing Center
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

C...........   MODULES for public variables
C...........   This module is the source inventory arrays
        USE MODSOURC

C.........  This module contains the information about the source category
        USE MODINFO
        
        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters

C...........   EXTERNAL FUNCTIONS 
        CHARACTER*2  CRLF

        EXTERNAL     CRLF

C...........   SUBROUTINE ARGUMENTS
        INTEGER,      INTENT    (IN) :: NSRCIN              ! no. sources
        INTEGER,      INTENT    (IN) :: NSTEPS              ! no. time steps
        REAL,         INTENT    (IN) :: MINTEMP             ! min of minimum vals
        REAL,         INTENT    (IN) :: MAXTEMP             ! max of maximum vals
        CHARACTER(*), INTENT    (IN) :: DESC                ! data description
        REAL,         INTENT(IN OUT) :: HOURBYSRC( NSRCIN,0:NSTEPS-1 ) ! hourly temp values

C...........   Other local variables
        INTEGER       H, L, S      ! counters and indices

        DOUBLE PRECISION :: MIN    ! tmp min value
        DOUBLE PRECISION :: MAX    ! tmp max value

        REAL          VAL          ! tmp value

        CHARACTER*300 BUFFER        ! formatted source info for messages
        CHARACTER*300 MESG          ! message buffer
        CHARACTER(LEN=SRCLEN3) CSRC ! tmp concat source characteristics

        CHARACTER*16 :: PROGNAME = 'ADJSHOUR' ! program name

C***********************************************************************
C   begin body of subroutine ADJSHOUR

C.........  Set double precision temperature variables for computing and
C           matching
        MIN = DBLE( MINTEMP )
        MAX = DBLE( MAXTEMP )
        
C.........  Loop through sources and check for minimum and maximum values
        DO S = 1, NSRCIN
            DO H = 0, NSTEPS - 1

                CSRC = CSOURC( S )

                VAL = DBLE( HOURBYSRC( S,H ) )

C.............  Screen for missing values
                IF( VAL < AMISS3 .OR. VAL == 0. ) CYCLE

                IF( VAL < MIN ) THEN

C..............  Round value up to minimum

                    CALL FMTCSRC( CSRC, NCHARS, BUFFER, L )
                    WRITE( MESG, 94020 )
     &                     'Increasing hourly '  // DESC // ' from',
     &                     VAL, 'to', MINTEMP, 'for source' //
     &                     CRLF() // BLANK10 // BUFFER( 1:L ) // '.'
                    CALL M3MESG( MESG )

                    HOURBYSRC( S,H ) = MINTEMP 

                ELSEIF( VAL > MAX ) THEN

C..............  Round value down to maximum

                    CALL FMTCSRC( CSRC, NCHARS, BUFFER, L )
                    WRITE( MESG, 94020 )
     &                     'Decreasing hourly '  // DESC // ' from',
     &                     VAL, ' to', MAXTEMP, 'for source' //
     &                     CRLF() // BLANK10 // BUFFER( 1:L ) // '.'
                    CALL M3MESG( MESG )

                    HOURBYSRC( S,H ) = MAXTEMP
                
                END IF

            END DO
        END DO

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I9, :, 1X ) )

94020   FORMAT( A, 4( 1X, F8.2, 1X, A ) )
 
        END SUBROUTINE ADJSHOUR
