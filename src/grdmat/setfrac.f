
        SUBROUTINE SETFRAC( FDEV, SRCID, SRGIDX, CELIDX, FIPIDX, NC, 
     &                      REPORT, CSRC, FRAC )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      This subroutine sets the surrogate fraction from the surrogates
C      module. It tests to make sure that the surrogate for the county is not
C      zero, and if it is, it applies the default surrogate. It also evaluates
C      the default surrogate to use from the environment.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Created by M. Houyoux 5/99
C
C**************************************************************************
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
C***************************************************************************

C...........   MODULES for public variables
C...........   This module contains the gridding surrogates tables
        USE MODSURG

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER*2     CRLF
        INTEGER         ENVINT
        INTEGER         FIND1

        EXTERNAL        CRLF, ENVINT, FIND1

C...........   SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN) :: FDEV          ! unit no. for report file
        INTEGER     , INTENT (IN) :: SRCID         ! source ID
        INTEGER     , INTENT (IN) :: SRGIDX        ! surrogate index
        INTEGER     , INTENT (IN) :: CELIDX        ! cell index
        INTEGER     , INTENT (IN) :: FIPIDX        ! FIPS code index
        INTEGER     , INTENT (IN) :: NC            ! no. src chars for msg
        LOGICAL     , INTENT (IN) :: REPORT        ! true: okay to report
        CHARACTER(*), INTENT (IN) :: CSRC          ! source chars
        REAL        , INTENT(OUT) :: FRAC          ! surrogate fraction

C...........   Local allocatable arrays...

        INTEGER          L2       !  indices and counters.
        INTEGER, SAVE :: ID2      !  tmp saved surrogate ID fallback
        INTEGER          IOS      !  i/o status

        INTEGER, SAVE :: DEFSRGID !  default surrogate ID
        INTEGER, SAVE :: ISDEF    !  default surrogate ID code index

        LOGICAL, SAVE :: FIRSTIME = .TRUE.  ! true: first time routine called

        CHARACTER*300   BUFFER    !  source fields buffer
        CHARACTER*300   MESG      !  message buffer 

        CHARACTER(LEN=SRCLEN3), SAVE :: LCSRC = ' ' ! prev call src chars string

        CHARACTER*16 :: PROGNAME = 'SETFRAC' ! program name

C***********************************************************************
C   begin body of subroutine SETFRAC

C.........  Determine default surrogate number from the environment
C.........  Default of 50 is population
        IF( FIRSTIME ) THEN

            DEFSRGID = ENVINT( 'SMK_DEFAULT_SRGID', MESG, 50, IOS )
            ISDEF = FIND1( DEFSRGID, NSRGS, SRGLIST )

            IF( ISDEF .LE. 0 ) THEN
        	WRITE( MESG,94010 ) 'WARNING: Fallback surrogate', 
     &             DEFSRGID, 'not found in surrogate list, resetting '//
     &             'it to ', SRGLIST( 1 )
        	CALL M3MSG2( MESG )
        	DEFSRGID = SRGLIST( 1 )
        	ISDEF = 1
            END IF

            FIRSTIME = .FALSE.

        END IF

C.........  Check if surrogate selected by cross-reference for this
C                   source is non-zero in the country/state/county code of
C                   interest.
        IF( SRGCSUM( SRGIDX,FIPIDX ) .EQ. 0. ) THEN

C.............  Write note about changing surrogate used for current
C               source if it has not yet been written
            IF( REPORT .AND. CSRC .NE. LCSRC ) THEN

                CALL FMTCSRC( CSRC, NC, BUFFER, L2 )

                WRITE( MESG,94010 ) 
     &                 'WARNING: Using fallback surrogate', DEFSRGID,
     &                 CRLF()// BLANK10 // 'to prevent zeroing ' //
     &                 'by original surrogate for:'
     &                 // CRLF()// BLANK10// BUFFER( 1:L2 ) //
     &                 ' Surrogate ID: ', SRGLIST( SRGIDX )
                CALL M3MESG( MESG )

C.................  Write warning for default fraction of zero
                IF( SRGCSUM( ISDEF,FIPIDX ) .EQ. 0. ) THEN

                    CALL FMTCSRC( CSRC, NC, BUFFER, L2 )
                    MESG = 'WARNING: Fallback surrogate data '//
     &                     'will cause zero emissions' // CRLF() //
     &                     BLANK10 // 'inside the grid for:'//
     &                     CRLF() // BLANK10 // BUFFER( 1:L2 )
                    CALL M3MESG( MESG )

                END IF

            END IF

C.............  Set surrogate fraction using default surrogate
            FRAC = SRGFRAC( ISDEF,CELIDX,FIPIDX )
            ID2  = DEFSRGID

C.........  Set surrogate fraction with cross-reference-selected
C           surrogate
        ELSE

            FRAC = SRGFRAC( SRGIDX, CELIDX, FIPIDX )
            ID2 = 0

        END IF

C.........  Write surrogate code used for each source
        IF( FDEV .GT. 0 .AND. CSRC .NE. LCSRC ) THEN

            WRITE( FDEV,93360 ) SRCID, SRGLIST( SRGIDX ), ID2

            LCSRC = CSRC  ! Store source info for next iteration

        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93360   FORMAT( I8, 1X, I4, 1X, I4 )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I9, :, 1X ) )
 
        END SUBROUTINE SETFRAC

