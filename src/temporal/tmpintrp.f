
        SUBROUTINE TMPINTRP( S, SRCTMPR, OSRC, T1, T2, PP, QQ )

C***********************************************************************
C  subroutine body starts at line 93
C
C  DESCRIPTION:
C 
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Created 10/99 by M. Houyoux
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

C...........   MODULES for public variables   
C...........   This module contains the source ararys
        USE MODSOURC

C...........   This module is the derived meteorology data for emission factors
        USE MODMET

C.........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constat parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER*2     CRLF
        INTEGER         ENVINT

        EXTERNAL        CRLF, ENVINT

C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT (IN) :: S       ! source ID
        REAL   , INTENT (IN) :: SRCTMPR ! temperature at source
        INTEGER, INTENT(OUT) :: OSRC    ! source ID
        INTEGER, INTENT(OUT) :: T1      ! valid tmpr index below SRCTMPR
        INTEGER, INTENT(OUT) :: T2      ! valid tmpr index above SRCTMPR 
        REAL   , INTENT(OUT) :: PP      ! factor 1
        REAL,    INTENT(OUT) :: QQ      ! factor 2

C...........   Other local variables

        INTEGER          L     ! counters and indices
        INTEGER, SAVE :: MAXCNT = 0     ! count of min warning
        INTEGER, SAVE :: MINCNT = 0     ! count of max warning
        INTEGER, SAVE :: MXWARN = 0     !  max warning messages to output

        REAL             SCR        ! scratch value
        REAL   , SAVE :: TINCINV    ! tmpr increment inverse

        LOGICAL, SAVE :: FIRSTIME = .TRUE.  ! true: first time routine called

        CHARACTER*300          BUFFER   ! source fields buffer
        CHARACTER*300          MESG     ! message buffer 
        CHARACTER(LEN=SRCLEN3) CSRC     ! tmp source chars string

        CHARACTER*16 :: PROGNAME = 'TMPINTRP' ! program name

C***********************************************************************
C   begin body of subroutine TMPINTRP

C.........  Compute the inverse temperature increment
        IF( FIRSTIME ) THEN
            TINCINV = 1./ TMMINVL
            MXWARN = ENVINT( WARNSET, ' ', 100, L )
            FIRSTIME = .FALSE.
        END IF

C.........  Set source characteristics string
        CSRC = CSOURC( S )

C.........   Get factors for interpolating non-diurnal emission factors
        IF( SRCTMPR .LT. AMISS3 ) THEN

            OSRC = OSRC + 1

            T1 = 1
            T2 = 1
            PP = 0.0
            QQ = 0.0

C.........  Trap source's temperature against minimum temperature available
C           in emission factors file
        ELSEIF ( SRCTMPR .LT. MINT_MIN ) THEN   

            IF( MINCNT .LE. MXWARN ) THEN

        	CALL FMTCSRC( CSRC, NCHARS, BUFFER, L )
        	WRITE( MESG,94020 )
     &             'WARNING: Bounding source temperature ', SRCTMPR, 
     &             'at minimum'// CRLF()// BLANK10 //
     &             'allowable value of ', MINT_MIN, ' for source:'//
     &             CRLF() // BLANK10 // BUFFER( 1:L )
        	CALL M3MESG( MESG )

                MINCNT = MINCNT + 1

            END IF

            PP = 1.0
            QQ = 0.0
            T1 = 1
            T2 = 2

C.........  Trap source's temperature against maximum temperature available
C           in emission factors file
        ELSE IF( SRCTMPR .GT. MAXT_MAX ) THEN

            IF( MAXCNT .LE. MXWARN ) THEN

        	CALL FMTCSRC( CSRC, NCHARS, BUFFER, L )
        	WRITE( MESG,94020 )
     &             'WARNING: Bounding source temperature', SRCTMPR, 
     &             'at maximum'// CRLF()// BLANK10 //
     &             'allowable value of', MAXT_MAX, 'for source:'//
     &             CRLF() // BLANK10 // BUFFER( 1:L )
        	CALL M3MESG( MESG )

                MAXCNT = MAXCNT + 1

            END IF

            PP = 0.0
            QQ = 1.0
            T1 = NTMPR - 1
            T2 = NTMPR

C.........  Otherwise, get valid temperature indices and set factors
        ELSE

            SCR = TINCINV * ( SRCTMPR - MINT_MIN )
            T1  = 1 + INT( SCR )
            T2  = T1 + 1
            QQ  = AMOD( SCR, 1.0 )
            PP  = 1.0 - QQ

        END IF		!  SRCTMPR trapped; coeffs & index found.

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94000   FORMAT( I2.2 )
 
94020   FORMAT( A, 2( 1X, F8.2, 1X, A ) )

        END SUBROUTINE TMPINTRP

