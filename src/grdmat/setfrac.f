
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
C      Modified by Gabe Cano 2/02 - deterministic/stochastic mode
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

C...........   This module contains the cross-reference tables
        USE MODXREF

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  i/o api constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER*2     CRLF
        INTEGER         ENVINT
        INTEGER         FIND1
        LOGICAL         ENVYN

        EXTERNAL        CRLF, ENVINT, FIND1, ENVYN

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

C...........   Local variables...

        INTEGER          L2       !  indices and counters.

        INTEGER, SAVE :: DEFSRGID !  default surrogate ID
        INTEGER, SAVE :: DEFIDX   !  default surrogate ID code index
        INTEGER, SAVE :: IDFOUND  !  surrogate ID from X-Ref
        INTEGER, SAVE :: IDUSED   !  tmp saved surrogate ID fallback
        INTEGER          IOS      !  i/o status
        INTEGER          IOU      !  i/o status for uncertainy variable
        INTEGER, SAVE :: LFIPIDX  !  last FIPs index used
        INTEGER, SAVE :: LSRCID   !  last SRCID used
        INTEGER          SRCXPOS  !  stores position in the X-Ref array
        INTEGER          TMPSRGID !  temporary surrogate ID from X-Ref
        INTEGER          TMPIDX   !  temporary surrogate ID code index
        INTEGER, SAVE :: USEIDX   !  surrogate ID code index used

        LOGICAL, SAVE :: FIRSTIME = .TRUE.  ! true: first time routine called
        LOGICAL          SRCMISS  !  true: a surrogate ID has not been utilized
        LOGICAL, SAVE :: UNCERT   !  true: uncertainty mode

        CHARACTER*300   BUFFER    !  source fields buffer
        CHARACTER*300   MESG      !  message buffer 

        CHARACTER(LEN=SRCLEN3), SAVE :: LCSRC = ' ' ! prev call src chars string

        CHARACTER*16 :: PROGNAME = 'SETFRAC' ! program name

            REAL,    SAVE :: URAND01  !  uniform quasi-random number

C***********************************************************************
C   begin body of subroutine SETFRAC

C.........  Determine default surrogate number from the environment
C.........  Default of 50 is population
        IF( FIRSTIME ) THEN

            CALL RANDOM_SEED

            MESG = 'Switch for uncertainty mode setting'
            UNCERT = ENVYN( 'SMK_UNCERT', MESG, .FALSE., IOU )

            DEFSRGID = ENVINT( 'SMK_DEFAULT_SRGID', MESG, 50, IOS )
            DEFIDX = FIND1( DEFSRGID, NSRGS, SRGLIST )

            IF( DEFIDX .LE. 0 ) THEN
        	WRITE( MESG,94010 ) 'WARNING: Fallback surrogate', 
     &             DEFSRGID, 'not found in surrogate list, resetting'//
     &             ' it to ', SRGLIST( 1 )
        	CALL M3MSG2( MESG )
        	DEFSRGID = SRGLIST( 1 )
        	DEFIDX = 1
            END IF

            FIRSTIME = .FALSE.

        END IF  !  end FIRSTIME if

        SRCXPOS = SGROWPOS( SRCID )

        IF ( SRGTOUSE( SRCID ) .EQ. IMISS3 ) THEN

            IF ( UNCERT ) THEN
                CALL STOCHASTIC
            ELSE
                CALL DETERMINISTIC
            END IF

C.............  Check that a surrogate ID was found
            IF ( IDFOUND .GT. 0 ) THEN
                USEIDX = TMPIDX
                SRGIDPOS( SRCID ) = TMPIDX
                SRGTOUSE( SRCID ) = IDFOUND

C.............  Assign NULL values for fallback retrieval
            ELSE
                USEIDX = 0
                SRGTOUSE( SRCID ) = 0

            END IF

        END IF

C.............  Retieve the information for the surrogate ID found 
        IF ( SRGTOUSE( SRCID ) .GT. 0 ) THEN 
            IDUSED = 0
            USEIDX = SRGIDPOS( SRCID )
            FRAC = SRGFRAC( USEIDX, CELIDX, FIPIDX )  

C.............  Otherwise use fallback surrogate and report zero fractions.
        ELSE
            IDUSED = DEFSRGID
            USEIDX = SRGIDPOS( SRCID )
            FRAC = SRGFRAC( DEFIDX, CELIDX, FIPIDX )  

C.............  Write note about changing surrogate used for current
C               source if it has not yet been written
            IF( REPORT .AND. CSRC .NE. LCSRC ) THEN

                CALL FMTCSRC( CSRC, NC, BUFFER, L2 )

                WRITE( MESG,94010 ) 
     &                 'WARNING: Using fallback surrogate', DEFSRGID,
     &                 CRLF()// BLANK10 // 'to prevent zeroing ' //
     &                 'by original surrogate for:'
     &                 // CRLF()// BLANK10// BUFFER( 1:L2 ) //
     &                 ' Surrogate ID: ', SRGLIST( USEIDX )
                CALL M3MESG( MESG )

C.................  Write warning for default fraction of zero
                IF( SRGCSUM( DEFIDX,FIPIDX ) .EQ. 0. ) THEN

                    CALL FMTCSRC( CSRC, NC, BUFFER, L2 )
                    MESG = 'WARNING: Fallback surrogate data '//
     &                     'will cause zero emissions' // CRLF() //
     &                     BLANK10 // 'inside the grid for:'//
     &                     CRLF() // BLANK10 // BUFFER( 1:L2 )
                    CALL M3MESG( MESG )

                END IF

            END IF  !  report check

        END IF  !  

C.........  Write surrogate code used for each source
        IF( FDEV .GT. 0 .AND. CSRC .NE. LCSRC ) THEN
            WRITE( FDEV,93360 ) SRCID, SRGLIST( USEIDX ), IDUSED
            LCSRC = CSRC  ! Store source info for next iteration
        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93360   FORMAT( I8, 1X, I4, 1X, I4 )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I9, :, 1X ) )

C******************  INTERNAL SUBPROGRAMS  *****************************

        CONTAINS

C.............  This internal subprogram stores the index of the surrogate 
C               codes from the surrogates file for each source.
            SUBROUTINE DETERMINISTIC

C.............  Local variables
            INTEGER          J        !  index variables

C----------------------------------------------------------------------

            IDFOUND = 0
            DO J = 1, INPAIRA( SRCXPOS )

C.................  Get next available surrogate code from X-Ref file
                TMPSRGID = ISRGCDA( SRCXPOS, J )
                TMPIDX = FIND1( TMPSRGID, NSRGS, SRGLIST )

                IF ( TMPIDX .LE. 0 ) CYCLE
                IF ( SRGCSUM( TMPIDX, FIPIDX ) .EQ. 0. ) CYCLE

C.................  Index of X-Ref array surrogate code to be used
                IDFOUND = J
                EXIT

            END DO  !  DO J

            RETURN

            END SUBROUTINE DETERMINISTIC

C***********************************************************************
C.............  Redistributes probability values to sum to 1 or uniformly 
C               assigns weighted uniform selection from the listed 
C               surrogate codes in the X-Ref file
            SUBROUTINE RENORMALIZE ( START, PAIRS, 
     &                              SRGCA, PROBA, SUMZERO )

C...........   SUBROUTINE ARGUMENTS
            INTEGER, INTENT (IN)     :: START          !  index starting point 
            INTEGER, INTENT (IN)     :: PAIRS          !  no. of SRG/prob. pairs
            INTEGER, INTENT (IN OUT) :: SRGCA( PAIRS ) !  array of surrogate codes
            REAL,    INTENT (IN OUT) :: PROBA( PAIRS ) !  array of prob values
            LOGICAL, INTENT (IN OUT) :: SUMZERO        !  TRUE: sum already reached 0

C.............  Local parameters
            REAL,    PARAMETER       :: NEARZERO=0.001 !  sum of prob values is 0.

C.............  Local variables
            INTEGER          I        !  index variables

            LOGICAL          ISZERO   !  TRUE: sum is zero on this iteration

            REAL             ADJFACT  !  adjustment factor
            REAL             PROBSUM  !  sum of prob values

C----------------------------------------------------------------------

            IF ( START .EQ. PAIRS ) THEN

                PROBA( PAIRS ) = 1.001

C.............  Begin renomalizing
            ELSE

C.............  Sum the probabilities to be used.
                PROBSUM = 0.0
                DO I = START, PAIRS
                    PROBSUM = PROBSUM + PROBA(I)
                END DO

                ISZERO = (PROBSUM .LE. NEARZERO )

                IF ( .NOT. ( SUMZERO ) ) THEN
                    IF ( ISZERO ) THEN
                        WRITE( MESG,94010 )
     &                     'WARNING: Sum of probabilities went to ' //
     &                     'zero.' // CRLF() // '         Remaining '
     &                     // 'surrogates will have an equal ' //
     &                     'probability of selection.'
                        CALL M3MESG( MESG )
                        SUMZERO = ISZERO
                    END IF
                END IF

C.............  If the sum of the probabilities is near 0 then evenly
C               distribute probability of selection for remaining
C               surrogates.  otherwise renormalize.
                IF ( PROBSUM .LE. NEARZERO ) THEN
                    DO I = START, PAIRS
                        PROBA(I) = 100.0 / REAL( PAIRS )
                    END DO
                ELSE
                    DO I = START, PAIRS
                        PROBA(I) = PROBA(I) * PROBSUM
                    END DO
                END IF

            END IF

            RETURN

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I9, :, 1X ) )

            END SUBROUTINE RENORMALIZE

C***********************************************************************
C.............  Performs a weighted uniform selection from the listed 
C               surrogate codes in the X-Ref file
            SUBROUTINE STOCHASTIC

C.............  Local variables
            INTEGER       :: ISEED = 27195 !  seed for random number

            INTEGER          IOS      !  I/O status variable
            INTEGER          ITMP     !  temporary space
            INTEGER          J        !  index variables
            INTEGER          PAIRS    !  number of surrogate/probability pairs
            INTEGER          START    !  starting point for selection
            INTEGER          TMPSRGID !  temporary surrogate ID

            LOGICAL          ISDEF    !  TRUE: surrogate code defined
            LOGICAL          ISZERO   !  TRUE: sum is 0

            REAL             PROBSUM  !  uniform quasi-random number
            REAL             RTMP     !  temporary space
            REAL             SELECSUM !  uniform quasi-random number
            REAL             TMPSUM   !  temporary sum for normalizing
            LOGICAL       :: PSUMZERO=.FALSE. !  sum of WRKPROB went to 0

C.............  Temporary workspace
            INTEGER, ALLOCATABLE :: WRKSRGC( : ) !  workspace for ISRGCDA
            REAL,    ALLOCATABLE :: WRKPROB( : ) !  workspace for RSPROBA

C----------------------------------------------------------------------

C.............  get the number of surrogate/probability pairs
            PAIRS = INPAIRA( SRCXPOS )

C.............  Allocate memory for workspace arrays
            ALLOCATE( WRKSRGC( PAIRS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'WRKSRGC', PROGNAME )
            ALLOCATE( WRKPROB( PAIRS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'WRKPROB', PROGNAME )

            DO J = 1, PAIRS
                WRKSRGC( J ) = J
                WRKPROB( J ) = RSPROBA( SRCXPOS, J )
            END DO

            IDFOUND = 0
            START = 1

            DO
C.................  Get quasi random number
                CALL RANDOM_NUMBER( URAND01 )

C.................  Check/ensure that probabilities add up to 1.0
                CALL RENORMALIZE( START, PAIRS, 
     &                           WRKSRGC, WRKPROB, PSUMZERO )
                SELECSUM = 0.0
                DO J = START, PAIRS

                    SELECSUM = SELECSUM + WRKPROB( J )
  
                    IF ( URAND01 .LE. SELECSUM ) THEN

                        TMPSRGID = ISRGCDA( SRCXPOS, WRKSRGC( J ) )
                        TMPIDX = FIND1( TMPSRGID, NSRGS, SRGLIST )

                        ISZERO = .FALSE.
                        ISDEF = ( TMPIDX .GT. 0 )
                        IF ( ISDEF ) 
     &                      ISZERO = (SRGCSUM( TMPIDX, FIPIDX) .EQ. 0.)

                        IF ( .NOT. (ISDEF) .OR. ISZERO ) THEN 

                            TMPSUM = 1.0 - WRKPROB( J )

                            RTMP = WRKPROB( J )
                            WRKPROB( J ) = WRKPROB( START )
                            WRKPROB( START ) = RTMP

                            ITMP = WRKSRGC( J )
                            WRKSRGC( J ) = WRKSRGC( START )
                            WRKSRGC( START ) = ITMP

                            EXIT
                        END IF

C.........................  A valid surrogate was found
                        IDFOUND = WRKSRGC( J )
                        EXIT

                    END IF  ! end random number check

                END DO  ! end J loop

                START = START + 1
                IF ( IDFOUND .GT. 0 .OR. START .GT. PAIRS ) EXIT

            END DO  ! end infinte DO

            DEALLOCATE( WRKSRGC, WRKPROB )

            RETURN

            END SUBROUTINE STOCHASTIC

        END SUBROUTINE SETFRAC
