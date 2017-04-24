
      SUBROUTINE XFRACTBL( NDIM, NXREF, INDX, CSRC, SPRF, FRAC )

C***********************************************************************
C  subroutine body starts at line  83
C
C  DESCRIPTION:
C       Construct data structures for in-XREF speciation-profile fractions.
C
C  PRECONDITIONS REQUIRED:
C           Read and SORTIC() the XREFS.
C           Must call before XREFTBL()
C           FRAC(K) < 0.0 for non-frac XREFs
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C       Created 8/2016 by Carlie J. Coats, Jr., UNC IE, to treat
C       xref fractional speciation-profiles
C
C****************************************************************************/
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C
C COPYRIGHT (C) 2016, Center Environmental Modeling for Policy Development
C UNC Institute for the Environment.  All Rights Reserved
C****************************************************************************/

C.........  MODULES for public variables

        USE MODSPRO, ONLY:
     &        CMBMAX, CMBCNT, CMBNP, CMBWGHT, CMBSPCD, CMBPRF, CMBDEX
        USE MODXREF, ONLY: CFIPTA, CSCCTA, ISPTA, CMACTA, CISICA, CSPRNA
        USE MODINFO, ONLY: EINAM, NCHARS

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   Local parameters

        INTEGER     , PARAMETER :: CSRCLEN = SSMLEN3 + POLLEN3
        CHARACTER(1), PARAMETER :: BLANK   = ' '

C...........   Subroutine arguments

        INTEGER            , INTENT(INOUT) :: NDIM          !  no. ungrouped x-ref entries
        INTEGER            , INTENT(INOUT) :: NXREF         !  no. ungrouped x-ref entries
        INTEGER            , INTENT(IN   ) :: INDX( NDIM )  !  index returned by SORTIC()
        CHARACTER(CSRCLEN) , INTENT(IN   ) :: CSRC( NDIM )  !  source chars (unsorted)
        CHARACTER(SPNLEN3) , INTENT(INOUT) :: SPRF( NDIM )  !  speciation profiles (unsorted)
        REAL               , INTENT(IN   ) :: FRAC( NDIM )  !  profile-fractions

C...........   EXTERNAL FUNCTIONS:

        LOGICAL, EXTERNAL :: ENVYN

C...........   Local variabless

        INTEGER     I, J, JJ, K, KK, L, L2, M, N, IOS
        INTEGER     NFREF
        LOGICAL     AFLAG, EFLAG, SUMCHECK
        REAL        WT, WSUM

        CHARACTER(SPNLEN3)  ASPRF
        CHARACTER(FIPLEN3)  AFIP
        CHARACTER(SCCLEN3)  ASCC
        CHARACTER(SICLEN3)  ASIC
        CHARACTER(MACLEN3)  AMCT
        CHARACTER(POLLEN3)  APOL
        CHARACTER(256)      MESG
        CHARACTER(300)      BUFFER

        CHARACTER(16), PARAMETER :: PROGNAME = 'XFRACTBL' ! subroutine name

C.........  Statement function:  "is not approximately one", with tolerance TOL=1%

        REAL,    PARAMETER :: TOL   = 1.0e-2
        REAL,    PARAMETER :: TOLSQ = TOL**2

        REAL    X
        LOGICAL NOT_ONE
        NOT_ONE( X ) = ( ( X - 1.0 )**2 .GT. TOLSQ )

C***********************************************************************
C   Begin body of subroutine COMBOTBL

        SUMCHECK = ENVYN( 'CMB_CHKFRACS',
     &         'Bad sum-of-weights for fractional-profiles is fatal?',
     &                    .TRUE., IOS )
        IF ( IOS .GT. 0 ) THEN
            MESG = 'Bad environment variable "CMB_CHKFRACS"'
            CALL M3EXIT( PROGNAME, 0,0, MESG, 2 )
        END IF

C.............  Count this-pollutant fractional (duplicate) references:

        EFLAG = .FALSE.
        NFREF = 0
        J     = 1
        DO              !!  looping on J by number L of duplicates...

            K = INDX( J )
            CALL FMTCSRC( CSRC(K), NCHARS, BUFFER, L2 )

            IF ( FRAC(K) .GE. 0.0 ) THEN
                L = 1
                DO JJ = J+1, NXREF
                    KK = INDX( JJ )
                    IF ( CSRC(K) .EQ. CSRC(KK) ) THEN
                        L = L + 1
                        IF ( FRAC(KK) .LT. 0.0 )  THEN
                            EFLAG = .TRUE.
                            MESG  = 'Inconsistent XREFs for '//
     &                               TRIM( BUFFER ) // BLANK // 
     &                               EINAM( ISPTA( K ) )
                            CALL M3MESG( MESG )
                        END IF
                    ELSE
                        EXIT
                    END IF
                END DO
                NFREF = NFREF + 1

            ELSE        !!  frac(k) < 0:  not a fractions xref
                L = 1
                DO JJ = J+1, NXREF
                    KK = INDX( JJ )
                    IF ( CSRC(K) .EQ. CSRC(KK) ) THEN
                        L = L + 1
                        IF ( FRAC(KK) .GE. 0.0 )  THEN
                            EFLAG = .TRUE.
                            MESG  = 'Inconsistent XREFs for '//
     &                               TRIM( BUFFER ) // BLANK //
     &                               EINAM( ISPTA( K ) )
                            CALL M3MESG( MESG )
                        END IF
                    ELSE
                        EXIT
                    END IF
                END DO
                IF ( L .GT. 1 ) THEN
                    WRITE( MESG, '( A, I3, A  )' )
     &                   'XREF has', L, ' duplicates for ' //
     &                   TRIM( BUFFER ) // BLANK // EINAM( ISPTA( K ) )
                    CALL M3MESG( MESG )
                END IF

            END IF

            J = J + L
            
            IF ( J .GT. NXREF ) EXIT
            
        END DO      !!  end loop on J

        IF ( EFLAG ) THEN
            MESG = 'ERROR:  mixed profile-fractions/non-fractions for the same XREF(s)'
            CALL M3EXIT( PROGNAME, 0,0, MESG, 2 )
        END IF

        IF ( NFREF .EQ. 0 )   RETURN

C.............  Allocate arrays for tables

        ALLOCATE(   CMBNP( NFREF ),
     &            CMBWGHT( NFREF,CMBMAX ),
     &            CMBSPCD( NFREF,CMBMAX ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CMBNP:CMBSPCD', PROGNAME )

        AFLAG = .FALSE.
        CMBNP = 0       !!  array
        M     = 0

        DO N = 1, NXREF

            K = INDX( N )
            IF ( SPRF( K )(1:1) .EQ. '_' )  CYCLE    !!  already-processed

            !!........  Loop:  count matches to CSRC(K)

            IF ( FRAC(K) .GE. 0.0 ) THEN

                L = 1
                DO J = N+1, NXREF
                    KK = INDX( J )
                    IF ( CSRC( K ) .EQ. CSRC( KK ) ) THEN
                        L = L + 1
                    ELSE
                        EXIT
                    END IF
                END DO

                !!........  Process matches to CSRC(K)

                M          = M + 1
                CMBNP( M ) = L
                ASPRF      = CMBPRF( M )    !! begins with underscore '_'
                WSUM       = 0.0
                I          = 0
                DO J = N, N+L-1
                    I              = I + 1      !!  I=1...L
                    KK             = INDX( J )
                    WT             = FRAC( KK )
                    WSUM           = WSUM + WT
                    CMBWGHT( M,I ) = WT
                    CMBSPCD( M,I ) = SPRF( KK )
                    SPRF( KK )     = ASPRF 
                END DO

                IF ( NOT_ONE( WSUM ) )  THEN
                    CALL FMTCSRC( CSRC(K), NCHARS, BUFFER, L2 )
                    APOL  = EINAM( ISPTA( K ) )
                    MESG = 'WARNING: XREF fracs do not add up to 1 for ' //
     &                     TRIM( BUFFER ) // BLANK // EINAM( ISPTA( K ) )
                    CALL M3MESG( MESG )
                    AFLAG = .TRUE.
                END IF

            END IF

        END DO

        CMBCNT = M

        IF ( AFLAG .AND. SUMCHECK ) THEN
            MESG = 'ERROR:  XREF profile-fractions do not add up to 1.0'
            CALL M3EXIT( PROGNAME, 0,0, MESG, 2 )
        END IF

        RETURN

      END SUBROUTINE XFRACTBL


