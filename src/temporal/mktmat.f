
        SUBROUTINE MKTMAT( NSRC, NPOL, JDATE, TZONE, NDAY, ZONES, 
     &                     TPF, MDEX, WDEX, DDEX, MONTH, DAYOW, 
     &                     DFAC, TMAT )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C       Construct temporal-coefficient-transform matrices for 
C       program TMPPOINT
C
C  PRECONDITIONS REQUIRED:
C       Temporal profile arrays for monthly, weekly, diurnal profiles.
C       MDEX, WDEX entries set to zero for month- or week-independent
C       source records.  Assumes that temporal profiles have already been
C       renormalized (if desired) and weighted by days in month.  Note that
C       the monthly and weekly profiles come in from the MODTPRO, but the
C       diurnal comes in as an argument because it may be the weekday or
C       weekend.
C
C  SUBROUTINES AND FUNCTIONS CALLED:  none
C
C  REVISION  HISTORY:
C       Copied from mktmat.F 2.6 by M Houyoux 1/99
C
C***********************************************************************
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
C***********************************************************************

C.........  MODULES for public variables
C.........  This module contains the temporal profile tables
        USE MODTPRO

C.........  This module contains data for day- and hour-specific data
        USE MODDAYHR

        IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   SUBROUTINE ARGUMENTS:

        INTEGER,INTENT (IN) :: NSRC                  ! number of sources
        INTEGER,INTENT (IN) :: NPOL                  ! number of pollutants
        INTEGER,INTENT (IN) :: JDATE                 ! date YYYYDDD
        INTEGER,INTENT (IN) :: TZONE                 ! time zone (5 for Eastern)
        INTEGER,INTENT (IN) :: NDAY                  ! act. no. of diurnal profs
        INTEGER,INTENT (IN) :: ZONES( NSRC )         ! src time zones
        INTEGER,INTENT (IN) :: TPF  ( NSRC )         ! src tmprl treatment flag
        INTEGER,INTENT (IN) :: MDEX ( NSRC, NPOL )   ! monthly profile codes
        INTEGER,INTENT (IN) :: WDEX ( NSRC, NPOL )   ! weekly  profile codes
        INTEGER,INTENT (IN) :: DDEX ( NSRC, NPOL )   ! diurnal profile codes
        INTEGER,INTENT (IN) :: MONTH( 24, 0:23 )     ! source time zone's 1...12
        INTEGER,INTENT (IN) :: DAYOW( 24, 0:23 )     ! source time zone's 1...7
        REAL   ,INTENT (IN) :: DFAC ( 24, NDAY )     ! diurnal profile
        REAL   ,INTENT(OUT) :: TMAT ( NSRC, NPOL, 24 )! temporal-profile coeffs

C...........   EXTERNAL FUNCTIONS:

        INTEGER         FIND1
        LOGICAL         ISDSTIME       !  true iff daylight savings time( date)
        REAL            YR2DAY

        EXTERNAL        FIND1, ISDSTIME, YR2DAY

C...........   Other Local variables:

        INTEGER         H, I, J, K, L, S, V  ! counters and indices

        INTEGER         DAY             !  day for source and hour pointer
        INTEGER         HCORR           !  daylight savings time correction
        INTEGER         MON             !  month for source and hour pointer

        REAL            FAC             !  partial matrix factor
        REAL            YRFAC           !  year to day factor

        CHARACTER*16  :: PROGNAME = 'MKTMAT'  ! program name

C***********************************************************************
C   begin body of subroutine  MKTMAT

C.......   Compute correct year-to-day conversion factor:

        YRFAC = YR2DAY( JDATE / 1000 )

C.......   Compute index correction (offset by 1 because of
C.......   1 + MOD(...) needed below

        HCORR = TZONE + 23

C.......   Compute TMAT for current group of pollutants

        DO V = 1, NPOL

            DO S = 1, NSRC

                L = DDEX( S,V )

C.................  Use day-specific data (no adjustments for month or weekday)
                J = FIND1( S, NDYSRC, INDXD )
                IF ( LDSPOA( V ) .AND. J .GT. 0 ) THEN

                    DO H = 1, 24

                        K = 1 + MOD( H + HCORR - ZONES( S ), 24 )
                        TMAT( S,V,H ) = DFAC( K,L )

                    END DO

C.................  Adjust for year-normal data
                ELSE IF ( MOD( TPF( S ), MTPRFAC ) .EQ. 0 ) THEN

                    DO H = 1, 24

                        MON = MONTH( H, ZONES( S ) )
                        DAY = DAYOW( H, ZONES( S ) )
                        FAC = MONFAC( MON,MDEX( S,V ) ) * 
     &                        WEKFAC( DAY,WDEX( S,V ) )
                        K = 1 + MOD( H + HCORR - ZONES( S ), 24 )
                        TMAT( S,V,H ) = FAC * DFAC( K,L )

                    END DO

C.................  Adjust for week-normal data assuming whole week normalizer
                ELSE IF ( MOD( TPF( S ), WTPRFAC ) .EQ. 0 ) THEN

                    DO H = 1, 24

                        DAY = DAYOW( H, ZONES( S ) )
                        FAC = YRFAC * WEKFAC( DAY,WDEX( S,V ) )
                        K = 1 + MOD( H + HCORR - ZONES( S ), 24 )
                        TMAT( S,V,H ) = FAC * DFAC( K,L )

                    END DO

C.................  Adjust for week-normal data assuming week-days normalizer
                ELSE IF ( MOD( TPF( S ), WDTPFAC ) .EQ. 0 ) THEN

                    DO H = 1, 24
 
                        DAY = DAYOW( H, ZONES( S ) )
                        FAC = YRFAC * XWKFAC( DAY,WDEX( S,V ) )
                        K = 1 + MOD( H + HCORR - ZONES( S ), 24 )
                        TMAT( S,V,H ) = FAC * DFAC( K,L )
 
                    END DO

                ELSE

                    DO H = 1, 24

                        K = 1 + MOD( H + HCORR - ZONES( S ), 24 )
                        TMAT( S,V,H ) = YRFAC * DFAC( K,L )

                    END DO

                END IF

            END DO  ! end loop on sources, S

        END DO      ! end loop on pollutants, V

        RETURN

        END SUBROUTINE MKTMAT

