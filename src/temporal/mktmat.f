
        SUBROUTINE MKTMAT( NSRC, NPOL, JDATE, TZONE, ZONES, 
     &                     TPF, MDEX, WDEX, DDEX, MONTH, DAYOW, 
     &                     TMAT )

C***********************************************************************
C  subroutine body starts at line 101
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
C       the monthly and weekly profiles come in from the MODTMPRL, but the
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
C COPYRIGHT (C) 2002, MCNC Environmental Modeling Center
C All Rights Reserved
C  
C See file COPYRIGHT for conditions of use.
C  
C Environmental Modeling Center
C MCNC
C P.O. Box 12889
C Research Triangle Park, NC  27709-2889
C  
C smoke@emc.mcnc.org
C  
C Pathname: $Source$
C Last updated: $Date$ 
C  
C***********************************************************************

C.........  MODULES for public variables
C.........  This module contains the temporal profile tables
        USE MODTMPRL, ONLY: MONFAC, WEKFAC, HRLFAC, XWKFAC

        IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   SUBROUTINE ARGUMENTS:

        INTEGER,INTENT (IN) :: NSRC                  ! number of sources
        INTEGER,INTENT (IN) :: NPOL                  ! number of pollutants
        INTEGER,INTENT (IN) :: JDATE                 ! date YYYYDDD
        INTEGER,INTENT (IN) :: TZONE                 ! time zone (5 for Eastern)
        INTEGER,INTENT (IN) :: ZONES( NSRC )         ! src time zones
        INTEGER,INTENT (IN) :: TPF  ( NSRC )         ! src tmprl treatment flag
        INTEGER,INTENT (IN) :: MDEX ( NSRC, NPOL )   ! monthly profile codes
        INTEGER,INTENT (IN) :: WDEX ( NSRC, NPOL )   ! weekly  profile codes
        INTEGER,INTENT (IN) :: DDEX ( NSRC, NPOL )   ! weekday diurnal profile codes
        INTEGER,INTENT (IN) :: MONTH( 24, 0:23 )     ! source time zone's 1...12
        INTEGER,INTENT (IN) :: DAYOW( 24, 0:23 )     ! source time zone's 1...7
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

C.............  Skip record in case of NGRP > 1 and all fields not used
            IF( DDEX( 1,V ) .LE. 0 ) CYCLE

            DO S = 1, NSRC

                L = DDEX( S,V )

C.................  Adjust for annual data, which should always use an
C                   average-day factor (if the emissions are an annual total,
C                   then don't want to base day-of-week adjustment on weekday
C                   assumption.  The reader routines should set 
                IF ( MOD( TPF( S ), MTPRFAC ) .EQ. 0 .AND.
     &                    MOD( TPF( S ), WTPRFAC ) .EQ. 0       ) THEN

                    DO H = 1, 24

                        MON = MONTH( H, ZONES( S ) )
                        DAY = DAYOW( H, ZONES( S ) )
                        FAC = MONFAC( MON,MDEX( S,V ) ) * 
     &                        WEKFAC( DAY,WDEX( S,V ) )
                        K   = 1 + MOD( H + HCORR - ZONES( S ), 24 )
                        TMAT( S,V,H ) = FAC * HRLFAC( K, L, DAY )

                    END DO

C.................  This is when annual-data field is used for storing average-
C                   day emissions by multiplying it by 365.  Have to undo that
C                   here.
C.................  Adjust for week-normal data assuming whole week normalizer
                ELSE IF ( MOD( TPF( S ), WTPRFAC ) .EQ. 0 ) THEN

                    DO H = 1, 24

                        DAY = DAYOW( H, ZONES( S ) )
                        FAC = YRFAC * WEKFAC( DAY,WDEX( S,V ) )
                        K   = 1 + MOD( H + HCORR - ZONES( S ), 24 )
                        TMAT( S,V,H ) = FAC * HRLFAC( K, L, DAY )

                    END DO

C.................  This is when annual-data field is used for storing average-
C                   weekday emissions by multiplying it by 365.  Have to undo 
C                   that here.
C.................  Adjust for week-normal data assuming weekdays normalizer
                ELSE IF ( MOD( TPF( S ), WDTPFAC ) .EQ. 0 ) THEN

                    DO H = 1, 24
 
                        DAY = DAYOW( H, ZONES( S ) )
                        FAC = YRFAC * XWKFAC( DAY,WDEX( S,V ) )
                        K   = 1 + MOD( H + HCORR - ZONES( S ), 24 )
                        TMAT( S,V,H ) = FAC * HRLFAC( K, L, DAY )
 
                    END DO

C.................  This is for day-specific data.
C.................  Adjust without day-of-week factors
                ELSE

                    DO H = 1, 24

                        DAY = DAYOW( H, ZONES( S ) )
                        K   = 1 + MOD( H + HCORR - ZONES( S ), 24 )
                        TMAT( S,V,H ) = YRFAC * HRLFAC( K, L, DAY )

                    END DO

                END IF

            END DO  ! end loop on sources, S

        END DO      ! end loop on pollutants, V

        RETURN

        END SUBROUTINE MKTMAT

