
        SUBROUTINE UPDTMAT( NSRC, NPOL, JDATE, TZONE, VIDX, HIDX, 
     &                      MONTH, DAYOW, TMAT )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C       This routine updates the temporal matrix for sources that have source-
C       specific hourly profiles.  Since the hourly profile fractions are read
C       in one variable and hour at a time, this routine is called for each
C       hour and each pollutant/activity.
C
C  PRECONDITIONS REQUIRED:
C       Source-specific hourly profiles read for VIDX and HIDX.  MONTH and
C       DAYOW arrays populated
C
C  SUBROUTINES AND FUNCTIONS CALLED:  none
C
C  REVISION  HISTORY:
C       Copied from UPDTMAT.F 2.6 by M Houyoux 1/99
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
C.........  This module contains the inventory arrays
        USE MODSOURC, ONLY: TPFLAG, TZONES

C...........   This module contains the cross-reference tables
        USE MODXREF, ONLY: MDEX, WDEX, DDEX

C.........  This module contains the temporal profile tables
        USE MODTMPRL, ONLY: MONFAC, WEKFAC, XWKFAC

C.........  This module contains data for day- and hour-specific data
        USE MODDAYHR, ONLY: INDXD, INDXH, NHRSRC, NDYSRC, EMACH,
     &                      LDSPOA

        IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters

C...........   SUBROUTINE ARGUMENTS:

        INTEGER,INTENT (IN) :: NSRC                  ! number of sources
        INTEGER,INTENT (IN) :: NPOL                  ! number of pollutants
        INTEGER,INTENT (IN) :: JDATE                 ! date YYYYDDD
        INTEGER,INTENT (IN) :: TZONE                 ! time zone (5 for Eastern)
        INTEGER,INTENT (IN) :: VIDX                  ! pol/act index
        INTEGER,INTENT (IN) :: HIDX                  ! hour index
        INTEGER,INTENT (IN) :: MONTH( 24, 0:23 )     ! source time zone's 1...12
        INTEGER,INTENT (IN) :: DAYOW( 24, 0:23 )     ! source time zone's 1...7
        REAL   ,INTENT(OUT) :: TMAT( NSRC, NPOL, 24 )! temporal-profile coeffs

C...........   EXTERNAL FUNCTIONS:

        INTEGER         FIND1
        REAL            YR2DAY

        EXTERNAL        FIND1, YR2DAY

C...........   Other Local variables:

        INTEGER         I, J, S  ! counters and indices

        INTEGER         DAY             !  day for source and hour pointer
        INTEGER         MON             !  month for source and hour pointer

        REAL            FAC             !  partial matrix factor
        REAL            YRFAC           !  year to day factor

        CHARACTER*16  :: PROGNAME = 'UPDTMAT'  ! program name

C***********************************************************************
C   begin body of subroutine  UPDTMAT

C.........  Compute correct year-to-day conversion factor:

        YRFAC = YR2DAY( JDATE / 1000 )

C.........  Compute TMAT for current group of pollutants

        DO I = 1, NHRSRC

            S = INDXH( I )

C.............  Skip if source index is 0 (past last source for current hour),
C               or if profile value is missine
            IF( S .EQ. 0 .OR. EMACH( I ) .LE. AMISS3 ) CYCLE

C.............  Use day-specific data (no adjustments for month or 
C               weekday)
            J = FIND1( S, NDYSRC, INDXD )
            IF ( LDSPOA( VIDX ) .AND. J .GT. 0 ) THEN

                TMAT( S,VIDX,HIDX ) = EMACH( I )

C.............  Adjust for year-normal data
            ELSE IF ( MOD( TPFLAG( S ), MTPRFAC ) .EQ. 0 ) THEN

                MON = MONTH( HIDX, TZONES( S ) )
                DAY = DAYOW( HIDX, TZONES( S ) )
                FAC = MONFAC( MON,MDEX( S,VIDX ) ) * 
     &                WEKFAC( DAY,WDEX( S,VIDX ) )
                TMAT( S,VIDX,HIDX ) = FAC * EMACH( I )

C.............  Adjust for week-normal data assuming whole week normalizer
            ELSE IF ( MOD( TPFLAG( S ), WTPRFAC ) .EQ. 0 ) THEN

                DAY = DAYOW( HIDX, TZONES( S ) )
                FAC = YRFAC * WEKFAC( DAY,WDEX( S,VIDX ) )
                TMAT( S,VIDX,HIDX ) = FAC * EMACH( I )


C.............  Adjust for week-normal data assuming week-days normalizer
            ELSE IF ( MOD( TPFLAG( S ), WDTPFAC ) .EQ. 0 ) THEN

                DAY = DAYOW( HIDX, TZONES( S ) )
                FAC = YRFAC * XWKFAC( DAY,WDEX( S,VIDX ) )
                TMAT( S,VIDX,HIDX ) = FAC * EMACH( I )

            ELSE

                TMAT( S,VIDX,HIDX ) = YRFAC * EMACH( I )

            END IF

        END DO          ! end loop hour-specific sources

        RETURN

        END SUBROUTINE UPDTMAT

