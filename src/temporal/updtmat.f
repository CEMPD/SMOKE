
        SUBROUTINE UPDTMAT( NSRC, IGRP, NPOL, JDATE, TZONE, VIDX, HIDX, 
     &                      MONTH, DAYOW, DAYOM )

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
C       Version 07/2014 by C.Coats for  new GENTPRO CSV profiles and cross-references
C
C***********************************************************************
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
C       Updated with USE M3UTILIO by Huy Tran UNC-IE on 2026-01
C***********************************************************************

C.........  MODULES for public variables
C.........  MODSOURC contains the inventory arrays
C.........  MODTMPRL contains the temporal profile tables
C.........  MODDAYHR contains data for day- and hour-specific data

        USE M3UTILIO

        USE MODSOURC, ONLY: TPFLAG, TZONES

        USE MODTMPRL, ONLY: MONFAC,  DOMFAC,  WEKFAC,  XWKFAC, IPOL2D, 
     &                      MTHPROF, DOMPROF, WEKPROF

        USE MODDAYHR, ONLY: INDXD, INDXH, NHRSRC, NDYSRC, EMACH, LDSPOA,
     &                      TMAT

        IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
C        INCLUDE 'PARMS3.EXT'    !  I/O API parameters

C...........   SUBROUTINE ARGUMENTS:

        INTEGER, INTENT(IN ) :: NSRC                    ! number of sources
        INTEGER, INTENT(IN ) :: IGRP                    ! pollutant group
        INTEGER, INTENT(IN ) :: NPOL                    ! number of pollutants
        INTEGER, INTENT(IN ) :: JDATE                   ! date YYYYDDD
        INTEGER, INTENT(IN ) :: TZONE                   ! time zone (5 for Eastern)
        INTEGER, INTENT(IN ) :: VIDX                    ! pol/act index
        INTEGER, INTENT(IN ) :: HIDX                    ! hour index
        INTEGER, INTENT(IN ) :: MONTH( 24, -23:23 )       ! source time zone's month-of-year 1...12
        INTEGER, INTENT(IN ) :: DAYOW( 24, -23:23 )       ! source time zone's day-of-week   1...7
        INTEGER, INTENT(IN ) :: DAYOM( 24, -23:23 )       ! source time zone's day-of-month  1...31

C...........   EXTERNAL FUNCTIONS:

C       INTEGER, EXTERNAL :: FIND1
C       REAL   , EXTERNAL :: YR2DAY

C...........   Other Local variables:

        INTEGER         I, J, S, V  ! counters and indices

        INTEGER         DAY, DOM        !  day for source and hour pointer
        INTEGER         MON             !  month for source and hour pointer
        INTEGER         IMON, IWEK, IDOM

        REAL            FAC             !  partial matrix factor
        REAL            YRFAC           !  year to day factor

        CHARACTER(16), PARAMETER :: PROGNAME = 'UPDTMAT'  ! program name

C***********************************************************************
C   begin body of subroutine  UPDTMAT

C.........  Compute correct year-to-day conversion factor:

        YRFAC = YR2DAY( JDATE / 1000 )            
        V     = IPOL2D( VIDX,IGRP )

C.........  Compute TMAT for current group of pollutants

        DO I = 1, NHRSRC

            S = INDXH( I )

            IMON = MTHPROF( S,V )
            IWEK = WEKPROF( S,V )
            IDOM = DOMPROF( S,V )

            MON = MONTH( HIDX, TZONES( S ) )
            DAY = DAYOW( HIDX, TZONES( S ) )
            DOM = DAYOM( HIDX, TZONES( S ) )

C.............  Skip if source index is 0 (past last source for current hour),
C               or if profile value is missing

            IF( S .EQ. 0 .OR. EMACH( I ) .LE. AMISS3 ) CYCLE

C.............  Use day-specific data (no adjustments for month or weekday)

            J = FIND1( S, NDYSRC, INDXD )
            IF ( LDSPOA( VIDX ) .AND. J .GT. 0 ) THEN

                TMAT( S,VIDX,HIDX ) = EMACH( I )

C.............  Adjust for year-normal data

            ELSE IF ( MOD( TPFLAG( S ), MTPRFAC ) .EQ. 0 ) THEN

                IF ( IMON .GT. 0 ) THEN
                    FAC = MONFAC( MON,IMON ) * DOMFAC( DOM,MON,IDOM )
                ELSE
                    FAC = MONFAC( MON,IMON ) * WEKFAC( DAY,IWEK )
                END IF
                TMAT( S,VIDX,HIDX ) = FAC * EMACH( I )

C.................  Adjust for week-normal data assuming whole week normalizer

            ELSE IF ( MOD( TPFLAG( S ), WTPRFAC ) .EQ. 0 ) THEN

                IF ( IMON .GT. 0 ) THEN
                    FAC = YRFAC * DOMFAC( DOM,MON,IDOM )
                ELSE
                    FAC = YRFAC * WEKFAC( DAY,IWEK )
                END IF
                TMAT( S,VIDX,HIDX ) = FAC * EMACH( I )


C.................  Adjust for week-normal data assuming week-days normalizer

            ELSE IF ( MOD( TPFLAG( S ), WDTPFAC ) .EQ. 0 ) THEN

                IF ( IMON .GT. 0 ) THEN
                    FAC = YRFAC * DOMFAC( DOM,MON,IDOM )
                ELSE
                    FAC = YRFAC * XWKFAC( DAY,IWEK )
                END IF
                TMAT( S,VIDX,HIDX ) = FAC * EMACH( I )

            ELSE

                TMAT( S,VIDX,HIDX ) = YRFAC * EMACH( I )

            END IF

        END DO          ! end loop hour-specific sources

        RETURN

        END SUBROUTINE UPDTMAT

