
        SUBROUTINE MKTMAT( NSRC, IGRP, NGRP, NPOL, JDATE, JTIME, TZONE, ZONES,
     &                     TPF, MONTH, DAYOW, DAYOM, PNAME, TMAT )

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
C***********************************************************************

C.........  MODULES for public variables
C.........  MODSOURC contains the inventory arrays
C.........  MODTMPRL contains the temporal profile tables

        USE MODSOURC, ONLY: IFIP

        USE MODTMPRL, ONLY: MONFAC, WEKFAC, HRLFAC, XWKFAC, DOMFAC, HOUR_TPROF,
     &                      METPROTYPE, METFACS, MONFAC_ORG, IPOL2D,
     &                      METPROF, MTHPROF, WEKPROF, DOMPROF, HRLPROF

        IMPLICIT NONE

C.........    INCLUDES:

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C.........    SUBROUTINE ARGUMENTS:

        INTEGER     , INTENT(IN   ) :: NSRC                 ! number of sources
        INTEGER     , INTENT(IN   ) :: IGRP                 ! pollutant group
        INTEGER     , INTENT(IN   ) :: NGRP                 ! #( pollutant groups )
        INTEGER     , INTENT(IN   ) :: NPOL                 ! number of pollutants in group
        INTEGER     , INTENT(IN   ) :: JDATE                ! date YYYYDDD
        INTEGER     , INTENT(IN   ) :: JTIME                ! time 10000
        INTEGER     , INTENT(IN   ) :: TZONE                ! time zone (5 for Eastern)
        INTEGER     , INTENT(IN   ) :: ZONES( NSRC )        ! src time zones
        INTEGER     , INTENT(IN   ) :: TPF  ( NSRC )        ! src tmprl treatment flag
        INTEGER     , INTENT(IN   ) :: MONTH( 24, 0:23 )    ! source time zone's month-of-year 1...12
        INTEGER     , INTENT(IN   ) :: DAYOW( 24, 0:23 )    ! source time zone's day-of-week   1...7
        INTEGER     , INTENT(IN   ) :: DAYOM( 24, 0:23 )    ! source time zone's day-of-month  1...31
        CHARACTER(*), INTENT(IN   ) :: PNAME                ! met-based temp profile file name
        REAL        , INTENT(  OUT) :: TMAT( NSRC,NPOL,24 ) ! temporal-profile coeffs

C.........  EXTERNAL FUNCTIONS:

        INTEGER, EXTERNAL :: FIND1, INDEX1, JULIAN
        LOGICAL, EXTERNAL :: ISDSTIME       !  true iff daylight savings time( date)
        REAL   , EXTERNAL :: YR2DAY

C.........  Other Local variables:

        INTEGER         H, I, J, K, L, S, V, VV  ! counters and indices
        INTEGER         IMET, IMTH, IWEK, IDOM, IHR

        INTEGER         DOM, DAY        !  day for source and hour pointer
        INTEGER         HCORR           !  daylight savings time correction
        INTEGER         MON             !  month for source and hour pointer
        INTEGER         TDATE, TTIME    !  tmp date and time

        REAL            FAC             !  partial matrix factor
        REAL            YRFAC           !  year to day factor

        CHARACTER(200)  MESG            !  line buffer

        CHARACTER(16), PARAMETER :: PROGNAME = 'MKTMAT'  ! program name

C***********************************************************************
C   begin body of subroutine  MKTMAT

C.........  Compute correct year-to-day conversion factor:

        YRFAC = YR2DAY( JDATE / 1000 )

C..........    Compute index correction (offset by 1 because of
C              1 + MOD(...) needed below

        HCORR = TZONE + 23

C.........    Compute TMAT for current group of pollutants

        DO VV = 1, NPOL

C.............  Skip record in case of NGRP > 1 and all fields not used
            
            V = IPOL2D( VV,IGRP )
            IF ( V .LE. 0 ) CYCLE

            DO S = 1, NSRC

                IMET = METPROF( S,V )
                IMTH = MTHPROF( S,V )
                IWEK = WEKPROF( S,V )
                IDOM = DOMPROF( S,V )

C.................  Adjust for annual data, which should always use an
C                   average-day factor (if the emissions are an annual total,
C                   then don't want to base day-of-week adjustment on weekday
C                   assumption.  The reader routines should set
C                   Also update TMAT using met-based temporal profiles

                IF ( MOD( TPF( S ), MTPRFAC ) .EQ. 0 .AND.
     &               MOD( TPF( S ), WTPRFAC ) .EQ. 0       ) THEN

C.....................  update TMAT with met-based profiles (METFACS)
                    
                    IF( IMET .GT. 0 ) THEN

                        CALL UPDATE_TMAT_METFACS( METPROTYPE, .TRUE. )

                    ELSE IF ( IDOM .GT. 0 ) THEN        !  use day-of-month profiles

                        DO H = 1, 24

                            MON  = MONTH( H, ZONES( S ) )
                            DOM  = DAYOM( H, ZONES( S ) )
                            DAY  = DAYOW( H, ZONES( S ) )
                            IHR = HRLPROF( S,DAY,V )
                            FAC = MONFAC_ORG( MON,IMTH ) *
     &                            DOMFAC( DOM,MON,IDOM )
                            K   = 1 + MOD( H + HCORR - ZONES( S ), 24 )
                            TMAT( S,VV,H ) = FAC * HRLFAC( K,IHR )

                        END DO

                    ELSE IF ( IWEK .GT. 0 ) THEN        !  use day-of-week profiles

                        DO H = 1, 24

                            MON = MONTH( H, ZONES( S ) )
                            DAY = DAYOW( H, ZONES( S ) )
                            IHR = HRLPROF( S,DAY,V )
                            FAC = MONFAC( MON,IMTH ) *
     &                            WEKFAC( DAY,IWEK )
                            K   = 1 + MOD( H + HCORR - ZONES( S ), 24 )
                            TMAT( S,VV,H ) = FAC * HRLFAC( K,IHR )

                        END DO

                    END IF

C..................  This is when annual-data field is used for storing 
C                    average-day emissions by multiplying it by 365.  
C                    Have to undo that here.
C                    Adjust for week-normal data assuming whole week normalizer
                    
                ELSE IF ( MOD( TPF( S ), WTPRFAC ) .EQ. 0 ) THEN
                
C.....................  update TMAT with met-based profiles (METFACS)
                    
                    IF( IMET .GT. 0 ) THEN

                        CALL UPDATE_TMAT_METFACS( METPROTYPE, .FALSE. )

                    ELSE IF ( IDOM .GT. 0 ) THEN        !  use day-of-month profiles

                        DO H = 1, 24

                            MON = MONTH( H, ZONES( S ) )
                            DAY = DAYOW( H, ZONES( S ) )
                            DOM = DAYOM( H, ZONES( S ) )
                            IHR = HRLPROF( S,DAY,V )
                            FAC = YRFAC * DOMFAC( DOM,MON,IDOM )
                            K   = 1 + MOD( H + HCORR - ZONES( S ), 24 )
                            TMAT( S,VV,H ) = FAC * HRLFAC( K,IHR )

                        END DO
                    
                    ELSE IF ( IWEK .GT. 0 ) THEN        !  use day-of-week profiles

                        DO H = 1, 24

                            DAY = DAYOW( H, ZONES( S ) )
                            IHR = HRLPROF( S,DAY,V )
                            FAC = YRFAC * WEKFAC( DAY,IWEK )
                            K   = 1 + MOD( H + HCORR - ZONES( S ), 24 )
                            TMAT( S,VV,H ) = FAC * HRLFAC( K,IHR )

                        END DO

                    END IF

C..................  This is when annual-data field is used for storing 
C                    average-weekday emissions by multiplying it by 365.  
C                    Have to undo that here.
C                    Adjust for week-normal data assuming weekdays normalizer
                    
                ELSE IF ( MOD( TPF( S ), WDTPFAC ) .EQ. 0 ) THEN

                    DO H = 1, 24

                        DAY = DAYOW( H, ZONES( S ) )
                        FAC = YRFAC * XWKFAC( DAY,IWEK )
                        IHR = HRLPROF( S,DAY,V )
                        K   = 1 + MOD( H + HCORR - ZONES( S ), 24 )
                        TMAT( S,VV,H ) = FAC * HRLFAC( K,IHR )

                    END DO

C..................  This is for day-specific data.
C                    Adjust without day-of-week factors
                    
                ELSE

                    DO H = 1, 24

                        DAY = DAYOW( H, ZONES( S ) )
                        IHR = HRLPROF( S,DAY,V )
                        K    = 1 + MOD( H + HCORR - ZONES( S ), 24 )
                        TMAT( S,VV,H ) = YRFAC * HRLFAC( K,IHR )

                    END DO

                END IF

            END DO  ! end loop on sources, S

        END DO      ! end loop on pollutants, V

        RETURN

C******************  INTERNAL SUBPROGRAMS  *****************************

        CONTAINS

C.............  This internal subprogram will update TMAT array with
C               Met-based temporal profiles generated by Gentpro

            SUBROUTINE UPDATE_TMAT_METFACS( PRFTYPE, IFLAG )


            CHARACTER(*), INTENT( IN ) :: PRFTYPE   !  met profile avg method
            LOGICAL,      INTENT( IN ) :: IFLAG     ! true: ann inv, false: avgday inv

C.............  Local variables

            INTEGER  HRDEX, H, I, J, K, L, IOS
            INTEGER  IYEAR, IDATE, MYEAR, MDAY, MHOUR

            REAL     METFAC, TOT, DIV, FAC

            CHARACTER(300) MESG
            CHARACTER(16)  TVARNAME

C----------------------------------------------------------------------

C............. Compute month and day of month using JDATE
            CALL DAYMON( JDATE, MON, MDAY )

C............. Search index of FIPS from METFACS

            SELECT CASE( PRFTYPE )

            CASE( 'MONTHLY' )       !  Renormalize METFACS using days of month
                TOT = 0.0
                FAC = 0.0
                DO I = 1, 12
                    TOT = TOT + METFACS( IMET, I, 1 )
                    FAC = FAC + MON_DAYS( I ) * METFACS( IMET, I, 1 )
                END DO

                IF( FAC > 0.0 ) DIV = TOT / FAC
                
                IF ( IDOM .GT. 0 ) THEN     !  use day-of-month profiles

                    DO H = 1, 24

                        MON = MONTH( H, ZONES( S ) )
                        DOM = DAYOM( H, ZONES( S ) )
                        DAY = DAYOW( H, ZONES( S ) )

                        METFAC = DIV * METFACS( IMET,MON,1 )
                        FAC    = METFAC * DOMFAC( DOM,MON,IDOM )
                        HRDEX  = HRLPROF( S,DAY,V )
                        K      = 1 + MOD( H + HCORR - ZONES( S ), 24 )
                        TMAT( S,VV,H ) = FAC * HRLFAC( K,HRDEX )

                    END DO

                ELSE IF ( IWEK .GT. 0 ) THEN     !  use day-of-week profiles

                    DO H = 1, 24

                        MON = MONTH( H, ZONES( S ) )
                        DAY = DAYOW( H, ZONES( S ) )

                        METFAC = DIV * METFACS( IMET,MON,1 )
                        FAC    = METFAC  * WEKFAC( DAY,IWEK )
                        K      = 1 + MOD( H + HCORR - ZONES( S ), 24 )
                        HRDEX  = HRLPROF( S,DAY,V )
                        TMAT( S,VV,H ) = FAC * HRLFAC( K,HRDEX )

                    END DO

                END IF          !  if IDOM>0; else if iwek>0

            CASE( 'DAILY' )

C.................  For Met-based daily profiles, renormalization is not needed.
C  CRITICAL::       Average day inventory (TPF=2[WTPRFAC]) is not applicable in this case.
C  NOTE::           it requires Gentpro updates to compute day of month temporal profiles
C                   to process monthly inventories through Temporal

                IF( .NOT. IFLAG ) THEN
                    MESG = 'ERROR: Can not apply day-of-year temporal profiles from ' //
     &                     'Genptro program to monthly inventory: NOT SUPPORTED YET'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF

                METFAC = METFACS( IMET,MON,MDAY )

                DO H = 1, 24

                    DAY    = DAYOW( H, ZONES( S ) )
                    K      = 1 + MOD( H + HCORR - ZONES( S ), 24 )
                    HRDEX  = HRLPROF( S,DAY,V )
                    TMAT( S,VV,H ) = METFAC * HRLFAC( K,HRDEX )

                END DO

            CASE( 'HOURLY' )

C.....................  Convert ann/avgday NH3 inventory to hourly NH3 before multiplying
C                       adjustment factor computed by Gentpro
                DO H = 1, 24

C.....................  Update temporal profile IDs to met-based temp profile IDs
C                       depending on whether inventory is monthly or annual
C                       APPLY hour-of-year profiles to annual inventory
C                       APPLY hour-of-month profiles to monthly inventory

                    IF( IFLAG ) THEN  !  annual total inventory

                        IF( HOUR_TPROF == 'YEAR' ) THEN
                            FAC         = 1.0
                            TVARNAME    = 'ANNTOT'

                        ELSE      ! convert annual to monthly before apply hour-of-month profiles from Gentpro
                            MON = MONTH( H, ZONES( S ) )
                            FAC = MONFAC_ORG( MON,IMTH )
                            TVARNAME    = 'MONTOT'

                        END IF

                    ELSE              !  monthly average day inventory
                        IF( HOUR_TPROF == 'YEAR' ) THEN
                            MESG = 'ERROR: Can not apply day-of-year '//
     &                          'temporal profiles from Genptro program to '//
     &                          'monthly inventory: Change HOURLY_TPROF_BASE to MONTH'
                            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                        END IF
                        MON = MONTH( H, ZONES( S ) )
                        FAC = YRFAC * MON_DAYS( MON )
                        TVARNAME    = 'MONTOT'

                    END IF

                    IF( .NOT. READ3( PNAME, TVARNAME, 1, TDATE, TTIME,
     &                        METFACS( :,1,1 ) ) ) THEN
                         MESG = 'Could not read ' // TRIM(TVARNAME) //
     &                          'from ' // TRIM( PNAME )
                         CALL M3EXIT( PROGNAME, TDATE, TTIME, MESG, 2 )
                    END IF
                    TOT = METFACS( IMET,1,1 )

                    IF( .NOT. READ3( PNAME, 'HRLSRC', 1, TDATE, TTIME,
     &                           METFACS( :,1,1 ) ) ) THEN
                        MESG = 'Could not read ' // TRIM( TVARNAME ) //
     &                         ' from ' // TRIM( PNAME )
                        CALL M3EXIT( PROGNAME, TDATE, TTIME, MESG, 2 )
                    END IF


                    TMAT( S,VV,H ) = FAC * METFACS( IMET,1,1 ) / TOT

                    CALL NEXTIME( TDATE, TTIME, 10000 )

                 END DO

            END SELECT

            END SUBROUTINE UPDATE_TMAT_METFACS

C-------------------------------------------------------------------------
C-------------------------------------------------------------------------

        END SUBROUTINE MKTMAT

