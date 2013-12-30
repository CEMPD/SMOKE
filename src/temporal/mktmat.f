
        SUBROUTINE MKTMAT( NSRC, NPOL, JDATE, JTIME, TZONE, ZONES, 
     &                     TPF, MDEX, WDEX, DDEX, MONTH, DAYOW, PNAME,
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
C.........  This module contains the temporal profile tables
        USE MODTMPRL, ONLY: MONFAC, WEKFAC, HRLFAC, XWKFAC, HOUR_TPROF,
     &                      METPROTYPE, METPROF, NMETPROF, METFACS,
     &                      MONFAC_ORG 
     
        USE MODSOURC, ONLY: CIFIP

        IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  i/o api parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.
        
C...........   SUBROUTINE ARGUMENTS:

        INTEGER,INTENT (IN) :: NSRC                  ! number of sources
        INTEGER,INTENT (IN) :: NPOL                  ! number of pollutants
        INTEGER,INTENT (IN) :: JDATE                 ! date YYYYDDD
        INTEGER,INTENT (IN) :: JTIME                 ! time 10000
        INTEGER,INTENT (IN) :: TZONE                 ! time zone (5 for Eastern)
        INTEGER,INTENT (IN) :: ZONES( NSRC )         ! src time zones
        INTEGER,INTENT (IN) :: TPF  ( NSRC )         ! src tmprl treatment flag
        INTEGER,INTENT (IN OUT ) :: MDEX ( NSRC, NPOL )   ! monthly profile codes
        INTEGER,INTENT (IN OUT ) :: WDEX ( NSRC, NPOL )   ! weekly  profile codes
        INTEGER,INTENT (IN OUT ) :: DDEX ( NSRC, NPOL )   ! weekday diurnal profile codes
        INTEGER,INTENT (IN) :: MONTH( 24, 0:23 )     ! source time zone's 1...12
        INTEGER,INTENT (IN) :: DAYOW( 24, 0:23 )     ! source time zone's 1...7
        CHARACTER(*),INTENT(IN):: PNAME              ! met-based temp profile file name
        REAL   ,INTENT(OUT) :: TMAT ( NSRC, NPOL, 24 )! temporal-profile coeffs

C...........   EXTERNAL FUNCTIONS:

        INTEGER         FINDC, INDEX1, JULIAN
        LOGICAL         ISDSTIME       !  true iff daylight savings time( date)
        REAL            YR2DAY

        EXTERNAL        FINDC, INDEX1, ISDSTIME, YR2DAY, JULIAN

C...........   Other Local variables:

        INTEGER         H, I, J, K, L, S, V  ! counters and indices

        INTEGER         DAY             !  day for source and hour pointer
        INTEGER         HCORR           !  daylight savings time correction
        INTEGER         MON             !  month for source and hour pointer
        INTEGER         TDATE, TTIME    !  tmp date and time

        REAL            FAC             !  partial matrix factor
        REAL            YRFAC           !  year to day factor
        
        CHARACTER(200)  MESG            !  line buffer

        CHARACTER(16) :: PROGNAME = 'MKTMAT'  ! program name

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
C                   Also update TMAT using met-based temporal profiles
                IF ( MOD( TPF( S ), MTPRFAC ) .EQ. 0 .AND.
     &                    MOD( TPF( S ), WTPRFAC ) .EQ. 0       ) THEN

C.....................  update TMAT with met-based profiles (METFACS)
                    IF( MDEX( S,V ) == 99999 .OR.
     &                  WDEX( S,V ) == 99999 .OR.
     &                  DDEX( S,V ) == 99999 ) THEN

                        CALL UPDATE_TMAT_METFACS( METPROTYPE,.TRUE. )

                    ELSE

                        DO H = 1, 24

                            MON = MONTH( H, ZONES( S ) )
                            DAY = DAYOW( H, ZONES( S ) )
                            FAC = MONFAC( MON,MDEX( S,V ) ) * 
     &                            WEKFAC( DAY,WDEX( S,V ) )
                            K   = 1 + MOD( H + HCORR - ZONES( S ), 24 )
                            TMAT( S,V,H ) = FAC * HRLFAC( K, L, DAY )

                        END DO

                    END IF

C.................  This is when annual-data field is used for storing average-
C                   day emissions by multiplying it by 365.  Have to undo that
C                   here.
C.................  Adjust for week-normal data assuming whole week normalizer
                ELSE IF ( MOD( TPF( S ), WTPRFAC ) .EQ. 0 ) THEN
C.....................  update TMAT with met-based profiles (METFACS)
                    IF( WDEX( S,V ) == 99999 .OR.
     &                  DDEX( S,V ) == 99999 ) THEN
 
                        CALL UPDATE_TMAT_METFACS( METPROTYPE, .FALSE. )

                    ELSE

                        IF( MDEX( S,V ) == 99999 ) THEN
                            MESG = 'WARNING: Skip applying MONTHLY ' //
     &                          'temporal profiles from Gentpro '//
     &                          'to averag day inventory'
                            CALL M3EXIT( PROGNAME,JDATE,JTIME,MESG,2 )
                        END IF

                        DO H = 1, 24

                            DAY = DAYOW( H, ZONES( S ) )
                            FAC = YRFAC * WEKFAC( DAY,WDEX( S,V ) )
                            K   = 1 + MOD( H + HCORR - ZONES( S ), 24 )
                            TMAT( S,V,H ) = FAC * HRLFAC( K, L, DAY )

                        END DO

                    END IF

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

C******************  INTERNAL SUBPROGRAMS  *****************************

        CONTAINS

C.............  This internal subprogram will update TMAT array with
C               Met-based temporal profiles generated by Gentpro 

            SUBROUTINE UPDATE_TMAT_METFACS( PRFTYPE, IFLAG )


            CHARACTER(*), INTENT( IN ) :: PRFTYPE    ! met profile avg method
            LOGICAL,      INTENT( IN ) :: IFLAG      ! true: ann inv, false: avgday inv
 
C.............  Local variables
            INTEGER  H, I, J, K, L, NP, IOS 
            INTEGER  IYEAR, IDATE, MYEAR, MDAY, MHOUR

            REAL     METFAC, TOT, DIV, FAC
    
            CHARACTER(300) MESG
            CHARACTER(16)  TVARNAME
            
C----------------------------------------------------------------------

C............. Compute month and day of month using JDATE
            CALL DAYMON( JDATE, MON, MDAY )

C............. Search index of FIPS from METFACS
            NP = FINDC( CIFIP( S ), NMETPROF, METPROF )

            SELECT CASE( PRFTYPE )

            CASE( 'MONTHLY' )

C.................  Renormalize METFACS using days of month
                TOT = 0.0
                FAC = 0.0
                DO I = 1, 12
                    TOT = TOT + METFACS( NP, I, 1 )
                    FAC = FAC + MON_DAYS( I ) * METFACS( NP, I, 1 )
                END DO

                IF( FAC > 0.0 ) DIV = TOT / FAC 

                DO H = 1, 24

                    MON = MONTH( H, ZONES( S ) )
                    DAY = DAYOW( H, ZONES( S ) )

                    METFAC = DIV * METFACS( NP,MON,1 )
                    FAC = METFAC  * WEKFAC( DAY,WDEX( S,V ) )
                    L = DDEX( S,V )
                    K   = 1 + MOD( H + HCORR - ZONES( S ), 24 )
                    TMAT( S,V,H ) = FAC * HRLFAC( K, L, DAY )

                END DO
            
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

                METFAC = METFACS( NP,MON,MDAY )
                DO H = 1, 24

                    MON = MONTH( H, ZONES( S ) )
                    DAY = DAYOW( H, ZONES( S ) )

                    FAC = METFAC
                    L = DDEX( S,V )
                    K   = 1 + MOD( H + HCORR - ZONES( S ), 24 )
                    TMAT( S,V,H ) = FAC * HRLFAC( K, L, DAY )

                END DO

            CASE( 'HOURLY' )

C.................  Get header description of hour-specific temporal profile input file
                IF( .NOT. DESC3( PNAME ) ) THEN
                    CALL M3EXIT( PROGNAME, 0, 0,
     &                   'Could not get description of file "'
     &                   // PNAME( 1:LEN_TRIM( PNAME ) ) // '"', 2 )
                END IF

C.................  Convert timezone (w daylight saving) to read 
C                   correct county-specific hr profiles
                TDATE = JDATE
                TTIME = JTIME
                IF( TTIME > 230000 ) TTIME = TTIME - 240000

C.................  Read profileIDs and hourly factors
                IF( .NOT. READ3( PNAME, 'COUNTIES', 1, TDATE, TTIME,
     &                           METPROF ) ) THEN
                    MESG = 'Could not read ' // 'COUNTIES' //
     &                     ' from ' // TRIM( PNAME )
                    CALL M3EXIT( PROGNAME, TDATE, TTIME, MESG, 2 )
                END IF

                NP = FINDC( CIFIP( S ), NROWS3D, METPROF )

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
                            MDEX( S,V ) = 99999
                            WDEX( S,V ) = 99999

                        ELSE      ! convert annual to monthly before apply hour-of-month profiles from Gentpro
                            MON = MONTH( H, ZONES( S ) )
                            FAC = MONFAC_ORG( MON,MDEX( S,V ) )
                            TVARNAME    = 'MONTOT'
                            WDEX( S,V ) = 99999
                           
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
                        WDEX( S,V ) = 99999

                    END IF

                    IF( .NOT. READ3( PNAME, TVARNAME, 1, TDATE, TTIME,
     &                        METFACS( :,1,1 ) ) ) THEN 
                         MESG = 'Could not read ' // TRIM(TVARNAME) //
     &                          'from ' // TRIM( PNAME )
                         CALL M3EXIT( PROGNAME, TDATE, TTIME, MESG, 2 )
                    END IF
                    TOT = METFACS( NP,1,1 )

                    IF( .NOT. READ3( PNAME, 'HRLSRC', 1, TDATE, TTIME, 
     &                           METFACS( :,1,1 ) ) ) THEN
                        MESG = 'Could not read ' // TRIM( TVARNAME ) //
     &                         ' from ' // TRIM( PNAME )
                        CALL M3EXIT( PROGNAME, TDATE, TTIME, MESG, 2 )
                    END IF

                    
                    TMAT( S,V,H ) = FAC * METFACS( NP,1,1 ) / TOT

                    CALL NEXTIME( TDATE, TTIME, 10000 ) 

                 END DO
                
            END SELECT

            END SUBROUTINE UPDATE_TMAT_METFACS

C-------------------------------------------------------------------------
C-------------------------------------------------------------------------

        END SUBROUTINE MKTMAT

