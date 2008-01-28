
        SUBROUTINE RDCOMBO( CDEV, ENAM )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C       This program 
C
C  REVISION  HISTORY:
C     Created 12/07 by M. Houyoux
C
C***************************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C 
C smoke@unc.edu
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C***************************************************************************

C...........   MODULES for public variables   
C...........   This module contains the speciation profile tables
        USE MODSPRO, ONLY: CMBNP, CMBSPCD, CMBWGHT, CMBMAX

C.........  This module contains the lists of unique source characteristics
        USE MODLISTS, ONLY: NINVIFIP, INVCFIP

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  i/o api constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER   CRLF
        INTEGER     ENVINT 
        INTEGER     FINDC 
        INTEGER     GETFLINE
        INTEGER     STR2INT
        LOGICAL     BLKORCMT
        REAL        STR2REAL  

        EXTERNAL    CRLF, ENVINT, FINDC, GETFLINE, STR2INT, STR2REAL
     &              BLKORCMT
 
C.........  SUBROUTINE ARGUMENTS
        INTEGER     , INTENT    (IN) :: CDEV    ! unit number of input file
        CHARACTER(*), INTENT    (IN) :: ENAM    ! pol/emis type name of interest

C.........  Local parameters
        INTEGER, PARAMETER :: STATETYP = 1
        INTEGER, PARAMETER :: CNTYTYP  = 2

C.........  Local allocatable arrays
        INTEGER, ALLOCATABLE :: CMBTYP( : )  ! type of last entry applied to county (state- or county-specific)

C.........  Local arrays
        INTEGER, SAVE :: ERRCNT( 6 )  
        INTEGER, SAVE :: WARNCNT( 2 )  
        REAL    :: CWEIGHT( CMBMAX )          ! tmp array for profile weights
        CHARACTER(POLLEN3) :: SEGMENT( 24 )  ! Segments of parsed lines

C.........  Other local variables
        INTEGER         F, I, N       !  counters and indices

        INTEGER         FMISS   !  tmp count of FIPs codes missing records in GSPRO_COMBO
        INTEGER         IOS     !  i/o status
        INTEGER         IREC    !  record counter
        INTEGER         NLINES  !  number of lines
        INTEGER         PERIOD  !  period to use to match to GSPRO_COMBO file
        INTEGER         PER     !  tmp period from input file
        INTEGER, SAVE :: MXERR                ! max no. errors
        INTEGER, SAVE :: MXWARN               ! max no. warnings
        INTEGER         NP      !  tmp number of profiles per combo
        INTEGER         NRECS   !  number of records read in current call

        LOGICAL       :: EFLAG    = .FALSE. ! true: error detected
        LOGICAL, SAVE :: FIRSTIME = .TRUE.

        CHARACTER(FIPLEN3) CFIP     !  tmp buffer for state/county FIPS code
        CHARACTER(STALEN3) CSTA     !  tmp buff for state FIPS code
        CHARACTER(STALEN3) PSTA     !  tmp buff for previous state FIPS code
        CHARACTER(POLLEN3) CPOL     !  tmp buffer for pollutant code
        CHARACTER(256)     MESG     !  message buffer
        CHARACTER(512)     LINE     !  line buffer

        CHARACTER(16) :: PROGNAME = 'RCOMBO' ! program name

C***********************************************************************
C   begin body of subroutine RDCOMBO

C.........  Perform one-time steps
        IF ( FIRSTIME ) THEN

            MXWARN = ENVINT( WARNSET , ' ', 100, I )
            MXERR  = ENVINT( ERRSET  , ' ', 100, I )

            ERRCNT  = 0  ! Array
            WARNCNT = 0  ! Array

            FIRSTIME = .FALSE.

        END IF

C.........   Evaluate environment variables and store in saved variable
        MESG = 'Period to read from GSPRO_COMBO file'
        PERIOD = ENVINT( 'SPCMAT_PERIOD', MESG, 1, IOS )

C.........  Allocate public arrays if not already allocated. This just needs to be done
C           once, since allocating for all FIPS codes (and not dependent on file size)
        IF ( .NOT. ALLOCATED( CMBNP ) ) THEN

            ALLOCATE( CMBNP( NINVIFIP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CMBNP', PROGNAME )

            ALLOCATE( CMBSPCD( NINVIFIP, CMBMAX ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CMBSPCD', PROGNAME )

            ALLOCATE( CMBWGHT( NINVIFIP, CMBMAX ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CMBWGHT', PROGNAME )

        END IF

C.........  Allocate local array if not already allocated
        IF ( .NOT. ALLOCATED( CMBTYP ) ) THEN
            ALLOCATE( CMBTYP( NINVIFIP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CMBTYP', PROGNAME )
        END IF

C.........  Initialize public and local arrays to 0
        CMBNP   = 0   ! array
        CMBSPCD = ' ' ! array
        CMBWGHT = 0.  ! array
        CMBTYP  = 0   ! array

C.........  Other initializations
        PSTA = '-9'

C.........  Get the number of lines in the input file
        NLINES = GETFLINE( CDEV, 'GSPRO Combos file' )

C.........  Loop through file and read until the end.
        IREC = 0
        NRECS = 0
        DO I = 1, NLINES

            READ ( CDEV, 93000, END=999, IOSTAT=IOS ) LINE
            IREC = IREC + 1

            IF ( IOS .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 )
     &              'I/O error', IOS,
     &              'reading profile combinations file at line', IREC
                CALL M3MESG( MESG )
                CYCLE
            END IF

C.............  Skip blank and comment lines
            IF( BLKORCMT( LINE ) ) CYCLE

C.............  Split out columns from line
            CALL PARSLINE( LINE, 24, SEGMENT )

C.............  Get period, pollutant, and number of profiles from line
            CPOL = SEGMENT( 1 )
            PER  = STR2INT( SEGMENT( 3 ) )
            NP   = STR2INT( SEGMENT( 4 ) )

C.............  If period or pollutant doesn't match, or if number
C               of profiles is <=0, then skip
            IF ( CPOL .NE. ENAM .OR. PER .NE. PERIOD .OR. NP <=0 ) CYCLE

C.............  Get FIPs code and pad with leading zeros, if needed
            CFIP = SEGMENT( 2 )
            CALL PADZERO( CFIP )

C.............  Get tmp value for number of profiles

C.............  If state-specific entry (county set to 000), then 
C               loop through all FIPs codes and apply information.
            IF ( CFIP(STALEN3+1:FIPLEN3) == '000' ) THEN

C.................  Extract state code from state/county FIPS code
                CSTA = CFIP(1:STALEN3)

                DO F = 1, NINVIFIP   ! Loop through inventory state/county FIPS codes

C.....................  If inventory state code matches COMBO record state,
C                       the store information for current county
                    IF( INVCFIP(F)(1:STALEN3) .EQ. CSTA ) THEN

C.........................  Give error if any county-specific entries already 
C                           applied in the input file
                        IF( CMBTYP( F ) == CNTYTYP .AND. 
     &                      ERRCNT(1) <=  MXERR          ) THEN

                            WRITE( MESG,94010 ) 'ERROR: State-specific '
     &                        //'record at line',IREC,'comes after '//
     &                        'county-specific record matching same '//
     &                        'source.'
                            CALL M3MESG( MESG )
                            ERRCNT(1) = ERRCNT(1) + 1
                            EFLAG = .TRUE.
                            CYCLE

                        END IF
            
C.........................  Give error if duplicate (non-zero) entry
                        IF ( CMBTYP( F ) == STATETYP .AND.
     &                       CSTA .NE. PSTA .AND.
     &                       ERRCNT(2) <= MXERR        ) THEN

                            WRITE( MESG,94010 )
     &                        'ERROR: Duplicate entry found at line', 
     &                        IREC,'for state "' // CSTA //
     &                        '" of GSPRO_COMBO file'
                            CALL M3MESG( MESG )
                            ERRCNT(2) = ERRCNT(2) + 1
                            EFLAG = .TRUE.
                            CYCLE

                        END IF

C.........................  Record  flag for state-specific entry
                        CMBTYP( F ) = STATETYP

C.........................  Store entry for current FIPs, pollutant, period
                        CMBNP( F ) = NP

C.........................  If on a new state (not just looping over all counties)
                        IF( CSTA .NE. PSTA ) THEN

C.............................  Convert fractions from strings to reals
                            DO N = 1,  NP
                                CWEIGHT( N ) = STR2REAL( SEGMENT( 4+N*2 ) )
                            END DO

C.............................  Check if profile fractions meet the +/- 0.001 criterion
C                               and renormalize if needed. Provide a warning if need  
C                               to renormalize.
                            CALL CHECK_AND_SET_FRACS

                        END IF

                        DO N = 1, NP

                            CMBSPCD( F,N )= SEGMENT( 3+N*2 )
                            CMBWGHT( F,N )= CWEIGHT( N )

                        END DO
                    
                    END IF

C.....................  Set previous state code for next iteration
                    PSTA = CSTA

                END DO

C.............  Otherwise, apply entry to a single county
            ELSE

C.................  Search for state/county FIPS code in list of these from
C                   the inventory. If not found, then skip this line.
                F = FINDC( CFIP, NINVIFIP, INVCFIP )
                IF ( F .LE. 0 ) CYCLE

C.................  Give error if duplicate (non-zero) entry
                IF ( CMBTYP( F ) == CNTYTYP .AND. 
     &               ERRCNT(3)   <= MXERR         ) THEN
                    WRITE( MESG,94010 )
     &                'ERROR: Duplicate entry found at line', IREC,
     &                'for county "'//CFIP//'" of GSPRO_COMBO file'
                    CALL M3MESG( MESG )
                    ERRCNT(3) = ERRCNT(3) + 1
                    EFLAG = .TRUE.
                    CYCLE
                END IF

C.................  Record  flag for county-specific entry
                CMBTYP( F ) = CNTYTYP

C.................  Store entry for current FIPs, pollutant, period
                CMBNP( F ) = NP

C.................  Convert fractions from strings to reals
                DO N = 1,  NP
                    CWEIGHT( N ) = STR2REAL( SEGMENT( 4+N*2 ) )
                END DO

C.................  Check if profile fractions meet the +/- 0.001 criterion
C                   and renormalize if needed. Provide a warning if need  
C                   to renormalize.
                CALL CHECK_AND_SET_FRACS

                DO N = 1, NP

                    CMBSPCD( F,N ) = SEGMENT( 3 + N*2 )
                    CMBWGHT( F,N ) = CWEIGHT( N )

                END DO

            END IF

            NRECS = NRECS + 1

        END DO

C.........  Check that at least one record present for the selected PERIOD
C.........  Note: want this always to be reported, regardless of MAXERROR
        IF( NRECS == 0 ) THEN

            WRITE( MESG, 94010 ) 'ERROR: No records found for period',
     &             PERIOD, 'and pollutant '//TRIM( ENAM ) // ' in '//
     &             'GSPRO_COMBO file.'
            CALL M3MESG( MESG )
            EFLAG = .TRUE.

C.........  Check that all counties were set with either the state or county
C           approach.
        ELSE

            FMISS = 0
            DO F = 1, NINVIFIP

C.................  Write error message for counties not found
                IF ( CMBTYP( F ) .EQ. 0 .AND. ERRCNT(5) <= MXERR ) THEN
                    FMISS = FMISS + 1
                    WRITE( MESG,94010 ) 'ERROR: State/County FIPS '//
     &                     INVCFIP( F )//' is missing data in the '//
     &                     'GSPRO_COMBO file'//CRLF()//BLANK10//
     &                     'for period', PERIOD, 'and pollutant '//
     &                     TRIM( ENAM )
                    CALL M3MESG( MESG )
                    ERRCNT(5) = ERRCNT(5) + 1
                    EFLAG = .TRUE.
                END IF

            END DO

C.............  Write additional message with total number of missing counties
C.............  Note: want this always to be reported, regardless of MAXERROR
            IF( FMISS > 0 ) THEN
                WRITE( MESG,94010 ) 'ERROR:',FMISS, 'State/County '//
     &                 'FIPS codes out of ',NINVIFIP, 'are missing'//
     &                 'entries in the GSPRO_COMBO file'//CRLF()//
     &                 BLANK10// 'for period', PERIOD,'and pollutant '//
     &                 TRIM( ENAM )
                CALL M3MESG( MESG )
            END IF

        END IF

C.........  Deallocate local memory
        IF( ALLOCATED( CMBTYP ) ) DEALLOCATE( CMBTYP )

C.........  If error found, abort program
        IF( EFLAG ) THEN
            MESG = 'Problem reading profile combos file GSPRO_COMBO'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF 

        RETURN

C.........  Error message for reaching the end of file too soon
999     MESG = 'End of file reached unexpectedly. ' //
     &         'Check format of GSPRO_COMBO file.'

        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

C******************  INTERNAL SUBPROGRAMS  *****************************

        CONTAINS

C.............  This internal subroutine writes the message when a default
C               speciation profile is unavailable for a given pollutant
            SUBROUTINE CHECK_AND_SET_FRACS

C................. Inherits the following variables:
C               CWEIGHT
C               NP
C               IREC

C.............  Local parameters
            REAL, PARAMETER :: TOLERANCE = 0.001

C.............  Local variables
            INTEGER   N
            REAL      SUM    ! sum of weights
            REAL      SUMINV ! inverse of SUM

C----------------------------------------------------------------------
            SUM = 0.0
            DO N = 1,  NP
                SUM = SUM + CWEIGHT( N )
            END DO

            IF( SUM > 1+TOLERANCE .OR. SUM < 1-TOLERANCE ) THEN

                SUMINV = 1. / SUM

C.................  Renormalize fractions
                DO N = 1, NP
                    CWEIGHT( N ) = CWEIGHT( N ) * SUMINV
                END DO

C.................  Send warning message that entry is being renormalized
                IF( SUM > 1+TOLERANCE .AND. WARNCNT(1) <= MXWARN ) THEN
                    WRITE( MESG,94010 ) 'WARNING: GSPRO_COMBO '//
     &                 'fractions summed > 1.001 at line',IREC,
     &                 'and were renormalized'
                    CALL M3MESG( MESG )
                    WARNCNT(1) = WARNCNT(1) + 1

                ELSE IF ( SUM < 1-TOLERANCE .AND. 
     &                    WARNCNT(2) <= MXWARN    ) THEN
                    WRITE( MESG,94010 ) 'WARNING: GSPRO_COMBO '//
     &                 'fractions summed < 0.999 at line',IREC,
     &                 'and were renormalized'
                    CALL M3MESG( MESG )
                    WARNCNT(2) = WARNCNT(2) + 1
                END IF

            END IF

C------------------- SUBPROGRAM FORMAT STATEMENTS ----------------------

C...............   Internal buffering formats............ 94xxx

94010       FORMAT( 10( A, :, I8, :, 1X ) )

            END SUBROUTINE CHECK_AND_SET_FRACS

        END SUBROUTINE RDCOMBO
