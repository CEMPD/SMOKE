
        SUBROUTINE RDEPSPT( FDEV, NRAWIN, WKSET, INY, NRAWOUT, IOS, 
     &                      IREC, ERFILDSC, EFLAG, NDROP, EDROP )

C***********************************************************************
C  subroutine body starts at line 175
C
C  DESCRIPTION:
C      This subroutine reads the EPS2.0 AFS format importing point source 
C      inventory data.
C      It can be called multiple times for multiple sets of files
C
C  PRECONDITIONS REQUIRED:
C      Files must be opened and their unit numbers stored in EDEV() in the
C      order listed in the description.
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C      Subroutines: I/O API subroutines, BLDCSRC, CHECKMEM 
C      Functions: I/O API functions, GETFLINE, YR2DAY
C
C  REVISION  HISTORY:
C      Copied from rdemspt.f by M. Houyoux (3/2000)
C
C****************************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2001, MCNC--North Carolina Supercomputing Center
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
C...........   This module is the inventory arrays
        USE MODSOURC

C.........  This module contains the lists of unique inventory information
        USE MODLISTS

C.........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'CONST3.EXT'    !  physical constants
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS and their descriptions:
        LOGICAL                CHKINT
        LOGICAL                CHKREAL
        CHARACTER*2            CRLF
        INTEGER                ENVINT
        LOGICAL                ENVYN
        INTEGER                FIND1
        INTEGER                FINDC   !  returns -1 for failure
        INTEGER                GETFLINE
        INTEGER                INDEX1
        INTEGER                JULIAN 
        INTEGER                SECSDIFF
        INTEGER                STR2INT
        REAL                   STR2REAL
        INTEGER                YEAR4  
        REAL                   YR2DAY  

        EXTERNAL CHKINT, CHKREAL, CRLF, ENVINT, ENVYN, FIND1, FINDC,
     &           GETFLINE, INDEX1, JULIAN, SECSDIFF, STR2INT, STR2REAL, 
     &           YEAR4, YR2DAY

C...........   SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN) :: FDEV      ! input file unit no.
        INTEGER     , INTENT (IN) :: NRAWIN    ! total raw record-count
        INTEGER     , INTENT (IN) :: WKSET     ! weekly profile interpretation
        INTEGER , INTENT (IN OUT) :: INY       ! inv year for this set of files
        INTEGER     , INTENT(OUT) :: NRAWOUT   ! valid raw record-count
        INTEGER     , INTENT(OUT) :: IOS       ! I/O status
        INTEGER     , INTENT(OUT) :: IREC      ! line number
        CHARACTER(*), INTENT(OUT) :: ERFILDSC  ! file desc of file in error
        LOGICAL     , INTENT(OUT) :: EFLAG     ! error flag 
        INTEGER     , INTENT(OUT) :: NDROP     ! number of records dropped
        REAL        , INTENT(OUT) :: EDROP( MXIDAT ) ! emis dropped per pol

C.........  Local allocatable arrays
        INTEGER, ALLOCATABLE, SAVE :: SRTDCOD( : )
        INTEGER, ALLOCATABLE, SAVE :: SRTDIDX( : )

C.........  Temporary variables for storing source characteristics.  These
C           variables must be the width of the fields for global source
C           characteristics definition for use in BLDCSRC.
        CHARACTER(LEN=PLTLEN3) FCID  ! tmp plant ID
        CHARACTER(LEN=CHRLEN3) PTID  ! tmp point ID
        CHARACTER(LEN=CHRLEN3) SKID  ! tmp stack ID
        CHARACTER(LEN=CHRLEN3) SGID  ! tmp segment ID

C...........   Other local variables

        REAL             TK, HT, DM, VL  ! Temporary stack parms
        REAL             CEFF    !  tmp control effectiveness
        REAL             DAY2YR  !  Local, leap-year-able, DAY to YEAR factor
        REAL             EMIS    !  tmp emission value
        REAL             FSAV    !  Flow saved
        REAL             LAT     !  tmp latitude  (y-dir)
        REAL             LON     !  tmp longitude (x-dir)
        REAL             OZEMIS  !  tmp ozone season emission value
        REAL             REFF    !  tmp rule effectiveness
        REAL             RPEN    !  tmp rule penetration   
        REAL             XLOC    !  Scratch X-coordinate (UTM X or lon)
        REAL             YLOC    !  Scratch Y-coordinate (UTM Y or lat)
        REAL             YREMIS  !  tmp annual emission value

        INTEGER          I, J, L             ! counters and indices

        INTEGER          COD     !  Temporary pollutant code number
        INTEGER          CPRI    !  tmp primary control code
        INTEGER          DD      !  tmp day 
        INTEGER          EDT     !  tmp end date
        INTEGER          ES      !  counter for emission file
        INTEGER          ETM     !  tmp end time
        INTEGER          FIP     !  tmp state/county code
        INTEGER          HH      !  tmp hour
        INTEGER          ICC     !  position of CNTRY in CTRYNAM
        INTEGER          IYY     !  tmp emissions year
        INTEGER          MM      !  tmp month
        INTEGER, SAVE :: MXWARN  !  maximum number of warnings
        INTEGER          NDAYS   !  tmp no days
        INTEGER, SAVE :: NWARN =0!  number of warnings in this routine
        INTEGER, SAVE :: NSRCPOL = 0 ! cumulative source x data variables
        INTEGER          SDT     !  tmp start date
        INTEGER          SEGMENT !  tmp integer segment number
        INTEGER          STACK   !  tmp integer stack ID
        INTEGER          STM     !  tmp start time
        INTEGER          SIC     !  tmp standard industrial code
        INTEGER          TPF     !  Temporary temporal ID
        INTEGER          ZONE    !  Temporary UTM zone

        LOGICAL, SAVE :: CFLAG    !  true: recalculate the velocity from flow
        LOGICAL, SAVE :: FIRSTIME = .TRUE.
        LOGICAL, SAVE :: FIRSTP   = .TRUE.
        LOGICAL, SAVE :: WFLAG    !  true: convert lat-lons to Western hemisphr

        CHARACTER*2            TMPAA !  tmp time period code
        CHARACTER*300          LINE  !  Input line from POINT file
        CHARACTER*300          MESG  !  Text for M3EXIT()
        CHARACTER(LEN=IOVLEN3) CPOL  !  Temporary pollutant code
        CHARACTER(LEN=FIPLEN3) CFIP  !  Character FIP code
        CHARACTER(LEN=POLLEN3) CCOD  !  Character pollutant index to INVDNAM
        CHARACTER(LEN=DSCLEN3) DESC  !  plant description
        CHARACTER(LEN=SCCLEN3) TSCC  !  Temporary character SCC
        CHARACTER(LEN=CHRLEN3) CHAR4 !  tmp 4th plant characteristic

        CHARACTER*16 :: PROGNAME = 'RDEPSPT' ! Program name

C***********************************************************************
C   begin body of subroutine RDEPSPT

C.........  Set up settings the first time the subroutine is called
        IF( FIRSTIME ) THEN

C.............  Get settings from the environment
            MESG = 'Flag for recalculating velocity'
            CFLAG = ENVYN( 'VELOC_RECALC', MESG, .FALSE., IOS )

            MESG = 'Western hemisphere flag'
            WFLAG = ENVYN( 'WEST_HSPHERE', MESG, .TRUE., IOS )

            MXWARN = ENVINT( WARNSET, ' ', 100, IOS )

C.............  Create sorted data codes
            ALLOCATE( SRTDIDX( MXIDAT ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SRTDIDX', PROGNAME )
            ALLOCATE( SRTDCOD( MXIDAT ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SRTDCOD', PROGNAME )

            DO I = 1, MXIDAT
                SRTDIDX( I ) = I
            END DO

            CALL SORTI1( MXIDAT, SRTDIDX, INVDCOD )

            DO I = 1, MXIDAT
                J = SRTDIDX( I )
                SRTDCOD( I ) = INVDCOD( J )
            END DO

            FIRSTIME = .FALSE.

        END IF

C.........  Reinitialize for multiple subroutine calls
        ICC    = -9
        FIRSTP = .TRUE.

C........................................................................
C.............  Head of the FDEV-read loop  .............................
C........................................................................

        ES   = NSRCPOL
        IREC = 0
        TPF  = MTPRFAC * WKSET
        ERFILDSC = 'emission'
        DO

C.............  Read a line of input file and check input status

            READ( FDEV, 93000, END=199, IOSTAT=IOS ) LINE
            IREC = IREC + 1

            IF ( IOS .NE. 0 ) THEN

                EFLAG = .TRUE.
                WRITE( MESG, 94010 )
     &              'I/O error', IOS, 
     &              'reading inventory file at line', IREC
                CALL M3MESG( MESG )
                CYCLE

            END IF

            L = LEN_TRIM( LINE )  ! store width of line and check

C.............  Skip blank lines
            IF( L .EQ. 0 ) CYCLE

C.............  Scan for header lines and check to ensure all are set
C               properly.  Note that data value names are not read.
            CALL GETHDR( MXIDAT, .TRUE., .TRUE., .FALSE., 
     &                   LINE, ICC, INY, I, IOS )

C.............  Interpret error status
            IF( IOS .GT. 0 ) THEN
                EFLAG = .TRUE.
            END IF                

C.............  If a header line was encountered, go to next line
            IF( IOS .GE. 0 ) CYCLE

C.............  Define day to year conversion factor and real type for integer 
C               missing value
            DAY2YR  = 1. / YR2DAY( INY )

C.............  Set pollutant name
            CPOL = ADJUSTL( LINE( 156:160 ) )
            L    = LEN_TRIM( CPOL )

C.............  Ensure that pollutant field is not blank
C.............  If it is blank, move on to next record
            IF( CPOL .EQ. ' ' ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: Pollutant field ' //
     &                 'is missing at line', IREC
                CALL M3MESG( MESG )
                CYCLE

C.............  Determine if pollutant ID is supplied as integer or string
C.............  If supplied as integer, look up name
            ELSE IF( CHKINT( CPOL ) ) THEN
                COD = STR2INT( CPOL )
                J = FIND1( COD, MXIDAT, SRTDCOD )

C.................  Error if pollutant ID is not in master list
                IF( J .LE. 0 ) THEN
                    WRITE( MESG,94010 ) 'ERROR: Source dropped: ' //
     &                     'pollutant number "' // CPOL( 1:L ) //
     &                     '" at line', IREC, CRLF() // BLANK10 //
     &                     'of emission file is not in master ' //
     &                     'pollutants list'
                    CALL M3MESG( MESG )
                    CYCLE

C.................  Store pollutant position in master list and name
                ELSE
                    COD  = SRTDIDX( J )            ! position in original list
                    CPOL = INVDNAM( COD )
                END IF

C.............  Otherwise, make sure name is in list of valid pollutants
            ELSE
                COD = INDEX1( CPOL, MXIDAT, INVDNAM )

C..................  Error if pollutant name is not in master list
                IF( COD .LE. 0 ) THEN
                    WRITE( MESG,94010 ) 'ERROR: Source dropped: ' //
     &                     'pollutant "' // CPOL( 1:L ) //
     &                     '" at line', IREC, CRLF() // BLANK10 //
     &                     'of emission file is not in master ' //
     &                     'pollutants list'
                    CALL M3MESG( MESG )
                    CYCLE

                END IF 

            END IF

C.............  Make sure all integer fields are integers...
 
C.............  Check state/county codes, error for missing
            IF( .NOT. CHKINT( LINE( 12:16 ) ) .OR.
     &          LINE( 12:16 ) .EQ. ' '             ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: State and county ' //
     &                 'code is non-integer or missing at line', IREC
                CALL M3MESG( MESG )
            END IF

C.............  Check SIC code, warning for missing
            IF( .NOT. CHKINT( LINE( 24:27 ) ) ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: SIC code is non-' //
     &                 'integer at line', IREC
                CALL M3MESG( MESG )

            ELSE IF ( NWARN .LT. MXWARN .AND. 
     &                LINE( 24:27 ) .EQ. ' ' ) THEN
                WRITE( MESG,94010 ) 'WARNING: Missing SIC code at ' //
     &                 'line', IREC, '. Default 0000 will be used.'
                CALL M3MESG( MESG )
                NWARN = NWARN + 1

            END IF

C.............  Check primary control code
            IF( .NOT. CHKINT( LINE( 75:77 ) ) ) THEN

                EFLAG = .TRUE.
                L = LEN_TRIM( CPOL )
                WRITE( MESG,94010 ) 'ERROR: Primary control code' //
     &                     'for ' //
     &                     CPOL( 1:L ) // ' are non-integer ' //
     &                     'at line', IREC
                CALL M3MSG2( MESG )

            END IF

C.............  Make sure that all of the needed real values are real...

C.............  Stack height, diam, exit temperature, flow, & velocity
            IF( .NOT. CHKREAL( LINE( 103:107 ) ) .OR.
     &          .NOT. CHKREAL( LINE( 109:113 ) ) .OR.
     &          .NOT. CHKREAL( LINE( 115:119 ) ) .OR.
     &          .NOT. CHKREAL( LINE( 121:125 ) )      ) THEN

                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: Stack parameters are ' //
     &                 'not numbers or have bad formatting' // CRLF() //
     &                 BLANK10 // 'at line', IREC
                CALL M3MSG2( MESG )

            END IF

C.............  Check stack coordinates, missing is an error
            IF( .NOT. CHKREAL( LINE( 78:87 ) ) .OR.
     &          .NOT. CHKREAL( LINE( 89:98 ) )      ) THEN

                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: latitude and/or ' //
     &                 'longitude are not numbers or have ' // CRLF() //
     &                 BLANK10 // 'bad formatting at line', IREC
                CALL M3MESG( MESG )

            ELSE IF( LINE( 78:87 ) .EQ. ' ' .OR.
     &               LINE( 89:98 ) .EQ. ' '      ) THEN
      
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: latitude and/or ' //
     &                 'longitude are missing at line', IREC
                CALL M3MESG( MESG )

            END IF

C.............  Check emissions values
            IF( .NOT. CHKREAL( LINE( 162:171 ) ) ) THEN
                EFLAG = .TRUE.
                L = LEN_TRIM( CPOL )
                WRITE( MESG,94010 ) 'ERROR: Emissions data value for '//
     &                 CPOL( 1:L ) // ' is not a number ' // CRLF()//
     &                 BLANK10 // 'or has bad formatting at line', IREC
                CALL M3MESG( MESG )

            ELSE IF( LINE( 162:171 ) .EQ. ' ' ) THEN
                L = LEN_TRIM( CPOL )
                WRITE( MESG,94010 ) 'WARNING: Emissions data for ' //
     &                 CPOL( 1:L ) // ' are missing at line', IREC
                CALL M3MESG( MESG )
                LINE( 162:171 ) = '0.'  ! to prevent STR2REAL warning

            END IF

C.............  Check emissions-associated values
            IF( .NOT. CHKREAL( LINE( 177:182 ) ) .OR.
     &          .NOT. CHKREAL( LINE( 184:189 ) ) .OR.
     &          .NOT. CHKREAL( LINE( 191:196 ) )      ) THEN

                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: control efficiency, ' //
     &                 'rule effectiveness, and/or rule penetration ' // 
     &                 CRLF() //BLANK10 // 
     &                 'have bad formatting at line', IREC
                CALL M3MESG( MESG )
            END IF

C.............  If there has been an error, do not try to store any of the
C               records.  Instead  go to next line of file.
            IF( EFLAG ) CYCLE
       
C.............  Now use the file format definition to parse the LINE into
C               the various data fields...

            FIP  = ICC * 100000 + STR2INT( LINE( 12:16 ) )
            FCID = ADJUSTL ( LINE(  40: 44 ) )  ! plant ID
            SKID = ADJUSTL ( LINE(  45: 48 ) )  ! stack ID
            PTID = ADJUSTL ( LINE(  50: 52 ) )  ! point ID
            SGID = ADJUSTL ( LINE(  54: 55 ) )  ! segment ID
            DESC = ADJUSTL ( LINE( 200:240 ) )  ! plant description
            TSCC = ADJUSTL ( LINE(  29: 38 ) )  ! SCC code

            SIC  = MAX( STR2INT ( LINE(  24: 27 ) ), 0 )
            CPRI = MAX( STR2INT ( LINE(  75: 77 ) ), 0 )

            HT   = STR2REAL( LINE( 103:107 ) )
            DM   = STR2REAL( LINE( 109:113 ) )
            TK   = STR2REAL( LINE( 115:119 ) )
            VL   = STR2REAL( LINE( 121:125 ) )
            ZONE = STR2INT ( LINE( 100:101 ) )
            YLOC = STR2REAL( LINE(  78: 87 ) )
            XLOC = STR2REAL( LINE(  89: 98 ) )

            CEFF  = STR2REAL( LINE( 177:182 ) )
            IF ( CEFF .LT. 0.0 )  CEFF = 0.0

            REFF  = STR2REAL( LINE( 184:189 ) )
            IF ( REFF .LT. 0.0 )  REFF = 100.0

            RPEN  = STR2REAL( LINE( 191:196 ) )
            IF ( RPEN .GE. 0.0 .AND. RPEN .NE. 100. )  THEN
                WRITE( MESG,94010 ) 'WARNING: Rule penetration ' //
     &                 'at line', IREC, 'for point sources is not'//
     &                 CRLF()// BLANK10// 'supported by SMOKE.'
                CALL M3MESG( MESG )
                
            END IF

            EMIS = STR2REAL( LINE( 162:171 ) )

            SDT  = STR2INT( LINE( 60:67 ) )
            EDT  = STR2INT( LINE( 69:76 ) )

C.............  Make adjustments to pad with zeros, if needed
            WRITE( CFIP,94120 ) FIP
            CALL PADZERO( CFIP )
            CALL PADZERO( TSCC )

C.............  Convert stack and segment to integer and then character to
C               remove leading zeros
            IF( CHKINT( SKID ) .AND. SKID .NE. ' ' ) THEN
                STACK   = STR2INT( SKID )
                WRITE( SKID, '(I4)' ) STACK
            END IF

            IF( CHKINT( SGID ) .AND. SGID .NE. ' ' ) THEN
                SEGMENT = STR2INT( SGID )
                WRITE( SGID, '(I2)' ) SEGMENT
            END IF

C.............  Store latitude and longitude, and make conversions from UTM,
C               if needed...

C.............  Data already provided in lat/lon
            IF( ZONE .LE. 0 ) THEN

C.................  Convert longitude if WEST_HSPHERE environment
C.................  variable is true, and longitude is negative.
                IF( WFLAG .AND. XLOC .GT. 0 ) THEN

                    XLOC = -XLOC

                END IF
 
        	LON = XLOC
        	LAT = YLOC

C.............  Convert to lat/lon from UTM
            ELSE
                XLOC = 1000.0 * XLOC                     ! km to meters
                YLOC = 1000.0 * YLOC                     ! km to meters
                CALL UTM2LL( XLOC, YLOC, ZONE, LON, LAT )

            END IF

C.............  Check and set time period type (Year/day/hourly)
            TMPAA = ADJUSTL( LINE( 57:58 ) )
            CALL UPCASE( TMPAA )

C.............  If interval indicator is blank, emissions are annual total
            IF ( TMPAA .EQ. '  ' ) THEN

                TPF = MTPRFAC * WKSET           !  use month, week profiles

C.............  Emissions are peak day 
C.............  NOTE- emissions will be stored as tons/day

            ELSE IF ( TMPAA( 1:1 ) .EQ. 'P' ) THEN

                HH  = MOD( SDT        , 100 )
                STM = 10000 * HH
                DD  = MOD( SDT / 100  , 100 )
                MM  = MOD( SDT / 10000, 100 )
                IYY  = YEAR4( SDT / 1000000 )
                SDT = IYY*1000 + JULIAN( IYY, MM, DD )  !  as Julian date

                HH  = MOD( EDT        , 100 )
                ETM = 10000 * HH
                DD  = MOD( EDT / 100  , 100 )
                MM  = MOD( EDT / 10000, 100 )
                IYY  = YEAR4( EDT / 1000000 )
                EDT = IYY*1000 + JULIAN( IYY, MM, DD )  !  as Julian date
                CALL NEXTIME( EDT, ETM, 10000 ) ! ETM means *through* hr HH

                NDAYS  = SECSDIFF( SDT, STM, EDT, ETM ) * SEC2DAY
                IF ( NDAYS .GT. 28 ) THEN
                    TPF = MTPRFAC * WKSET       !  use month, week profiles
                ELSE IF ( NDAYS .GT. 1 ) THEN
                    TPF = WKSET                 !  use only week profiles
                ELSE
                    TPF = 1                     !  use only hourly profiles
                END IF                                

                IF( FIRSTP .AND. INY .NE. IYY ) THEN
                    FIRSTP = .FALSE.
                    WRITE( MESG,94010 ) 'WARNING: Year', INY,
     &                     'given in header, but year', IYY,
     &                     'appears in inventory input file'
                    CALL M3MESG( MESG )
                END IF

C.............  Emissions are over a special interval
            ELSE IF ( TMPAA( 1:1 ) .EQ. 'S' ) THEN     !  special interval

                HH  = MOD( SDT        , 100 )
                STM = 10000 * HH
                DD  = MOD( SDT / 100  , 100 )
                MM  = MOD( SDT / 10000, 100 )
                IYY  = YEAR4( SDT / 1000000 )
                SDT = IYY*1000 + JULIAN( IYY, MM, DD )  !  as Julian date

                HH  = MOD( EDT        , 100 )
                ETM = 10000 * HH
                DD  = MOD( EDT / 100  , 100 )
                MM  = MOD( EDT / 10000, 100 )
                IYY  = YEAR4( EDT / 1000000 )
                EDT = IYY*1000 + JULIAN( IYY, MM, DD )  !  as Julian date
                CALL NEXTIME( EDT, ETM, 10000 ) ! ETM means *through* hr HH

                NDAYS  = SECSDIFF( SDT, STM, EDT, ETM ) * SEC2DAY
                IF ( NDAYS .GT. 28 ) THEN
                    TPF = MTPRFAC * WKSET       !  use month, week profiles
                ELSE IF ( NDAYS .GT. 1 ) THEN
                    TPF = WKSET                 !  use only week profiles
                ELSE
                    TPF = 1                     !  use only hourly profiles
                END IF

                EMIS = 0.0      !  note ZERO for day-specific "SP" records

                IF( NWARN .LT. MXWARN ) THEN
                    WRITE( MESG, 94010 ) 'WARNING: Period-specific ' //
     &                 'emissions being set to 0.0 in inventory ' //
     &                 'at line'// CRLF() // BLANK10, IREC,
     &                 '. Actual emissions will be used later ' //
     &                 'by the Temporal program.'
                    CALL M3MESG( MESG )
                    NWARN = NWARN + 1
                END IF

                IF( FIRSTP .AND. INY .NE. IYY ) THEN
                    FIRSTP = .FALSE.
                    WRITE( MESG,94010 ) 'WARNING: Year', INY,
     &                         'given in header, but year', IYY,
     &                         'appears in inventory input file'
                    CALL M3MESG( MESG )
                END IF

            ELSE                                        !  unrecognized type

                WRITE( MESG,94010 ) 
     &                  'ERROR: Unsupported time period type "' // 
     &                  TMPAA // '" at line', IREC, 
     &                  ' in inventory input file' 
                CALL M3MESG( MESG )
                EFLAG = .TRUE.
                CYCLE          !  to head of PDEV-read loop

            END IF          !  tests on record type line( 57:58 )

C.............  Set annual or ozone-season value, depending on type of data 
C               available.
            IF( TMPAA .EQ. 'PO' ) THEN
                YREMIS = 0.
                OZEMIS = EMIS
            ELSE
                YREMIS = EMIS
                OZEMIS = 0.
            END IF
  
C.............  Time to store data in unsorted lists if we've made it this far
            ES = ES + 1

            IF ( ES .LE. NRAWIN ) THEN

                IFIPA  ( ES ) = FIP
                ISICA  ( ES ) = SIC
                TPFLGA ( ES ) = TPF
                INVYRA ( ES ) = INY
                STKHTA ( ES ) = HT
                STKDMA ( ES ) = DM
                STKTKA ( ES ) = TK
                STKVEA ( ES ) = VL
                XLOCAA ( ES ) = LON
                YLOCAA ( ES ) = LAT
                CSCCA  ( ES ) = TSCC
                CORISA ( ES ) = ORSBLNK3
                CBLRIDA( ES ) = BLRBLNK3
                CPDESCA( ES ) = DESC
                POLVLA ( ES,NEM ) = YREMIS
                POLVLA ( ES,NOZ ) = OZEMIS
                POLVLA ( ES,NCE ) = CEFF
                POLVLA ( ES,NRE ) = REFF
                POLVLA ( ES,NC1 ) = CPRI

                WRITE( CCOD,94125 ) COD
 
                CHAR4 = TSCC
                CALL BLDCSRC( CFIP, FCID, SKID, PTID, SGID, 
     &                        CHAR4, CHRBLNK3, CCOD, CSOURCA( ES ) )

            END IF          !  if S in range

        END DO           !  to head of FDEV-read loop

199     CONTINUE        !  end of the FDEV-read loop

        CLOSE( FDEV )

        WRITE( MESG,94010 ) 
     &         'EMISSION FILE processed:'  // CRLF() // BLANK5 //
     &         '   This-file  EPS2.0 SOURCE  record-count', ES-NSRCPOL,
     &         CRLF() // BLANK5 //
     &         '   Cumulative EPS2.0 SOURCE  record-count', ES

        CALL M3MSG2( MESG )

        NSRCPOL = ES        !  cumulative emissions lines

C........................................................................

        IF( NSRCPOL .GT. NRAWIN ) THEN

            EFLAG = .TRUE.
            MESG = 'Memory allocation insufficient for EPS2.0 inventory'
            CALL M3MSG2( MESG )

        ELSE
            NRAWOUT = NSRCPOL

        END IF		!  if overflow or if errors

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94120   FORMAT( I6.6 )

94125   FORMAT( I5 )

94300   FORMAT( A, I2.2, A, I2.2, A )

        END SUBROUTINE RDEPSPT
