
        SUBROUTINE RDEMSPT( EDEV, INY, NRAWIN, WKSET, NRAWOUT, IOS, 
     &                      IREC, ERFILDSC, EFLAG, NDROP, EDROP )

C***********************************************************************
C  subroutine body starts at line 232
C
C  DESCRIPTION:
C      This subroutine reads the EMS-95 format for one set of 5 files
C      (device.pt, emission.pt, facility.pt, stack.pt, process.pt)
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
C      Copied from emspoint.F by M. Houyoux (10/98)
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
        CHARACTER*2            CRLF
        LOGICAL                ENVYN
        INTEGER                FINDC   !  returns -1 for failure
        INTEGER                GETFLINE
        INTEGER                INDEX1
        INTEGER                LBLANK
        INTEGER                STR2INT
        REAL                   STR2REAL
        REAL                   YR2DAY  

        EXTERNAL CRLF, ENVYN, FINDC, GETFLINE, INDEX1, LBLANK, 
     &           STR2INT, STR2REAL, YR2DAY

C...........   SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN) :: EDEV( 5 ) !  unit no.: dv, em, fc, sk, pr
        INTEGER     , INTENT (IN) :: INY       !  inv year for this set of files
        INTEGER     , INTENT (IN) :: NRAWIN    !  total raw record-count
        INTEGER     , INTENT (IN) :: WKSET     !  weekly profile interpretation
        INTEGER     , INTENT(OUT) :: NRAWOUT   ! valid raw record-count
        INTEGER     , INTENT(OUT) :: IOS       ! I/O status
        INTEGER     , INTENT(OUT) :: IREC      ! line number
        CHARACTER(*), INTENT(OUT) :: ERFILDSC  ! file desc of file in error
        LOGICAL     , INTENT(OUT) :: EFLAG     ! error flag 
        INTEGER    , INTENT(INOUT):: NDROP     ! number of records dropped
        REAL       , INTENT(INOUT):: EDROP( MXIDAT ) ! emis dropped per pol

C...........   Local parameters, indpendent
        INTEGER, PARAMETER :: DVIDLEN = 12
        INTEGER, PARAMETER :: EDSCLEN = 40
        INTEGER, PARAMETER :: FCIDLEN = 15
        INTEGER, PARAMETER :: SKIDLEN = 12
        INTEGER, PARAMETER :: PRIDLEN = 12

C...........   Local parameters, dependent
        INTEGER, PARAMETER :: FCKYLEN = FIPLEN3 + FCIDLEN
        INTEGER, PARAMETER :: SKKYLEN = FCKYLEN + SKIDLEN
        INTEGER, PARAMETER :: DVKYLEN = SKKYLEN + DVIDLEN
        INTEGER, PARAMETER :: PRKYLEN = DVKYLEN + PRIDLEN

C...........   Point sources table for facility.pt files (unsorted; sorted)

        INTEGER                FS           ! Ctr & number of fc records
        INTEGER, ALLOCATABLE:: INDXFCA( : ) ! subscript array for sort
        INTEGER, ALLOCATABLE:: FCZONEA( : ) ! UTM zone from facility files
        REAL,    ALLOCATABLE:: FCCRDXA( : ) ! Facility X coordinates
        REAL,    ALLOCATABLE:: FCCRDYA( : ) ! Facility Y coordinates
        CHARACTER(LEN=FCKYLEN), ALLOCATABLE:: FCKEYA ( : ) ! FIP // facility ID
        CHARACTER(LEN=EDSCLEN), ALLOCATABLE:: FCDESCA( : ) ! facility name 

        INTEGER,                ALLOCATABLE:: FCIDX( : ) ! index to sort lists
        CHARACTER(LEN=FCKYLEN), ALLOCATABLE:: FCKEY( : ) ! FIP // facility ID
        CHARACTER(LEN=FCKYLEN)                FKEY       ! Tmp facility key
        CHARACTER(LEN=FCIDLEN)                FCID       ! Tmp facility ID

C...........   Point sources table for stack.pt files  (unsorted; sorted)

        INTEGER                SS           ! Ctr & number of sk records
        INTEGER, ALLOCATABLE:: INDXSKA( : ) ! subscript table for sort
        INTEGER, ALLOCATABLE:: IFCKEYA( : ) ! substript for facility arrays
        REAL,    ALLOCATABLE:: SKDIAMA( : ) ! stack diam from sk files [m]
        REAL,    ALLOCATABLE:: SKHEITA( : ) ! stack height from stack files [m]
        REAL,    ALLOCATABLE:: SKTEMPA( : ) ! stack exit temp from sk files [C]
        REAL,    ALLOCATABLE:: SKCRDXA( : ) ! stack UTM easting coord [m]
        REAL,    ALLOCATABLE:: SKCRDYA( : ) ! stack UTM northing coord [m]
        REAL,    ALLOCATABLE:: SKVELOA( : ) ! stack exit velocity [m/sec]
        REAL,    ALLOCATABLE:: SKFLOWA( : ) ! sk flow when not match v [m^3/sec]
        CHARACTER(LEN=SKKYLEN), ALLOCATABLE:: SKKEYA( : ) ! FIP // FCID // SKID

        INTEGER,                ALLOCATABLE:: SKIDX( : ) ! index to sort list
        CHARACTER(LEN=SKKYLEN), ALLOCATABLE:: SKKEY( : ) ! FIP // FCID // SKID
        CHARACTER(LEN=SKKYLEN)                SKEY       ! Tmp stack key
        CHARACTER(LEN=SKIDLEN)                SKID       ! Tmp stack ID

C...........   Point sources table for process.pt files  (unsorted; sorted)

        INTEGER                PS           ! Ctr & number of pr records
        INTEGER, ALLOCATABLE:: INDXPA( : )  ! subscript table
        INTEGER, ALLOCATABLE:: PSSCCA( : )  ! SCC code from process files
        CHARACTER(LEN=PRKYLEN), ALLOCATABLE:: PSKEYA( : )! FIP//FCID//SKID//DVID//PRID

        INTEGER,                ALLOCATABLE:: PSSCC ( : )! SCC from pt files
        CHARACTER(LEN=PRKYLEN), ALLOCATABLE:: PSKEY( : ) ! FIP//FCID//SKID//DVID//PRID
        CHARACTER(LEN=PRKYLEN)                PKEY       ! Tmp process key
        CHARACTER(LEN=PRIDLEN)                PRID       ! Tmp process ID

C...........   Point sources table for device.pt files   (unsorted; sorted)

        INTEGER                DS ! Counter & number of device records
        INTEGER, ALLOCATABLE:: INDXDVA( : ) ! subscript table for sort
        INTEGER, ALLOCATABLE:: DVSICA ( : ) ! SIC code from device files
        INTEGER, ALLOCATABLE:: DVIWEKA( : ) ! Weekly prof code from device files
        INTEGER, ALLOCATABLE:: DVIDIUA( : ) ! Hourly prof code...
        CHARACTER(LEN=DVKYLEN), ALLOCATABLE:: DVKEYA ( : ) ! FIP//FCID//SKID//DVID

        INTEGER,                ALLOCATABLE:: DVIDX( : ) ! index to sort lists
        CHARACTER(LEN=DVKYLEN), ALLOCATABLE:: DVKEY( : ) ! FIP//FCID//SKID//DVID
        CHARACTER(LEN=DVKYLEN)                DKEY       ! Tmp device key
        CHARACTER(LEN=DVIDLEN)                DVID       ! Tmp device ID

C.........  Temporary variables for building string of source chars.  These
C           variables must be the width of the fields for global source
C           characteristics definition for use in BLDCSRC.
        CHARACTER(LEN=PLTLEN3) FCIDOUT  ! tmp plant ID for output
        CHARACTER(LEN=CHRLEN3) SKIDOUT  ! tmp stack ID for output
        CHARACTER(LEN=CHRLEN3) DVIDOUT  ! tmp device ID for output
        CHARACTER(LEN=CHRLEN3) PRIDOUT  ! tmp process ID for output

C...........   File units and logical/physical names
        INTEGER         DDEV    !  Unit number for device file
        INTEGER         FDEV    !  Unit number for facility file
        INTEGER         MDEV    !  Unit number for emission file
        INTEGER         SDEV    !  Unit number for stack file
        INTEGER         PDEV    !  Unit number for process file

C...........   Other local variables

        REAL            STKF, STKH, STKD, STKT, STKV  ! Temporary stack parms
        REAL            CEFF    !  Temporary control effectiveness
        REAL            DAY2YR  !  Local, leap-year-able, DAY to YEAR factor
        REAL            EMIS    !  Temporary emission value
        REAL            FSAV    !  Flow saved
        REAL            REFF    !  Temporary rule effectiveness
        REAL            XLOC    !  Scratch X-coordinate (UTM X or lon)
        REAL            YLOC    !  Scratch Y-coordinate (UTM Y or lat)
        REAL            XX      !  Temporary X-coordinate (lon)
        REAL            YY      !  Temporary Y-coordinate (lat)

        INTEGER         I, J                ! counters and indices
        INTEGER         K1, K2, K3, K4      ! counters and indices
        INTEGER         L2                  ! counters and indices

        INTEGER          CSS     !  Start of non-blank character string
        INTEGER          COD     !  Temporary pollutant code number
        INTEGER          DLINE, FLINE, SLINE, PLINE ! Number of lines in files
        INTEGER          ES      !  counter for emission file
        INTEGER          FIP, SCC, SIC  ! Temporary fip, scc, sic
        INTEGER          TDIU    !  Temporary hourly profile code
        INTEGER          TWEK    !  Temporary weekly profile code
        INTEGER          LDC, LFC, LSC, LPC  ! Lengths of DVID,FCID,SKID,PRID
        INTEGER, SAVE :: NSRCSAV = 0 ! Cumulative source count
        INTEGER          TPF     !  Temporary temporal ID
        INTEGER          ZONE    !  Temporary UTM zone

        LOGICAL, SAVE :: CFLAG    !  true: recalculate the velocity from flow
        LOGICAL, SAVE :: FIRSTIME = .TRUE.
        LOGICAL          RULFLAG  !  Rule effective file(s): TRUE if exist
        LOGICAL, SAVE :: WFLAG    !  true: convert lat-lons to Western hemisphr

        CHARACTER*2            TMPAA !  tmp time period code
        CHARACTER*10, SAVE  :: FIPFMT! formt to write co/st/cy to string
        CHARACTER*10, SAVE  :: FMTSCC!  format for writing integer SCC to char
        CHARACTER*300          LINE  !  Input line from POINT file
        CHARACTER*300          MESG  !  Text for M3EXIT()
        CHARACTER(LEN=IOVLEN3) CPOL  !  Temporary pollutant code
        CHARACTER(LEN=FIPLEN3) CFIP  !  Character FIP code
        CHARACTER(LEN=POLLEN3) CCOD  !  Character pollutant index to INVPNAM
        CHARACTER(LEN=SCCLEN3) TSCC  !  Temporary character SCC

        CHARACTER*16 :: PROGNAME = 'RDEMSPT' ! Program name

C***********************************************************************
C   begin body of subroutine RDEMSPT

C.........  Set up settings the first time the subroutine is called
        IF( FIRSTIME ) THEN

            FIRSTIME = .FALSE.

C.............  Get settings from the environment
            MESG = 'Flag for recalculating velocity'
            CFLAG = ENVYN( 'VELOC_RECALC', MESG, .FALSE., IOS )

            MESG = 'Western hemisphere flag'
            WFLAG = ENVYN( 'WEST_HSPHERE', MESG, .TRUE., IOS )

C.............  Create format for writing SCC to character
            WRITE( FMTSCC, 94300 ) '(I', SCCLEN3, '.', SCCLEN3, ')'
            WRITE( FIPFMT, '("(I",I2.2,".",I2.2,")")' ) FIPLEN3, FIPLEN3

        ENDIF

C........................................................................
C.......  Get number of lines in the 5 files and allocate memory ........
C........................................................................

        DLINE = GETFLINE( EDEV( 1 ), 'EMS-95 device file' )
        FLINE = GETFLINE( EDEV( 3 ), 'EMS-95 facility file' )
        PLINE = GETFLINE( EDEV( 4 ), 'EMS-95 process file' )
        SLINE = GETFLINE( EDEV( 5 ), 'EMS-95 stack file' ) 

        ALLOCATE( INDXDVA( DLINE ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INDXDVA', PROGNAME )
        ALLOCATE( DVSICA( DLINE ), STAT=IOS )
        CALL CHECKMEM( IOS, 'DVSICA', PROGNAME )
        ALLOCATE( DVIWEKA( DLINE ), STAT=IOS )
        CALL CHECKMEM( IOS, 'DVIWEKA', PROGNAME )
        ALLOCATE( DVIDIUA( DLINE ), STAT=IOS )
        CALL CHECKMEM( IOS, 'DVIDIUA', PROGNAME )
        ALLOCATE( DVKEYA( DLINE ), STAT=IOS )
        CALL CHECKMEM( IOS, 'DVKEYA', PROGNAME )
        ALLOCATE( DVIDX( DLINE ), STAT=IOS )
        CALL CHECKMEM( IOS, 'DVIDX', PROGNAME )
        ALLOCATE( DVKEY( DLINE ), STAT=IOS )
        CALL CHECKMEM( IOS, 'DVKEY', PROGNAME )

        ALLOCATE( INDXFCA( FLINE ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INDXFCA', PROGNAME )
        ALLOCATE( FCZONEA( FLINE ), STAT=IOS )
        CALL CHECKMEM( IOS, 'FCZONEA', PROGNAME )
        ALLOCATE( FCCRDXA( FLINE ), STAT=IOS )
        CALL CHECKMEM( IOS, 'FCCRDXA', PROGNAME )
        ALLOCATE( FCCRDYA( FLINE ), STAT=IOS )
        CALL CHECKMEM( IOS, 'FCCRDYA', PROGNAME )
        ALLOCATE( FCKEYA( FLINE ), STAT=IOS )
        CALL CHECKMEM( IOS, 'FCKEYA', PROGNAME )
        ALLOCATE( FCDESCA( FLINE ), STAT=IOS )
        CALL CHECKMEM( IOS, 'FCDESCA', PROGNAME )
        ALLOCATE( FCIDX( FLINE ), STAT=IOS )
        CALL CHECKMEM( IOS, 'FCIDX', PROGNAME )
        ALLOCATE( FCKEY( FLINE ), STAT=IOS )
        CALL CHECKMEM( IOS, 'FCKEY', PROGNAME )

        ALLOCATE( INDXSKA( SLINE ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INDXSKA', PROGNAME )

        ALLOCATE( IFCKEYA( SLINE ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IFCKEYA', PROGNAME )

        ALLOCATE( SKDIAMA( SLINE ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SKDIAMA', PROGNAME )
        ALLOCATE( SKHEITA( SLINE ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SKHEITA', PROGNAME )
        ALLOCATE( SKTEMPA( SLINE ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SKTEMPA', PROGNAME )
        ALLOCATE( SKCRDXA( SLINE ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SKCRDXA', PROGNAME )
        ALLOCATE( SKCRDYA( SLINE ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SKCRDYA', PROGNAME )
        ALLOCATE( SKVELOA( SLINE ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SKVELOA', PROGNAME )
        ALLOCATE( SKFLOWA( SLINE ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SKFLOWA', PROGNAME )
        ALLOCATE( SKKEYA( SLINE ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SKKEYA', PROGNAME )
        ALLOCATE( SKIDX( SLINE ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SKIDX', PROGNAME )
        ALLOCATE( SKKEY( SLINE ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SKKEY', PROGNAME )

        ALLOCATE( INDXPA( PLINE ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INDXPA', PROGNAME )
        ALLOCATE( PSSCCA( PLINE ), STAT=IOS )
        CALL CHECKMEM( IOS, 'PSSCCA', PROGNAME )
        ALLOCATE( PSKEYA( PLINE ), STAT=IOS )
        CALL CHECKMEM( IOS, 'PSKEYA', PROGNAME )
        ALLOCATE( PSSCC( PLINE ), STAT=IOS )
        CALL CHECKMEM( IOS, 'PSSCC', PROGNAME )
        ALLOCATE( PSKEY( PLINE ), STAT=IOS )
        CALL CHECKMEM( IOS, 'PSKEY', PROGNAME )
        
C.........  Set local variables with unit numbers
        DDEV = EDEV( 1 )
        MDEV = EDEV( 2 )
        FDEV = EDEV( 3 )
        PDEV = EDEV( 4 )
        SDEV = EDEV( 5 )

C.........  Define day to year conversion factor and real type for integer 
C           missing value
        DAY2YR  = 1. / YR2DAY( INY )
C........................................................................
C.............  Head of the FDEV-read loop  .............................
C........................................................................

        FS   = 0     ! For facility.pt's
        IREC = 0
        ERFILDSC = 'facility'
        DO

C.............  Read a line of facility.pt file; Check input status
            LINE = BLANK5
            READ( FDEV, 93000, END=112, IOSTAT=IOS ) LINE
            IREC = IREC + 1

            IF ( IOS .GT. 0 ) GO TO 999  ! to end of subroutine

            ZONE = STR2INT( LINE( 43:44 ) )

C.............  Check zone conversion to integer
            IF( ZONE .EQ. IMISS3 ) THEN

                WRITE( MESG,94010 )
     &                 'UTM zone blank at line', IREC, 'in facility ' //
     &                 'file, so will assume lat/lon coordinates.'
                CALL M3MESG( MESG )
                CYCLE  ! Go on to next iteration

            END IF

            FS = FS + 1

            IF ( FS .LE. FLINE ) THEN

                WRITE ( CFIP,FIPFMT ) 
     &                   1000 * STR2INT( LINE( 1:2 ) ) +
     &                          STR2INT( LINE( 3:5 ) )

                CSS  = LBLANK( LINE( 6:20 ) )
                FCID = LINE( MIN(CSS+6,20):20 )
                IF( FCID .EQ. ' ' ) FCID = EMCMISS3

                INDXFCA( FS ) = FS
                FCKEYA ( FS ) = CFIP // FCID
                FCDESCA( FS ) = LINE( 45:84 )
                FCCRDXA( FS ) = STR2REAL( LINE( 25:33 ) )
                FCCRDYA( FS ) = STR2REAL( LINE( 34:42 ) )
                FCZONEA( FS ) = ZONE

            END IF          !  if fs in bounds

        ENDDO      !  to head of FDEV-read loop
112     CONTINUE   !  end of the FDEV-read loop

        CLOSE( FDEV )

        WRITE( MESG,94010 ) 
     &         'FACILITY FILE processed:' // CRLF() // BLANK5 //
     &         '   Actual FACILITY record count ', FS

        CALL M3MSG2( MESG )

        IF ( FS .GT. FLINE ) THEN

            EFLAG = .TRUE.
            MESG = 'Memory allocation insufficient for FACILITY file.'
            CALL M3MSG2( MESG )

        ELSE	!  else sort input facility table:

            CALL SORTIC( FS, INDXFCA, FCKEYA )
            DO I = 1, FS
                J = INDXFCA( I )
                FCKEY( I ) = FCKEYA( J )
                FCIDX( I ) = J
            ENDDO

        END IF		!  if facility table overflow or not

C........................................................................
C.............  Head of the SDEV-read loop  .............................
C........................................................................

        SS   = 0     ! For stack.pt's
        IREC = 0
        ERFILDSC = 'stack'
        DO

C.............  Read a line of stack.pt file and check input status
            READ( SDEV, 93000, END=123, IOSTAT=IOS ) LINE
            IREC = IREC + 1
            IF ( IOS .GT. 0 ) GO TO 999  ! to end of subroutine

C.............  Get lookup into facilities table
            WRITE( CFIP,FIPFMT ) 1000*STR2INT( LINE( 1:2 ) ) +
     &                               STR2INT( LINE( 3:5 ) )

C.............  Find source match in facilities table
            CSS  = LBLANK( LINE( 6:20 ) )
            FCID = LINE( MIN(CSS+6,20):20 )
            IF( FCID .EQ. ' ' ) FCID = EMCMISS3
            FKEY = CFIP // FCID
            K1   = FINDC( FKEY, FS, FCKEY )
            LFC  = LEN_TRIM( FCID )
                    
            CSS  = LBLANK( LINE( 21:32 ) )
            SKID = LINE( MIN(CSS+21,32):32 )
            IF( SKID .EQ. ' ' ) SKID = EMCMISS3

            IF ( K1 .LE. 0 ) THEN
                WRITE( MESG,94010 )
     &                'Stack not in plant recs: FIP=' // CFIP // 
     &                '; Facility=' // FCID( 1:LFC ) //
     &                '; Stack='    // SKID //
     &                '; line= ', IREC
                CALL M3MESG( MESG )
                CYCLE
            END IF

C.............  Convert units for stack parameters

            STKD = STR2REAL( LINE( 33:40 ) )
            STKH = STR2REAL( LINE( 41:47 ) )
            STKT = STR2REAL( LINE( 48:54 ) ) 
            STKV = STR2REAL( LINE( 55:61 ) )
            STKF = STR2REAL( LINE( 62:71 ) )

            IF ( STKD .GT. 0.0 ) STKD= STKD * FT2M     ! Diam ft to m
            IF ( STKH .GT. 0.0 ) STKH= STKH * FT2M     ! Ht,  ft to m
            IF ( STKT .GT. 0.0 ) STKT= (STKT-32.)*FTOC + CTOK ! Temp, F to K
            IF ( STKV .GT. 0.0 ) STKV= STKV * FT2M     ! Veloc ft/s to m/s
            IF ( STKF .GT. 0.0 ) STKF= STKF * FLWE2M   ! Flow ft^3/min m^3/s

C.............  Initialize FSAV using velocity and diameter
            FSAV = STKV*0.25*PI*STKD*STKD

C.............  Calculate velocity from flow and diameter                
            IF ( CFLAG           .OR.
     &         ( STKV .LE. 0.0   .AND.
     &           STKD .GT. 0.0   .AND.
     &           STKF .GT. 0.0 )       ) THEN

                STKV = STKF / ( 0.25 * PI * STKD * STKD )

C.............  Compare flow to velocity and diameter. Set to exact flow input
C               value only if it is consistent with velocity and diameter
            ELSEIF( STKF. GT. 0 ) THEN 

                IF( ( STKF - FSAV ) / STKF .GT. 0.001 ) THEN
                    FSAV = STKF
                END IF

            END IF

C.............   Resolve coordinates for this stack:

            XLOC = STR2REAL( LINE( 72:80 ) ) ! supposed to be UTM
            YLOC = STR2REAL( LINE( 81:89 ) ) ! supposed to be UTM
                
            IF ( XLOC .LE. 0.0 .OR. YLOC .LE. 0.0 ) THEN
                    
                XLOC = FCCRDXA( FCIDX( K1 ) )
                YLOC = FCCRDYA( FCIDX( K1 ) )

                IF( XLOC .LE. 0.0 .OR. YLOC .LE. 0.0 ) THEN  ! Warning

                    WRITE( MESG,94010 ) 'Source dropped because ' //
     &                     '(X,Y) less than or equal to 0.0 found ' //
     &                     'in stack and facility files for:' //
     &                     CRLF() // BLANK5 //
     &                     'FIP=' // CFIP // '; Facility=' // 
     &                     FCID( 1:LFC ) // '; and Stack=' // SKID
                    CALL M3MESG( MESG )

                    CYCLE       ! to end of do-loop

                END IF              !  if xloc, yloc missing

            END IF

C.............   Create XX and YY for use at end of loop
C.............   Convert UTM coordinates to lat-lon using double precision
C.............   If zone is < 0, then assume already in lat-lon

            ZONE = FCZONEA( FCIDX( K1 ) )
            IF( ZONE .GT. 0 ) THEN
                CALL UTM2LL( XLOC, YLOC, ZONE, XX, YY )

            ELSE
                XX = XLOC
                YY = YLOC
                IF( WFLAG .AND. XX .GT. 0 ) XX = -XX  ! To western hemisphere
            ENDIF

C.............  Store stack file information

            SS = SS + 1

            IF ( SS .LE. SLINE ) THEN

                SKKEYA ( SS ) = CFIP // FCID( 1:LFC ) // SKID
                INDXSKA( SS ) = SS
                IFCKEYA( SS ) = K1
                SKDIAMA( SS ) = STKD
                SKHEITA( SS ) = STKH
                SKTEMPA( SS ) = STKT
                SKVELOA( SS ) = STKV
                SKFLOWA( SS ) = FSAV
                SKCRDXA( SS ) = XX
                SKCRDYA( SS ) = YY

            END IF          !  if ss in bounds

        ENDDO           !  to head of SDEV-read loop
123     CONTINUE        !  end of the SDEV-read loop

        CLOSE( SDEV )

        WRITE( MESG,94010 ) 
     &         'STACK FILE processed:' // CRLF() // BLANK5 //
     &         '   Actual  STACK   record-count ', SS
        CALL M3MSG2( MESG )

        IF ( SS .GT. SLINE ) THEN	!  overflow:

            EFLAG = .TRUE.
            MESG = 'Memory allocation insufficient for STACK file.'
            CALL M3MSG2( MESG )

        ELSE	!  else sort input stack parameter table:

            CALL SORTIC( SS, INDXSKA, SKKEYA )
            DO I = 1, SS
                J = INDXSKA( I )
                SKKEY( I ) = SKKEYA( J )
                SKIDX( I ) = J
            ENDDO

        END IF		!  if stack parameter table overflow or not

C........................................................................
C.............  Head of the PDEV-read loop  .............................
C........................................................................

        PS   = 0     ! For process.pt's
        IREC = 0
        ERFILDSC = 'process'
        DO

C.............  Read a line of process.pt file, and check input status

            READ( PDEV, 93000, END=134, IOSTAT=IOS ) LINE
            IREC = IREC + 1
            IF ( IOS .GT. 0 ) GO TO 999

C.............  Convert and check SCC value

            SCC = STR2INT( LINE( 57:64 ) )
            IF ( SCC .EQ. IMISS3 ) THEN  ! SCC not there

                WRITE( MESG,94010 ) 'WARNING: Source dropped.' //
     &                 CRLF() // BLANK5 // 'Missing or alpha-' //
     &                 'numeric SCC in process file at line', IREC
                CALL M3MESG( MESG )
                CYCLE

            ELSEIF ( SCC .LE. 9999999 ) THEN  ! SCC must be 8 digits

                WRITE( MESG,94010 ) 'WARNING: Source dropped. ' //
     &                 'Invalid SCC "' // LINE( 57:64 ) //
     &                 '" in process file at line', IREC
                CALL M3MESG( MESG )
                CYCLE

            END IF

            PS = PS + 1

            IF ( PS .LE. PLINE ) THEN

                WRITE ( CFIP,FIPFMT ) 
     &               1000 * STR2INT( LINE( 1:2 ) ) +
     &                      STR2INT( LINE( 3:5 ) )

                CSS  = LBLANK( LINE( 6:20 ) )
                FCID = LINE( MIN(CSS+6,20):20 )
                IF( FCID .EQ. ' ' ) FCID = EMCMISS3
                LFC  = LEN_TRIM( FCID )

                CSS  = LBLANK( LINE( 21:32 ) )
                SKID = LINE( MIN(CSS+21,32):32 )
                IF( SKID .EQ. ' ' ) SKID = EMCMISS3
                LSC  = LEN_TRIM( SKID )

                CSS  = LBLANK( LINE( 33:44 ) )
                DVID = LINE( MIN(CSS+33,44):44 )
                IF( DVID .EQ. ' ' ) DVID = EMCMISS3
                LDC  = LEN_TRIM( DVID )

                CSS  = LBLANK( LINE( 45:56 ) )
                PRID = LINE( MIN(CSS+45,56):56 )
                IF( PRID .EQ. ' ' ) PRID = EMCMISS3
                LPC  = LEN_TRIM( PRID )

                PSKEYA( PS ) = CFIP // 
     &                         FCID( 1:LFC ) // SKID( 1:LSC ) //
     &                         DVID( 1:LDC ) // PRID( 1:LPC )

                INDXPA( PS ) = PS
                PSSCCA( PS ) = SCC

            END IF          !  if ps in bounds

        ENDDO           !  to head of PDEV-read loop
134     CONTINUE        !  end of the PDEV-read loop

        CLOSE( PDEV )

        WRITE( MESG,94010 ) 
     &         'PROCESS FILE processed' // CRLF() // BLANK5 //
     &         '   Actual PROCESS  record-count ', PS
        CALL M3MSG2( MESG )

        IF ( PS .GT. PLINE ) THEN  !  overflow

            EFLAG = .TRUE.
            MESG = 'Memory allocation insufficient for PROCESS file.'
            CALL M3MSG2( MESG )

        ELSE	!  sort input process parameter table

            CALL SORTIC( PS, INDXPA,  PSKEYA )
            DO I = 1, PS
                J = INDXPA( I )
                PSSCC( I ) = PSSCCA( J )
                PSKEY( I ) = PSKEYA( J )
            ENDDO

        END IF	!  if input process parameter table overflow or not


C........................................................................
C.............  Head of the DDEV-read loop  .............................
C........................................................................

        DS   = 0     ! For device.pt's
        IREC = 0
        ERFILDSC = 'device'
        DO

            READ( DDEV, 93000, END=145, IOSTAT=IOS ) LINE
            IREC = IREC + 1
            IF ( IOS .GT. 0 ) GO TO 999

C.............  Convert and check SIC value

            SIC = STR2INT( LINE( 45:48 ) )
            IF ( SIC .EQ. IMISS3 ) THEN  ! SIC not there

                WRITE( MESG,94010 ) 
     &                 'Missing SIC in device file at line', IREC,
     &                 CRLF() // BLANK5 // '   Setting to 0000'
                CALL M3MESG( MESG )
                SIC = 0

            ELSEIF ( SIC .EQ. 0 ) THEN

                WRITE( MESG,94010 )
     &                 'Default SIC "0000" in device file at line', IREC
                CALL M3MESG( MESG )

            ELSE IF ( SIC .LT. 111 ) THEN  ! valid SIC codes from 0111 to 9999

                WRITE( MESG,94010 )
     &                 'Invalid SIC "' // LINE( 45:48 ) //
     &                 '" in device file at line', IREC, 
     &                 CRLF() // BLANK5 // '   Setting to 0000'
                CALL M3MESG( MESG )
                SIC = 0

            END IF          !  if sic zero or not

C.................  Convert and check temporal profile numbers

C temp      TMON = STR2INT( LINE( ) )
            TWEK = STR2INT( LINE( 123:124 ) )
            TDIU = STR2INT( LINE( 121:122 ) )

            IF( TWEK .LT. 0 ) TWEK = 0  ! Treat missing as default
            IF( TDIU .LT. 0 ) TDIU = 0  ! Treat missing as default

            IF( TWEK .EQ. 0 ) THEN
                WRITE( MESG,94010 ) 
     &                'Default weekly profile', TWEK, 
     &                'in device file at line', IREC
                CALL M3MESG( MESG )
            END IF

            IF( TDIU .LE. 0 ) THEN
                WRITE( MESG,94010 ) 
     &                 'Default hourly profile', TDIU, 
     &                 'in device file at line', IREC
                CALL M3MESG( MESG )
            END IF

            DS = DS + 1

            IF ( DS .LE. DLINE ) THEN

                WRITE ( CFIP,FIPFMT ) 
     &                   1000 * STR2INT( LINE( 1:2 ) ) + 
     &                          STR2INT( LINE( 3:5 ) )

                CSS  = LBLANK( LINE( 6:20 ) )
                FCID = LINE( MIN(CSS+6,20):20 )
                IF( FCID .EQ. ' ' ) FCID = EMCMISS3
                LFC  = LEN_TRIM( FCID )

                CSS  = LBLANK( LINE( 21:32 ) )
                SKID = LINE( MIN(CSS+21,32):32 )
                IF( SKID .EQ. ' ' ) SKID = EMCMISS3
                LSC  = LEN_TRIM( SKID )

                CSS  = LBLANK( LINE( 33:44 ) )
                DVID = LINE( MIN(CSS+33,44):44 )
                IF( DVID .EQ. ' ' ) DVID = EMCMISS3
                LDC  = LEN_TRIM( DVID )

                INDXDVA( DS ) = DS
                DVKEYA ( DS ) = CFIP          // FCID( 1:LFC ) //
     &                          SKID( 1:LSC ) // DVID( 1:LDC )

                DVSICA ( DS ) = SIC
C temp          DVIMONA( DS ) = TMON
                DVIWEKA( DS ) = TWEK
                DVIDIUA( DS ) = TDIU

            END IF          !  if ds in bounds

        ENDDO           !  to head of DDEV-read loop
145     CONTINUE        !  end of the DDEV-read loop

        CLOSE( DDEV )

        WRITE( MESG,94010 ) 
     &         'DEVICE FILE processed' // CRLF() // BLANK5 //
     &         '   Actual  DEVICE  record-count ', DS
        CALL M3MSG2( MESG )

        IF ( DS .GT. DLINE ) THEN	!  overflow

            EFLAG = .TRUE.
            MESG = 'Memory allocation insufficient for DEVICE file.'
            CALL M3MSG2( MESG )

        ELSE	!  sort input device parameter table:

            CALL SORTIC( DS, INDXDVA, DVKEYA )
            DO I = 1, DS
                J = INDXDVA( I )
                DVIDX ( I ) = J
                DVKEY ( I ) = DVKEYA ( J )
            ENDDO

        END IF	!  if evice parameter table overflow or not


C........  Read rule effectiveness file
        RULFLAG = .FALSE.   ! No rule effectiveness files
        REFF    = 100.0

C........................................................................
C.............  Head of the MDEV-read loop  .............................
C........................................................................

        ES   = NSRCSAV
        IREC = 0
        ERFILDSC = 'emission'
        DO

C.............  Read a line of emission.pt file and check input status

            READ( MDEV, 93000, END=199, IOSTAT=IOS ) LINE
            IREC = IREC + 1
            IF ( IOS .GT. 0 ) GO TO 999

C.............  Find pollutant name in master list to set index COD
C.............  NOTE- Pollutant names here and in INVPNAM converted to uppercase
            CSS  = LBLANK ( LINE( 57:61 ) )
            CPOL = LINE   ( MIN(CSS+57,61):61 )
            CALL UPCASE( CPOL )
            COD  = INDEX1( CPOL, MXIDAT, INVDNAM )

            IF( COD .LE. 0 ) THEN
                L2 = LEN_TRIM( CPOL )
                WRITE( MESG,94010 )  'Source dropped: ' //
     &                 'pollutant name "' // CPOL( 1:L2 ) // 
     &                 '" in emission file at line', IREC,
     &                 CRLF() // BLANK5 // 
     &                 'is not in master pollutants list'
                CALL M3MESG( MESG )
                CYCLE      !  to head of loop
            END IF

C.............  Check and set emissions value

            EMIS = STR2REAL( LINE( 88:100 ) )
            IF ( EMIS .LT. 0.0 )  THEN
                WRITE( MESG,94010 )  'Source dropped: ' //
     &                 'bad emissions value "' // LINE( 88:100 ) //
     &                 '" in emission file at line', IREC
                CALL M3MESG( MESG )
                CYCLE
            END IF

C.............  Set emission.pt file source arrays and temporary arrays
C.............  Assume that all sources will appear in emission.pt files
C.............  and give warning later if this is not the case.

            FIP  = 1000 * STR2INT( LINE( 1:2 ) ) +
     &                    STR2INT( LINE( 3:5 ) )
            WRITE( CFIP,FIPFMT ) FIP
            CSS  = LBLANK( LINE( 6:20 ) )

            FCID = LINE( MIN(CSS+6,20):20 )
            IF( FCID .EQ. ' ' ) FCID = EMCMISS3 
            LFC  = LEN_TRIM( FCID )
            FKEY = CFIP // FCID

            CSS  = LBLANK( LINE( 21:32 ) )
            SKID = LINE( MIN(CSS+21,32):32 )
            IF( SKID .EQ. ' ' ) SKID = EMCMISS3 
            LSC  = LEN_TRIM( SKID )
            SKEY = FKEY( 1:LEN_TRIM( FKEY ) ) // SKID
            K2   = FINDC( SKEY, SS, SKKEY )

            IF( K2 .LE. 0 ) THEN   ! Key not found
                WRITE( MESG,94010 ) 'Emission file record dropped ' //
     &                 'because not in stack file.' //
     &                 CRLF() // BLANK5 // 'FIP=' // CFIP // 
     &                 '; Facility=' // FCID( 1:LFC ) //
     &                 '; Stack='    // SKID
                CALL M3MESG( MESG )

                NDROP = NDROP + 1
                EDROP( COD ) = EDROP( COD ) + EMIS

                CYCLE  ! To next emissions source line

            END IF          !  if stack key not found

            CSS  = LBLANK( LINE( 33:44 ) )
            DVID = LINE( MIN(CSS+33,44):44 )
            IF( DVID .EQ. ' ' ) DVID = EMCMISS3
            LDC  = LEN_TRIM( DVID )
            DKEY = SKEY( 1:LEN_TRIM( SKEY ) ) // DVID
            K3   = FINDC( DKEY, DS, DVKEY )

            IF( K3 .LE. 0 ) THEN    ! Key not found
                WRITE( MESG,94010 ) 'Emission file record dropped ' //
     &                 'because not in device file.' //
     &                 CRLF() // BLANK5 // 'FIP=' // CFIP // 
     &                 '; Facility=' // FCID( 1:LFC ) //
     &                 '; Stack='    // SKID( 1:LSC ) //
     &                 '; Device='   // DVID
                CALL M3MESG( MESG )

                NDROP = NDROP + 1
                EDROP( COD ) = EDROP( COD ) + EMIS

                CYCLE  ! To next source line

            END IF          !  if device key not found

            CSS  = LBLANK( LINE( 45:56 ) )
            PRID = LINE( MIN(CSS+45,56):56 )
            IF( PRID .EQ. ' ' ) PRID = EMCMISS3
            LPC  = LEN_TRIM( PRID )
            PKEY = DKEY( 1:LEN_TRIM( DKEY ) ) // PRID
            K4   = FINDC( PKEY, PS, PSKEY )

            IF( K4 .LE. 0 ) THEN    ! Key not found
                WRITE( MESG,94010 ) 'Emission file record dropped ' //
     &                 ' because not in process file.' //
     &                 CRLF() // BLANK5 // ' FIP=' // CFIP // 
     &                 '; Facility=' // FCID( 1:LFC ) //
     &                 '; Stack='    // SKID( 1:LSC ) //
     &                 '; Device='   // DVID( 1:LDC ) //
     &                 '; Process='  // PRID
                CALL M3MESG( MESG )

                NDROP = NDROP + 1
                EDROP( COD ) = EDROP( COD ) + EMIS

                CYCLE  ! To next emissions source line

            END IF          !   if process key not found

            WRITE( TSCC, FMTSCC ) PSSCC( K4 )

C.............  Check and set time period type (Year/day/hourly)
            TMPAA = LINE( 114:115 )
            CALL UPCASE( TMPAA )
            IF ( TMPAA .EQ. 'AA' ) THEN 

                TPF = MTPRFAC * WTPRFAC       !  use month, week profiles

            ELSE IF ( TMPAA .EQ. 'AD' ) THEN 

                TPF  = WKSET                !  use week profiles
                EMIS = DAY2YR * EMIS

            ELSE IF ( TMPAA .EQ. 'DS' ) THEN

                TPF = 1                     !  use only hourly profiles
                EMIS = DAY2YR * EMIS

            ELSE                            !  unrecognized type

                NDROP = NDROP + 1
                EDROP( COD ) = EDROP( COD ) + EMIS
                WRITE( MESG,94010 )  'Source dropped: unsupported ' //
     &                 'time period type "' // TMPAA //
     &                 '" in emission file at line', IREC
                CALL M3MESG( MESG )
                CYCLE          !  to head of MDEV-read loop

            END IF          !  tests on record type line( 57:58 )

C.............  Set source control parameters (will get rule effectiveness
C.............  from another file.

            CEFF  = STR2REAL( LINE( 126:132 ) )
            CEFF = CEFF * 100.0

C.............  Copy source characteristics from strings with local 
C               length to strings with global lengths
            FCIDOUT = FCID
            SKIDOUT = SKID
            DVIDOUT = DVID
            PRIDOUT = PRID

C.............  Time to store data in unsorted lists if we've made it this far
            ES = ES + 1

            IF ( ES .LE. NRAWIN ) THEN

                IFIPA  ( ES ) = FIP
                ISICA  ( ES ) = DVSICA( DVIDX( K3 ) )
                TPFLGA ( ES ) = TPF
                INVYRA ( ES ) = INY
                IWEKA  ( ES ) = DVIWEKA( DVIDX( K3 ) )
                IDIUA  ( ES ) = DVIDIUA( DVIDX( K3 ) )
                STKHTA ( ES ) = SKHEITA( SKIDX( K2 ) )
                STKDMA ( ES ) = SKDIAMA( SKIDX( K2 ) )
                STKTKA ( ES ) = SKTEMPA( SKIDX( K2 ) )
                STKVEA ( ES ) = SKVELOA( SKIDX( K2 ) ) 
                XLOCAA ( ES ) = SKCRDXA( SKIDX( K2 ) )
                YLOCAA ( ES ) = SKCRDYA( SKIDX( K2 ) )
                POLVLA ( ES,NEM ) = EMIS
                POLVLA ( ES,NCE ) = CEFF
                POLVLA ( ES,NRE ) = REFF
                K1            = IFCKEYA( SKIDX( K2 ) )
                CSCCA  ( ES ) = TSCC
                CORISA ( ES ) = ORSBLNK3
                CBLRIDA( ES ) = BLRBLNK3
                CPDESCA( ES ) = FCDESCA( FCIDX( K1 ) )

                WRITE( CCOD,94125 ) COD
 
                CALL BLDCSRC( CFIP, FCIDOUT, SKIDOUT, DVIDOUT, PRIDOUT, 
     &                        CHRBLNK3, CHRBLNK3, CCOD, CSOURCA( ES ) )

            END IF          !  if S in range

        ENDDO           !  to head of MDEV-read loop

199     CONTINUE        !  end of the MDEV-read loop

        CLOSE( MDEV )

        WRITE( MESG,94010 ) 
     &         'EMISSION FILE processed:'  // CRLF() // BLANK5 //
     &         '   This-file  EMS-95 SOURCE  record-count', ES-NSRCSAV,
     &         CRLF() // BLANK5 //
     &         '   Cumulative EMS-95 SOURCE  record-count', ES

        CALL M3MSG2( MESG )

        NSRCSAV = ES        !  cumulative emissions lines

C........................................................................

        IF( NSRCSAV .GT. NRAWIN ) THEN

            EFLAG = .TRUE.
            MESG = 'Memory allocation insufficient for EMS-95 inventory'
            CALL M3MSG2( MESG )

        ELSE
            NRAWOUT = NSRCSAV

        END IF		!  if overflow or if errors

C.........  Return from subroutine (but deallocate memory first!)

999     CONTINUE

C.........  Deallocate memory before returning from subroutine

        DEALLOCATE( INDXDVA )
        DEALLOCATE( DVSICA )
        DEALLOCATE( DVIWEKA )
        DEALLOCATE( DVIDIUA )
        DEALLOCATE( DVKEYA )
        DEALLOCATE( DVIDX )
        DEALLOCATE( DVKEY )

        DEALLOCATE( INDXFCA )
        DEALLOCATE( IFCKEYA )
        DEALLOCATE( FCZONEA )
        DEALLOCATE( FCCRDXA )
        DEALLOCATE( FCCRDYA )
        DEALLOCATE( FCKEYA )
        DEALLOCATE( FCDESCA )
        DEALLOCATE( FCIDX )
        DEALLOCATE( FCKEY )

        DEALLOCATE( INDXSKA )
        DEALLOCATE( SKDIAMA )
        DEALLOCATE( SKHEITA )
        DEALLOCATE( SKTEMPA )
        DEALLOCATE( SKCRDXA )
        DEALLOCATE( SKCRDYA )
        DEALLOCATE( SKVELOA )
        DEALLOCATE( SKFLOWA )
        DEALLOCATE( SKKEYA )
        DEALLOCATE( SKIDX )
        DEALLOCATE( SKKEY )

        DEALLOCATE( INDXPA )
        DEALLOCATE( PSSCCA )
        DEALLOCATE( PSKEYA )
        DEALLOCATE( PSSCC )
        DEALLOCATE( PSKEY )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94125   FORMAT( I5 )

94300   FORMAT( A, I2.2, A, I2.2, A )

        END SUBROUTINE RDEMSPT
