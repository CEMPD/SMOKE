 
C Version @(#)$Id$ $Source$ $Date$ 

         PROGRAM EMSPOINT

C***********************************************************************
C  program body starts at line 242
C
C  DESCRIPTION:
C       Sort data contained in EMS-95 raw point source files and 
C       construct Models3/EDSS point source file
C
C  PRECONDITIONS REQUIRED:
C       EMS-95 Input data for point sources, concatonated into single
C       emission.pt, facility.pt, stack.pt, process.pt, and device.pt 
C       files.
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       Models-3 I/O
C       FIND1, PROMPTFFILE, PROMPTMFILE, GETYN, TRIMLEN
C
C  REVISION  HISTORY:
C       Prototype  5/96 by MR Houyoux (based on V1.3 rawpoint.F)
C***********************************************************************

      IMPLICIT NONE

C...........   INCLUDES:
        
        INCLUDE 'PTDIMS3.EXT'   !  point-source dimensioning parameters
        INCLUDE 'GRDIMS3.EXT'   !  grid dimensioning parameters
        INCLUDE 'CHDIMS3.EXT'   !  emis chem info (both inventory and model)
        INCLUDE 'TMDIMS3.EXT'   !  time related parameters
        INCLUDE 'CONST3.EXT'    !  physical constants
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.


C...........   PARAMETER

        REAL            FTOC    ! Farenheit to Celcius ratio
        REAL            FLWE2M  ! Flow in cubic ft per min to cubic m per sec
        REAL            FT2M    ! Feet to meters conversion

        PARAMETER     ( FTOC   = 5./9., 
     &                  FLWE2M = 4.71947E-4, 
     &                  FT2M   = 0.3048      )

C...........   EXTERNAL FUNCTIONS and their descriptions:
        
        LOGICAL         DSCGRID
        LOGICAL         ENVYN
        INTEGER         FIND1   !  returns -1 for failure
        INTEGER         GETNUM
        LOGICAL         GETYN
        INTEGER         INDEX1
        INTEGER         JULIAN
        LOGICAL         LAMBERT
        INTEGER         LBLANK
        LOGICAL         LL2LAM
        LOGICAL         UTM2LAM
        INTEGER         PROMPTFFILE
        CHARACTER*16    PROMPTMFILE
        INTEGER         SECSDIFF
        LOGICAL         SETLAM
        INTEGER         STR2INT
        REAL            STR2REAL
        INTEGER         TRIMLEN
        
        EXTERNAL DSCGRID, ENVYN, FIND1, GETNUM, GETYN, INDEX1, JULIAN,
     &           LAMBERT, LBLANK, LL2LAM, UTM2LAM, PROMPTFFILE, 
     &           PROMPTMFILE, SECSDIFF, SETLAM, STR2INT, STR2REAL, 
     &           TRIMLEN


C...........   LOCAL VARIABLES and their descriptions:
C...........   NOTE that SSC (point-source-category) ID's are 10-digit
C...........   unsigned integers which may be treated as a leading 7-digit
C...........   field, and a trailing 3-digit field.  *7 and *3 arrays below
C...........   follow this scheme with parallel arrays

        INTEGER     UZONE       !  UTM zone for the file as a whole
        INTEGER     TZONE       !  time zone
        INTEGER     IWEK        !  Temporary weekly profile code
        INTEGER     IDIU        !  Temporary hourly profile code
        INTEGER     COD         !  Temporary pollutant code number
        INTEGER     FIP
        INTEGER     LFIP, LSCC, LPLT, LSTK
        INTEGER     NDROP       !  Count of number of dropped emis records
        INTEGER     SCC
        INTEGER     SCCL        !  Previous SCC in loop
        INTEGER     SIC
        INTEGER     PLANT
        INTEGER     STACK
        INTEGER     ZONE
        INTEGER     TPF
        INTEGER     INY
        REAL        DAY2YR       !  Local, leap-year-able, DAY to YEAR factor
        REAL        CEFF
        REAL        REFF
        REAL        RPEN
        REAL        XLOC
        REAL        YLOC
        REAL        EMIS
        REAL        STKDA        !  Temporary stack diameter
        REAL        STKF         !  Temporary stack flow
        REAL        STKH         !  Temporary stack height
        REAL        STKD         !  Temporary stack diameter
        REAL        STKT         !  Temporary stack exit temperature
        REAL        STKV         !  Temporary stack exit velocity
        REAL        EDROP( NIPOL ) ! Sum of dropped emissions
        CHARACTER*5  CPOL        !  Temporary pollutant code
        CHARACTER*12 DVID        !  Temporary device ID
        CHARACTER*12 PRID        !  Temporary process ID
        CHARACTER*12 SKID        !  Temporary stack ID
        CHARACTER*15 FCID        !  Temporary facility ID
        
C...........   Point sources table for emission.pt files (some are in unsorted)
        CHARACTER*15 EMFCID( NPEMAX )    ! facility id from emission files
        CHARACTER*12 EMDVID( NPEMAX )    ! device id from emission files
        CHARACTER*12 EMSKID( NPEMAX )    ! stack id from emission files
        CHARACTER*12 EMPRID( NPEMAX )    ! process id from emission files

C...........   Point sources table for facility.pt files
        INTEGER      FS                  ! Counter & number of facility records
        INTEGER      FCZONE( NPFMAX )    ! UTM zone from facility files
        REAL         FCCRDX( NPFMAX )    ! Facility X coordinates
        REAL         FCCRDY( NPFMAX )    ! Facility Y coordinates
        CHARACTER*20 FCKEY ( NPFMAX )    ! FIP // facility ID
        CHARACTER*40 FCDESC( NPFMAX )    ! facility name from facility files

C...........   Point sources table for stack.pt files
        INTEGER      SS                  ! Counter & number of stack records
        REAL         SKDIAM( NPSMAX )    ! stack diameter from stack files [m]
        REAL         SKHEIT( NPSMAX )    ! stack height from stack files [m]
        REAL         SKTEMP( NPSMAX )    ! stack exit temp from stack files [C]
        REAL         SKCRDX( NPSMAX )    ! stack UTM easting coord [m]
        REAL         SKCRDY( NPSMAX )    ! stack UTM northing coord [m]
        REAL         SKVELO( NPSMAX )    ! stack exit velocity [m^3/sec]
        CHARACTER*32 SKKEY ( NPSMAX )    ! FIP // facility ID // stack ID

C...........   Point sources table for process.pt files
        INTEGER      PS                  ! Counter & number of process records
        INTEGER      PSSCC ( NPPMAX )    ! SCC code from process files
        CHARACTER*56 PSKEY ( NPPMAX )    ! FIP // FCID // SKID // DVID // PRID

C...........   Point sources table for device.pt files
        INTEGER      DS                  ! Counter & number of device records
        INTEGER      DVSIC ( NPDMAX )    ! SIC code from device files
        INTEGER      DVIWEK( NPDMAX )    ! Weekly profile code from device files
        INTEGER      DVIDIU( NPDMAX )    ! Hourly profile code...
        CHARACTER*44 DVKEY ( NPDMAX )    ! FIP // FCID // SKID // DVID

C...........   Point Sources Table input unsorted copy (*A); 
C...........   sorted final version (image EMISREC of data record)
        
        INTEGER     NSRC             !  current source-count
        INTEGER     INDEXA( NPEMAX ) !  subscript table for SORTI4()
        INTEGER     IDXSCC( NPEMAX ) !  subscript table for SORTI1() for SCC
        INTEGER     IFIPA ( NPEMAX ) !  source FIPS (county) ID
        INTEGER     ISCCA ( NPEMAX ) !  source SCC
        INTEGER     ISICA ( NPEMAX ) !  source SIC
        INTEGER     IPLANA( NPEMAX ) !  Plant ID
        INTEGER     ISTACA( NPEMAX ) !  Stack ID
        INTEGER     ICODEA( NPEMAX ) !  spcecies-subscript
        INTEGER     TPFLGA( NPEMAX ) !  applicability of temporal profile types
        INTEGER     INVYRA( NPEMAX ) !  inventory year
        INTEGER     IDIUA ( NPEMAX ) !  Hourly profile code for each source
        INTEGER     IWEKA ( NPEMAX ) !  Weekly profile code for each source
        REAL        XLOCAA( NPEMAX ) !  UTM X-location (m)
        REAL        YLOCAA( NPEMAX ) !  UTM Y-location (m)
        REAL        STKHTA( NPEMAX ) !  stack height   (m)
        REAL        STKDMA( NPEMAX ) !  stack diameter (m)
        REAL        STKTKA( NPEMAX ) !  exhaust temperature (deg K)
        REAL        STKVEA( NPEMAX ) !  exhaust velocity    (m/s)
        REAL        RULPEA( NPEMAX ) !  rule penetration   fraction
        REAL        RULEFA( NPEMAX ) !  rule effectiveness fraction
        REAL        CTLEFA( NPEMAX ) !  control efficiency fraction
        REAL        EMISVA( NPEMAX ) !  emissions values (tons/yr)

C.......   Common EMISREC holds an entire output record.  
C.......   Order of arrays in EMISREC _must_ match order of
C.......   variables in the output file.

        INTEGER     NPOINT           !  current source-count
        INTEGER     IFIP  ( NPSRC )  !  source FIPS (county) ID
        INTEGER     ISCC  ( NPSRC )  !  source SCC
        INTEGER     ISIC  ( NPSRC )  !  source SIC
        INTEGER     IPLANT( NPSRC )  !  Plant ID
        INTEGER     ISTACK( NPSRC )  !  Stack ID
        INTEGER     TZONES( NPSRC )  !  time zones
        INTEGER     TPFLAG( NPSRC )  !  applicability of temporal profile types
        INTEGER     INVYR ( NPSRC )  !  inventory year for this record
        REAL        XLOCA ( NPSRC )  !  UTM X-location (m)
        REAL        YLOCA ( NPSRC )  !  UTM Y-location (m)
        REAL        STKHT ( NPSRC )  !  stack height   (m)
        REAL        STKDM ( NPSRC )  !  stack diameter (m)
        REAL        STKTK ( NPSRC )  !  exhaust temperature (deg K)
        REAL        STKVE ( NPSRC )  !  exhaust velocity    (m/s)
        REAL        RULPEN( NPSRC, NIPOL )  !  rule penetration   fraction
        REAL        RULEFF( NPSRC, NIPOL )  !  rule effectiveness fraction
        REAL        CTLEFF( NPSRC, NIPOL )  !  control efficiency fraction
        REAL        EMISV ( NPSRC, NIPOL )  !  emissions values (tons/yr)

        COMMON /EMISREC / IFIP, ISIC, ISCC, IPLANT, ISTACK, TZONES,
     &                    TPFLAG, INVYR, XLOCA, YLOCA, 
     &                    STKHT, STKDM, STKTK, STKVE, 
     &                    CTLEFF, RULEFF, RULPEN, EMISV

        CHARACTER*16  ENAME
        
C...........   Other local variables

        INTEGER         ADEV    !  Unit number for output actual SCCs
        INTEGER         CDEV    !  Unit number for EMS95/SMOKE facility/stack ID
        INTEGER         CSS     !  Start of non-blank character string
        INTEGER         S, I, J, V         !  counters and indices
        INTEGER         K, K1, K2, K3, K4  !  counters and indices
        INTEGER         EDEV    !  Unit number for emission.pt file
        INTEGER         DDEV    !  Unit number for device.pt
        INTEGER         FDEV    !  Unit number for facility.pt file
        INTEGER         PDEV    !  Unit number for process.pt
        INTEGER         SDEV    !  Unit number for stack.pt file
        INTEGER         TDEV    !  Unit number for output temporal profile #s
        INTEGER         ZDEV    !  for time zone file
        INTEGER         LDEV    !  log-device

        LOGICAL         EFLAG   !  input verification:  TRUE iff ERROR
        LOGICAL         RULFLAG !  rule effective file(s): TRUE if exist
        LOGICAL         SFLAG   !  input verification:  report missing species

        REAL            XX, YY  !  scratch location variables
                                
        INTEGER         IOS     !  I/O status
        INTEGER         IREC    !  input line (record) number
        INTEGER         IFIL    !  input file counter
        CHARACTER*5     CFIP    !  Character FIP code
        CHARACTER*16    SCRBUF  !  scratch buffer
        CHARACTER*16    COORDN  !  coordinate system name
        CHARACTER*16    DVNAME  !  device.pt   file logical name buffer 
        CHARACTER*16    EMNAME  !  emission.pt file logical name buffer 
        CHARACTER*16    FCNAME  !  facility.pt file logical name buffer 
        CHARACTER*16    PSNAME  !  process.pt  file logical name buffer 
        CHARACTER*16    SKNAME  !  stack.pt    file logical name buffer 
        CHARACTER*60    KEY     !  Temporary matching KEY
        CHARACTER*60    LKEY    !  Temporary last matching KEY
        CHARACTER*256   MESG    !  text for M3EXIT()
        CHARACTER*240   LINE    !  input line from POINT file


C***********************************************************************
C   begin body of program EMSPOINT
        
        LDEV = INIT3()

        WRITE( *,92000 ) 
     &  ' ',
     &  'Program EMSPOINT to take the concatonated EMS-95 point',
     &  'source foundation files and produce the POINT SOURCE',
     &  'EMISSIONS VECTOR file. ',
     &  ' ',
     &  'You will need to enter the logical names for the input and',
     &  'output files (and to have set them prior to program launch,',
     &  'using "setenv <logicalname> <pathname>").',
     &  'Optional checking that all species are reported for each ',
     &  'source may be turned on via "setenv RAW_SRC_CHECK Y".',
     &  ' ',
     &  'You may use END_OF-FILE (control-D) to quit the program',
     &  'during logical-name entry. Default responses are given in',
     &  'brackets [LIKE THIS] and can be accepted by hitting the',
     &  '<RETURN> key.',
     &  ' '
        
        IF ( .NOT. GETYN( 'Continue with program?', .TRUE. ) ) THEN
            WRITE( *,92000 ) ' ', 'Ending program EMSPOINT', ' '
            CALL EXIT( 0 )
        END IF

        SFLAG = ENVYN( 'RAW_SRC_CHECK', 
     &                 'EMSPOINT check for missing species-records',
     &                 .FALSE.,
     &                 IOS )

C.......   Get the coordinate system name and parameters
C.......   to put into file header:

        IF ( .NOT. DSCGRID( GRDNM, COORDN, GDTYP3D, 
     &              P_ALP3D, P_BET3D, P_GAM3D, XCENT3D, YCENT3D,
     &              XORIG3D, YORIG3D, XCELL3D, YCELL3D, 
     &              NCOLS3D, NROWS3D, NTHIK3D ) ) THEN

            SCRBUF = GRDNM
            MESG   = '"' // SCRBUF( 1:TRIMLEN( SCRBUF ) ) //
     &               '" not found in GRIDDESC file'
            CALL M3EXIT( 'EMSPOINT', 0, 0, MESG, 2 )

        END IF

        IF      ( GDTYP3D .EQ. UTMGRD3 ) THEN
            UZONE = NINT( P_ALP3D )
        ELSE IF ( GDTYP3D .EQ. LAMGRD3 ) THEN
            IF ( .NOT. LAMBERT( COORDN, P_ALP3D, P_BET3D, P_GAM3D,
     &                          XCENT3D, YCENT3D ) ) THEN
                MESG = 'Bad grid: "' // GRDNM // '"'
                CALL M3EXIT( 'EMSPOINT', 0, 0, MESG, 2 )
            END IF
        ELSE IF ( GDTYP3D .EQ. LATGRD3 ) THEN
            !  do nothing -- no setup necessary
        ELSE
            WRITE( MESG, 94010 ) 
     &          'Unsupported coordinate system type', GDTYP3D
            CALL M3EXIT( 'EMSPOINT', 0, 0, MESG, 2 )
        END IF

C.......   Get file name; open input time zone file

        ZDEV = PROMPTFFILE( 
     &         'Enter logical name for TIME ZONE file',
     &         .TRUE., .TRUE., 'ZONES', 'EMSPOINT' )

C.......   Read the point source files:
C.......   Get file name; open input raw point source file

        WRITE( *,92000 ) 
     &  ' ', 
     &  'Enter the logical file names for ',
     &  'EMS-95 RAW POINT SOURCE files...', 
     &  ' '
        IFIL =    0
        S    =    0     ! For sources (emission.pt)
        FS   =    0     ! For facility.pt's
        SS   =    0     ! For stack.pt's
        PS   =    0     ! For process.pt's
        DS   =    0     ! For device.pt's
        INY  = 1994
        EFLAG = .FALSE.

101     CONTINUE

            IFIL = IFIL + 1

            WRITE( EMNAME, 94040 ) 'EMRAW', IFIL
            WRITE( FCNAME, 94040 ) 'FCRAW', IFIL
            WRITE( SKNAME, 94040 ) 'SKRAW', IFIL
            WRITE( PSNAME, 94040 ) 'PSRAW', IFIL
            WRITE( DVNAME, 94040 ) 'DVRAW', IFIL

            EDEV = PROMPTFFILE( 
     &         'Enter the next EMISSION.PT FILE, or "NONE" to end list',
     &          .TRUE., .TRUE., EMNAME, 'EMSPOINT' )
            IF ( EDEV .LT. 0 ) GO TO  123    !  end of loop readin input files

            FDEV = PROMPTFFILE( 
     &              'Enter the next FACILITY.PT FILE',
     &              .TRUE., .TRUE., FCNAME, 'EMSPOINT' )

            SDEV = PROMPTFFILE( 
     &              'Enter the next STACK.PT FILE',
     &              .TRUE., .TRUE., SKNAME, 'EMSPOINT' )

            PDEV = PROMPTFFILE( 
     &              'Enter the next PROCESS.PT FILE',
     &              .TRUE., .TRUE., PSNAME, 'EMSPOINT' )

            DDEV = PROMPTFFILE( 
     &              'Enter the next DEVICE.PT FILE',
     &              .TRUE., .TRUE., DVNAME, 'EMSPOINT' )

            INY = GETNUM( IMISS3, -IMISS3, INY + 1,
     &                    'Enter INVENTORY YEAR for these files' )

            WRITE( *,92000 ) ' '

C.............   Calculate DAY2YR factor
            K = INY / 1000        !  year-number
            IF ( MOD( K,4 ) .NE. 0 ) THEN           !  nonleap years
                DAY2YR = 365.0
            ELSE IF ( MOD( K,100 ) .NE. 0 ) THEN    !  noncentury leap years
                DAY2YR = 366.0
            ELSE IF ( MOD( K,400 ) .NE. 0 ) THEN    !  century nonleap years
                DAY2YR = 365.0
            ELSE                                    !  leap centuries
                DAY2YR = 366.0
            END IF

C........................................................................
C.............  Head of the EDEV-read loop  .............................
C........................................................................
            IREC  =  0
111         CONTINUE

                IREC = IREC + 1

C.................  Read a line of emission.pt file
                READ( EDEV, 93000, END=112, IOSTAT=IOS ) LINE

C.................  Check input status
                IF ( IOS .NE. 0 ) THEN
                    WRITE( MESG, 94010 ) 
     &                  'Error ', IOS,  'reading "' // 
     &                  EMNAME( 1:TRIMLEN( EMNAME ) ) //
     &                  '" at line', IREC
                    CALL M3EXIT( 'EMSPOINT', 0, 0, 2 )
                END IF

C.................  Check pollutant code and set index I
                CSS  = LBLANK( LINE( 57:61 ) )
                CPOL = LINE( CSS+57:61 )
                I = INDEX1( CPOL( 1:TRIMLEN( CPOL ) ), NIPOL, EINAM )
                IF ( I  .LE. 0 ) THEN
                    WRITE( MESG,94010 ) 
     &                  'Bad line', IREC, 
     &                  'Pollutant code "' // LINE(57:61) // '" in "' //
     &                  EMNAME( 1:TRIMLEN( EMNAME ) ) // '"'
                    CALL M3MESG( MESG )
                    GO TO  111      !  to head of loop
                END IF

C.................  Check and set emissions value
                EMIS  = STR2REAL( LINE( 88:100 ) )
                IF ( EMIS .LT. 0.0 )  THEN
                    WRITE( MESG,94010 ) 
     &                  'Bad line', IREC, '" in "' //
     &                  EMNAME( 1:TRIMLEN( EMNAME ) ) // 
     &                  '"--emis value "' // LINE( 88:100 ) // '"'
                    CALL M3MESG( MESG )
                    GO TO  111
                END IF

C.................  Check and set time period type (Year/day/hourly)
                IF ( LINE( 114:115 ) .EQ. 'AA' ) THEN 

                    TPF = MTPRFAC * WTPRFAC     !  use month, week profiles

                ELSE IF ( LINE( 114:115 ) .EQ. 'AD' ) THEN 

                    TPF  = WDTPFAC              !  use week profiles
                    EMIS = DAY2YR * EMIS

                ELSE IF ( LINE( 114:115 ) .EQ. 'AS' ) THEN

                    TPF = 1                     !  use only hourly profiles

                ELSE                            !  unrecognized type
                    WRITE( MESG,94010 ) 
     &                  'Bad line', IREC, '" in "' //
     &                  EMNAME( 1:TRIMLEN( EMNAME ) ) // 
     &                  '"--Unsupported time period type "' // 
     &                  LINE( 57:58 ) // '"'
                    CALL M3MESG( MESG )
                    GO TO  111          !  to head of EDEV-read loop

                END IF          !  tests on record type line( 57:58 )

C.................  Set source control parameters (will get rule effectiveness
C.................  from another file.
                CEFF  = STR2REAL( LINE( 126:132 ) )
                IF ( CEFF .LT. 0.0 )  THEN 
                    CEFF = 0.0
                ELSE
                    CEFF = CEFF * 100.0
                ENDIF

                RPEN = 100.0

C.................  Set emission.pt file source arrays and temporary arrays
C.................  Assume that all sources will appear in emission.pt files
C.................  and give warning later if this is not the case.

                S = S + 1

                IF ( S .LE. NPEMAX ) THEN

                    INDEXA( S ) = S         !  index-table for use in SORTI4()

                    IFIPA ( S ) = STR2INT( LINE( 1:2 ) )*1000 + 
     &                            STR2INT( LINE( 3:5 ) )
                    ICODEA( S ) = I         !  subscript into EINAM(*)
                    TPFLGA( S ) = TPF
                    INVYRA( S ) = INY
                    CTLEFA( S ) = CEFF
                    RULPEA( S ) = RPEN
                    EMISVA( S ) = EMIS

                    CSS  = LBLANK( LINE( 6:20 ) )
                    FCID = LINE( CSS+6:20 )

                    CSS  = LBLANK( LINE( 21:32 ) )
                    SKID = LINE( CSS+21:32 )
                    IF( SKID .EQ. ' ' ) SKID = '*' 
                    
                    CSS  = LBLANK( LINE( 33:44 ) )
                    DVID = LINE( CSS+33:44 )
                    IF( DVID .EQ. ' ' ) DVID = '*' 

                    CSS  = LBLANK( LINE( 45:56 ) )
                    PRID = LINE( CSS+45:56 )
                    IF( PRID .EQ. ' ' ) PRID = '*' 

                    EMFCID( S ) = FCID
                    EMSKID( S ) = SKID
                    EMDVID( S ) = DVID
                    EMPRID( S ) = PRID

                END IF

            GO TO  111      !  to head of EDEV-read loop
112         CONTINUE        !  end of the EDEV-read loop

C........................................................................
C.............  Head of the FDEV-read loop  .............................
C........................................................................
            IREC = 0
113         CONTINUE

                IREC = IREC + 1

C.................  Read a line of facility.pt file
                READ( FDEV, 93000, END=114, IOSTAT=IOS ) LINE

C.................  Check input status
                IF ( IOS .NE. 0 ) THEN
                    WRITE( MESG, 94010 ) 
     &                  'Error ', IOS,  'reading "' // 
     &                  FCNAME( 1:TRIMLEN( FCNAME ) ) //
     &                  '" at line', IREC
                    CALL M3EXIT( 'EMSPOINT', 0, 0, 2 )
                END IF

                ZONE = STR2INT( LINE( 43:44 ) )

C.................  Check zone conversion to integer
                IF( ZONE .EQ. IMISS3 ) THEN


                    WRITE( MESG,94010 )
     &                  'Zone blank at line', IREC,
     &                  'in "' // FCNAME( 1:TRIMLEN( FCNAME ) ) // 
     &                  '", so will assume LAT/LON coordinates.'
                    CALL M3MESG( MESG )
                    GO TO 113

                ENDIF

                FS = FS + 1

                IF ( FS .LE. NPFMAX ) THEN

                    WRITE ( CFIP,93010 ) STR2INT( LINE( 1:2 ) ),
     &                                   STR2INT( LINE( 3:5 ) )

                    CSS  = LBLANK( LINE( 6:20 ) )
                    FCID = LINE( CSS+6:20 )
                    IF( FCID .EQ. ' ' ) FCID = '*' 
                    
                    FCKEY ( FS ) = CFIP // FCID( 1:TRIMLEN( FCID ) )
                    FCDESC( FS ) = LINE( 45:84 )
                    FCCRDX( FS ) = STR2REAL( LINE( 25:33 ) )
                    FCCRDY( FS ) = STR2REAL( LINE( 34:42 ) )
                    FCZONE( FS ) = ZONE

                END IF

            GO TO  113      !  to head of FDEV-read loop
114         CONTINUE        !  end of the FDEV-read loop

C........................................................................
C.............  Head of the SDEV-read loop  .............................
C........................................................................
            IREC = 0
115         CONTINUE

                IREC = IREC + 1

C.................  Read a line of stack.pt file
                READ( SDEV, 93000, END=116, IOSTAT=IOS ) LINE

C.................  Check input status
                IF ( IOS .NE. 0 ) THEN
                    WRITE( MESG, 94010 ) 
     &                  'Error ', IOS,  'reading "' // 
     &                  SKNAME( 1:TRIMLEN( SKNAME ) ) //
     &                  '" at line', IREC
                    CALL M3EXIT( 'EMSPOINT', 0, 0, MESG, 2 )
                END IF

C.................  Convert units for stack parameters
                STKD = STR2REAL( LINE( 33:40 ) ) * FT2M    ! Diam ft to m
                STKH = STR2REAL( LINE( 41:47 ) ) * FT2M    ! Heit ft to m
                STKT = CTOK + FTOC *                       ! Temp F to K 
     &               ( STR2REAL( LINE( 48:54 ) ) - 32 )    
                STKV = STR2REAL( LINE( 55:61 ) ) * FT2M    ! Veloc ft/s to m/s
                STKF = STR2REAL( LINE( 62:71 ) ) * FLWE2M  ! Flow ft^3/min m^3/s

C.................  If velocity blank or zero...
                IF( STKV .LE. 0. ) THEN
C.....................  Calculate from flow and diameter
                    IF( STKD .GT. 0. .AND. 
     &                  STKF .GT. 0.       ) THEN

                        STKDA = PI * STKD * STKD * 0.25      ! 1/area
                        STKV  = STKF * STKDA
C.....................  Apply default
                    ELSE 
                        STKV  = DEFSTKV
                        WRITE( MESG,94050 )
     &                     'Applying default velocity for stack record', 
     &                      IREC, ': VELOC=', STKV, 'm/s'
                        CALL M3MESG( MESG )
                    ENDIF
                ENDIF

C.................  If diameter blank or zero, apply default
                IF( STKD .LE. 0. ) THEN
                    STKD  = DEFSTKD
                    WRITE( MESG,94050 )
     &                 'Applying default diameter for stack record', 
     &                  IREC, ': DIAM=', STKD, 'm'
                    CALL M3MESG( MESG )
                ENDIF

C.................  If temperature blank or zero, apply default
                IF( STKT .LE. 0. ) THEN
                    STKT  = DEFSTKT
                    WRITE( MESG,94050 )
     &                 'Applying default temperature for stack record', 
     &                  IREC, ': TEMP=', STKT, 'deg K'
                    CALL M3MESG( MESG )
                ENDIF

C.................  Store stack file information
                SS = SS + 1

                IF ( SS .LE. NPSMAX ) THEN

                    WRITE ( CFIP,93010 ) STR2INT( LINE( 1:2 ) ),
     &                                   STR2INT( LINE( 3:5 ) )

                    CSS  = LBLANK( LINE( 6:20 ) )
                    FCID = LINE( CSS+6:20 )
                    IF( FCID .EQ. ' ' ) FCID = '*' 

                    CSS  = LBLANK( LINE( 21:32 ) )
                    SKID = LINE( CSS+21:32 )
                    IF( SKID .EQ. ' ' ) SKID = '*' 
                    
                    SKKEY ( SS ) = CFIP // 
     &                             FCID( 1:TRIMLEN( FCID ) ) //
     &                             SKID( 1:TRIMLEN( SKID ) )

                    SKDIAM( SS ) = STKD
                    SKHEIT( SS ) = STKH
                    SKTEMP( SS ) = STKT
                    SKVELO( SS ) = STKV
                    SKCRDX( SS ) = STR2REAL( LINE( 72:80 ) )! Supossed to be UTM
                    SKCRDY( SS ) = STR2REAL( LINE( 81:89 ) )! Supossed to be UTM

                END IF

            GO TO  115      !  to head of SDEV-read loop
116         CONTINUE        !  end of the SDEV-read loop

C........................................................................
C.............  Head of the PDEV-read loop  .............................
C........................................................................
            IREC = 0
117         CONTINUE

                IREC = IREC + 1

C.................  Read a line of process.pt file
                READ( PDEV, 93000, END=118, IOSTAT=IOS ) LINE

C.................  Check input status
                IF ( IOS .NE. 0 ) THEN
                    WRITE( MESG, 94010 ) 
     &                  'Error ', IOS,  'reading "' // 
     &                  PSNAME( 1:TRIMLEN( PSNAME ) ) //
     &                  '" at line', IREC
                    CALL M3EXIT( 'EMSPOINT', 0, 0, MESG, 2 )
                END IF

C.................  Convert and check SCC value
                SCC = STR2INT( LINE( 57:64 ) )
                IF ( SCC .LE. 0 ) THEN
                    WRITE( MESG,94010 )
     &                  'Bad SCC:', SCC, 'at line', IREC,
     &                  'in "' // PSNAME( 1:TRIMLEN( PSNAME ) ) // '"'
                    CALL M3MESG( MESG )
                    GO TO 117

                END IF

                PS = PS + 1

                IF ( PS .LE. NPPMAX ) THEN

                    WRITE ( CFIP,93010 ) STR2INT( LINE( 1:2 ) ),
     &                                   STR2INT( LINE( 3:5 ) )

                    CSS  = LBLANK( LINE( 6:20 ) )
                    FCID = LINE( CSS+6:20 )
                    IF( FCID .EQ. ' ' ) FCID = '*' 

                    CSS  = LBLANK( LINE( 21:32 ) )
                    SKID = LINE( CSS+21:32 )
                    IF( SKID .EQ. ' ' ) SKID = '*' 
                    
                    CSS  = LBLANK( LINE( 33:44 ) )
                    DVID = LINE( CSS+33:44 )
                    IF( DVID .EQ. ' ' ) DVID = '*' 

                    CSS  = LBLANK( LINE( 45:56 ) )
                    PRID = LINE( CSS+45:56 )
                    IF( PRID .EQ. ' ' ) PRID = '*' 

                    PSKEY ( PS ) = CFIP // 
     &                             FCID( 1:TRIMLEN( FCID ) ) //
     &                             SKID( 1:TRIMLEN( SKID ) ) //
     &                             DVID( 1:TRIMLEN( DVID ) ) //
     &                             PRID( 1:TRIMLEN( PRID ) )

                    PSSCC ( PS ) = SCC

                END IF

            GO TO  117      !  to head of PDEV-read loop
118         CONTINUE        !  end of the PDEV-read loop

C........................................................................
C.............  Head of the DDEV-read loop  .............................
C........................................................................
            IREC = 0
119         CONTINUE

                IREC = IREC + 1

C.................  Read a line of device.pt file
                READ( DDEV, 93000, END=120, IOSTAT=IOS ) LINE

C.................  Check input status
                IF ( IOS .NE. 0 ) THEN
                    WRITE( MESG, 94010 ) 
     &                  'Error ', IOS,  'reading "' // 
     &                  DVNAME( 1:TRIMLEN( DVNAME ) ) //
     &                  '" at line', IREC
                    CALL M3EXIT( 'EMSPOINT', 0, 0, MESG, 2 )
                END IF

C.................  Convert and check SIC value
                SIC = STR2INT( LINE( 45:48 ) )
                IF ( SIC .EQ. 0 ) THEN
                    WRITE( MESG,94010 )
     &                  'Default SIC:', SIC, 'at line', IREC,
     &                  'in "' // DVNAME( 1:TRIMLEN( DVNAME ) ) // '"'
                    CALL M3MESG( MESG )
                ELSEIF ( SIC .LT. 0 ) THEN
                    WRITE( MESG,94010 )
     &                  'Missing SIC:', SIC, 'at line', IREC,
     &                  'in "' // DVNAME( 1:TRIMLEN( DVNAME ) ) // '"'
                    CALL M3MESG( MESG )

                    WRITE( MESG,94010 )
     &                  '     Setting to default of 0'
                    CALL M3MESG( MESG )
                    SIC = 0
                END IF

C.................  Convert and check temporal profile numbers
C temp          IMON = STR2INT( LINE( ) )
                IWEK = STR2INT( LINE( 123:124 ) )
                IDIU = STR2INT( LINE( 121:122 ) )

                IF( IWEK .LE. 0 ) THEN
                    WRITE( MESG,94010 ) 
     &                  'Default weekly profile', IWEK, 'at line', IREC,
     &                  'in "' // DVNAME( 1:TRIMLEN( DVNAME ) ) // '"'
                    CALL M3MESG( MESG )
                ENDIF

                IF( IDIU .LE. 0 ) THEN
                    WRITE( MESG,94010 ) 
     &                  'Default hourly profile', IDIU, 'at line', IREC,
     &                  'in "' // DVNAME( 1:TRIMLEN( DVNAME ) ) // '"'
                    CALL M3MESG( MESG )
                ENDIF

                DS = DS + 1

                IF ( DS .LE. NPDMAX ) THEN

                    WRITE ( CFIP,93010 ) STR2INT( LINE( 1:2 ) ),
     &                                   STR2INT( LINE( 3:5 ) )

                    CSS  = LBLANK( LINE( 6:20 ) )
                    FCID = LINE( CSS+6:20 )
                    IF( FCID .EQ. ' ' ) FCID = '*' 

                    CSS  = LBLANK( LINE( 21:32 ) )
                    SKID = LINE( CSS+21:32 )
                    IF( SKID .EQ. ' ' ) SKID = '*' 
                    
                    CSS  = LBLANK( LINE( 33:44 ) )
                    DVID = LINE( CSS+33:44 )
                    IF( DVID .EQ. ' ' ) DVID = '*' 

                    DVKEY ( DS ) = CFIP // 
     &                             FCID( 1:TRIMLEN( FCID ) ) //
     &                             SKID( 1:TRIMLEN( SKID ) ) //
     &                             DVID( 1:TRIMLEN( DVID ) )

                    DVSIC ( DS ) = SIC
C temp              DVIMON( DS ) = IMON
                    DVIWEK( DS ) = IWEK
                    DVIDIU( DS ) = IDIU

                END IF

            GO TO  119      !  to head of DDEV-read loop
120         CONTINUE        !  end of the DDEV-read loop

C............  Read rule effectiveness file
C................. NOTE: Will be added later.  Not needed for SESARM
            RULFLAG = .FALSE.   ! No rule effectiveness files

        GO TO 101   !  to head of loop on input files
        
C........................................................................
123     CONTINUE        !  end loop on input files

        WRITE( *,92000 ) ' '

        NSRC = S
        IF ( NSRC .GT. NPEMAX ) THEN
            EFLAG = .TRUE.
            WRITE( *,92010 ) 
     &          'Record count                ', NSRC,
     &          'Max dimensioned record-count', NPEMAX
            CALL M3MESG( 'Max record-count exceeded in EMISSION files.')
            WRITE( *,92000 ) ' '
           
        END IF
        
        IF ( FS .GT. NPFMAX ) THEN
            EFLAG = .TRUE.
            WRITE( *,92010 ) 
     &          'Record count                ', FS,
     &          'Max dimensioned record-count', NPFMAX
            CALL M3MESG('Max record-count exceeded in FACILITY files.')
            WRITE( *,92000 ) ' '
           
        END IF
        
        IF ( SS .GT. NPSMAX ) THEN
            EFLAG = .TRUE.
            WRITE( *,92010 ) 
     &          'Record count                ', SS,
     &          'Max dimensioned record-count', NPSMAX
            CALL M3MESG('Max record-count exceeded in STACK files.')
            WRITE( *,92000 ) ' '
           
        END IF
        
        IF ( PS .GT. NPPMAX ) THEN
            EFLAG = .TRUE.
            WRITE( *,92010 ) 
     &          'Record count                ', PS,
     &          'Max dimensioned record-count', NPPMAX
            CALL M3MESG('Max record-count exceeded in PROCESS files.')
            WRITE( *,92000 ) ' '
           
        END IF
        
        IF ( DS .GT. NPDMAX ) THEN
            EFLAG = .TRUE.
            WRITE( *,92010 ) 
     &          'Record count                ', DS,
     &          'Max dimensioned record-count', NPDMAX
            CALL M3MESG('Max record-count exceeded in DEVICE files.')
            WRITE( *,92000 ) ' '
           
        END IF

C...........   Exit if any errors
        IF ( EFLAG ) THEN
            CALL M3EXIT( 'EMSPOINT', 0, 0, 
     &           'Error reading POINT SOURCE FOUNDATION files.', 2 )
        END IF

        ADEV = PROMPTFFILE( 
     &          'Enter the name of the ACTUAL SCC output file',
     &          .FALSE., .TRUE., 'PSCC', 'EMSPOINT' )

        TDEV = PROMPTFFILE( 
     &          'Enter the name of the TEMPORAL XREF output file',
     &          .FALSE., .TRUE., 'PTREF', 'EMSPOINT' )

        CDEV = PROMPTFFILE( 
     &          'Enter the name of the Facility/Stack name output file',
     &          .FALSE., .TRUE., 'PSRC', 'EMSPOINT' )

C...........   Initialize variables for keeping track of dropped emissions
        NDROP = 0
        DO 124 COD = 1, NIPOL
            EDROP( COD ) = 0.
  124   CONTINUE

C...........   Loop through emissions records and fill in information
C...........   from other tables

        WRITE( *,92000 ) ' ', 'Processing RAW POINT SOURCE file...', ' '

        J = 0
        DO 125 S = 1, NSRC

C.............  Set pollutant type code for use in keeping track of how
C               much emissions we are dropping
            COD = ICODEA( S )

C.............  Get UTM zone and default coordinates from facility records
            WRITE( CFIP,93020 ) IFIPA( S )

            KEY  = CFIP // EMFCID( S )
            LKEY = KEY

            K1 = INDEX1( KEY( 1:TRIMLEN( KEY ) ), FS, FCKEY )

            IF( K1 .LE. 0 ) THEN    ! Key not found
                WRITE( MESG,94010 )
     &                'Source dropped: not in facility recs: FIP=' //
     &                 CFIP // '; Facility =' // EMFCID( S )
                CALL M3MESG( MESG )

                NDROP = NDROP + 1
                EDROP( COD ) = EDROP( COD ) + EMISVA( S )

                GO TO 125  ! To end of do-loop

            ENDIF

            ZONE        = FCZONE( K1 )
            XLOC        = FCCRDX( K1 )
            YLOC        = FCCRDY( K1 )

C.............  Get stack parameters from stack records

            KEY  = LKEY( 1:TRIMLEN( LKEY ) ) // EMSKID( S )
            LKEY = KEY

            K2 = INDEX1( KEY( 1:TRIMLEN( KEY ) ), SS, SKKEY )

            IF( K2 .LE. 0 ) THEN   ! Key not found
                WRITE( MESG,94010 )
     &                'Source dropped: not in stack recs: FIP=' //
     &                 CFIP // '; Facility=' // 
     &                 EMFCID( S )( 1:TRIMLEN( EMFCID( S ) ) ) //
     &                '; Stack=' // EMSKID( S )
                CALL M3MESG( MESG )

                NDROP = NDROP + 1
                EDROP( COD ) = EDROP( COD ) + EMISVA( S )

                GO TO 125  ! To end of do-loop

            ENDIF

C.............   Decide which coordinates to use (facility or stack)
            IF( SKCRDX( K2 ) .GT. 0. .AND. SKCRDY( K2 ) .GT. 0. ) THEN
                XLOC = SKCRDX( K2 )  ! Use stack X value
                YLOC = SKCRDY( K2 )  ! Use stack Y value

            ELSE                    ! Check facility coordinates

                IF( XLOC .LE. 0. .OR. YLOC .LE. 0. ) THEN  ! Warning
                    WRITE( MESG,94010 )
     &                'Source dropped: bad coords for stack rec', K2,
     &                 ': FIP='// CFIP // '; Facility=' // 
     &                 EMFCID( S )( 1:TRIMLEN( EMFCID( S ) ) ) //
     &                '; and Stack=' // EMSKID( S )
                    CALL M3MESG( MESG )

                    NDROP = NDROP + 1
                    EDROP( COD ) = EDROP( COD ) + EMISVA( S )

                    GO TO 125            ! to end of do-loop
                ENDIF
            ENDIF

C.............   Create XLOC and YLOC for use at end of loop
C.............   For UTM Grid...
            IF ( GDTYP3D .EQ. UTMGRD3 ) THEN   !  UTM output
 
                IF ( ZONE .EQ. IMISS3 ) THEN        !  lat-lon input
 
                    CALL LL2UTM( XLOC, YLOC, UZONE, XX, YY )
                    XLOC = XX
                    YLOC = YY
 
                ELSE IF ( ZONE .NE. UZONE ) THEN    !  UTM input w/ diff zone
 
                    CALL UTM2LL( XLOC, YLOC, ZONE, XX,   YY )
                    CALL LL2UTM( XX,   YY,  UZONE, XLOC, YLOC )

                END IF  !  if zone missing (lat-lon), or different

C.............   Set XLOC and YLOC for Lambert Grid ...
            ELSE IF ( GDTYP3D .EQ. LAMGRD3 ) THEN
 
                IF ( ZONE .EQ. IMISS3 ) THEN        !  lat-lon input

                    IF ( .NOT. LL2LAM( XLOC, YLOC, XX, YY ) ) THEN
                        WRITE( MESG,94010 )
     &                    'Bad coordinates: LL2LAM() failed for FIP=' //
     &                     CFIP // '; Facility=' // EMFCID( S ) //
     &                    '; and Stack=' // EMSKID( S )
                        CALL M3MESG( MESG )

                        WRITE( MESG,94060 )
     &                    '    Setting coordinates to ( ', 
     &                         XORIG3D, ',', YORIG3D, ')'
                        CALL M3MESG( MESG )

                        XX = XORIG3D
                        YY = YORIG3D

                    END IF

                    XLOC = XX
                    YLOC = YY

                ELSE                                !  utm input
 
                    XX = XLOC
                    YY = YLOC

                    IF ( .NOT. UTM2LAM( XX, YY, ZONE,
     &                                  XLOC, YLOC ) ) THEN
                        WRITE( MESG,94010 )
     &                    'Bad coordinates: UTM2LAM() failed for FIP='//
     &                     CFIP // '; Facility=' // EMFCID( S ) //
     &                    '; and Stack=' // EMSKID( S )
                        CALL M3MESG( MESG )

                        WRITE( MESG,94060 )
     &                    '    Setting coordinates to ( ', 
     &                         XORIG3D, ',', YORIG3D, ')'
                        CALL M3MESG( MESG )

                        XLOC = XORIG3D
                        YLOC = YORIG3D
                    END IF
 
                END IF      !  if zone missing (lat-lon) or not (utm)

C.............   Set XLOC and YLOC for LAT-LON Grid ...
            ELSE IF ( GDTYP3D .EQ. LATGRD3 ) THEN
 
                IF ( ZONE .NE. IMISS3 ) THEN          ! utm input
 
                    XX = XLOC
                    YY = YLOC
                    CALL UTM2LL( XX, YY, ZONE, XLOC, YLOC )
 
                END IF
 
            END IF          !  if coord type UTM or Lambert or Lat-Lon

C.............  Get SIC from device records
            KEY  = LKEY( 1:TRIMLEN( LKEY ) ) // EMDVID( S )
            LKEY = KEY

            K3 = INDEX1( KEY( 1:TRIMLEN( KEY ) ), DS, DVKEY )

            IF( K3 .LE. 0 ) THEN    ! Key not found
                WRITE( MESG,94010 ) 
     &                'Source dropped: not in device recs: FIP=' //
     &                 CFIP // '; Facility=' // 
     &                 EMFCID( S )( 1:TRIMLEN( EMFCID( S ) ) ) //
     &                '; Stack=' // 
     &                 EMSKID( S )( 1:TRIMLEN( EMSKID( S ) ) ) //
     &                '; Device=' // 
     &                 EMDVID( S )( 1:TRIMLEN( EMDVID( S ) ) )
                CALL M3MESG( MESG )

                NDROP = NDROP + 1
                EDROP( COD ) = EDROP( COD ) + EMISVA( S )

                GO TO 125  ! To end of do-loop

            ENDIF

C.............  Get SCC from process records
            KEY  = LKEY( 1:TRIMLEN( LKEY ) ) // EMPRID( S )
            LKEY = KEY

            K4 = INDEX1( KEY( 1:TRIMLEN( KEY ) ), PS, PSKEY )

            IF( K4 .LE. 0 ) THEN    ! Key not found
                WRITE( MESG,94010 )
     &                'Source dropped: not in process recs: FIP=' //
     &                 CFIP // '; Facility=' // 
     &                 EMFCID( S )( 1:TRIMLEN( EMFCID( S ) ) ) //
     &                '; Stack=' // 
     &                 EMSKID( S )( 1:TRIMLEN( EMSKID( S ) ) ) //
     &                '; Device=' // 
     &                 EMDVID( S )( 1:TRIMLEN( EMDVID( S ) ) ) //
     &                '; Process=' // 
     &                 EMPRID( S )( 1:TRIMLEN( EMPRID( S ) ) )
                CALL M3MESG( MESG )

                NDROP = NDROP + 1
                EDROP( COD ) = EDROP( COD ) + EMISVA( S )

                GO TO 125  ! To end of do-loop

            ENDIF

C.............  Get rule effectiveness from prulefac records
            REFF = 100.
            IF( RULFLAG ) THEN   ! If rule effectiveness files exists
                REFF = 100.  ! After add read, will do match on ASCII keys here
            ENDIF

C.............  If make it here (haven't skipped to 125 ) then we've
C.............  decided that source is valid, so index by one and store.

            J = J + 1

            INDEXA( J ) = J         !  index-table for use in SORTI4()
            IDXSCC( J ) = J         !  index-table for use in SORTI1()

            IFIPA ( J ) = IFIPA ( S )
            IPLANA( J ) = K1
            ISTACA( J ) = K2
            ICODEA( J ) = ICODEA( S )
            TPFLGA( J ) = TPFLGA( S )
            INVYRA( J ) = INVYRA( S )
            STKHTA( J ) = SKHEIT( K2 )
            STKDMA( J ) = SKDIAM( K2 )
            STKTKA( J ) = SKTEMP( K2 )
            STKVEA( J ) = SKVELO( K2 )
            XLOCAA( J ) = XLOC
            YLOCAA( J ) = YLOC
            ISICA ( J ) = DVSIC ( K3 )
            IWEKA ( J ) = DVIWEK( K3 )
            IDIUA ( J ) = DVIDIU( K3 )
            ISCCA ( J ) = PSSCC ( K4 )
            CTLEFA( J ) = CTLEFA( S )
            RULPEA( J ) = RULPEA( S )
            RULEFA( J ) = REFF
            EMISVA( J ) = EMISVA( S )

  125   CONTINUE   ! End of source processing loop

C.........  Set new source count
        NSRC = J

C.........  Report how many sources were dropped
        IF( NDROP .GT. 0 ) THEN
            WRITE( *,92000 ) ' '
            WRITE( MESG,94010 ) 'NOTE:', NDROP, 
     &           'sources dropped. New source count is', NSRC
            CALL M3MESG( MESG )

C.........  Report dropped emissions totals
            DO 126 I = 1, NIPOL

                WRITE( MESG,94060 ) 
     &             '     Dropped ' 
     &             // EINAM( I )( 1:TRIMLEN(EINAM( I )) ) //
     &             ' emissions: ', EDROP( I ), 'tons/year'
                CALL M3MESG( MESG )

  126       CONTINUE
            WRITE( *,92000 ) ' '

        ENDIF

C.......   Use SORTI4() to perform an indirect sort by FIPS,SCC,PLANT,STACK
C.......   then permute the records according to the result:

        CALL SORTI4( NSRC, INDEXA, IFIPA, ISCCA, IPLANA, ISTACA )

        DO  132  I = 1, NIPOL	!  initializations
        DO  131  S = 1, NPSRC
            EMISV ( S,I ) = BADVAL3
            RULPEN( S,I ) = BADVAL3
            RULEFF( S,I ) = BADVAL3
            CTLEFF( S,I ) = BADVAL3
131     CONTINUE
132     CONTINUE		!  end initializations

        LFIP  = -1
        LSCC  = -1
        LPLT  = -1
        LSTK  = -1
        J     =  0

        DO  133  S = 1, NSRC

            I = INDEXA( S )

            FIP   = IFIPA ( I )
            SCC   = ISCCA ( I )
            PLANT = IPLANA( I )
            STACK = ISTACA( I )

C.............  If first time got to this source...
            IF ( FIP   .NE. LFIP  .OR.
     &           SCC   .NE. LSCC  .OR.
     &           PLANT .NE. LPLT  .OR.
     &           STACK .NE. LSTK       ) THEN

                J = J + 1

                IF( J .LE. NPSRC ) THEN

                    LFIP = FIP
                    LSCC = SCC
                    LPLT = PLANT
                    LSTK = STACK

C.....................  Set final arrays
                    IFIP  ( J ) = FIP
                    ISCC  ( J ) = SCC
                    ISIC  ( J ) = ISICA ( I )
                    IPLANT( J ) = PLANT
                    ISTACK( J ) = STACK
                    TPFLAG( J ) = TPFLGA( I )
                    INVYR ( J ) = INVYRA( I )
                    XLOCA ( J ) = XLOCAA( I )
                    YLOCA ( J ) = YLOCAA( I )
                    STKHT ( J ) = STKHTA( I )
                    STKDM ( J ) = STKDMA( I )
                    STKTK ( J ) = STKTKA( I )
                    STKVE ( J ) = STKVEA( I )

C.....................  Write out ASCII temporal profiles file for source
C temp (using hard-coded monthly profile )

                    WRITE( TDEV, 93040 ) 1, IWEKA( I ), IDIUA( I )

                ELSE

                    WRITE( MESG,94010 ) 
     &                'NOTE: Skipping source', J, 'to avoid overflow.'
                    CALL M3MESG( MESG )

                    GO TO 133  ! To end of loop

                ENDIF

            END IF      !  if first encounter with this source

            V    = ICODEA( I )
            EMIS = EMISVA( I )
            IF ( EMISV( J,V ) .LE. AMISS3 ) THEN

                EMISV ( J,V ) = EMIS
                CTLEFF( J,V ) = CTLEFA( I )
                RULEFF( J,V ) = RULEFA( I )
                RULPEN( J,V ) = RULPEA( I )

            ELSE IF ( EMIS .GT. 0.0 ) THEN        !  note "SP" ~~~> 0
                          
                EMISV ( J,V ) = EMISV ( J,V ) + EMIS

            END IF

133     CONTINUE

        NPOINT = J

        IF ( NPOINT .NE. NPSRC ) THEN

            WRITE( *,92000 ) ' '
            WRITE( MESG,94010 ) 
     &          'Actual      number of sources', NPOINT,
     &          'Dimensioned number of sources', NPSRC,       
     &          'Do not match!'
            CALL M3MESG ( MESG )
            WRITE( *,92000 ) ' '

            IF( NPOINT. GT. NPSRC ) NPOINT = NPSRC

        END IF

        IF ( SFLAG ) THEN       !  if missing records are fatal:

            EFLAG = .FALSE.

            DO  146  V = 1, NIPOL
            DO  145  S = 1, NPSRC

               IF ( EMISV( S,V ) .LT. AMISS3 ) THEN

                    EFLAG = .TRUE.
                    WRITE( MESG,94020 )
     &                  'Missing record:  FIP:', IFIP( S ),
     &                  'SCC:',     ISCC( S ), 
     &                  'Plant:',   IPLANT( S ),
     &                  'Stack:',   ISTACK( S ),
     &                  'Species: ' // EINAM( V )
                    CALL M3MESG( MESG )

                END IF

145         CONTINUE
146         CONTINUE

            IF ( EFLAG ) THEN
                CALL M3EXIT( 'EMSPOINT', 0, 0,
     &                       'Missing species recs in input file', 2 )
            END IF

        ELSE            !  "missing" ~~~> 0

            DO  148  V = 1, NIPOL
            DO  147  S = 1, NPSRC
                IF ( EMISV( S,V ) .LT. AMISS3 ) EMISV( S,V ) = 0.0
147         CONTINUE
148         CONTINUE

        END IF          !  sflag or not:  check for non-set EMIS values

C.......   Write out EMS-95 / SMOKE facility-stack table

        DO 150 S = 1, NPOINT

            FIP = IFIP  ( S )
            SCC = ISCC  ( S )
            J   = IPLANT( S )
            K   = ISTACK( S )

C.............  Retrive stack ID from SKKEY by taking out FCKEY
            I    = TRIMLEN( FCKEY( J ) )
            FCID = FCKEY( J )( 6:I )
            SKID = SKKEY( K )( I+1:32 )

            WRITE( CDEV,93050 ) FIP, SCC, J, K, FCID, SKID, FCDESC( J )

150     CONTINUE

C.......   Process the time zone file

        WRITE( *,92000 ) ' ', 'Processing TIME ZONE file...', ' '
        IREC = 0

155     CONTINUE        !  head of time-zone loop

            READ( ZDEV,*,END=188,IOSTAT=IOS ) FIP, TZONE 
            IREC = IREC + 1

            IF ( IOS .NE. 0 ) THEN
                WRITE( *,92010 ) 
     &              'Unit number', ZDEV, 
     &              'I/O Status ', IOS, 
     &              'Line number', IREC
                    CALL M3EXIT( 'EMSPOINT', 0, 0, 
     &              'Error reading TIME ZONE file.', 2 )
            END IF

            IF ( FIP .EQ. 0 ) THEN              !  fallback -- all sources

                DO  156  S = 1, NPSRC
                    TZONES( S ) = TZONE
156             CONTINUE

            ELSE IF ( MOD( FIP, 100 ) .EQ. 0 ) THEN     !  state-specific zone

                FIP = FIP / 1000
                DO  166  S = 1, NPSRC
                    IF ( FIP .EQ. IFIP( S ) / 1000 ) THEN
                        TZONES( S ) = TZONE
                    END IF
166             CONTINUE

            ELSE                                        !  county-specific zone

                DO  177  S = 1, NPSRC
                    IF ( FIP .EQ. IFIP( S ) ) THEN
                        TZONES( S ) = TZONE
                    END IF
177             CONTINUE

            END IF      !  if fip zero, or nn000, or not.

            GO TO  155

188     CONTINUE        !  end of time-zone loop


C.......   Get file name; open output point sources file

        FTYPE3D = GRDDED3
        P_ALP3D = DBLE( UZONE )
        P_BET3D = DBLE( AMISS3 )
        P_GAM3D = DBLE( AMISS3 )
        XCENT3D = 0.0D0
        YCENT3D = 0.0D0
        XORIG3D = DBLE( AMISS3 )
        YORIG3D = DBLE( AMISS3 )
        SDATE3D = 0 !  n/a
        STIME3D = 0 !  n/a
        TSTEP3D = 0             !  time independent
        NVARS3D = 4 * NIPOL + 14
        NCOLS3D = 1
        NROWS3D = NPSRC    !  number of rows = # of point sources.
        NLAYS3D = 1
        NTHIK3D = 1
        GDTYP3D = GDTYP3D
        VGTYP3D = IMISS3
        VGTOP3D = AMISS3
        GDNAM3D = COORDN

        FDESC3D( 1 ) = 'NC Point Source emissions values.'
        FDESC3D( 2 ) = 'Inventory  type "' // LINE(  1: 1 ) // '"'
        FDESC3D( 3 ) = 'Projection year "' // LINE(  3: 4 ) // '"'
        FDESC3D( 4 ) = 'Base year       "' // LINE(  6: 7 ) // '"'
        FDESC3D( 5 ) = 'Emission   type "' // LINE(  9:10 ) // '"'
        FDESC3D( 7 ) = 'Beginning  Date "' // LINE( 60:67 ) // '"'
        FDESC3D( 8 ) = 'Ending     Date "' // LINE( 69:76 ) // '"'

        IF ( LINE( 57:58 ) .EQ. 'PO' ) THEN
            FDESC3D( 6 ) = 'Period     type "' // LINE( 57:58 ) //
     &                     '" (peak ozone day)'
        ELSE IF ( LINE( 57:58 ) .EQ. 'PC' ) THEN
            FDESC3D( 6 ) = 'Period     type "' // LINE( 57:58 ) // 
     &                     '" (peak CO day)'
        ELSE IF ( LINE( 57:58 ) .EQ. '  ' ) THEN
            FDESC3D( 6 ) = 'Period     type "' // LINE( 57:58 ) // 
     &                     '" (annual total)'
        ELSE
            FDESC3D( 6 ) = 'Period     type "' // LINE( 57:58 ) // '"'
        END IF

        DO  211  I = 9, MXDESC3
            FDESC3D( I ) = ' '
211     CONTINUE

        J = 1
        VNAME3D( J ) = 'IFIP'
        VTYPE3D( J ) = M3INT
        UNITS3D( J ) = 'n/a'
        VDESC3D( J ) = 'FIP code for counties'
        J = J + 1

        VNAME3D( J ) = 'ISIC'
        VTYPE3D( J ) = M3INT
        UNITS3D( J ) = 'n/a'
        VDESC3D( J ) = 'Source Industrial Code'
        J = J + 1

        VNAME3D( J ) = 'ISCC'
        VTYPE3D( J ) = M3INT
        UNITS3D( J ) = 'n/a'
        VDESC3D( J ) = 'Point Source Category Code'
        J = J + 1

        VNAME3D( J ) = 'IPLANT'
        VTYPE3D( J ) = M3INT
        UNITS3D( J ) = 'n/a'
        VDESC3D( J ) = 'Point Source plant-ID'
        J = J + 1

        VNAME3D( J ) = 'ISTACK'
        VTYPE3D( J ) = M3INT
        UNITS3D( J ) = 'n/a'
        VDESC3D( J ) = 'Point Source stack-ID'
        J = J + 1

        VNAME3D( J ) = 'TZONES'
        VTYPE3D( J ) = M3INT
        UNITS3D( J ) = 'n/a'
        VDESC3D( J ) = 'Time zone for site'
        J = J + 1

        VNAME3D( J ) = 'TPFLAG'
        VTYPE3D( J ) = M3INT
        UNITS3D( J ) = 'T|2? T|3?'
        VDESC3D( J ) = 'Use week(2), month(3) temporal profiles or not'
        J = J + 1

        VNAME3D( J ) = 'INVYR'
        VTYPE3D( J ) = M3INT
        UNITS3D( J ) = 'year AD'
        VDESC3D( J ) = 'Year of inventory for this record'
        J = J + 1

        VNAME3D( J ) = 'XLOCA'
        VTYPE3D( J ) = M3REAL
        IF( GDTYP3D .EQ. LATGRD3 ) THEN
            UNITS3D( J ) = 'deg LON'
        ELSE
            UNITS3D( J ) = 'meters'
        ENDIF
        VDESC3D( J ) = 'X coordinate for site'
        J = J + 1

        VNAME3D( J ) = 'YLOCA'
        VTYPE3D( J ) = M3REAL
        IF( GDTYP3D .EQ. LATGRD3 ) THEN
            UNITS3D( J ) = 'deg LAT'
        ELSE
            UNITS3D( J ) = 'meters'
        ENDIF
        VDESC3D( J ) = 'Y coordinate for site'
        J = J + 1

        VNAME3D( J ) = 'STKHT'
        VTYPE3D( J ) = M3REAL
        UNITS3D( J ) = 'm'
        VDESC3D( J ) = 'Stack height'
        J = J + 1

        VNAME3D( J ) = 'STKDM'
        VTYPE3D( J ) = M3REAL
        UNITS3D( J ) = 'm'
        VDESC3D( J ) = 'Stack diameter'
        J = J + 1

        VNAME3D( J ) = 'STKTK'
        VTYPE3D( J ) = M3REAL
        UNITS3D( J ) = 'deg K'
        VDESC3D( J ) = 'Stack exhaust temperature'
        J = J + 1

        VNAME3D( J ) = 'STKVE'
        VTYPE3D( J ) = M3REAL
        UNITS3D( J ) = 'm/s'
        VDESC3D( J ) = 'Stack exhaust velocity'
        J = J + 1

        DO  222  V = 1 , NIPOL

            VNAME3D( J ) = 'CTLEFF_' // EINAM( V )
            VTYPE3D( J ) = M3REAL
            UNITS3D( J ) = 'n/a'
            VDESC3D( J ) = 
     &      'Control efficiency (in [0,100], or "MISSING": < -9.0E36)'
            J = J + 1

222     CONTINUE        !  end loop on inventory pollutants I

        DO  233  V = 1, NIPOL

            VNAME3D( J ) = 'RULEFF_' // EINAM( V )
            VTYPE3D( J ) = M3REAL
            UNITS3D( J ) = 'n/a'
            VDESC3D( J ) = 
     &      'Rule Effectiveness  (in [0,100], or "MISSING": < -9.0E36)'
            J = J + 1

233     CONTINUE        !  end loop on inventory pollutants I

        DO  244  V = 1 , NIPOL

            VNAME3D( J ) = 'RULPEN_' // EINAM( V )
            VTYPE3D( J ) = M3REAL
            UNITS3D( J ) = 'n/a'
            VDESC3D( J ) = 
     &      'Rule penetration (in [0,100], or "MISSING": < -9.0E36)'
            J = J + 1

244     CONTINUE        !  end loop on inventory pollutants I

        DO  255  V = 1 , NIPOL

            VNAME3D( J ) = EINAM( V )
            VTYPE3D( J ) = M3REAL
            UNITS3D( J ) = 'tons/year'
            VDESC3D( J ) = 
     &        LINE( 1:1 ) // LINE( 57:58 )  // ' emissions totals'
            J = J + 1

255     CONTINUE        !  end loop on inventory pollutants I

        ENAME = PROMPTMFILE( 
     &          'Enter logical name for POINTS output file',
     &          FSUNKN3, 'PNTS', 'EMSPOINT' )

C.......   Write out the point source emissions values:

        WRITE( *,92000 ) ' ', 'Writing out POINT output file...', ' '

        IF ( .NOT. WRITE3( ENAME, ALLVAR3, 0, 0,  IFIP ) ) THEN
            CALL M3EXIT( 'EMSPOINT', 0, 0, 
     &                   'Error writing output file "' //
     &                   ENAME( 1:TRIMLEN( ENAME ) ) // '"',  2 )
        END IF

C.......   Sort and Write out the point source emissions values:

        CALL SORTI1( NSRC, IDXSCC, ISCCA )

        SCCL = IMISS3
        DO 299 I = 1, NSRC

            SCC = ISCCA( IDXSCC( I ) )

            IF( SCC .NE. SCCL ) THEN
                WRITE( ADEV, 93030 ) SCC
            ENDIF

            SCCL = SCC
        
299     CONTINUE


999     CONTINUE          !  exit program:  normal completion

        MESG = 'Successful completion of Program EMSPOINT'
        CALL M3EXIT( 'EMSPOINT', 0, 0, MESG, 0 )


C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )

92010   FORMAT( 5X, A, :, I10 )


C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

93010   FORMAT( I2.2, I3.3 )

93020   FORMAT( I5.5 )

93030   FORMAT( I8.8 )

93040   FORMAT( I5, ',' ,I5, ',' ,I5 )

93050   FORMAT( I5.5, ',', I8.8, ',', 2 (I6, ','), 
     &          A15 , ',', A12 , ',', A40 )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I7, :, 1X ) )

94020   FORMAT( A, I7, X, A, I10, X, A, A, X, A, I4, X, A, I6 )

94030   FORMAT( 'Emis set to zero for source ', 
     &           I5, ':',I8, ':', A, ':', I3.3, ' at line ', I6 )

94040   FORMAT( A, I2.2 )

94050   FORMAT( A, I7, 10( A, :, E10.3, :, 1X ) )

94060   FORMAT( 10( A, :, E10.3, :, 1X ) )


        END

