
C Version "@(#)$Id$ $Source$ $Date$ 

         PROGRAM RAWAREA

C***********************************************************************
C  program body starts at line 146
C
C  DESCRIPTION:
C       Construct Models3/EDSS area source file from data contained in
C       EPS/AMS-style area source file.
C
C  PRECONDITIONS REQUIRED:
C       Sorted, cut-down input data for area sources.
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       Models-3 I/O
C       PROMPTFFILE, PROMPTMFILE, STR2INT, STR2REAL
C
C  REVISION  HISTORY:
C       Prototype  2/95 by CJC.
C       Version   11/95 by CJC reads, combines, and sorts multiple input 
C       files ARAW01, ARAW02, ..., each with its own inventory-year.
C       Supports both UTM and Lambert-based coordinate systems.
C       Version 8/96 by CJC does bounds-checking in loop 144; better
C       treatment of "missing".
C***********************************************************************

      IMPLICIT NONE


        INCLUDE 'ARDIMS3.EXT'   !  area-source dimensioning parameters
        INCLUDE 'CHDIMS3.EXT'   !  emis chem parms (inventory + model)
        INCLUDE 'GRDIMS3.EXT'   !  grid parameters
        INCLUDE 'TMDIMS3.EXT'   !  time parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.


C...........   EXTERNAL FUNCTIONS and their descriptions:

        LOGICAL         ENVYN
        INTEGER         FIND1
        INTEGER         GETNUM
        LOGICAL         GETYN
        INTEGER         JULIAN
        LOGICAL         DSCGRID
        INTEGER         PROMPTFFILE
        CHARACTER*16    PROMPTMFILE
        INTEGER         SECSDIFF
        INTEGER         STR2INT
        REAL            STR2REAL
        INTEGER         TRIMLEN

        EXTERNAL DSCGRID, ENVYN, FIND1, GETNUM, GETYN, JULIAN, 
     &           PROMPTFFILE, PROMPTMFILE, SECSDIFF, STR2INT, STR2REAL,
     &           TRIMLEN


C...........   PARAMETER:

        INTEGER         NSMAX
        PARAMETER     ( NSMAX = NASRC * NIPOL )


C...........   LOCAL VARIABLES and their descriptions:
C...........   NOTE that ASC (area-source-category) ID's are 10-digit
C...........   unsigned integers which may be treated as a leading 7-digit
C...........   field, and a trailing 3-digit field.  *7 and *3 arrays below
C...........   follow this scheme with parallel arrays

        INTEGER     TZONE       !  time zone
        INTEGER     FIP
        INTEGER     ID7
        INTEGER     ID3
        INTEGER     IPOL
        INTEGER     TPF
        INTEGER     INY
        INTEGER     SDT, EDT, STM, ETM
        INTEGER     MM, DD, HH, YY, DAYS
        REAL        CEFF
        REAL        REFF
        REAL        RPEN
        REAL        EMIS
        REAL        VMISS
        INTEGER     LFIP, LID7, LID3
                                 
C...........   Area Sources Table input unsorted copy (*A); 
C...........   sorted final version (image EMISREC of data record)

        INTEGER     NSRC             !  record-count
        INTEGER     INDEXA( NSMAX )  !  subscript table
        INTEGER     IFIPSA( NSMAX )  !  source FIPS (county) ID
        INTEGER     ISCC7A( NSMAX )  !  leading-7  digits of source ASC
        INTEGER     ISCC3A( NSMAX )  !  trailing-3 digits of source ASC
        INTEGER     ICODEA( NSMAX )  !  inventory species
        INTEGER     TPFLGA( NSMAX )  !  applicability of temporal profile types
        INTEGER     INVYRA( NSMAX )  !  year inventory taken
        REAL        CTLEFA( NSMAX )  !  control-factor (efficiency)
        REAL        RULEFA( NSMAX )  !  control-rule effectiveness
        REAL        RULPEA( NSMAX )  !  control-rule pentration
        REAL        EMISVA( NSMAX )  !  emissions values.

C.......   Common EMISREC holds an entire output record.  
C.......   Order of arrays in EMISREC _must_ match order of
C.......   variables in the output file.

        INTEGER     NAREA            !  current source-count
        INTEGER     IFIPS ( NASRC )  !  source FIPS (county) ID
        INTEGER     ISCC7 ( NASRC )  !  leading-7  digits of source ASC
        INTEGER     ISCC3 ( NASRC )  !  trailing-3 digits of source ASC
        INTEGER     ZONES ( NASRC )  !  time zones
        INTEGER     TPFLAG( NASRC )  !  applicability of temporal profile types
        INTEGER     INVYR ( NASRC )  !  year inventory taken
        REAL        CTLEFF( NASRC, NIPOL )  !  control-factor (efficiency)
        REAL        RULEFF( NASRC, NIPOL )  !  control-rule effectiveness
        REAL        RULPEN( NASRC, NIPOL )  !  control-rule pentration
        REAL        EMISV ( NASRC, NIPOL )  !  emissions values.

        COMMON  / EMISREC / IFIPS,  ISCC7,  ISCC3,  ZONES, INVYR, 
     &                      TPFLAG, CTLEFF, RULEFF, RULPEN, EMISV

C...........   Other local variables

        INTEGER         S, I, J, V      !  loop counters.
        INTEGER         ADEV            !  for area-source file
        INTEGER         ZDEV            !  for time-zone file
        INTEGER         LDEV            !  log-device
        INTEGER         IOS             !  I/O status
        INTEGER         IREC            !  input line (record) number
        INTEGER         IFIL            !  input file counter
        LOGICAL         VFLAG   !  input verification:  TRUE iff ERROR
        LOGICAL         SFLAG   !  input verification:  report missing species
        CHARACTER*16    SCRBUF  !  scratch name buffer
        CHARACTER*16    COORDN  !  coordinate system name
        CHARACTER*16    RNAME   !  logical name for emission source input file
        CHARACTER*16    ENAME   !  logical name for emission source output file
        CHARACTER*240   LINE    !  input line from AREA file
        CHARACTER*256   MESG    !  for M3EXIT() output


C***********************************************************************
C   begin body of program RAWAREA

        LDEV = INIT3()

        WRITE( *,92000 ) 
     &  ' ',
     &  'Program RAWAREA to take multiple EPS/AMS area source files',
     &  'and produce the AREA SOURCE EMISSIONS VECTOR file. ',
     &  ' ',
     &  'You will need to enter the logical names for the input and',
     &  'output files (and to have set them prior to program launch,',
     &  'using "setenv <logicalname> <pathname>").',
     &  'Optional checking that all species are reported for each ',
     &  'source may be turned on via "setenv RAW_SRC_CHECK Y".',
     &  ' ',
     &  'You may use END_OF-FILE (control-D) to quit the program',
     &  'during logical-name entry. Default responses are indicated',
     &  'in brackets [LIKE THIS].',
     &  ' '

        IF ( .NOT. GETYN( 'Continue with program?', .TRUE. ) ) THEN
            WRITE( *,92000 ) ' ', 'Ending program RAWAREA', ' '
            CALL EXIT( 0 )
        END IF

        SFLAG = ENVYN( 'RAW_SRC_CHECK', 
     &                 'RAWAREA check for missing species-records',
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
            CALL M3EXIT( 'RAWAREA', 0, 0, MESG, 2 )

        END IF


C.......   Get file name; open input time zone file

        ZDEV = PROMPTFFILE( 
     &           'Enter logical name for TIME ZONE file',
     &           .TRUE., .TRUE., 'ZONES', 'RAWAREA' )


C.......   Get file name; open input raw area source file

        WRITE( *,92000 ) 
     &  ' ', 
     &  'Now enter the list of logical file names for ',
     &  'RAW AREA SOURCE files...', 
     &  ' '
        IFIL =    0
        S    =    0
        INY  = 1989

101     CONTINUE

            IFIL = IFIL + 1
            WRITE( RNAME, 94040 ) 'ARAW', IFIL

            ADEV = PROMPTFFILE( 
     &          'Enter the NEXT FILE, or "NONE" to end list',
     &          .TRUE., .TRUE., RNAME, 'RAWAREA' )

            IF ( ADEV .LT. 0 ) GO TO 141        !  end of input files.

            INY = GETNUM( IMISS3, -IMISS3, INY+1,
     &                    'Enter INVENTORY YEAR for this file' )

C.......   Process this area source file

            IREC = 0

111         CONTINUE        !  head of the ADEV-read loop

                READ( ADEV, 93000, END=133, IOSTAT=IOS )LINE

                IF ( IOS .NE. 0 ) THEN
                    WRITE( *,92010 ) 
     &                  'Unit number', ADEV, 
     &                  'I/O Status ', IOS, 
     &                  'Line number', IREC
                    CALL M3EXIT( 'RAWAREA', 0, 0, 
     &                  'Error reading "' //
     &                  RNAME( 1:TRIMLEN( RNAME ) ) // '"', 2 )
                END IF

                IREC = IREC + 1

                IPOL = STR2INT ( LINE( 58:62 ) )
                IF ( IPOL .EQ. 42602 )  IPOL = 42603  !  NO2 ~~> NOX correction
                I    = FIND1(IPOL, NIPOL, PCODE )
                IF ( I .LE. 0 ) THEN
                    WRITE( MESG,94020 ) 
     &                  'Bad line', IREC, 
     &                  'Pollutant code "' // LINE( 162:171 ) // 
     &                  '" in "' //
     &                  RNAME( 1:TRIMLEN( RNAME ) ) // '"'
                    CALL M3MESG( MESG )
                    GO TO  111
                END IF

                EMIS = STR2REAL( LINE( 64:73 ) )
                IF ( EMIS .LT. 0.0 )  THEN
                    WRITE( MESG,94020 ) 
     &                  'Bad line', IREC, 
     &                  'Emis value "' // LINE( 162:171 ) // 
     &                  '" in "' //
     &                  RNAME( 1:TRIMLEN( RNAME ) ) // '"'
                    CALL M3MESG( MESG )
                    VFLAG = .TRUE.
                    GO TO  111
                END IF

                FIP  = STR2INT( LINE(  9:13 ) )
                ID7  = STR2INT( LINE( 26:32 ) )
                ID3  = STR2INT( LINE( 33:35 ) )

                CEFF = STR2REAL( LINE( 79:84 ) )
                IF ( CEFF .LT. 0 ) CEFF = 0.0

                REFF = STR2REAL( LINE( 86:91 ) )
                IF ( REFF .LT. 0 ) REFF = 100.0

                RPEN = STR2REAL( LINE( 93:97 ) )
                IF ( RPEN .LT. 0 ) RPEN = 100.0

                SDT  = STR2INT( LINE( 40:47 ) )
                EDT  = STR2INT( LINE( 49:56 ) )

                IF ( LINE( 37:38 ) .EQ. '  ' ) THEN             !  annual total

                    TPF = MTPRFAC * WTPRFAC             !  use month, week profiles

                ELSE IF ( LINE( 37:37 ) .EQ. 'P' ) THEN !  peak O3 or CO day

                    HH  = MOD( SDT        , 100 )
                    STM = 10000 * HH
                    DD  = MOD( SDT / 100  , 100 )
                    MM  = MOD( SDT / 10000, 100 )
                    YY  = 1900 + SDT / 1000000
                    SDT = 1900 + JULIAN( YY, MM, DD )   !  as Julian date

                    HH  = MOD( EDT        , 100 )
                    ETM = 10000 * HH
                    DD  = MOD( EDT / 100  , 100 )
                    MM  = MOD( EDT / 10000, 100 )
                    YY  = 1900 + EDT / 1000000 
                    EDT = 1900 + JULIAN( YY, MM, DD )        !  as Julian date
                    CALL NEXTIME( EDT, ETM, 10000 )     ! ETM means *through* hr HH

                    DAYS  = SECSDIFF( SDT, STM, EDT, ETM ) / SEC2DAY
                    IF ( DAYS .GT. 28 ) THEN
                        TPF = MTPRFAC * WTPRFAC     !  use month, week profiles
                    ELSE IF ( DAYS .GT. 1 ) THEN
                        TPF = WTPRFAC           !  use only week profiles
                    ELSE
                        TPF = 1                 !  use only hourly profiles
                    END IF                                

                    EMIS = DAY2YR * EMIS

                ELSE IF ( LINE( 37:37 ) .EQ. 'S' ) THEN !  special intervalC

                    HH  = MOD( SDT        , 100 )
                    STM = 10000 * HH
                    DD  = MOD( SDT / 100  , 100 )
                    MM  = MOD( SDT / 10000, 100 )
                    YY  = 1900 + SDT / 1000000
                    SDT = 1900 + JULIAN( YY, MM, DD )   !  as Julian date

                    HH  = MOD( EDT        , 100 )
                    ETM = 10000 * HH
                    DD  = MOD( EDT / 100  , 100 )
                    MM  = MOD( EDT / 10000, 100 )
                    YY  = 1900 + EDT / 1000000 
                    EDT = 1900 + JULIAN( YY, MM, DD )        !  as Julian date
                    CALL NEXTIME( EDT, ETM, 10000 )     ! ETM means *through* hr HH

                    DAYS  = SECSDIFF( SDT, STM, EDT, ETM ) / SEC2DAY

C                   IF ( DAYS .GT. 28 ) THEN
C                        TPF = MTPRFAC * WTPRFAC     !  use month, week profiles
C                   ELSE IF ( DAYS .GT. 1 ) THEN
C                        TPF = WTPRFAC          !  use only week profiles
C                   ELSE
C                        TPF = 1                        !  use only hourly profiles
C                   END IF

                    IF ( DAYS .GT. 1 ) THEN
                        WRITE( MESG,94020 ) 
     &                  'Bad line', IREC,  
     &                  'SPECIAL time period type for', 
     &                  DAYS, ' days not supported'
                        CALL M3MESG( MESG )
                        VFLAG = .TRUE.
                        GO TO  111
                    END IF
                    EMIS = DAY2YR * EMIS

                ELSE                                    !  unrecognized type

                    WRITE( MESG,94020 ) 
     &                  'Bad line', IREC, 
     &                  'Unsupported time period type "' // 
     &                  LINE( 37:38 ) // '"'
                    CALL M3MESG( MESG )
                    VFLAG = .TRUE.
                    GO TO  111

                END IF  !  if annual, peak, special, or unrecognized

                S = S + 1
                IF ( S .LE. NSMAX ) THEN

                    INDEXA( S ) = S !  index table for later SORTI3()
                    IFIPSA( S ) = FIP  
                    ISCC7A( S ) = ID7 
                    ISCC3A( S ) = ID3
                    ICODEA( S ) = I
                    TPFLGA( S ) = TPF
                    INVYRA( S ) = INY

                    EMISVA( S ) = EMIS
                    CTLEFA( S ) = CEFF
                    RULEFA( S ) = REFF
                    RULPEA( S ) = RPEN

                    LFIP = FIP  
                    LID7 = ID7 
                    LID3 = ID3

                END IF              !  if S in bounds

                GO TO  111          !  to head of ADEV-read loop

133         CONTINUE        !  end of the ADEV-read loop


            IF ( VFLAG ) THEN
                CALL M3EXIT( 'RAWAREA', 0, 0, 
     &                       'Error reading AREA SOURCE file.', 2 )
            END IF      !  if VFLAG

            GO TO 101   !  to process next input file

141     CONTINUE        !  all input files now processed

        NSRC = S
        IF ( S .GT. NSMAX ) THEN
            WRITE( *,92010 ) 
     &              'Actual record count         ', NSRC,
     &              'Max dimensioned record-count', NSMAX
            CALL M3EXIT( 'RAWAREA', 0, 0, 
     &        'Max record-count exceeded in AREA SOURCE file.', 2 )

        END IF


C.......   Use SORTI3() to perform an indirect sort by FIPS,SCC7,SCC3;
C.......   then permute the records according to the result:

        WRITE( *,92000 ) ' ', 'Sorting RAW AREA SOURCE file...', ' '

        CALL SORTI3( NSRC, INDEXA, IFIPSA, ISCC7A, ISCC3A )

        IF ( SFLAG ) THEN
            VMISS = BADVAL3
        ELSE
            VMISS = 0.0
        END IF

        VFLAG = .FALSE.
        LFIP  = IMISS3
        LID7  = IMISS3
        LID3  = IMISS3
        J     = 0

        DO  144  S = 1, NSRC

            I   = INDEXA( S )
            FIP = IFIPSA( I )
            ID7 = ISCC7A( I )
            ID3 = ISCC3A( I )

            IF ( FIP .NE. LFIP  .OR.
     &           ID7 .NE. LID7  .OR.
     &           ID3 .NE. LID3 ) THEN

                J = J + 1
                LFIP = FIP
                LID7 = ID7
                LID3 = ID3

                IF ( J .LE. NASRC ) THEN
                    IFIPS ( J ) = FIP
                    ISCC7 ( J ) = ID7
                    ISCC3 ( J ) = ID3
                    TPFLAG( J ) = TPFLGA( I )
                    INVYR ( J ) = INVYRA( I )
                END IF      !  if J in bounds

                DO  143  I = 1, NIPOL
                    EMISV ( J,I ) = VMISS
                    CTLEFF( J,I ) = BADVAL3
                    RULPEN( J,I ) = BADVAL3
                    RULEFF( J,I ) = BADVAL3
143             CONTINUE

            ELSE IF ( J .LE. NASRC ) THEN

                V = ICODEA( I )

                IF ( EMISV ( J,V ) .LE. 0.0 ) THEN 

                    EMISV ( J,V ) = EMISVA( I )
                    CTLEFF( J,V ) = CTLEFA( I )
                    RULEFF( J,V ) = RULEFA( I )
                    RULPEN( J,V ) = RULPEA( I )

                ELSE IF ( EMISVA( I ) .GT. 0.0 ) THEN

                    VFLAG = .TRUE.
                    WRITE( MESG,94020 )
     &                  'Duplicate emissions record:  FIP:', FIP,
     &                  'ASC:', ID7, ID3, 
     &                   'Species:', PCODE( V )
                    CALL M3MESG( MESG )

                END IF  !  if new pollutant for this source, or not

            END IF      !  if new source encountered; else if J in bounds

144     CONTINUE

        NAREA = J
        IF ( J .NE. NASRC ) THEN
            WRITE( MESG,94010 ) 
     &          'Actual=', NAREA, ', Dimensioned=', NASRC,       
     &          'numbers of sources do not match!'
            CALL M3EXIT( 'RAWAREA', 0, 0, MESG, 2 )            
        END IF

        IF ( VFLAG ) THEN
            CALL M3EXIT( 'RAWAREA', 0, 0, 
     &                   'Duplicate records found in input file', 2 )
        END IF

        IF ( SFLAG ) THEN

            VFLAG = .FALSE.

            DO  146  V = 1, NIPOL
            DO  145  S = 1, NASRC
 
               IF ( EMISV( S,V ) .LT. AMISS3 ) THEN

                    VFLAG = .TRUE.
                    WRITE( MESG,94020 )
     &                  'Missing emissions record:  FIP:', FIP,
     &                  'ASC:', ID7, ID3, 'Species:', PCODE( V )
                    CALL M3MESG( MESG )

                END IF

145         CONTINUE
146         CONTINUE

            IF ( VFLAG ) THEN
                CALL M3EXIT( 'RAWAREA', 0, 0,
     &                       'Missing species recs in input file', 2 )
            END IF

        END IF          !  sflag or not:  check for non-set EMIS values


C.......   Process the time zone file

        WRITE( *,92000 ) ' ', 'Processing TIME ZONE file...', ' '
        IREC = 0

155     CONTINUE        !  head of time-zone loop

            READ( ZDEV,*,END=199,IOSTAT=IOS ) FIP, TZONE 
            IREC = IREC + 1

            IF ( IOS .NE. 0 ) THEN
                WRITE( *,92010 ) 
     &              'Unit number', ZDEV, 
     &              'I/O Status ', IOS, 
     &              'Line number', IREC
                    CALL M3EXIT( 'RAWAREA', 0, 0, 
     &              'Error reading TIME ZONE file.', 2 )
            END IF
            
            IF ( FIP .EQ. 0 ) THEN              !  fallback -- all sources

                DO  166  S = 1, NASRC
                    ZONES( S ) = TZONE
166             CONTINUE

            ELSE IF ( MOD( FIP, 100 ) .EQ. 0 ) THEN     !  state-specific zone

                FIP = FIP / 1000
                DO  177  S = 1, NASRC
                    IF ( FIP .EQ. IFIPS( S ) / 1000 ) THEN
                        ZONES( S ) = TZONE
                    END IF
177             CONTINUE

            ELSE                                        !  county-specific zone

                DO  188  S = 1, NASRC
                    IF ( FIP .EQ. IFIPS( S ) ) THEN
                        ZONES( S ) = TZONE
                    END IF
188             CONTINUE

            END IF      !  if fip zero, or nn000, or not.

            GO TO  155

199     CONTINUE        !  end of time-zone loop


C.......   Get file name; open output area sources file
C.......   Note that coordinate system definition was put into 
C.......   FDESC3.EXT data structures by DSCGRID().

        FTYPE3D = GRDDED3
        SDATE3D = 0 !  n/a
        STIME3D = 0 !  n/a
        TSTEP3D = 0             !  time independent
        NVARS3D = 4 * NIPOL + 6
        NCOLS3D = 1
        NROWS3D = NASRC     !  number of rows = # of area sources.
        NLAYS3D = 1
        NTHIK3D = 1
        VGTYP3D = IMISS3
        VGTOP3D = AMISS3
        GDNAM3D = COORDN
        FDESC3D( 1 ) = 'NC Area Source emissions values.'
        FDESC3D( 2 ) = 'Inventory  type "' // LINE(  1: 1 ) // '"'
        FDESC3D( 3 ) = 'Projection year "' // LINE(  2: 3 ) // '"'
        FDESC3D( 4 ) = 'Base year       "' // LINE(  5: 6 ) // '"'
        FDESC3D( 5 ) = 'Emission   type "' // LINE(  7: 8 ) // '"'
        FDESC3D( 6 ) = 'Period     type "' // LINE( 37:38 ) // '"'
        FDESC3D( 7 ) = 'Beginning  Date "' // LINE( 40:47 ) // '"'
        FDESC3D( 8 ) = 'Ending     Date "' // LINE( 49:56 ) // '"'
           
        J = 1
        VNAME3D( J ) = 'FIP'
        UNITS3D( J ) = 'n/a'
        VDESC3D( J ) = 'FIP code for counties'
        VTYPE3D( J ) = M3INT
        J = J + 1
           
        VNAME3D( J ) = 'ASC7'
        UNITS3D( J ) = 'n/a'
        VDESC3D( J ) = 'Area Source Category code digits 1-7'
        VTYPE3D( J ) = M3INT
        J = J + 1
           
        VNAME3D( J ) = 'ASC3'
        UNITS3D( J ) = 'n/a'
        VDESC3D( J ) = 'Area Source Category code digits 8-10'
        VTYPE3D( J ) = M3INT
        J = J + 1

        VNAME3D( J ) = 'ZONES'
        UNITS3D( J ) = 'hours from GMT'
        VDESC3D( J ) = 'Std. time zones (0 for GMT, 5 for Eastern)'
        VTYPE3D( J ) = M3INT
        J = J + 1

        VNAME3D( J ) = 'INVYR'
        UNITS3D( J ) = 'year AD'
        VDESC3D( J ) = 'Year of inventory for this record'
        VTYPE3D( J ) = M3INT
        J = J + 1

        VNAME3D( J ) = 'TPFLAG'
        UNITS3D( J ) = 'T|2? T|3?'
        VDESC3D( J ) = 'Use week(2), month(3) temporal profiles or not'
        VTYPE3D( J ) = M3INT
        J = J + 1

        DO  211  I = 1, NIPOL
            VNAME3D( J ) = 'CTLEFF_' // EINAM( I )
            UNITS3D( J ) = 'n/a'
            VDESC3D( J ) = 
     &      'control efficiency (in [0,100], or "MISSING": < -9.0E36)'
            VTYPE3D( J ) = M3REAL
            J = J + 1
211     CONTINUE
           
        DO  212  I = 1, NIPOL
            VNAME3D( J ) = 'RULEFF_' // EINAM( I )
            UNITS3D( J ) = 'n/a'
            VDESC3D( J ) = 
     &      'Rule effectiveness (in [0,100], or "MISSING": < -9.0E36)'
            VTYPE3D( J ) = M3REAL
            J = J + 1
212     CONTINUE
           
        DO  213  I = 1, NIPOL
            VNAME3D( J ) = 'RULPEN_' // EINAM( I )
            UNITS3D( J ) = 'n/a'
            VDESC3D( J ) = 
     &      'Rule penetration (in [0,100], or "MISSING": < -9.0E3J)'
            VTYPE3D( J ) = M3REAL
            J = J + 1
213     CONTINUE
           
        DO  214  I = 1, NIPOL
            VNAME3D( J ) = EINAM( I )
            UNITS3D( J ) = 'tons/year'
            VDESC3D( J ) = LINE( 1:1 ) // LINE( 37:38 ) 
     &                                 // ' emissions totals'
            VTYPE3D( J ) = M3REAL
            J = J + 1
214     CONTINUE

        ENAME = PROMPTMFILE( 
     &          'Enter logical name for AREA SOURCE output file',
     &          FSUNKN3, 'AREA', 'RAWAREA' )


C.......   Write out the area source emissions values:

        WRITE( *,92000 ) ' ', 'Writing out AREA SOURCES file...', ' '

        IF ( .NOT. WRITE3( ENAME, 'ALL', 0, 0,  IFIPS ) ) THEN
            CALL M3EXIT( 'RAWAREA', 0, 0,  
     &                   'Error writing AREA OUTPUT file' , 2 )
        END IF


999   CONTINUE          !  exit program

      CALL M3EXIT( 'RAWAREA', 0, 0,
     &             'Normal completion of Program RAWAREA', 0 )


C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )

92010   FORMAT( 5X, A, :, I10 )


C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )


C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 3( A, :, I7, :, 1X ) )

94020   FORMAT( A, 1X, I5, 1X, A, I7.7, I3.3, 1X, A, I5 )

94040   FORMAT( A, I2.2 )

        END

