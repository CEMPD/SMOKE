C ****************************************************************************
C *
C * Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C *                System
C * File: @(#)$Id$
C *
C * COPYRIGHT (C) 1996, MCNC--North Carolina Supercomputing Center
C * All Rights Reserved
C *
C * See file COPYRIGHT for conditions of use.
C *
C * Environmental Programs Group
C * MCNC--North Carolina Supercomputing Center
C * P.O. Box 12889
C * Research Triangle Park, NC  27709-2889
C *
C * env_progs@mcnc.org
C *
C * Pathname: $Source$
C * Last updated: $Date$ 
C *
C ****************************************************************************/

         PROGRAM RAWPOINT

C***********************************************************************
C  program body starts at line 206
C
C  DESCRIPTION:
C       Sort data contained in EPS/AFS-style point source file and 
C       construct Models3/EDSS point source file
C
C  PRECONDITIONS REQUIRED:
C       Input data for point sources.
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       Models-3 I/O
C       FIND1, PROMPTFFILE, PROMPTMFILE, GETYN, TRIMLEN
C
C  REVISION  HISTORY:
C       Prototype  6/95 by CJC.
C       Version  10/95 does sorting internally on FIPS,SCC,PLANT,STACK.
C       Version  11/95 reads, combines, and sorts multiple input files
C       PRAW01, PRAW02, ..., each with its own inventory-year.
C       Supports both UTM, LAT-LON, and Lambert-based coordinate systems.
C***********************************************************************

      IMPLICIT NONE

C...........   INCLUDES:
        
        INCLUDE 'PTDIMS3.EXT'   !  point-source dimensioning parameters
        INCLUDE 'GRDIMS3.EXT'   !  grid dimensioning parameters
        INCLUDE 'CHDIMS3.EXT'   !  emis chem info (both inventory and model)
        INCLUDE 'TMDIMS3.EXT'   !  time related parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.


C...........   PARAMETER

        INTEGER         NSMAX
        PARAMETER     ( NSMAX  = NPSRC * NIPOL )


C...........   EXTERNAL FUNCTIONS and their descriptions:
        
        LOGICAL         DSCGRID
        LOGICAL         ENVYN
        INTEGER         FIND1   !  returns -1 for failure
        INTEGER         GETNUM
        LOGICAL         GETYN
        INTEGER         JULIAN
        LOGICAL         SETLAM
        LOGICAL         LL2LAM
        LOGICAL         UTM2LAM
        INTEGER         PROMPTFFILE
        CHARACTER*16    PROMPTMFILE
        INTEGER         SECSDIFF
        INTEGER         STR2INT
        REAL            STR2REAL
        INTEGER         TRIMLEN
        
        EXTERNAL DSCGRID, ENVYN, FIND1, GETNUM, GETYN, JULIAN,
     &           SETLAM, LL2LAM, UTM2LAM, PROMPTFFILE, PROMPTMFILE,
     &           SECSDIFF, STR2INT, STR2REAL, TRIMLEN


C...........   LOCAL VARIABLES and their descriptions:
C...........   NOTE that SSC (point-source-category) ID's are 10-digit
C...........   unsigned integers which may be treated as a leading 7-digit
C...........   field, and a trailing 3-digit field.  *7 and *3 arrays below
C...........   follow this scheme with parallel arrays

        INTEGER     DAY2YR      !  Conversion factor for day to year
        INTEGER     UZONE       !  UTM zone for the file as a whole
        INTEGER     TZONE       !  time zone
        INTEGER     FIP
        INTEGER     SCC
        INTEGER     SIC
        INTEGER     PLANT
        INTEGER     STACK
        INTEGER     ZONE
        INTEGER     TPF
        INTEGER     INY
        INTEGER     SDT, EDT, STM, ETM
        INTEGER     MM, DD, HH, DAYS
        REAL        CEFF
        REAL        REFF
        REAL        RPEN
        REAL        SLOC
        REAL        XLOC
        REAL        YLOC
        REAL        EMIS
        REAL        STKH
        REAL        STKD
        REAL        STKT
        REAL        STKV
        INTEGER     IPOL
        INTEGER     LFIP, LSCC, LPLT, LSTK
                                 
C...........   Point Sources Table input unsorted copy (*A); 
C...........   sorted final version (image EMISREC of data record)
        
        INTEGER     NSRC             !  current source-count
        INTEGER     INDEXA( NSMAX )  !  subscript table for SORTI4()
        INTEGER     IDXSCC( NSMAX )  !  subscript table for SORTI1() for SCC
        INTEGER     IFIPA ( NSMAX )  !  source FIPS (county) ID
        INTEGER     ISCCA ( NSMAX )  !  source SCC
        INTEGER     ISICA ( NSMAX )  !  source SIC
        INTEGER     IPLANA( NSMAX )  !  Plant ID
        INTEGER     ISTACA( NSMAX )  !  Stack ID
        INTEGER     ICODEA( NSMAX )  !  spcecies-subscript
        INTEGER     TPFLGA( NSMAX )  !  applicability of temporal profile types
        INTEGER     INVYRA( NSMAX )  !  inventory year
        REAL        XLOCAA( NSMAX )  !  UTM X-location (m)
        REAL        YLOCAA( NSMAX )  !  UTM Y-location (m)
        REAL        STKHTA( NSMAX )  !  stack height   (m)
        REAL        STKDMA( NSMAX )  !  stack diameter (m)
        REAL        STKTKA( NSMAX )  !  exhaust temperature (deg K)
        REAL        STKVEA( NSMAX )  !  exhaust velocity    (m/s)
        REAL        RULPEA( NSMAX )  !  rule penetration   fraction
        REAL        RULEFA( NSMAX )  !  rule effectiveness fraction
        REAL        CTLEFA( NSMAX )  !  control efficiency fraction
        REAL        EMISVA( NSMAX )  !  emissions values (tons/yr)

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

        INTEGER         S, I, J, V !  loop counters.
        INTEGER         ADEV    !  Unit number for output actual SCCs
        INTEGER         PDEV    !  for point-source file
        INTEGER         ZDEV    !  for time zone file
        INTEGER         LDEV    !  log-device

        LOGICAL         SFLAG   !  input verification:  report missing species
        LOGICAL         VFLAG   !  input verification:  TRUE iff ERROR
        LOGICAL         WFLAG   !  input verification:  convert LON to west hemis

        REAL            XX, YY  !  scratch location variables
                                
        INTEGER         IOS     !  I/O status
        INTEGER         IREC    !  input line (record) number
        INTEGER         IFIL    !  input file counter
        CHARACTER*16    SCRBUF  !  scratch buffer
        CHARACTER*16    COORDN  !  coordinate system name
        CHARACTER*16    RNAME   !  input file name buffer
        CHARACTER*256   MESG    !  text for M3EXIT()
        CHARACTER*240   LINE    !  input line from POINT file


C***********************************************************************
C   begin body of program RAWPOINT
        
        LDEV = INIT3()

        WRITE( *,92000 ) 
     &  ' ',
     &  'Program RAWPOINT to take the sorted EPS/AFS point source ',
     &  'file, and produce the POINT SOURCE EMISSIONS VECTOR file. ',
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
            WRITE( *,92000 ) ' ', 'Ending program RAWPOINT', ' '
            CALL EXIT( 0 )
        END IF

        SFLAG = ENVYN( 'RAW_SRC_CHECK', 
     &                 'RAWPOINT check for missing species-records',
     &                 .FALSE.,
     &                 IOS )
     
        WFLAG = ENVYN( 'WEST_HSPHERE', 
     &                 'Western hemisphere flag', .TRUE., IOS )

        IF( WFLAG ) THEN
          CALL M3MESG( 'NOTE: All non-negative longitude values will' )
          CALL M3MESG( '      be converted to Western Hemisphere.' )
        ENDIF

C.......   Get the coordinate system name and parameters
C.......   to put into file header:

        IF ( .NOT. DSCGRID( GRDNM, COORDN, GDTYP3D, 
     &              P_ALP3D, P_BET3D, P_GAM3D, XCENT3D, YCENT3D,
     &              XORIG3D, YORIG3D, XCELL3D, YCELL3D, 
     &              NCOLS3D, NROWS3D, NTHIK3D ) ) THEN

            SCRBUF = GRDNM
            MESG   = '"' // SCRBUF( 1:TRIMLEN( SCRBUF ) ) //
     &               '" not found in GRIDDESC file'
            CALL M3EXIT( 'RAWPOINT', 0, 0, MESG, 2 )

        END IF

        IF      ( GDTYP3D .EQ. UTMGRD3 ) THEN
            UZONE = NINT( P_ALP3D )
        ELSE IF ( GDTYP3D .EQ. LAMGRD3 ) THEN
            IF ( .NOT. SETLAM( P_ALP3D, P_BET3D, P_GAM3D, 
     &                         XCENT3D, YCENT3D ) ) THEN
                MESG = 'Bad grid ' // GDNAM3D
                CALL M3EXIT( 'RAWPOINT', 0, 0, MESG, 2 )
            END IF
        ELSE IF ( GDTYP3D .EQ. LATGRD3 ) THEN
            !  do nothing -- no setup necessary
        ELSE
            WRITE( MESG, 94010 ) 
     &          'Unsupported coordinate system type', GDTYP3D
            CALL M3EXIT( 'RAWPOINT', 0, 0, MESG, 2 )
        END IF


C.......   Get file name; open input time zone file

        ZDEV = PROMPTFFILE( 
     &         'Enter logical name for TIME ZONE file',
     &         .TRUE., .TRUE., 'ZONES', 'RAWPOINT' )

C.......   Read the point source files:
C.......   Get file name; open input raw point source file

        WRITE( *,92000 ) 
     &  ' ', 
     &  'Now enter the list of logical file names for ',
     &  'RAW POINT SOURCE files...', 
     &  ' '
        IFIL =    0
        S    =    0
        INY  = 1989

101     CONTINUE

            IFIL = IFIL + 1
            WRITE( RNAME, 94040 ) 'PRAW', IFIL

            PDEV = PROMPTFFILE( 
     &              'Enter the NEXT FILE, or "NONE" to end list',
     &              .TRUE., .TRUE., RNAME, 'RAWPOINT' )
            IF ( PDEV .LT. 0 ) GO TO  123       !  end of loop readin input files

            INY = GETNUM( IMISS3, -IMISS3, INY + 1,
     &                    'Enter INVENTORY YEAR for this file' )

            IREC  =  0
            VFLAG = .FALSE.

C.............   Calculate DAY2YR factor
            I = INY / 1000        !  year-number
            IF ( MOD( I,4 ) .NE. 0 ) THEN           !  nonleap years
                DAY2YR = 365.0
            ELSE IF ( MOD( I,100 ) .NE. 0 ) THEN    !  noncentury leap years
                DAY2YR = 366.0
            ELSE IF ( MOD( I,400 ) .NE. 0 ) THEN    !  century nonleap years
                DAY2YR = 365.0
            ELSE                                    !  leap centuries
                DAY2YR = 366.0
            END IF

111         CONTINUE        !  head of the PDEV-read loop

                READ( PDEV, 93000, END=122, IOSTAT=IOS )LINE
                IREC = IREC + 1

                IF ( IOS .NE. 0 ) THEN
                    WRITE( MESG, 94010 ) 
     &                  'Error ', IOS,  'reading "' // 
     &                  RNAME( 1:TRIMLEN( RNAME ) ) //
     &                  '" at line', IREC
                    CALL M3EXIT( 'RAWPOINT', 0, 0, MESG, 2 )
                END IF

                IPOL  = STR2INT ( LINE( 156:160 ) )
                IF ( IPOL .EQ. 42602 ) IPOL = 42603         !  fix NO2 ~~> NOx
                I = FIND1( IPOL, NIPOL, PCODE )
                IF ( I  .LT. 0 ) THEN
                    WRITE( MESG,94010 ) 
     &                  'Bad line', IREC, 
     &                  'Pollutant code "' // LINE( 156:160 ) // 
     &                  '" in "' //
     &                  RNAME( 1:TRIMLEN( RNAME ) ) // '"'
                    CALL M3MESG( MESG )
                    GO TO  111      !  to head of loop
                END IF

                EMIS  = STR2REAL( LINE( 162:171 ) )
                IF ( EMIS .LT. 0.0 )  THEN
                    WRITE( MESG,94010 ) 
     &                  'Bad line', IREC, 
     &                  '" in "' //
     &                  RNAME( 1:TRIMLEN( RNAME ) ) // 
     &                  '"--emis value "' // LINE( 162:171 ) // '"'
                    CALL M3MESG( MESG )
                    VFLAG = .TRUE.
                    GO TO  111
                END IF

                FIP   = STR2INT ( LINE(  12: 16 ) )
                SIC   = STR2INT ( LINE(  24: 27 ) )
                SCC   = STR2INT ( LINE(  29: 38 ) )
                PLANT = STR2INT ( LINE(  40: 44 ) )
                STACK = STR2INT ( LINE(  46: 48 ) )
                XLOC  = STR2REAL( LINE(  78: 87 ) )
                YLOC  = STR2REAL( LINE(  89: 98 ) )
                ZONE  = STR2INT ( LINE( 100:101 ) )
                STKH  = STR2REAL( LINE( 103:107 ) )
                STKD  = STR2REAL( LINE( 109:113 ) )
                STKT  = STR2REAL( LINE( 115:119 ) )
                STKV  = STR2REAL( LINE( 121:125 ) )

                CEFF  = STR2REAL( LINE( 177:182 ) )
                IF ( CEFF .LT. 0.0 )  CEFF = 0.0

                REFF  = STR2REAL( LINE( 184:189 ) )
                IF ( REFF .LT. 0.0 )  REFF = 100.0

                RPEN  = STR2REAL( LINE( 191:196 ) )
                IF ( RPEN .LT. 0.0 )  RPEN = 100.0

                SDT  = STR2INT( LINE( 60:67 ) )
                EDT  = STR2INT( LINE( 69:76 ) )

                IF ( LINE( 57:58 ) .EQ. '  ' ) THEN         !  annual total

                    TPF = MTPRFAC * WTPRFAC         !  use month, week profiles

                ELSE IF ( LINE( 57:57 ) .EQ. 'P' ) THEN     !  peak O3 or CO day

                    HH  = MOD( SDT        , 100 )
                    STM = 10000 * HH
                    DD  = MOD( SDT / 100  , 100 )
                    MM  = MOD( SDT / 10000, 100 )
                    YY  = 1900 + SDT / 1000000
                    SDT = 1900 + JULIAN( YY, MM, DD )       !  as Julian date

                    HH  = MOD( EDT        , 100 )
                    ETM = 10000 * HH
                    DD  = MOD( EDT / 100  , 100 )
                    MM  = MOD( EDT / 10000, 100 )
                    YY  = 1900 + EDT / 1000000 
                    EDT = 1900 + JULIAN( YY, MM, DD )        !  as Julian date
                    CALL NEXTIME( EDT, ETM, 10000 ) ! ETM means *through* hr HH

                    DAYS  = SECSDIFF( SDT, STM, EDT, ETM ) / SEC2DAY
                    IF ( DAYS .GT. 28 ) THEN
                        TPF = MTPRFAC * WTPRFAC     !  use month, week profiles
                    ELSE IF ( DAYS .GT. 1 ) THEN
                        TPF = WTPRFAC               !  use only week profiles
                    ELSE
                        TPF = 1                     !  use only hourly profiles
                    END IF                                

                    EMIS = DAY2YR * EMIS

                ELSE IF ( LINE( 57:57 ) .EQ. 'S' ) THEN     !  special interval

                    HH  = MOD( SDT        , 100 )
                    STM = 10000 * HH
                    DD  = MOD( SDT / 100  , 100 )
                    MM  = MOD( SDT / 10000, 100 )
                    YY  = 1900 + SDT / 1000000
                    SDT = 1900 + JULIAN( YY, MM, DD )       !  as Julian date

                    HH  = MOD( EDT        , 100 )
                    ETM = 10000 * HH
                    DD  = MOD( EDT / 100  , 100 )
                    MM  = MOD( EDT / 10000, 100 )
                    YY  = 1900 + EDT / 1000000 
                    EDT = 1900 + JULIAN( YY, MM, DD )        !  as Julian date
                    CALL NEXTIME( EDT, ETM, 10000 ) ! ETM means *through* hr HH

                    DAYS  = SECSDIFF( SDT, STM, EDT, ETM ) / SEC2DAY
                    IF ( DAYS .GT. 28 ) THEN
                        TPF = MTPRFAC * WTPRFAC     !  use month, week profiles
                    ELSE IF ( DAYS .GT. 1 ) THEN
                    ELSE
                        TPF = 1                     !  use only hourly profiles
                    END IF

                    EMIS = 0.0      !  note ZERO for day-specific "SP" records

                    WRITE( MESG, 94030 ) FIP, SCC, PLANT, STACK, IREC
                    CALL M3MESG( MESG )

                ELSE                                        !  unrecognized type

                    WRITE( MESG,94010 ) 
     &                  'Bad line', IREC, 
     &                  '" in "' //
     &                  RNAME( 1:TRIMLEN( RNAME ) ) // 
     &                  '"--Unsupported time period type "' // 
     &                  LINE( 57:58 ) // '"'
                    CALL M3MESG( MESG )
                    VFLAG = .TRUE.
                    GO TO  111          !  to head of PDEV-read loop

                END IF          !  tests on record type line( 57:58 )

C................  Make changes in XLOC and YLOC if they are provided in LAT-LON

               IF( ZONE .EQ. IMISS3 ) THEN

C....................  Trade XLOC/YLOC because of EPS format
                   SLOC = YLOC
                   YLOC = XLOC
                   XLOC = SLOC

C....................  Convert longitude if WEST_HSPHERE environment
C....................  variable is true, and lonitude is negative.
                   IF( WFLAG .AND. XLOC .GT. 0 ) THEN

                       XLOC = -XLOC

                   ENDIF
               ENDIF

C...........   Normalize the UTM zone and coordinates
C...........   NOTE that AFS files have the temerity to use KM units
C...........   and still call what they are doing "UTM":

                IF ( GDTYP3D .EQ. UTMGRD3 ) THEN

                    IF ( ZONE .EQ. IMISS3 ) THEN	!  lat-lon input

                        CALL LL2UTM( XLOC, YLOC, UZONE, XX, YY )
                        XLOC = XX
                        YLOC = YY

                    ELSE IF ( ZONE .NE. UZONE ) THEN

                        XLOC = 1000.0 * XLOC    ! conversion to meters
                        YLOC = 1000.0 * YLOC    ! conversion to meters
                        CALL UTM2LL( XLOC, YLOC, ZONE, XX,   YY )
                        CALL LL2UTM( XX,   YY,  UZONE, XLOC, YLOC )
                        ZONE = UZONE

                    ELSE

                        XLOC = 1000.0 * XLOC    ! conversion to meters
                        YLOC = 1000.0 * YLOC    ! conversion to meters

                    END IF  !  if zone missing (lat-lon), or different, or not

                ELSE IF ( GDTYP3D .EQ. LAMGRD3 ) THEN

                    IF ( ZONE .EQ. IMISS3 ) THEN        !  lat-lon input

                        IF ( .NOT. LL2LAM( XLOC, YLOC, XX, YY ) ) THEN
                            VFLAG = .TRUE.
                            WRITE( MESG,94010 ) 
     &                          'Bad line', IREC, 
     &                          '" in "' //
     &                          RNAME( 1:TRIMLEN( RNAME ) ) // 
     &                          'LL2LAM() failed for "' // 
     &                          LINE( 57:58 ) // '"'
                            CALL M3MESG( MESG )
                            VFLAG = .TRUE.
                            GO TO  111          !  to head of PDEV-read loop
                        END IF

                        XLOC = XX
                        YLOC = YY

                    ELSE				!  utm input

                        XX = 1000.0 * XLOC    ! conversion to meters
                        YY = 1000.0 * YLOC    ! conversion to meters

                        IF ( .NOT. UTM2LAM( XX, YY, ZONE, 
     &                                      XLOC, YLOC ) ) THEN
                            VFLAG = .TRUE.
                            WRITE( MESG,94010 ) 
     &                          'Bad line', IREC, 
     &                          '" in "' //
     &                          RNAME( 1:TRIMLEN( RNAME ) ) // 
     &                          'LL2LAM() failed for "' // 
     &                          LINE( 57:58 ) // '"'
                            CALL M3MESG( MESG )
                            VFLAG = .TRUE.
                            GO TO  111          !  to head of PDEV-read loop

                        END IF

                    END IF	!  if zone missing (lat-lon) or not (utm)

                ELSE IF ( GDTYP3D .EQ. LATGRD3 ) THEN

                    IF ( ZONE .NE. IMISS3 ) THEN

                        XX = 1000.0 * XLOC    ! conversion to meters
                        YY = 1000.0 * YLOC    ! conversion to meters
                        CALL UTM2LL( XX, YY, ZONE, XLOC, YLOC )

                    END IF	! if UTM input location

                END IF		!  if coord type UTM or Lambert or Lat-Lon

                S = S + 1

                IF ( S .LE. NSMAX ) THEN

                    INDEXA( S ) = S         !  index-table for use in SORTI4()
                    IDXSCC( S ) = S         !  index-table for use in SORTI1()

                    IFIPA ( S ) = FIP
                    ISICA ( S ) = SIC
                    ISCCA ( S ) = SCC
                    IPLANA( S ) = PLANT
                    ISTACA( S ) = STACK
                    ICODEA( S ) = I         !  subscript into ACODE(*)
                    INVYRA( S ) = INY
                    TPFLGA( S ) = TPF
                    XLOCAA( S ) = XLOC
                    YLOCAA( S ) = YLOC
                    STKHTA( S ) = STKH
                    STKDMA( S ) = STKD
                    STKTKA( S ) = STKT
                    STKVEA( S ) = STKV

                    CTLEFA( S ) = CEFF
                    RULEFA( S ) = REFF
                    RULPEA( S ) = RPEN
                    EMISVA( S ) = EMIS

                END IF

                GO TO  111          !  to head of PDEV-read loop


122         CONTINUE        !  end of the PDEV-read loop

            IF ( VFLAG ) THEN
                CALL M3EXIT( 'RAWPOINT', 0, 0, 
     &                       'Error reading POINT SOURCE file.', 2 )
            END IF

            GO TO 101   !  to head of loop on input files
        
123     CONTINUE        !  end loop on input files

        NSRC = S
        IF ( NSRC .GT. NSMAX ) THEN
            WRITE( *,92010 ) 
     &          'Record count                ', IREC,
     &          'Max dimensioned record-count', NSMAX
            CALL M3EXIT( 'RAWPOINT', 0, 0, 
     &        'Max record-count exceeded in POINT SOURCE file.', 2 )
           
        END IF

        ADEV = PROMPTFFILE(
     &          'Enter the name of the ACTUAL SCC output file',
     &          .FALSE., .TRUE., 'PSCC', 'EMSPOINT' )
        
C.......   Use SORTI4() to perform an indirect sort by FIPS,SCC,PLANT,STACK
C.......   then permute the records according to the result:

        WRITE( *,92000 ) ' ', 'Processing RAW POINT SOURCE file...', ' '

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

            IF ( FIP   .NE. LFIP  .OR.
     &           SCC   .NE. LSCC  .OR.
     &           PLANT .NE. LPLT  .OR.
     &           STACK .NE. LSTK  ) THEN

                J = J + 1

                LFIP = FIP
                LSCC = SCC
                LPLT = PLANT
                LSTK = STACK

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

            WRITE( MESG,94010 ) 
     &          'Actual      number of sources', NPOINT,
     &          'Dimensioned number of sources', NPSRC,       
     &          'Do not match!'
            CALL M3EXIT( 'RAWPOINT', 0, 0, MESG, 2 )

        END IF

        IF ( SFLAG ) THEN       !  if missing records are fatal:

            VFLAG = .FALSE.

            DO  146  V = 1, NIPOL
            DO  145  S = 1, NPSRC

               IF ( EMISV( S,V ) .LT. AMISS3 ) THEN

                    VFLAG = .TRUE.
                    WRITE( MESG,94020 )
     &                  'Missing record:  FIP:', IFIP( S ),
     &                  'SCC:',     ISCC( S ), 
     &                  'Plant:',   IPLANT( S ),
     &                  'Stack:',   ISTACK( S ),
     &                  'Species:', PCODE( V )
                    CALL M3MESG( MESG )

                END IF

145         CONTINUE
146         CONTINUE

            IF ( VFLAG ) THEN
                CALL M3EXIT( 'RAWPOINT', 0, 0,
     &                       'Missing species recs in input file', 2 )
            END IF

        ELSE            !  "missing" ~~~> 0

            DO  148  V = 1, NIPOL
            DO  147  S = 1, NPSRC
                IF ( EMISV( S,V ) .LT. AMISS3 ) EMISV( S,V ) = 0.0
147         CONTINUE
148         CONTINUE

        END IF          !  sflag or not:  check for non-set EMIS values


C.......   Process the time zone file

        WRITE( *,92000 ) ' ', 'Processing TIME ZONE file...', ' '
        IREC = 0

144     CONTINUE        !  head of time-zone loop

            READ( ZDEV,*,END=188,IOSTAT=IOS ) FIP, TZONE 
            IREC = IREC + 1

            IF ( IOS .NE. 0 ) THEN
                WRITE( *,92010 ) 
     &              'Unit number', ZDEV, 
     &              'I/O Status ', IOS, 
     &              'Line number', IREC
                    CALL M3EXIT( 'RAWPOINT', 0, 0, 
     &              'Error reading TIME ZONE file.', 2 )
            END IF

            IF ( FIP .EQ. 0 ) THEN              !  fallback -- all sources

                DO  155  S = 1, NPSRC
                    TZONES( S ) = TZONE
155             CONTINUE

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

            GO TO  144

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
        NROWS3D = NPSRC     !  number of rows = # of point sources.
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
     &          FSUNKN3, 'PNTS', 'RAWPOINT' )


C.......   Write out the point source emissions values:

        WRITE( *,92000 ) ' ', 'Writing out POINT output file...', ' '

        IF ( .NOT. WRITE3( ENAME, ALLVAR3, 0, 0,  IFIP ) ) THEN
            CALL M3EXIT( 'RAWPOINT', 0, 0, 
     &                   'Error writing output file "' //
     &                   ENAME( 1:TRIMLEN( ENAME ) ) // '"',  2 )
        END IF

C.......   Sort and Write out the point source SCCs:

        CALL SORTI1( NSRC, IDXSCC, ISCCA )

        LSCC = IMISS3
        DO 299 I = 1, NSRC

            SCC = ISCCA( IDXSCC( I ) )

            IF( SCC .NE. LSCC ) THEN
                WRITE( ADEV, 93030 ) SCC
            ENDIF

            LSCC = SCC

299     CONTINUE

999     CONTINUE          !  exit program:  normal completion

        MESG = 'Successful completion of Program RAWPOINT'
        CALL M3EXIT( 'RAWPOINT', 0, 0, MESG, 0 )


C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )

92010   FORMAT( 5X, A, :, I10 )


C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

93010   FORMAT( A16 )

93030   FORMAT( I8.8 )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 3( A, :, I7, :, 1X ) )

94020   FORMAT( A, I7, X, A, I10, X, A, I8, X, A, I4, X, A, I6 )

94030   FORMAT( 'Emis set to zero for source ', 
     &           I5, ':',I8, ':', I8, ':', I3.3, ' at line ', I6 )

94040   FORMAT( A, I2.2 )


        END

