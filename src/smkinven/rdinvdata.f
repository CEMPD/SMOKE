
        SUBROUTINE RDINVDATA( FDEV, FNAME, NRAWBP, NONPOINT )

C***********************************************************************
C  subroutine body starts at line 133
C
C  DESCRIPTION:
C      This subroutine controls reading an ASCII inventory file for any source
C      category from one of many formats.  It determines the format and
C      calls the appropriate reader subroutines. It controls the looping
C      through multiple files when a list-formatted file is used as input.
C      This routine only reads the data (emissions and activities) from the
C      inventories.
C
C  PRECONDITIONS REQUIRED:
C      Input file unit FDEV opened
C      Inventory pollutant list created: MXIDAT, INVDCOD, and INVDNAM
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C      Subroutines:
C      Functions:
C
C  REVISION  HISTORY:
C       Created 1/03 by C. Seppanen (based on rdinven.f)
C
C       Version June 2016 by Carlie Coats:  add fugitive-emissions properties
C**************************************************************************
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
C***************************************************************************

C...........   MODULES for public variables
C...........   This module is the inventory arrays
        USE MODSOURC, ONLY: NSTRECS, SRCSBYREC, RECIDX,         ! file numbers and records
     &                      INDEXA, IPOSCODA, ICASCODA,         ! unsorted source characteristics
     &                      INRECA, SRCIDA, POLVLA,
     &                      INVYR, IDIU, IWEK, TPFLAG,          ! sorted integer characteristics
     &                      XLOCA, YLOCA, XLOC1, YLOC1, XLOC2,  ! sorted real characteristics
     &                      YLOC2, STKHT, STKDM, STKTK, STKVE,
     &                      CORIS, CBLRID, CPDESC, CERPTYP,     ! sorted character characteristics
     &                      CMACT, CNAICS, CSRCTYP, CNEIUID, CEXTORL,
     &                      CINTGR, CISIC, CSHAPE,
     &                      FUGHGT, FUGWID, FUGLEN, FUGANG

C.........  This module contains the information about the source category
        USE MODINFO, ONLY: CATEGORY, NEM, NDY, NEF, NCE, NRE, NRP,
     &                     NC1, NC2, NPPOL, NSRC, NPACT, INV_MON

C.........  This module contains the lists of unique inventory information
        USE MODLISTS, ONLY: FILFMT, LSTSTR, MXIDAT, INVDCNV, INVDNAM,
     &                      NUNIQCAS, UNIQCAS, UCASNKEP, ITNAMA,
     &                      SCASIDX, UCASIDX, UCASNPOL, ITKEEPA, ITFACA,
     &                      EMISBYCAS, RECSBYCAS, EMISBYPOL, INVSTAT,
     &                      NINVTBL, ITCASA, FIREFLAG, FIREFF10

C.........  This module is for mobile-specific data
        USE MODMOBIL, ONLY: SCCMAPFLAG, SCCMAPLIST,
     &                      NSCCMAP, EXCLSCCFLAG

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'CONST3.EXT'    !  physical and mathematical constants
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:

        LOGICAL         BLKORCMT
        LOGICAL         CHKINT
        LOGICAL         CHKREAL
        CHARACTER(2)    CRLF
        INTEGER         ENVINT
        LOGICAL         ENVYN
        INTEGER         FINDC
        INTEGER         GETINVYR
        INTEGER         GETFLINE
        INTEGER         INDEX1
        INTEGER         STR2INT
        REAL            STR2REAL
        REAL*8          STR2DBLE
        REAL            YR2DAY
        INTEGER         FIND1FIRST

        EXTERNAL        CHKINT, CHKREAL, CRLF, ENVINT, ENVYN, FINDC,
     &                  GETINVYR, INDEX1, STR2INT, STR2REAL, STR3DBLE,
     &                  YR2DAY, GETFLINE, BLKORCMT, FIND1FIRST

C...........   SUBROUTINE ARGUMENTS
        INTEGER,      INTENT (IN) :: FDEV         ! unit no. of inv file
        CHARACTER(*), INTENT (IN) :: FNAME        ! logical name of file
        INTEGER,      INTENT (IN) :: NRAWBP       ! no. sources with pollutants
        LOGICAL,      INTENT(OUT) :: NONPOINT     ! true: processing nonpoint inventory

C...........   Local parameters
        INTEGER, PARAMETER :: DATALEN3 = 25      ! length of data field

C...........   Dropped emissions
        INTEGER         NDROP             !  number of records dropped
        REAL            EDROP  ( MXIDAT ) !  total dropped for each pol/activity

C...........   File units and logical/physical names
        INTEGER         TDEV        !  file listed in list formatted input file

C...........   Output from individual reader routines
        CHARACTER(DATALEN3), ALLOCATABLE :: READDATA( :,: )  ! data values
        CHARACTER(IOVLEN3),  ALLOCATABLE :: READPOL ( : )    ! pollutant names

C...........   Other local variables
        INTEGER         I, J, JJ, K, KK, K1, L, NP, SP !  counters and indices
        INTEGER         L1, L2, L3, L4, L5, L6, L7, L8, L9

        INTEGER         CURFIL      !  current file from list formatted inventory
        INTEGER         CURFMT      !  format of current inventory file
        INTEGER         CURSRC      !  current source number
        INTEGER      :: INY = 0     !  tmp inventory year
        INTEGER         IOS         !  i/o status
        INTEGER         INVYEAR     !  inventory year from inventory file
        INTEGER         IREC        !  no. of records read
        INTEGER         ISTREC      !  no. of records stored
        INTEGER         IZONE       !  UTM zone
        INTEGER      :: LSTYR = 0   !  inventory year from list file
        INTEGER         NSCC        !  tmp no of reference SCCs
        INTEGER         MXWARN      !  maximum number of warnings
        INTEGER         NLINE       !  no. of lines in list file
        INTEGER         NPOLPERCAS  !  no. of pollutants per CAS number
        INTEGER      :: NPOLPERLN = 1  !  no. of pollutants per line of inventory file
        INTEGER      :: NWARN = 0   !  current number of warnings
        INTEGER         POLCOD      !  pollutant code
        INTEGER         TMPSTREC    !  temporary stored record counter
        INTEGER         TPF         !  temporal adjustments setting
        INTEGER         SCASPOS     !  position of CAS number in sorted array
        INTEGER         UCASPOS     !  position of CAS number in unique array
        INTEGER         WKSET       !  setting for wkly profile TPFLAG component
        INTEGER      :: FWCOUNT = 0 !  Filling warning count
        INTEGER      :: SAVNVAR = 1 !  Saved number of variables on previous line

        REAL            CEFF        !  tmp control effectiveness
        REAL            EANN        !  annual-ave emission value
        REAL            EMFC        !  emission factor
        REAL            EDAY        !  average day emission value
        REAL            REFF        !  rule effectiveness
        REAL            RPEN        !  rule penetration
        REAL            CPRI        !  primary control code
        REAL            CSEC        !  secondary control code
        REAL            DAY2YR      !  factor to convert from daily data to annual
        REAL            YEAR2DAY    !  factor to convert from annual to daily
        REAL            POLFAC      !  factor for current pollutant
        REAL            POLANN      !  annual emissions for current pollutant
        REAL            RBUF        !  tmp real value
        REAL            REALFL      !  tmp exit flow rate
        REAL            XLOCA1      !  x-dir link coord 1
        REAL            YLOCA1      !  y-dir link coord 1
        REAL            XLOCA2      !  x-dir link coord 2
        REAL            YLOCA2      !  y-dir link coord 2
        REAL            XLOC        !  tmp x coord
        REAL            YLOC        !  tmp y coord
        REAL            FUGSCR

        LOGICAL      :: ACTFLAG = .FALSE. ! true: current pollutant is activity
        LOGICAL      :: CFLAG             ! true: recalc vel w/ flow & diam
        LOGICAL      :: EFLAG   = .FALSE. ! true: error occured
        LOGICAL      :: FFLAG   = .FALSE. ! true: fill annual data with average day data
        LOGICAL      :: HDRFLAG           ! true: current line is part of header
        LOGICAL      :: AVEFLAG = .FALSE. ! true: Aveday inv is processed
        LOGICAL      :: PMCFLG  = .FALSE. ! true: add additional PMC var
        LOGICAL      :: LNKFLAG = .FALSE. ! true: current line has link information
        LOGICAL      :: LSTFLG  = .FALSE. ! true: using list-fmt inventory file
        LOGICAL      :: LSTTIME = .FALSE. ! true: last time through
        LOGICAL      :: NOPOLFLG= .FALSE. ! true: no pollutants stored for this line
        LOGICAL      :: WFLAG   = .FALSE. ! true: all lat-lons to western hemi
        LOGICAL      :: FIRSTIME = .TRUE. ! true: first time through

        CHARACTER(25)       X1        ! x-dir link coord 1
        CHARACTER(25)       Y1        ! y-dir link coord 1
        CHARACTER(25)       X2        ! x-dir link coord 2
        CHARACTER(25)       Y2        ! y-dir link coord 2
        CHARACTER(2)     :: ZONE = ' ' ! UTM zone
        CHARACTER(FIPLEN3)  CFIP      ! fips code
        CHARACTER(RWTLEN3)  CROAD     ! road class no.
        CHARACTER(LNKLEN3)  CLNK      ! link ID

        CHARACTER(NEILEN3) :: NEID = ' ' ! NEI unique ID
        CHARACTER(ORSLEN3) :: CORS = ' ' ! DOE plant ID
        CHARACTER(BLRLEN3) :: BLID = ' ' ! boiler ID
        CHARACTER(EXTLEN3) :: EXTORL = ' ' ! Extended ORL vars
        CHARACTER(40)       DESC      ! plant description
        CHARACTER(ERPLEN3) :: ERPTYP = ' ' ! emissions release point type
        CHARACTER(16)       HT        ! stack height
        CHARACTER(16)       DM        ! stack diameter
        CHARACTER(16)       TK        ! exit temperature
        CHARACTER(16)       FL        ! flow rate
        CHARACTER(16)       VL        ! exit velocity
        CHARACTER(SICLEN3)  SIC       ! SIC
        CHARACTER(SHPLEN3)  SHAPE     ! SHAPE_ID
        CHARACTER(MACLEN3)  MACT      ! MACT code
        CHARACTER(NAILEN3) :: NAICS = ' '  ! NAICS code
        CHARACTER(STPLEN3) :: SRCTYP = ' ' ! source type code
        CHARACTER(SCCLEN3)  TSCC      ! tmp SCC
        CHARACTER           CTYPE     ! coordinate type
        CHARACTER(16)       LAT       ! stack latitude
        CHARACTER(16)       LON       ! stack longitude
        CHARACTER(16)       FUGHT     ! fugitive emissions release height
        CHARACTER(16)       FUGWD     !  " width  (Y dim)
        CHARACTER(16)       FUGLN     !  " length (X dim)
        CHARACTER(16)       FUGAR     !  " area
        CHARACTER(16)       FUGAN     !  " angle

        CHARACTER(IOVLEN3) POLNAM     !  tmp pollutant name
        CHARACTER(300)     INFILE     !  input file line buffer
        CHARACTER(3000)    LINE       !  input file line buffer
        CHARACTER(300)     MESG       !  message buffer

        CHARACTER(16) :: PROGNAME =  'RDINVDATA' ! program name

C***********************************************************************
C   begin body of subroutine RDINVDATA

C.........  Check if inventory file is list format
        IF( ALLOCATED( LSTSTR ) ) THEN
            LSTFLG = .TRUE.

            NLINE = SIZE( LSTSTR )

            DO I = 1, SIZE( LSTSTR )
                IF( LSTSTR( I ) == ' ' ) THEN
                    NLINE = I - 1
                    EXIT
                END IF
            END DO
        ELSE
            NLINE = 1
        END IF

C.........   Initialize variables for keeping track of dropped emissions
        NDROP = 0
        EDROP = 0.  ! array

C.........  Initialize nonpoint flag to false
        NONPOINT = .FALSE.

C.........  Check if any files are ORL nonpoint format
        DO I = 1, NLINE
            IF( FILFMT( I ) == ORLNPFMT ) THEN
                NONPOINT = .TRUE.
                EXIT
            END IF
        END DO

C.........  Set weekly profile interpretation flag...
        WKSET = WTPRFAC
        MESG = 'NOTE: Setting inventory to use full-week '//
     &         'normalizer for weekly profiles'


C.........  Write message
        CALL M3MSG2( MESG )

C.........  Get annual data setting from environment
        MESG = 'Fill in 0. annual data based on average day data.'
        FFLAG = ENVYN( 'FILL_ANNUAL', MESG, .FALSE., IOS )

C.........  Get point specific settings
        IF( CATEGORY == 'POINT' ) THEN
            MESG = 'Flag for recalculating velocity'
            CFLAG = ENVYN( 'VELOC_RECALC', MESG, .FALSE., IOS )
        END IF

        IF( CATEGORY == 'POINT' .OR. CATEGORY == 'MOBILE' ) THEN
            MESG = 'Western hemisphere flag'
            WFLAG = ENVYN( 'WEST_HSPHERE', MESG, .FALSE., IOS )
        END IF

C.........  Get maximum number of warnings
        MXWARN = ENVINT( WARNSET, ' ', 100, IOS )

C.........  Set default inventory characteristics (declared in MODINFO) including NPPOL
        DO I = 1, NLINE
            IF( FILFMT( I ) > 0 ) THEN
                CALL INITINFO( FILFMT( 1 ) )
                EXIT
            END IF
        END DO

C.........  Allocate memory for storing inventory data
        ALLOCATE( INDEXA( NRAWBP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INDEXA', PROGNAME )
        ALLOCATE( POLVLA( NRAWBP,NPPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'POLVLA', PROGNAME )
        ALLOCATE( INRECA( NRAWBP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INRECA', PROGNAME )
        ALLOCATE( IPOSCODA( NRAWBP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IPOSCODA', PROGNAME )
        ALLOCATE( ICASCODA( NRAWBP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICASCODA', PROGNAME )

        ALLOCATE( TPFLAG( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TPFLAG', PROGNAME )
        ALLOCATE( INVYR( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INVYR', PROGNAME )

        ALLOCATE( CSRCTYP( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSRCTYP', PROGNAME )
        ALLOCATE( CINTGR ( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CINTGR', PROGNAME )
        ALLOCATE( CEXTORL ( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CEXTORL', PROGNAME )

        CSRCTYP = ' '       ! array
        CINTGR  = ' '       ! array
        CEXTORL = ' '       ! array

        IF( CATEGORY == 'AREA' ) THEN
            ALLOCATE( CISIC ( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CISIC', PROGNAME )
            ALLOCATE( CSHAPE ( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CSHAPE', PROGNAME )

            CISIC  = REPEAT( '0', SICLEN3 )         ! array
            CSHAPE = ' '
        END IF

        IF( CATEGORY == 'MOBILE' ) THEN
            ALLOCATE( XLOC1( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'XLOC1', PROGNAME )
            ALLOCATE( YLOC1( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'YLOC1', PROGNAME )
            ALLOCATE( XLOC2( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'XLOC2', PROGNAME )
            ALLOCATE( YLOC2( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'YLOC2', PROGNAME )

            XLOC1 = BADVAL3  ! array
            YLOC1 = BADVAL3  ! array
            XLOC2 = BADVAL3  ! array
            YLOC2 = BADVAL3  ! array
        END IF

        IF( CATEGORY == 'POINT' .OR. NONPOINT ) THEN
            ALLOCATE( CMACT( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CMACT', PROGNAME )
            ALLOCATE( CNAICS( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CNAICS', PROGNAME )

            CMACT   = ' '       ! array
            CNAICS  = ' '       ! array
        END IF

        IF( CATEGORY == 'POINT' ) THEN
            ALLOCATE( CISIC ( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CISIC', PROGNAME )
            ALLOCATE( IDIU  ( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'IDIU', PROGNAME )
            ALLOCATE( IWEK  ( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'IWEK', PROGNAME )
            ALLOCATE( STKHT ( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'STKHT', PROGNAME )
            ALLOCATE( STKDM ( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'STKDM', PROGNAME )
            ALLOCATE( STKTK ( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'STKTK', PROGNAME )
            ALLOCATE( STKVE ( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'STKVE', PROGNAME )
            ALLOCATE( XLOCA ( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'XLOCA', PROGNAME )
            ALLOCATE( YLOCA ( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'YLOCA', PROGNAME )
            ALLOCATE( CNEIUID ( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CNEIUID', PROGNAME )
            ALLOCATE( CORIS ( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CORIS', PROGNAME )
            ALLOCATE( CBLRID( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CBLRID', PROGNAME )
            ALLOCATE( CPDESC( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CPDESC', PROGNAME )
            ALLOCATE( CERPTYP( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CERPTYP', PROGNAME )
            ALLOCATE( FUGHGT( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'FUGHGT', PROGNAME )
            ALLOCATE( FUGWID( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'FUGWID', PROGNAME )
            ALLOCATE( FUGLEN( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'FUGLEN', PROGNAME )
            ALLOCATE( FUGANG( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'FUGANG', PROGNAME )

            CISIC    = REPEAT( '0', SICLEN3 )
            IDIU     = 0         ! array
            IWEK     = 0         ! array
            STKHT    = 0.        ! array
            STKDM    = 0.        ! array
            STKTK    = 0.        ! array
            STKVE    = 0.        ! array
            XLOCA    = IMISS3    ! array
            YLOCA    = IMISS3    ! array
            CNEIUID  = ' '       ! array
            CORIS    = ORSBLNK3  ! array
            CBLRID   = BLRBLNK3  ! array
            CPDESC   = ' '       ! array
            CERPTYP  = ' '       ! array
            FUGHGT   = 0.        ! array
            FUGWID   = 0.        ! array
            FUGLEN   = 0.        ! array
            FUGANG   = 0.        ! array
        END IF

C.........  Initialize pollutant-specific values as missing
        POLVLA  = BADVAL3  ! array

C.........  If inventory is list format, open first file for reading
        CURFIL = 1

        IF( LSTFLG ) THEN
            LINE = LSTSTR( CURFIL )

C.............  Skip #LIST lines (must be first)
            IF( INDEX( LINE, 'LIST' ) > 0 ) THEN
                CURFIL = CURFIL + 1
                LINE = LSTSTR( CURFIL )
            END IF

C.............  Check for inventory year packet
            IF( GETINVYR( LINE ) > 0 ) THEN
                CURFIL = CURFIL + 1  ! move to next file in list
            END IF

C.............  Store path of file name
            INFILE = LSTSTR( CURFIL )

C.............  Open current file
            OPEN( FDEV, FILE=INFILE, STATUS='OLD', IOSTAT=IOS )

C.............  Check for errors while opening file
            IF( IOS /= 0 ) THEN

                WRITE( MESG,94010 ) 'Problem at line ', CURFIL, 'of ' //
     &             TRIM( FNAME ) // '.' // ' Could not open file:' //
     &             CRLF() // BLANK5 // TRIM( INFILE )
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            ELSE
                WRITE( MESG,94010 ) 'Successful OPEN for ' //
     &             'inventory file:' // CRLF() // BLANK5 //
     &             TRIM( INFILE )
                CALL M3MSG2( MESG )

            END IF

C.............  Set default inventory characteristics that depend on file format
            CALL INITINFO( FILFMT( CURFIL ) )

C.........  Otherwise, rewind individual file
        ELSE
            REWIND( FDEV )

        END IF

C.........  Allocate memory to store emissions and pollutant from a single line
C.........  For now, set number of pollutants per line to 1
        CURFMT = FILFMT( CURFIL )
        NPOLPERLN = 1
        IF( CURFMT == MEDSFMT ) NPOLPERLN = 6
        ALLOCATE( READDATA( NPOLPERLN,NPPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'READDATA', PROGNAME )
        ALLOCATE( READPOL( NPOLPERLN ), STAT=IOS )
        CALL CHECKMEM( IOS, 'READPOL', PROGNAME )

        CURFMT = FILFMT( CURFIL )

        IREC = 0    ! current record number
        ISTREC = 0  ! current stored record
        SP = 0      ! current source with pollutant index

C.........  Loop through inventory files and read data
        DO

          READ( FDEV, 93000, IOSTAT=IOS ) LINE

          IREC = IREC + 1

          IF( IOS > 0 ) THEN
              EFLAG = .TRUE.
              WRITE( MESG,94010 ) 'I/O error', IOS,
     &           'reading inventory file at line', IREC
              CALL M3MESG( MESG )
              CYCLE
          END IF

C.............  Check if we've reached the end of the file
            IF( IOS < 0 ) THEN

C.................  If list format, try to open next file
                IF( LSTFLG ) THEN

C.....................  Close current file and reset counter
                    CLOSE( FDEV )
                    IREC = 0

C.....................  Advance to next file
                    CURFIL = CURFIL + 1

C.....................  Check if there are more files to read
                    IF( CURFIL <= NLINE ) THEN
                        LINE = LSTSTR( CURFIL )

C.........................  Check for #LIST line
                        IF( INDEX( LINE, 'LIST' ) > 0 ) THEN
                            CURFIL = CURFIL + 1  ! move to next file in list
                            LINE = LSTSTR( CURFIL )
                        END IF

C.........................  Check for INVYEAR packet
                        IF( GETINVYR( LINE ) > 0 ) THEN
                            LSTYR = GETINVYR( LINE )
                            CURFIL = CURFIL + 1  ! more to next file in list
                        END IF

C.........................  Make sure there are still files to read
                        IF( CURFIL > NLINE ) THEN
                            LSTTIME = .TRUE.
                            EXIT
                        END IF

                      INFILE = LSTSTR( CURFIL )

                      OPEN( FDEV, FILE=INFILE, STATUS='OLD',
     &                      IOSTAT=IOS )

C.......................  Check for errors while opening file
                      IF( IOS /= 0 ) THEN
                          WRITE( MESG,94010 ) 'Problem at line ',
     &                       CURFIL, 'of ' // TRIM( FNAME ) //
     &                       '.' // ' Could not open file:' //
     &                       CRLF() // BLANK5 // TRIM( INFILE )
                          CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

                      ELSE
                          WRITE( MESG,94010 )
     &                      'Successful OPEN for ' //
     &                      'inventory file(s):' // CRLF() //
     &                      BLANK5 // TRIM( INFILE )
                          CALL M3MSG2( MESG )

                      END IF

C.......................  Set default inventory characteristics that depend on file format
                      CALL INITINFO( FILFMT( CURFIL ) )
                      CURFMT = FILFMT( CURFIL )

C.......................  Reallocate memory to store emissions from a single line
                      NPOLPERLN = 1
                      IF( CURFMT == MEDSFMT ) NPOLPERLN = 6
                      DEALLOCATE( READDATA, READPOL )
                      ALLOCATE( READDATA( NPOLPERLN,NPPOL ), STAT=IOS )
                      CALL CHECKMEM( IOS, 'READDATA', PROGNAME )
                      ALLOCATE( READPOL( NPOLPERLN ), STAT=IOS )
                      CALL CHECKMEM( IOS, 'READPOL', PROGNAME )
                      SAVNVAR = 1

C.......................  Skip back to the beginning of the loop
                      CYCLE

C...................  Otherwise, no more files to read, so exit
                  ELSE
                      LSTTIME = .TRUE.
                      EXIT
                  END IF

C...............  Otherwise, not a list file, so exit
              ELSE
                  LSTTIME = .TRUE.
                  EXIT
              END IF

          END IF   ! end check for end of file

C...........  Skip blank lines
          IF( LINE == ' ' ) CYCLE
          EXTORL = ' '
          SIC    = ' '
          SHAPE  = ' '

C...........  Process line depending on file format and source category
          SELECT CASE( CURFMT )

            CASE( FF10FMT )
                SELECT CASE( CATEGORY )
                CASE( 'AREA' )
                    CALL RDDATAFF10AR( LINE, READDATA, READPOL, INVYEAR,
     &                                 SRCTYP, SIC, SHAPE, EXTORL,
     &                                 HDRFLAG, AVEFLAG, EFLAG )
                    NPOLPERLN = 1   ! have to set fake value to force reporting

                CASE( 'MOBILE' )
                    CALL RDDATAFF10MB( LINE, READDATA, READPOL, INVYEAR,
     &                        SRCTYP, TSCC, EXTORL, HDRFLAG, AVEFLAG, EFLAG )
                    NPOLPERLN = 1
                    LNKFLAG = .FALSE.

                CASE( 'POINT' )
                    CALL RDDATAFF10PT( LINE, READDATA, READPOL,
     &                                INVYEAR, DESC, ERPTYP, SRCTYP,
     &                                HT, DM, TK, FL, VL, SIC,
     &                                MACT, NAICS, CTYPE, LAT, LON, ZONE,
     &                                NEID, CORS, BLID,
     &                                FUGHT, FUGWD, FUGLN, FUGAN,
     &                                EXTORL, HDRFLAG,
     &                                AVEFLAG, EFLAG )
                    NPOLPERLN = 1
                    IF( READPOL( 1 ) == 'ACRESBURNED' ) FIREFF10 = .TRUE.  ! fire in FF10 format
                END SELECT

            CASE( MEDSFMT )
                SELECT CASE( CATEGORY )
                CASE( 'POINT' )
                    CALL RDDATAMEDSPT( LINE, READDATA, READPOL,
     &                                NPOLPERLN, INVYEAR, CORS, BLID,
     &                                DESC, HT, DM, TK, FL, VL, SIC,
     &                                LAT, LON, HDRFLAG )
                END SELECT

            CASE( ORLFMT )
                SELECT CASE( CATEGORY )
                CASE( 'AREA' )
                    CALL RDDATAORLAR( LINE, READDATA, READPOL, INVYEAR,
     &                                SRCTYP, EXTORL, HDRFLAG, EFLAG )
                    NPOLPERLN = 1   ! have to set fake value to force reporting

                CASE( 'MOBILE' )
                    CALL RDDATAORLMB( LINE, READDATA, READPOL, INVYEAR,
     &                                SRCTYP, EXTORL, HDRFLAG, EFLAG )
                    NPOLPERLN = 1
                    LNKFLAG = .FALSE.

                CASE( 'POINT' )
                    CALL RDDATAORLPT( LINE, READDATA, READPOL,
     &                                INVYEAR, DESC, ERPTYP, SRCTYP,
     &                                HT, DM, TK, FL, VL, SIC, MACT,
     &                                NAICS, CTYPE, LAT, LON, ZONE,
     &                                NEID, CORS, BLID, FUGHT, FUGAR,
     &                                EXTORL, HDRFLAG, EFLAG )
                    NPOLPERLN = 1
                END SELECT

            CASE( ORLNPFMT )
                CALL RDDATAORLNP( LINE, READDATA, READPOL, INVYEAR,
     &                            SIC, MACT, SRCTYP, NAICS, EXTORL,
     &                            HDRFLAG, EFLAG )

            CASE( ORLFIREFMT )
                CALL RDDATAORLFR( LINE, READDATA, READPOL,
     &                            NPOLPERLN, INVYEAR, DESC, SIC, MACT,
     &                            CTYPE, LAT, LON, HDRFLAG, EFLAG)
            END SELECT

C...........  Check for header lines
          IF( HDRFLAG ) THEN

C.................  Reallocate emissions memory with
C                   proper number of pollutants per line
              IF( ( CURFMT == MEDSFMT .OR.
     &              CURFMT == ORLFIREFMT )
     &              .AND. NPOLPERLN .NE. SAVNVAR ) THEN
                  DEALLOCATE( READDATA, READPOL )
                  ALLOCATE( READDATA( NPOLPERLN,NPPOL ), STAT=IOS )
                  CALL CHECKMEM( IOS, 'READDATA', PROGNAME )
                  ALLOCATE( READPOL( NPOLPERLN ), STAT=IOS )
                  CALL CHECKMEM( IOS, 'READPOL', PROGNAME )
                  SAVNVAR = NPOLPERLN
              END IF

C...............  Calculate day to year conversion factor
              IF( INVYEAR /= 0 ) THEN
                  IF( LSTYR > 0 .AND. INVYEAR /= LSTYR ) THEN
                      WRITE( MESG,94010 ) 'NOTE: Using year', LSTYR,
     &                       'from list file, and not year', INVYEAR,
     &                       'from inventory file.'
                      CALL M3MSG2( MESG )

                      INVYEAR = LSTYR
                  END IF

                  YEAR2DAY = YR2DAY( INVYEAR )
                  DAY2YR = 1. / YEAR2DAY
              END IF

              CYCLE
          END IF

C...........  Set inventory year in case there are no header lines
          IF( INVYEAR == 0 .OR.
     &      ( LSTYR > 0 .AND. INVYEAR /= LSTYR ) ) THEN
              INVYEAR = LSTYR

              YEAR2DAY = YR2DAY( INVYEAR )
              DAY2YR = 1. / YEAR2DAY
          END IF

C...........  Make sure some emissions are kept for this source
          IF( NPOLPERLN == 0 ) THEN
              CYCLE
          END IF

C...........  SCC mapping loop : Mobile activity inventory use only.
          KK = 0
          NSCC = 0
          IF( CATEGORY == 'MOBILE' ) THEN
          IF( SCCMAPFLAG ) THEN
              CALL PADZERO( TSCC )
              KK   = INDEX1( TSCC, NSCCMAP, SCCMAPLIST( :,1 ) )
              IF( KK > 0 ) THEN
                  NSCC = STR2INT( SCCMAPLIST( KK,3 ) )
              ELSE
                  IF( EXCLSCCFLAG ) CYCLE    ! drop SCCs not listed in SCCXREF file
              END IF
          END IF
          END IF

C.............  loop over mapped SCC
          DO JJ = 0, NSCC

            IF( JJ > 0 .AND. KK > 0 ) IREC = IREC + 1     ! increment no of records by reference SCCs
            IF( SCCMAPFLAG .AND. KK > 0 ) TSCC = SCCMAPLIST( KK+JJ,2 )

C.............  Check that mobile link info is correct
            IF( CATEGORY == 'MOBILE' .AND. LNKFLAG ) THEN
                IF( X1 == ' ' .OR. Y1 == ' ' .OR.
     &              X2 == ' ' .OR. Y2 == ' ' .OR. ZONE == ' ' ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 'ERROR: Missing link ' //
     &                     'coordinates or UTM zone at line', IREC
                    CALL M3MSG2( MESG )
                END IF

                IF( .NOT. CHKREAL( X1 ) .OR.
     &              .NOT. CHKREAL( Y1 ) .OR.
     &              .NOT. CHKREAL( X2 ) .OR.
     &              .NOT. CHKREAL( Y2 )      ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 'ERROR: Link coordinates ' //
     &                     'are not numbers or have bad formatting' //
     &                     CRLF() // BLANK10 // 'at line', IREC
                    CALL M3MSG2( MESG )
                END IF

                IF( .NOT. CHKINT( ZONE ) ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 'ERROR: UTM zone is not a ' //
     &                     'number or is badly formatted at line', IREC
                    CALL M3MSG2( MESG )
                END IF
            END IF

C.............  Check that point source information is correct
            IF( CATEGORY == 'POINT' .AND.
     &          CURFMT /= ORLFIREFMT ) THEN
                IF( .NOT. CHKREAL( HT ) .OR.
     &              .NOT. CHKREAL( DM ) .OR.
     &              .NOT. CHKREAL( TK ) .OR.
     &              .NOT. CHKREAL( FL ) .OR.
     &              .NOT. CHKREAL( VL )     )THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 'ERROR: Stack parameters ' //
     &                     'are not numbers or have bad formatting' //
     &                     CRLF() // BLANK10 // 'at line', IREC
                    CALL M3MSG2( MESG )
                END IF

                IF ( CURFMT == FF10FMT  ) THEN
                    IF( .NOT. CHKREAL( FUGHT ) .OR.
     &                  .NOT. CHKREAL( FUGWD ) .OR.
     &                  .NOT. CHKREAL( FUGLN ) .OR.
     &                  .NOT. CHKREAL( FUGAN ) ) THEN
                        EFLAG = .TRUE.
                        WRITE( MESG,94010 ) 'ERROR: FUG parameters ' //
     &                     'are not numbers or have bad formatting' //
     &                     CRLF() // BLANK10 // 'at line', IREC
                        CALL M3MSG2( MESG )
                    END IF
                ELSE IF ( CURFMT == ORLFMT  ) THEN
                    IF( .NOT. CHKREAL( FUGHT ) .OR.
     &                  .NOT. CHKREAL( FUGAR ) ) THEN
                        EFLAG = .TRUE.
                        WRITE( MESG,94010 ) 'ERROR: FUG parameters ' //
     &                     'are not numbers or have bad formatting' //
     &                     CRLF() // BLANK10 // 'at line', IREC
                        CALL M3MSG2( MESG )
                    END IF                
                END IF
                

C.................  Check stack height, diameter, and temperature values are
C                   great than zero. reset it to a missing value '-9.0' if not.
                IF( STR2REAL( DM ) == 0.0 ) THEN
                    IF( NWARN < MXWARN ) THEN
                    WRITE( MESG,94010 ) 'WARNING: Stack diameter ' //
     &                 'is equal to 0 at line', IREC ,': reset it to -9'
                    CALL M3MESG( MESG )
                    DM = '  -9.0'        ! reset it to a missing val
                    NWARN = NWARN + 1
                    END IF

                ELSE IF( STR2REAL( TK ) == 0.0 ) THEN
                    IF( NWARN < MXWARN ) THEN
                    WRITE( MESG,94010 ) 'WARNING: Stack temperature ' //
     &                 'is equal to 0 at line', IREC ,': reset it to -9'
                    CALL M3MESG( MESG )
                    TK = '-9.0'           ! reset it to a missing val
                    NWARN = NWARN + 1
                    END IF

                ELSE IF( STR2REAL( VL ) == 0.0 ) THEN
                    IF( NWARN < MXWARN ) THEN
                    WRITE( MESG,94010 ) 'WARNING: Stack velocity ' //
     &                 'is equal to 0 at line', IREC ,': reset it to -9'
                    VL = '     -9.0'      ! reset it to a missing val
                    NWARN = NWARN + 1
                    END IF

                END IF

                IF( .NOT. CHKREAL( LAT ) .OR.
     &              .NOT. CHKREAL( LON )      ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 'ERROR: Latitude and/or ' //
     &                     'longitude are not numbers or have bad ' //
     &                     'formatting' // CRLF() // BLANK10 //
     &                     'at line', IREC
                    CALL M3MSG2( MESG )
                ELSE IF( LAT == ' ' .OR. LON == ' ' ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 'ERROR: Latitude and/or ' //
     &                     'longitude are missing at line', IREC
                    CALL M3MSG2( MESG )
                END IF
            END IF

            IF( ( CATEGORY == 'POINT' .AND. CURFMT /= ORLFIREFMT ) .OR. CURFMT == ORLNPFMT ) THEN

C.................  Check ORL specific values
                IF( CURFMT == ORLFMT ) THEN

                    IF( CTYPE /= 'U' .AND. CTYPE /= 'L' ) THEN
                        EFLAG = .TRUE.
                        WRITE( MESG,94010 ) 'ERROR: Invalid ' //
     &                         'coordinate type at line', IREC,
     &                         CRLF() // BLANK10 // 'Valid ' //
     &                         'values are "U" or "L"'
                        CALL M3MESG( MESG )
                    END IF

                    IF( CTYPE == 'U' .AND.
     &                  ( ZONE == ' ' .OR. ZONE == '0' ) ) THEN
                        EFLAG = .TRUE.
                        WRITE( MESG,94010 ) 'ERROR: Missing ' //
     &                         'or invalid UTM zone at line', IREC
                        CALL M3MESG( MESG )
                    END IF

                END IF

            END IF

C.................  Check ORL FIRE specific values
            IF( CATEGORY == 'POINT' .AND. CURFMT == ORLFIREFMT ) THEN

                IF( .NOT. CHKREAL( LAT ) .OR.
     &              .NOT. CHKREAL( LON )      ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 'ERROR: Latitude and/or ' //
     &                     'longitude are not numbers or have bad ' //
     &                     'formatting' // CRLF() // BLANK10 //
     &                     'at line', IREC
                    CALL M3MSG2( MESG )
                ELSE IF( LAT == ' ' .OR. LON == ' ' ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 'ERROR: Latitude and/or ' //
     &                     'longitude are missing at line', IREC
                    CALL M3MSG2( MESG )
                END IF

                IF( MACT == ' ' ) THEN
                    IF( NWARN < MXWARN ) THEN
                        WRITE( MESG,94010 ) 'WARNING: Missing NFDRS ' //
     &                         'code at line', IREC, '. Default ' //
     &                         'UNKNWN will be used.'
                        CALL M3MESG( MESG )
                        NWARN = NWARN + 1
                    END IF
                    MACT = 'UNKNWN'
                END IF

            END IF

C.............  Get current CAS number position and check that it is valid
            IF( CURFMT == ORLFMT .OR. CURFMT == ORLNPFMT .OR.
     &          CURFMT == FF10FMT .OR. CURFMT == ORLFIREFMT ) THEN
                POLNAM = READPOL( 1 )
                UCASPOS = FINDC( POLNAM, NUNIQCAS, UNIQCAS )
                IF( UCASPOS < 1 ) THEN
                    IF( NWARN < MXWARN ) THEN
                        WRITE( MESG,94010 ) 'WARNING: ' //
     &                      'CAS number ' // TRIM( POLNAM ) //
     &                      ' at line', IREC, ' is not in the ' //
     &                      'inventory pollutants list;' // CRLF() //
     &                      BLANK5 // 'the source will be dropped'
                        CALL M3MESG( MESG )
                        NWARN = NWARN + 1
                    END IF
                    CYCLE

C.................  Check if any part of the CAS number is kept;
C                   could skip rest of loop since no emissions will
C                   by stored, but need values for reporting
                ELSE
                    IF( UCASNKEP( UCASPOS ) == 0 ) THEN
                        NOPOLFLG = .TRUE.
                    ELSE
                        NOPOLFLG = .FALSE.
                    END IF
                END IF

C.................  For non-ORL sources and ORL fires, set UCASPOS to use in
C                   ICASCODA
            ELSE
                UCASPOS = 0
                NOPOLFLG = .FALSE.
            END IF

C.............  Increment number of stored records and double check that we are
C               where we're supposed to be
            IF( .NOT. NOPOLFLG ) THEN
                ISTREC = ISTREC + 1

C.................  Make sure ISTREC isn't more than NSTRECS (shouldn't ever happen)
                IF( ISTREC > NSTRECS ) THEN
                    MESG = 'INTERNAL ERROR: Reached end of records ' //
     &                     'while file still being read'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF

                IF( SRCSBYREC( RECIDX( ISTREC ),1 ) /= CURFIL .OR.
     &              SRCSBYREC( RECIDX( ISTREC ),2 ) /= IREC       ) THEN
                    MESG = 'INTERNAL ERROR: Current record does ' //
     &                 'not match expected record'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF

C.................  Get source ID for current source
                CURSRC = SRCIDA( SRCSBYREC( RECIDX( ISTREC ),3 ) )

            ELSE
                CURSRC = 0

            END IF

C.............  Loop through all pollutants for current line
            DO I = 1, NPOLPERLN

                POLNAM = READPOL( I )
                UCASPOS = FINDC( POLNAM, NUNIQCAS, UNIQCAS )

C.................  Loop through values for this pollutant
C                   Technically, activities will not have NPPOL values,
C                   but we put in blanks when reading the data to avoid problems
C                   Could check if pollutant is an activity, but that would
C                   require an extra search to get the pollutant code

C.................  No need to check these for fires format, since we know
C                   that most values will not be populated intentionally.
                IF( .NOT. FIREFLAG ) THEN
                    DO J = 1, NPPOL
                        IF( .NOT. CHKREAL( READDATA( I,J ) ) ) THEN
                            EFLAG = .TRUE.
                            IF( NWARN < MXWARN ) THEN
                                WRITE( MESG,94010 ) 'ERROR: Emission data, ' //
     &                             'control percentages, and/or emission ' //
     &                             CRLF() // BLANK10 // 'factor for ' //
     &                             TRIM( POLNAM ) // ' are not a number ' //
     &                             'or have bad formatting at line', IREC
                                CALL M3MESG( MESG )
                            END IF
                            EXIT
                        END IF
                    END DO  ! end loop over data values

                    IF( READDATA( I,1 ) == ' ' .AND.
     &                  READDATA( I,2 ) == ' '       ) THEN
                        IF( NWARN < MXWARN ) THEN
                            WRITE( MESG,94010 ) 'WARNING: Missing annual' //
     &                         ' AND average day emissions for ' //
     &                         TRIM( POLNAM ) // ' at line', IREC
                            CALL M3MESG( MESG )
                            NWARN = NWARN + 1
                        END IF
                        READDATA( I,1 ) = '0.'
                        READDATA( I,2 ) = '0.'
                    END IF

                END IF   ! not ORL fires format

C.................  Skip rest of loop if an error has occured
                IF( EFLAG ) CYCLE

C.................  If format is not ORL or ORL Fires, find code corresponding
C                   to current pollutant
                ACTFLAG = .FALSE.
                IF( CURFMT /= ORLFMT .AND. CURFMT /= ORLNPFMT .AND.
     &              CURFMT /= ORLFIREFMT .AND. CURFMT /= FF10FMT ) THEN
                    POLCOD = INDEX1( POLNAM, MXIDAT, INVDNAM )

                    IF( POLCOD == 0 ) THEN
                        WRITE( MESG,94010 ) 'ERROR: Unknown  ' //
     &                      'pollutant ' // TRIM( POLNAM ) //
     &                      ' at line', IREC
                        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                    END IF

                    IF( INVSTAT( POLCOD ) < 0 ) THEN
                        ACTFLAG = .TRUE.
                    END IF

C.................  For ORL fires, #DATA line will have used CAS numbers, so
C                   need to translate these into pollutant names.
                ELSE IF ( CURFMT == ORLFIREFMT ) THEN
                    POLCOD = INDEX1( POLNAM, NINVTBL, ITCASA )
                    POLCOD = INDEX1( ITNAMA(POLCOD), MXIDAT, INVDNAM )
                    IF( POLCOD == 0 ) THEN
                        WRITE( MESG,94010 ) 'ERROR: Unknown  ' //
     &                      'pollutant ' // TRIM( POLNAM ) //
     &                      ' at line', IREC
                        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                    END IF

                END IF

C.................  Convert data to real numbers and check for missing values
                EANN = STR2REAL( READDATA( I,NEM ) )

                IF( EANN < AMISS3 .OR. EANN == -9 ) THEN
                    IF( NWARN < MXWARN .AND. .NOT. FIREFLAG ) THEN
                        IF( ACTFLAG ) THEN
                            WRITE( MESG,94010 ) 'WARNING: Missing ' //
     &                          'inventory data for ' //
     &                          TRIM( POLNAM ) // ' at line', IREC
                        ELSE
                            WRITE( MESG,94010 ) 'WARNING: Missing ' //
     &                          'annual emissions for ' //
     &                          TRIM( POLNAM ) // ' at line', IREC
                        END IF
                        CALL M3MESG( MESG )
                        NWARN = NWARN + 1
                    END IF

                    EANN = BADVAL3
                END IF

                EDAY = STR2REAL( READDATA( I,NDY ) )

                IF( EDAY < AMISS3 .OR. EDAY == -9 ) THEN
                    IF( NWARN < MXWARN .AND. .NOT. FIREFLAG ) THEN
                        IF( ACTFLAG ) THEN
                            WRITE( MESG,94010 ) 'WARNING: Missing ' //
     &                          'average day inventory data for ' //
     &                          TRIM( POLNAM ) // ' at line', IREC
                        ELSE
                            WRITE( MESG,94010 ) 'WARNING: Missing ' //
     &                          'average day emissions for ' //
     &                          TRIM( POLNAM ) // ' at line', IREC
                        END IF
                        CALL M3MESG( MESG )
                        NWARN = NWARN + 1
                    END IF

                    EDAY = BADVAL3
                END IF

C.................  For area and point, convert control percentages
                IF( CATEGORY == 'AREA' .OR. CATEGORY == 'POINT' ) THEN
                    EMFC = STR2REAL( READDATA( I,NEF ) )

C.....................  MRH removed emission factors warning on 6/6/07 because
C                       SMOKE doesn't need emission factor, so why create warning?
                    IF( EMFC < AMISS3 .OR. EMFC == -9 ) THEN
                        EMFC = BADVAL3
                    END IF

                    CEFF = STR2REAL( READDATA( I,NCE ) )

                    IF( CEFF < AMISS3 .OR. CEFF == -9 ) THEN
                        IF( NWARN < MXWARN .AND. .NOT. FIREFLAG ) THEN
                            WRITE( MESG,94010 ) 'WARNING: Missing ' //
     &                         'control efficiency for ' //
     &                         TRIM( POLNAM ) // ' at line', IREC
                            CALL M3MESG( MESG )
                            NWARN = NWARN + 1
                        END IF

                        CEFF = BADVAL3
                    END IF

                    REFF = STR2REAL( READDATA( I,NRE ) )

                    IF( REFF < AMISS3 .OR. REFF == -9 ) THEN
                        IF( NWARN < MXWARN .AND. .NOT. FIREFLAG ) THEN
                            WRITE( MESG,94010 ) 'WARNING: Missing ' //
     &                         'rule effectiveness for ' //
     &                         TRIM( POLNAM ) // ' at line', IREC
                            CALL M3MESG( MESG )
                            NWARN = NWARN + 1
                        END IF

                        REFF = BADVAL3
                    END IF
                END IF

                IF( CATEGORY == 'AREA' ) THEN
                    RPEN = STR2REAL( READDATA( I,NRP ) )

                    IF( RPEN < AMISS3 .OR. RPEN == -9 ) THEN
                        IF( NWARN < MXWARN ) THEN
                            WRITE( MESG,94010 ) 'WARNING: Missing ' //
     &                         'rule penetration for ' //
     &                         TRIM( POLNAM ) // ' at line', IREC
                            CALL M3MESG( MESG )
                            NWARN = NWARN + 1
                        END IF

                        RPEN = BADVAL3
                    END IF
                END IF

                IF( CATEGORY == 'POINT' ) THEN
                    CPRI = STR2REAL( READDATA( I,NC1 ) )

C.....................  MRH removed primary control code warning on 6/6/07 because
C                       SMOKE doesn't need code, so why create warning?
                    IF( CPRI < AMISS3 .OR. CPRI == -9 ) THEN
                        CPRI = BADVAL3
                    END IF

                    CSEC = STR2INT( READDATA( I,NC2 ) )

C.....................  MRH removed secondary control code warning on 6/6/07 because
C                       SMOKE doesn't need code, so why create warning?
                    IF( CSEC < AMISS3 .OR. CSEC == -9 ) THEN
                        CSEC = BADVAL3
                    END IF
                END IF

C.................  Set the default temporal resolution of the data
                TPF = MTPRFAC * WKSET

C.................  Replace annual data with average day information
C                   Only do this if current pollutant is not an activity, flag is set,
C                   annual data is less than or equal to zero and average day
C                   data is greater than zero
                IF( .NOT. ACTFLAG .AND.
     &                 EANN <= 0. .AND.
     &                 EDAY >  0.       ) THEN

C.....................  Fill annual inventory using average day inventory
                    IF( FFLAG ) THEN
                        EANN = EDAY * DAY2YR  ! fill annual inv (aveday*365)

                        IF( FWCOUNT < MXWARN .AND. CURFMT /= FF10FMT ) THEN
                           WRITE(MESG,94010) 'WARNING: Using average day '//
     &                       'emissions to fill in annual emissions' //
     &                       CRLF()// BLANK10// 'for ' //TRIM( POLNAM ) //
     &                       ' at line', IREC
                           CALL M3MESG( MESG )
                           FWCOUNT = FWCOUNT + 1
                        END IF
                    END IF

C.....................  Remove monthly factors for this source
                    TPF = WKSET

                END IF

C.................  Set if averday inv is processed in FF10 fomrat (no-activity inv)
                IF( .NOT. ACTFLAG .AND. AVEFLAG ) TPF = WKSET

C.................  Calculate average day emissions from annual data if needed
C                   Only do this if current pollutant is not an activity,
C                   annual data is greater than zero, and average day data is
C                   zero or negative
                IF( .NOT. ACTFLAG .AND.
     &                 EANN >=  0. .AND.
     &                 EDAY <  0.       ) THEN
                     EDAY = EANN * YEAR2DAY
                END IF

C.................  If current format is ORL, check if current CAS number
C                   is split
                IF( CURFMT == ORLFMT .OR. CURFMT == ORLNPFMT .OR.
     &              CURFMT == FF10FMT .OR. CURFMT == ORLFIREFMT ) THEN
                    NPOLPERCAS = UCASNPOL( UCASPOS )

C.....................  Store emissions by CAS number for reporting
                    IF( EANN > 0. ) THEN
                      EMISBYCAS( UCASPOS ) = EMISBYCAS( UCASPOS ) + EANN
                    END IF

                    RECSBYCAS( UCASPOS ) = RECSBYCAS( UCASPOS ) + 1

                ELSE
                    NPOLPERCAS = 1
                    POLFAC = 1.
                END IF

                DO J = 0, NPOLPERCAS - 1

C.....................  If ORL format, set current pollutant
                    IF( CURFMT == ORLFMT .OR. CURFMT == ORLNPFMT .OR.
     &                  CURFMT == FF10FMT .OR. CURFMT == ORLFIREFMT ) THEN

                        SCASPOS = UCASIDX( UCASPOS ) + J

C.........................  Set factor for this CAS number
                        POLFAC = ITFACA( SCASIDX( SCASPOS ))

C.........................  Multiply annual emissions by factor
                        IF( EANN > 0. ) THEN
                            POLANN = EANN * POLFAC
                        ELSE
                            POLANN = EANN
                        END IF

C.........................  Store emissions by pollutant for reporting
                        IF( POLANN > 0. ) THEN
                            EMISBYPOL( SCASPOS ) =
     &                          EMISBYPOL( SCASPOS ) + POLANN
                        END IF

C.........................  Make sure current pollutant is kept
                        IF( ITKEEPA( SCASIDX( SCASPOS ) ) ) THEN
                            POLNAM = ITNAMA( SCASIDX( SCASPOS ) )
                        ELSE
                            CYCLE
                        END IF

C.........................  Find code corresponding to current pollutant
                        POLCOD = INDEX1( POLNAM, MXIDAT, INVDNAM )
                        IF( POLCOD == 0 ) THEN
                            WRITE( MESG,94010 ) 'ERROR: Unknown  ' //
     &                          'pollutant ' // TRIM( POLNAM ) //
     &                          ' at line', IREC
                            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                        END IF

                    ELSE
                        POLANN = EANN
                    END IF

C.....................  Store data in unsorted order
                    SP = SP + 1

                    IF( SP <= NRAWBP ) THEN

                        INRECA  ( SP ) = CURSRC    ! map to sorted source number
                        INDEXA  ( SP ) = SP        ! index for sorting POLVLA
                        IPOSCODA( SP ) = POLCOD    ! pollutant code
                        ICASCODA( SP ) = UCASPOS   ! CAS number (set to 0 for non-ORL sources)

                        IF( POLANN > 0. ) THEN
                          POLVLA( SP,NEM ) = INVDCNV( POLCOD ) *
     &                                       POLANN
                        ELSE
                          POLVLA( SP,NEM ) = POLANN
                        END IF

                        IF( EDAY > 0. ) THEN
                            POLVLA( SP,NDY ) = EDAY * POLFAC
                        ELSE
                            POLVLA( SP,NDY ) = EDAY
                        END IF

                        IF( FIREFLAG ) POLVLA( SP,NDY )=BADVAL3

                        IF( CATEGORY == 'AREA' .OR.
     &                      CATEGORY == 'POINT'    ) THEN
                            POLVLA( SP,NEF ) = EMFC
                            POLVLA( SP,NCE ) = CEFF
                            POLVLA( SP,NRE ) = REFF
                        END IF

                        IF( CATEGORY == 'AREA' ) THEN
                            POLVLA( SP,NRP ) = RPEN
                        END IF

                        IF( CATEGORY == 'POINT' ) THEN
                            POLVLA( SP,NC1 ) = CPRI
                            POLVLA( SP,NC2 ) = CSEC
                        END IF
                    END IF

                END DO  ! end loop through pols per CAS number

            END DO  ! end loop through pols per line

C.............  Skip rest of loop if no pollutants are kept
            IF( NOPOLFLG ) CYCLE

            INVYR ( CURSRC ) = INVYEAR
            TPFLAG( CURSRC ) = TPF

C.............  Store SIC, SHAPE if not blank
            IF ( ASSOCIATED( CISIC ) .AND. SIC /= ' ' ) THEN
                CALL PADZERO( SIC )
                CISIC( CURSRC ) = SIC
            END IF
            IF ( ASSOCIATED( CSHAPE ) .AND. SHAPE /= ' ' ) THEN
                CSHAPE( CURSRC ) = SHAPE
            END IF

            IF( CATEGORY == 'POINT' .OR.
     &          CURFMT == ORLNPFMT .OR. CURFMT == ORLFMT .OR.
     &          CURFMT == FF10FMT  .OR. CURFMT /= MEDSFMT ) THEN

                 IF( ( CATEGORY == 'POINT' .AND. CURFMT == ORLFMT )
     &               .OR. (CATEGORY == 'POINT' .AND. CURFMT == FF10FMT)
     &               .OR. CURFMT == ORLNPFMT               ) THEN
                     IF( MACT == '-9' ) MACT = ' '
                     IF( NAICS == '-9' ) NAICS = ' '
                     CALL PADZERO( MACT )
                     CALL PADZERO( NAICS )
                     CMACT  ( CURSRC ) = MACT
                     CNAICS ( CURSRC ) = NAICS
                     CEXTORL( CURSRC ) = ADJUSTL( EXTORL )
                 END IF

                 IF( CURFMT /= MEDSFMT ) THEN
                     CALL PADZERO( SRCTYP )
                     CSRCTYP( CURSRC ) = SRCTYP
                     CEXTORL( CURSRC ) = ADJUSTL( EXTORL )
                 END IF

                IF( CATEGORY == 'POINT' .AND. CURFMT /= MEDSFMT ) THEN
                    CALL PADZERO( ERPTYP )
                    CERPTYP( CURSRC ) = ERPTYP

C.....................  Convert UTM values to lat-lon
                    IF( CTYPE == 'U' ) THEN
                        IZONE = STR2INT( ZONE )
                        XLOCA1 = STR2REAL( LON )
                        YLOCA1 = STR2REAL( LAT )
                        CALL UTM2LL( XLOCA1, YLOCA1, IZONE, XLOC, YLOC )
                        WRITE( LON,'(F9.5)' ) XLOC
                        WRITE( LAT,'(F9.5)' ) YLOC
                    END IF
                END IF
            END IF

            IF( CATEGORY == 'POINT' ) THEN
                IF( CURFMT /= ORLFIREFMT ) THEN
                    STKHT   ( CURSRC ) = STR2REAL( HT )
                    STKDM   ( CURSRC ) = STR2REAL( DM )
                    STKTK   ( CURSRC ) = STR2REAL( TK )
                    STKVE   ( CURSRC ) = STR2REAL( VL )
                    XLOCA   ( CURSRC ) = STR2DBLE( LON )
                    YLOCA   ( CURSRC ) = STR2DBLE( LAT )
                    CPDESC  ( CURSRC ) = DESC
                    CNEIUID ( CURSRC ) = ADJUSTR( NEID )
                    CORIS   ( CURSRC ) = ADJUSTR( CORS )
                    CBLRID  ( CURSRC ) = ADJUSTR( BLID )
                    CEXTORL ( CURSRC ) = ADJUSTL( EXTORL )

C.....................  Convert units on values
                    IF( STKHT( CURSRC ) < 0. ) THEN
                        STKHT( CURSRC ) = 0.
                        IF( NWARN < MXWARN ) THEN
                        WRITE( MESG,94010 ) 'WARNING: Negative stack ' //
     &                      'height :: Reset to zero at line', IREC
                        CALL M3MESG( MESG )
                        NWARN = NWARN +1
                        END IF
                    END IF
                    STKHT( CURSRC ) = STKHT( CURSRC ) * FT2M   ! ft to m

                    IF( STKDM( CURSRC ) < 0. ) THEN
                        STKDM( CURSRC ) = 0.
                        IF( NWARN < MXWARN ) THEN
                        WRITE( MESG,94010 ) 'WARNING: Negative stack ' //
     &                      'diameter :: Reset to zero at line', IREC
                        CALL M3MESG( MESG )
                        NWARN = NWARN +1
                        END IF
                    END IF
                    STKDM( CURSRC ) = STKDM( CURSRC ) * FT2M   ! ft to m

                    IF( STKVE( CURSRC ) < 0. ) THEN
                        STKVE( CURSRC ) = 0.
                        IF( NWARN < MXWARN ) THEN
                        WRITE( MESG,94010 ) 'WARNING: Negative stack ' //
     &                      'exit velocity :: Reset to zero at line', IREC
                        CALL M3MESG( MESG )
                        NWARN = NWARN +1
                        END IF
                    END IF
                    STKVE( CURSRC ) = STKVE( CURSRC ) * FT2M  ! ft/s to m/s

                    IF( STKTK( CURSRC ) < 0. ) THEN
                        STKTK( CURSRC ) = 0.
                        IF( NWARN < MXWARN ) THEN
                        WRITE( MESG,94010 ) 'WARNING: Negative stack ' //
     &                      'exit temperature :: Reset to zero at line', IREC
                        CALL M3MESG( MESG )
                        NWARN = NWARN +1
                        END IF
                    ELSE
                        STKTK( CURSRC ) = ( STKTK( CURSRC ) - 32.0 ) *   ! F to K
     &                                    FTOC + CTOK
                    END IF

                END IF

C.................  Recompute velocity if that option has been set
                IF( CFLAG ) THEN
                    RBUF = 0.25 * PI *
     &                     STKDM( CURSRC ) * STKDM( CURSRC )

                    REALFL = STR2REAL( FL )
                    IF( REALFL < 0. ) THEN
                        REALFL = 0.
                        WRITE( MESG,94010 ) 'WARNING: Negative stack flowrate ::'//
     &                      ' Computed stack exit velocity will be zero at line', IREC
                        CALL M3MESG( MESG )
                    END IF
                    REALFL = REALFL * FT2M3                 ! ft^3/s to m^3/s

                    IF( RBUF > 0 ) THEN
                        STKVE( CURSRC ) = REALFL / RBUF
                    ELSE
                        STKVE( CURSRC ) = 0.0
                        WRITE( MESG,94010 ) 'WARNING: Failed to compute stack ' //
     &                      'exit velocity due to zero diameter :: '//
     &                      'Reset to zero at line', IREC
                        CALL M3MESG( MESG )
                    END IF

                END IF

                IF ( CURFMT == FF10FMT  ) THEN
                    FUGHGT( CURSRC ) = MAX( 0.0, FT2M * STR2REAL( FUGHT ) )
                    FUGWID( CURSRC ) = MAX( 0.0, FT2M * STR2REAL( FUGWD ) )
                    FUGLEN( CURSRC ) = MAX( 0.0, FT2M * STR2REAL( FUGLN ) )
                    FUGANG( CURSRC ) = MAX( 0.0, STR2REAL( FUGAN ) )

                ELSE IF ( CURFMT == ORLFMT  ) THEN
                    FUGSCR = FT2M * SQRT( MAX( 0.0, STR2REAL( FUGAR ) ) )
                    FUGHGT( CURSRC ) = MAX( 0.0, FT2M * STR2REAL( FUGHT ) )
                    FUGWID( CURSRC ) = FUGSCR
                    FUGLEN( CURSRC ) = FUGSCR
                    FUGANG( CURSRC ) = 0.0
                END IF

C.................  Correct hemisphere for stack longitude
                IF( WFLAG .AND. XLOCA( CURSRC ) > 0. ) THEN
                    XLOCA( CURSRC ) = -XLOCA( CURSRC )
                END IF

            END IF

            IF( FIREFF10 ) THEN     ! reset stack parameters when processing fire in FF10 format
                STKHT  ( CURSRC ) = BADVAL3
                STKDM  ( CURSRC ) = BADVAL3
                STKTK  ( CURSRC ) = BADVAL3
                STKVE  ( CURSRC ) = BADVAL3
            END IF

            IF( CURFMT == ORLFIREFMT )THEN
                CALL PADZERO( SRCTYP )
                CALL PADZERO( ERPTYP )
                CALL PADZERO( MACT )
                CALL PADZERO( NAICS )

                CSRCTYP( CURSRC ) = SRCTYP
                CERPTYP( CURSRC ) = ERPTYP
                CMACT  ( CURSRC ) = MACT
                CNAICS ( CURSRC ) = NAICS
                XLOCA  ( CURSRC ) = STR2DBLE( LON )
                YLOCA  ( CURSRC ) = STR2DBLE( LAT )
                STKHT  ( CURSRC ) = BADVAL3
                STKDM  ( CURSRC ) = BADVAL3
                STKTK  ( CURSRC ) = BADVAL3
                STKVE  ( CURSRC ) = BADVAL3
                CPDESC ( CURSRC ) = DESC
                CORIS  ( CURSRC ) = ADJUSTR( CORS )
                CBLRID ( CURSRC ) = ADJUSTR( BLID )

C.................  Correct hemisphere for stack longitude
                IF( WFLAG .AND. XLOCA( CURSRC ) > 0. ) THEN
                    XLOCA( CURSRC ) = -XLOCA( CURSRC )
                END IF

            END IF

          END DO  ! end loop through reference SCC if applicable

        END DO  ! end loop through records array

C.........  Deallocate local memory, if its allocated
        IF( ALLOCATED( READDATA ) ) DEALLOCATE( READDATA )
        IF( ALLOCATED( READPOL  ) ) DEALLOCATE( READPOL )

C.........  Abort if there was an error
        IF( EFLAG ) THEN
            MESG = 'Error reading data from inventory file ' //
     &              TRIM( FNAME )
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Sort inventory data by source
        CALL M3MESG( 'Sorting inventory data by source ' //
     &               'and pollutant...' )

        CALL SORTI2( SP, INDEXA, INRECA, IPOSCODA )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

93000   FORMAT( A )

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94060   FORMAT( 10( A, :, E10.3, :, 1X ) )

        END SUBROUTINE RDINVDATA
