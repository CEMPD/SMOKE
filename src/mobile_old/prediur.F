 
C Version "@(#)$Id$ $Source$ $Date$ 
 
        PROGRAM PREDIUR
 
C***********************************************************************
C  program body starts at line
C
C  DESCRIPTION:
C     This program inputs gridded, time-dependent temperature data, a mobile 
C     list file, and an ungridding matrix and
C     determines the min/max temperature combinations for each mobile source 
C     and for each emission factor set.
C
C  PRECONDITIONS REQUIRED:
C     PREDIUR handles modeling period in multiple of a whole day. Input for
C     "Output duration" has to be a multiple of 24.
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Prototype 11/95 by MRH
C
C***********************************************************************

        IMPLICIT NONE
 
C...........   INCLUDES:

        INCLUDE 'GRDIMS3.EXT'    !  Gridding include parameters
        INCLUDE 'MBDIMS3.EXT'    !  Mobile include paramaters
        INCLUDE 'PARMS3.EXT'     !  I/O API parameters
        INCLUDE 'IODECL3.EXT'    !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'     !  I/O API file description data structures.
        INCLUDE 'CONST3.EXT'     !  physical and mathematical constants

C...........   EXTERNAL FUNCTIONS and their descriptions:

        INTEGER      FIND1
        INTEGER      FIND3
        INTEGER      GETNUM
        REAL         GETREAL
        LOGICAL      GETYN
        INTEGER      LOCAT1
        INTEGER      PROMPTFFILE
        CHARACTER*16 PROMPTMFILE
        INTEGER      RDMPLIST
        INTEGER      TRIMLEN

        EXTERNAL     FIND1, FIND3, GETNUM, GETREAL, GETYN, LOCAT1,
     &               PROMPTFFILE, PROMPTMFILE, RDMPLIST, TRIMLEN

C...........   LOCAL PARAMETERS and their descriptions:

        REAL            CTOF   !  Celcius to Fahrnht conversion constant
        REAL            FTOC   !  Fahrnht to Celcius conversion constant
        INTEGER         MXTFIL !  maximum number of input temperature files

        PARAMETER     (   CTOF   = 9.0 / 5.0 
     &                  , FTOC   = 1.0 / CTOF
     &                  , MXTFIL = 20 
     &                )

C...........   LOCAL VARIABLES and their descriptions:

C...........   Minimum and maximum source temperature variables
        INTEGER   TMMI ( NMSRC )   !  min/max combo index
        REAL      TKMAX( NMSRC )   !  source-maximum daily temperature
        REAL      TKMIN( NMSRC )   !  source-minimum daily temperature

        COMMON / MINMAX / TMMI, TKMIN, TKMAX

C...........   Temperature index per emission factor PSI variables
        INTEGER   NPSI   (        NMSRC )  ! Intrmdt number of PSIs per source
        INTEGER   PSI    ( MXPPS, NMSRC )  ! Intrmdt storage of PSIs per source

        INTEGER   NTIDX  (        MXPSI )  !
        INTEGER   TIDXPSI( MXTMMI, MXPSI )  !

C...........   Mobile Sources input and output arrays
 
        INTEGER   IFIP  ( NMSRC )         !  FIP codes for sources
        INTEGER   IRCLAS( NMSRC )         !  EPA road classification code
        INTEGER   ILINK ( NMSRC )         !  link number or 0 for non-link

C...........   MPLIST scenario-specific emissions factor XREF tables:
 
        INTEGER   NMPLIST                 !  number of XREF entries
        INTEGER   MPFIP( MXPLST )         !  FIP codes in XREF entries
        INTEGER   MPRDT( MXPLST )         !  road class codes in XREF entries
        INTEGER   MPLNK( MXPLST )         !  link IDs in XREF entries
        INTEGER   MPREF( 0:23, MXPLST )   !  hour&source-specific xrefs
        INTEGER   REFCNT                  !  count of total unique PSIs
        INTEGER   REFLIST( MXPSI )        !  list of unique PSIs

C...........   Gridded Temperature
        REAL      TA( NGRID )             !  one-layer slice of temperature

C...........   Ungridding matrix:
 
        INTEGER   NU( NMSRC )
        INTEGER   IU( NMATX )
        REAL      CU( NMATX )
 
        COMMON  / UGRIDMAT / NU, IU, CU

C...........   File units and logical names:

        INTEGER         EDEV   !  unit number for output TMMI per PSI
        INTEGER         LDEV   !  unit number for IO/API
        INTEGER         NLDEV  !  unit number for emissions factors xref file
        CHARACTER*16    ANAME  !  logical name for mobile-source input file
        CHARACTER*16    ONAME  !  logical name for output ungridded min/max temp
        CHARACTER*16    UNAME  !  logical name for ungridding-matrix input file

        CHARACTER*16    TNAME( MXTFIL ) !  logical names for surface temp files

C...........   Other local variables:

        INTEGER    F, H, I, J, JJ, K, L, S, T  ! Counters and pointers
        INTEGER    FIP     !  temporary FIP code
        INTEGER    IFIL    !  counter for temperature files
        INTEGER    IZONE   !  zone of time stamps in input temperature files
        INTEGER    JDATE   !  output date counter (YYYYDDD) in GMT
        INTEGER    JTIME   !  output time counter (HHMMSS)  in GMT
        INTEGER    KU      !  counter of grid cells per source
        INTEGER    LNK     !  temporary FIP code
        INTEGER    MAXTL   !  intermediate temperature index value
        INTEGER    MDATE   !  temperature time zone date counter (YYYYDDD)
        INTEGER    MTIME   !  temperature time zone time counter (HHMMSS)
        INTEGER    NSTEPS  !  number of time steps to process temperature data
        INTEGER    NTFIL   !  total number of input temperature files
        INTEGER    OSRC    !  number of sources outside grid
        INTEGER    RDT     !  temporary FIP code
        INTEGER    TMMIT(4)!  temporary TMMI indices for all four min/max
        INTEGER    TSTEP   !  time step of input temperature data (HHMMSS)
        INTEGER    TZONE   !  zone to determine output days
        INTEGER    FIRSTIME!  the starting time for modeling

        REAL       TMAX   !  deg F temporary source max temperature on interval
        REAL       TMIN   !  deg F temporary source min temperature on interval
        REAL       MSAV   !  temporary saved minimum temperature
        REAL       TASRC  !  temporary temperature, interpolated to sources
        REAL       TF     !  deg F temporary temperature storage

        LOGICAL         EFLAG   !  Error flag 
        LOGICAL         OFLAG   !  sources outside grid

        CHARACTER*16    BUFFER  !  Temporary file name buffer 
        CHARACTER*16    TKDAT   !  logical file name buffer for temperature 
        CHARACTER*256   MESG    !  message buffer

        
C***********************************************************************
C   begin body of program PREDIUR

        LDEV = INIT3()
 
        WRITE( *,92000 )
     &  ' ',
     &  'This program inputs gridded, time-dependent temperature',
     &  'data, a mobile list file, and an ungridding matrix and',
     &  'determines the min/max temperature combinations for each',
     &  'mobile source and for each emission factor set.',
     &  ' ',
     &  'You will need to enter the logical names for the input and',
     &  'output files (and to have set them prior to program launch,',
     &  'using "setenv <logicalname> <pathname>").',
     &  ' ',
     &  'You may use END_OF-FILE (control-D) to quit the program',
     &  'during logical-name entry. Default responses are given in',
     &  'brackets [LIKE THIS] and can be accepted by hitting the',
     &  '<RETURN> key.',
     &  ' '

        IF ( .NOT. GETYN( 'Continue with program?', .TRUE. ) ) THEN
            WRITE( *,92000 ) ' ', 'Ending program PREDIUR', ' '
            CALL EXIT( 0 )
        END IF

C.........  Prompt for SMOKE mobile inventory file
        ANAME = PROMPTMFILE(
     &          'Enter logical name for VMT INVENTORY input file',
     &          FSREAD3, 'MOBL', 'PREDIUR' )

C.........  Prompt for Ungridding matrix
        UNAME = PROMPTMFILE(
     &          'Enter logical name for UNGRIDDING MATRIX input file',
     &          FSREAD3, 'MUGMAT', 'PREDIUR' )

C.........  Prompt for MPLIST (scenario specific parameter scheme list)
        NLDEV = PROMPTFFILE(
     &          'Enter logical name for PARAMETER LIST input file',
     &         .TRUE., .TRUE., 'MPLIST', 'PREDIUR' )

C.........  Read FIP, road class, and Link IDs from inventory file

        WRITE( *,92000 ) ' ', 'Reading VMT INVENTORY file ...'

        IF ( .NOT. READ3( ANAME, 'IFIP', ALLAYS3, 0,0, IFIP ) ) THEN
            CALL M3EXIT( 'PREDIUR', 0, 0,
     &                   'Error reading IFIP from file' // ANAME, 2)
        END IF
 
        IF ( .NOT. READ3( ANAME, 'IRCLAS', ALLAYS3, 0,0, IRCLAS ) ) THEN
            CALL M3EXIT( 'PREDIUR', 0, 0,
     &                   'Error reading IRCLAS from file' // ANAME, 2 )
        END IF
 
        IF ( .NOT. READ3( ANAME, 'ILINK', ALLAYS3, 0,0, ILINK ) ) THEN
            CALL M3EXIT( 'PREDIUR', 0, 0,
     &                   'Error reading ILINK from file' // ANAME, 2 )
        END IF

C.........  Read Ungridding Matrix

        WRITE( *,92000 ) ' ', 'Reading UNGRIDDING MATRIX ...'

        IF ( .NOT. READ3( UNAME, ALLVAR3, ALLAYS3, 0, 0, NU ) ) THEN
            MESG = 'Could not read ungridding matrix from ' // UNAME
            CALL M3EXIT( 'PREDIUR', 0, 0, MESG, 2 )
        END IF

C.........  Prompt for time zone of temperature data.
        IZONE  = 0
        IZONE  = GETNUM( -12, 12, IZONE,
     &             'Enter time zone of temperature files '//
     &             '(0 for GMT, 5 for EST, etc.)' )

C.........  Loop for prompting for temperature data files:
        WRITE( *,92000 )
     &  ' ',
     &  'Now enter the list of logical file names for ',
     &  'surface temperature input files',
     &  ' '

        IFIL =    0
   66   CONTINUE

            IFIL = IFIL + 1
            WRITE( TKDAT, 94040 ) 'TKDAT', IFIL

C.............  Prompt for next temperature file (Met data file)
            BUFFER = PROMPTMFILE(
     &               'Enter the NEXT FILE, or "NONE" to end list',
     &               FSREAD3, TKDAT, 'PREDIUR' )

            IF ( BUFFER .EQ. 'NONE' ) THEN
                GO TO 70

            ELSEIF( IFIL .GT. MXTFIL ) THEN

                WRITE( MESG, 94010 ) 
     &                 'Maximum number of temperature files' //
     &                 'exceeded: MXTFIL = ', MXTFIL
                CALL M3EXIT( 'PREDIUR', 0, 0, MESG, 2 )

            ELSE
                TNAME( IFIL ) = BUFFER

            ENDIF

        GO TO 66

   70   CONTINUE  ! End of loop for temperature file names
        NTFIL = IFIL - 1

C.........  Read header of first temperature file
        IF ( .NOT. DESC3( TNAME( 1 ) ) ) THEN
            CALL M3EXIT( 'PREDIUR', 0, 0,
     &                  'Could not get description of file "'
     &                  // TNAME( 1 )( 1:TRIMLEN( TNAME( 1 ) ) )
     &                  // '"', 2 )
        END IF

C.........  Save time step for processing temperatures
        TSTEP = TSTEP3D

C.........  Prompt for time zone for output file (to get local time to use
C.........     for beginning of the day)
        TZONE  = 0  ! Keep data in MINMAXT in GMT for increased portability
C        TZONE  = GETNUM( -12, 12, TZONE,
C     &             'Enter time zone for output min/max '//
C     &             'per day (0 for GMT)' )
        JDATE  = GETNUM( 0, 9999999, SDATE3D,
     &                   'Enter starting date in GMT (YYYYDDD)' )
        JTIME  = GETNUM( 0, 999999, STIME3D,
     &                   'Enter starting time in GMT (HHMMSS)' )
        FIRSTIME = JTIME
75      NSTEPS = GETNUM( 1, 999999, 24,
     &                   'Enter output duration (hours)' )
        IF ( MOD( NSTEPS, 24 ) .NE. 0 ) THEN
           WRITE( *,92000 ) 'Output duration has to be a multiple of 24'
           GO TO 75
        END IF 
 
C.........  Prompt for temperature increment to use in output file
cmh     TINTV = GETREAL( 0.0001, 50., 2.,
cmh  &                   'Enter temperature increment [deg F]' )

C.........  Set date/time for temperature based on start date/time
        MDATE = JDATE
        MTIME = JTIME
        CALL NEXTIME( MDATE, MTIME, (TZONE - IZONE)*10000 )

C.........  Abort if start time is < than first time for temperatures
        IF( ( MDATE .EQ. SDATE3D .AND. MTIME .LT. STIME3D ) .OR.
     &        MDATE .LT. SDATE3D ) THEN

            WRITE( *,92010 )
     &      ' ',
     &      'Requested start date/time:', MDATE, MTIME, IZONE,
     &      'First temperature data data/time:',SDATE3D,STIME3D,IZONE,
     &      ' '
            CALL M3EXIT( 'PREDIUR', 0, 0,
     &                  'Requested start is earlier than start '
     &                  // 'of first temperature file', 2 )

        ENDIF

C.........  Initialize source-based output file NetCDF parameters
        FTYPE3D = GRDDED3
        P_ALP3D = DBLE( BADVAL3 )
        P_BET3D = DBLE( BADVAL3 )
        P_GAM3D = DBLE( BADVAL3 )
        P_GAM3D = DBLE( BADVAL3 )
        XCENT3D = 0.0D0
        YCENT3D = 0.0D0
        XORIG3D = DBLE( BADVAL3 )
        YORIG3D = DBLE( BADVAL3 )
        SDATE3D = JDATE
        STIME3D = JTIME
        TSTEP3D = 240000
        NVARS3D = 3
        NCOLS3D = 1
        NROWS3D = NMSRC     !  number of rows = # of mobile sources.
        NLAYS3D = 1
        NTHIK3D = 1
        GDTYP3D = UTMGRD3
        VGTYP3D = IMISS3
        VGTOP3D = BADVAL3
 
        J = 1
        VNAME3D( J ) = 'TMMI'
        UNITS3D( J ) = 'n/a'
        VDESC3D( J ) = 'Temperature combination index'
        VTYPE3D( J ) = M3INT
        J = J + 1

        VNAME3D( J ) = 'TKMIN' 
        UNITS3D( J ) = 'deg K' 
        VDESC3D( J ) = 'Minimum source temperature for day' 
        VTYPE3D( J ) = M3REAL
        J = J + 1

        VNAME3D( J ) = 'TKMAX'
        UNITS3D( J ) = 'deg K'
        VDESC3D( J ) = 'Maximum source temperature for day'
        VTYPE3D( J ) = M3REAL

        FDESC3D( 1 ) = 'Minimum and maximum temperatures per source'
        FDESC3D( 2 ) = 'per day.  Also has an index which is the'
        FDESC3D( 3 ) = 'row number for the temperature combination'
        FDESC3D( 4 ) = 'in a sorted list of all min/max temperature'
        FDESC3D( 5 ) = 'combinations (including illegal ones), with'
        FDESC3D( 6 ) = 'max cycling fastest, given:'
        WRITE( FDESC3D( 7  ), 94030 ) 'MINT_MIN = ', MINT_MIN, 'deg F'
        WRITE( FDESC3D( 8  ), 94030 ) 'MINT_MAX = ', MINT_MAX, 'deg F'
        WRITE( FDESC3D( 9  ), 94030 ) 'MAXT_MIN = ', MAXT_MIN, 'deg F'
        WRITE( FDESC3D( 10 ), 94030 ) 'MAXT_MAX = ', MAXT_MAX, 'deg F'
        WRITE( FDESC3D( 11 ), 94030 ) 'INCREMENT= ', TINTV, 'deg F'

C.........  Prompt for Source-based Output file
        ONAME = PROMPTMFILE(
     &      'Enter logical name for UNGRIDDED MIN/MAX TEMPERATURE file',
     &      FSUNKN3, 'MINMAXT', 'PREDIUR' )

C.........  Prompt for EF-based output file of temp combo indices
        EDEV  = PROMPTFFILE(
     &          'Enter logical name for output EF-REF/TEMP file',
     &          .FALSE., .TRUE., 'MEFTEMP', 'PREDIUR' )

C...........   Read and process emissions factors cross reference file MPLIST
 
        WRITE( *,92000 )
     &         ' ', 'Reading EMISSION FACTORS XREF file ...', ' '
 
        NMPLIST = RDMPLIST( NLDEV, MPFIP, MPRDT, MPLNK,
     &                      MPREF, REFCNT, REFLIST )

        IF ( NMPLIST .LE. 0 ) THEN
            CALL M3EXIT( 'PREDIUR', 0,0, 'Error reading MPLIST', 2 )
        END IF

C...........   Map emissions factor indexes onto sources, and store each
C...........   unique, pure emission factor index:
 
        WRITE( *,92000 )
     &         ' ', 'Processing EMISSION FACTORS XREF ...', ' '

        EFLAG = .FALSE.
        DO  93 S = 1, NMSRC
 
            FIP = IFIP  ( S )
            RDT = IRCLAS( S )
            LNK = ILINK ( S )
 
C............   Link-specific matches:
 
            F = FIND3( FIP, RDT, LNK, NMPLIST, MPFIP, MPRDT, MPLNK )
            IF ( F .GT. 0 ) GO TO 77 
 
C............   FIP-RDT matches:
 
            F = FIND3( FIP, RDT, 0, NMPLIST, MPFIP, MPRDT, MPLNK )
            IF ( F .GT. 0 ) GO TO 77
 
C............   state--RDT matches:
 
            F = FIND3( 1000 * ( FIP / 1000 ), RDT, 0, 
     &                 NMPLIST, MPFIP, MPRDT, MPLNK )
            IF ( F .GT. 0 ) GO TO 77
 
C............   State-only matches:
            F = FIND3( FIP, 0, 0, NMPLIST, MPFIP, MPRDT, MPLNK )
            IF ( F .GT. 0 ) GO TO 77
 
C............   RDT-only matches:
 
            F = FIND3( 0, RDT, 0, NMPLIST, MPFIP, MPRDT, MPLNK )
            IF ( F .GT. 0 ) GO TO 77
 
C............   fallback

            F = FIND3( 0, 0, 0, NMPLIST, MPFIP, MPRDT, MPLNK )
            IF ( F .GT. 0 ) GO TO 77
 
C...........   If you get to here:  error
 
            EFLAG = .TRUE.
            WRITE( MESG,94010 )
     &            'No emission factor xref for FIP', FIP,
     &            'RCLASS', RDT, 'Link', LNK
            CALL M3MESG( MESG )
 
  77        CONTINUE

            J = 0

C.............  Loop through all time steps for this source's PSIs
            DO 91 T = 0, 23

                K = LOCAT1( MPREF( T,F ), J, PSI( 1,S ) )

                IF( J .EQ. 0 ) THEN

                    J = 1
                    PSI( 1,S ) = MPREF( T,F ) 

                ELSEIF( K .GT. 0 ) THEN

C.....................  Check for overflow...
                    IF( J .EQ. MXPPS ) THEN
                        WRITE( MESG, 94010 )
     &                         'Maximum dimension MXPPS=', MXPPS,
     &                         'was exceeded'
                        CALL M3EXIT( 'PREDIUR', 0, 0, MESG, 2 )
 
C.....................  If no overflow...
                    ELSE

C.........................  Add one to count of PSIs for source
                        J = J + 1
 
C.........................  Shift stored PSIs down one for MPREF(T,F)
                        DO 88  JJ = J, K+1, -1
                            PSI( JJ,S ) = PSI( JJ-1,S )
   88                   CONTINUE
 
C.........................  Insert new MPREF for source (S)
                        PSI ( K,S ) = MPREF(T,F)

                    ENDIF
                ENDIF
   91       CONTINUE

            NPSI( S ) = J

   93   CONTINUE

        IF ( EFLAG ) THEN
            CALL M3EXIT( 'PREDIUR', 0, 0,
     &                   'Bad emission factor matchup', 2 )
        END IF


C.........  Initialize NTIDX and TIDXPSI 

        DO 99 K = 1, REFCNT
            NTIDX( K ) = 0

            DO 96 J = 1, MXTMMI
                TIDXPSI( J,K ) = 0
   96       CONTINUE

   99   CONTINUE

C.........  Initialize constants for loop
        MAXTL = INT( ( MAXT_MAX - MAXT_MIN ) / TINTV ) + 1

        WRITE( *,92000 ) 
     &        ' ', 'Processing TEMPERATURE data ...', ' '

C.........  Loop through days/hours of temperature files

        IFIL = 1
        DO 499 H = 1, NSTEPS

C.................................................................
C.............  Process last hour of a day or last hour
C.................................................................

C           IF( ( JTIME/10000 .EQ. 23 .OR. H .EQ. NSTEPS ) 
C    &          .AND. H .NE. 1 ) THEN
            IF( MOD( H, 24 ) .EQ. 23 ) THEN

C.................  Try to read current temperature file
                IF ( .NOT. READ3( TNAME( IFIL ),
     &               'TA', 1, MDATE, MTIME, TA ) ) THEN

C.....................  If last read failed, check if can go to next file
                    IF( IFIL .LT. NTFIL ) THEN
                        IFIL = IFIL + 1

C.........................  Read next file if possible
                        IF ( .NOT. READ3( TNAME( IFIL ),
     &                       'TA', 1, MDATE, MTIME, TA ) ) THEN
                            MESG = 'Could not read TA from ' // 
     &                             TNAME( IFIL )
                            CALL M3EXIT 
     &                           ( 'PREDIUR', MDATE, MTIME, MESG, 2 )
                        ENDIF

                    ELSE
                        MESG = 'Could not read TA from '//TNAME( IFIL )
                        CALL M3EXIT ( 'PREDIUR', MDATE, MTIME, MESG, 2 )
                    ENDIF
                END IF

                KU = 0
                DO 199 S = 1, NMSRC

C.....................  Apply ungridding matrix 
                    TASRC = 0.0
                    DO 101 J = 1, NU( S )
                        KU = KU + 1
                        TASRC = TASRC + TA( IU( KU ) ) * CU( KU )
  101               CONTINUE

C.....................  Update minimum and maximum per source
                    IF( NU( S ) .EQ. 0 ) THEN    !  source outside grid
 
                        OFLAG = .TRUE.
                        OSRC = OSRC + 1
                        TKMIN( S ) = 0.0
                        TKMAX( S ) = 0.0

                    ELSEIF( TASRC .LT. TKMIN( S ) ) THEN
                        TKMIN( S ) = TASRC

                    ELSEIF( TASRC .GT. TKMAX( S ) ) THEN
                        TKMAX( S ) = TASRC

                    ENDIF

C.....................  Round min temperature DOWN to nearest on interval

                    TF   = CTOF * ( TKMIN( S ) - CTOK ) + 32.

                    IF    ( TF .LT. MINT_MIN ) THEN
                        WRITE( MESG, 94020 ) 
     &                         'Setting TKMIN for source', S, 'from',
     &                         TF, '[deg F] to', MINT_MIN
                        CALL M3MESG( MESG )

                        TKMIN( S ) = CTOK + 
     &                               ( FTOC * ( MINT_MIN - 32.0 ) )
                        TMIN = MINT_MIN

                    ELSEIF( TF .GT. MINT_MAX ) THEN
                        WRITE( MESG, 94020 )
     &                         'Setting TKMIN for source', S, 'from',
     &                         TF, '[deg F] to', MINT_MAX
                        CALL M3MESG( MESG )

                        TKMIN( S ) = CTOK + 
     &                               ( FTOC * ( MINT_MAX - 32.0 ) )
                        TMIN = MINT_MAX

                    ELSE  ! Note: round DOWN
                        TMIN = MINT_MIN + TINTV *
     &                         INT( ( TF - MINT_MIN )/TINTV )

                    ENDIF

C.....................  Round max temperature to nearest on interval

                    TF   = CTOF * ( TKMAX( S ) - CTOK ) + 32.

                    IF    ( TF .LT. MAXT_MIN ) THEN
                        WRITE( MESG, 94020 )
     &                         'Setting TKMAX for source',
     &                         S, 'from', TF, '[deg F] to', MAXT_MIN
                        CALL M3MESG( MESG )

                        TKMAX( S ) = CTOK + 
     &                               ( FTOC * ( MAXT_MIN - 32.0 ) )
                        TMAX = MAXT_MIN

                    ELSEIF( TF .GT. MAXT_MAX ) THEN
                        WRITE( MESG, 94020 )
     &                         'Setting TKMAX for source',
     &                         S, 'from', TF, '[deg F] to', MAXT_MAX
                        CALL M3MESG( MESG )

                        TKMAX( S ) = CTOK + 
     &                               ( FTOC * ( MAXT_MAX - 32.0 ) )
                        TMAX = MAXT_MAX

                    ELSE  ! Note: round DOWN
                        TMAX = MAXT_MIN + TINTV *
     &                         INT( ( TF - MAXT_MIN )/TINTV )
                    ENDIF

C.....................  Ensure relationship between TMIN and TMAX is valid
C.........................  When difference greater than forty, 
C                  .......  set TMIN from TMAX, and adjust down on interval
                    IF( (TMAX - TMIN) .GE. 40 ) THEN
                        MSAV= TMIN
                        TMIN = TMAX - 39.
                        TMIN = MINT_MIN + TINTV *
     &                        INT( ( TMIN - MINT_MIN )/TINTV )


                        WRITE( MESG, 94020 ) 'Temperature: TMAX - TMIN'
     &                         // ' is > 40 for source', S,
     &                         '.  Changing TMIN from', MSAV, 'to', TMIN
                        CALL M3MESG( MESG )

                    ENDIF

C.....................  If TMIN and TMAX differ by less than delta T, need
C.....................  to avoid bad TMMI values.  So so set TMMI to zero,
C.....................  which means no diurnal emissions!
                    IF( TMAX - TMIN .LE. MAX( TINTV, 1.0 ) ) THEN
               
                        TMMI ( S )= 0
                        TMMIT( 1 )= 0
                        TMMIT( 2 )= 0
                        TMMIT( 3 )= 0
                        TMMIT( 4 )= 0

                        WRITE( MESG, 94050 ) 'Temperature: TMAX - TMIN'
     &                         // ' is <= ', MAX( TINTV, 1.0 ), 
     &                         'for source', S,
     &                         '.  No diurnal emissions will exist.'
                        CALL M3MESG( MESG )

                    ELSE

C.........................  Create vector for temperature combos index
                        TMMI( S ) = MAXTL * 
     &                          NINT( ( TMIN - MINT_MIN ) / TINTV ) + 
     &                          NINT( ( TMAX - MAXT_MIN ) / TINTV ) + 1

C.........................  Store the four indices needed for this source's 
C.........................  temperatures on sorted order
                        TMMIT( 1 )= TMMI( S )
                        TMMIT( 2 )= TMMIT( 1 ) + 1
                        TMMIT( 3 )= MIN( TMMIT( 1 ) + MAXTL, NMINMAX )
                        TMMIT( 4 )= MIN( TMMIT( 3 ) + 1    , NMINMAX )

                    ENDIF

C.....................  Store min/max temp index in all PSI(s)
C               ......  Loop through all PSI(s) for source
                    DO 151 I = 1, NPSI( S )

C.........................  Get index into EF reference list
                        K = FIND1( PSI( I,S ), REFCNT, REFLIST )

C.........................  Loop through the four TMMIs for source S
                        DO 141 J = 1, 4

C.............................  Locat TMMIT in temp combo list for this PSI
                            L = LOCAT1( 
     &                          TMMIT( J ), NTIDX( K ), TIDXPSI( 1,K ) )

                            IF( NTIDX( K ) .EQ. 0 ) THEN

                                NTIDX( K ) = 1
                                TIDXPSI( 1,K ) = TMMIT( J )

                            ELSEIF( L .GT. 0 ) THEN

C.................................  Check for overflow...
                                IF( NTIDX( K ) .EQ. MXTMMI ) THEN
                                    WRITE( MESG, 94010 )
     &                                   'Maximum dimension MXTMMI=',
     &                                   MXTMMI, 'was exceeded'
                                    CALL M3EXIT( 'PREDIUR', JDATE,
     &                                           JTIME, MESG, 2    )

C.................................  If no overflow...
                                ELSE

C.....................................  Add one to count of TMMIs for PSI (K)
                                    NTIDX( K ) = NTIDX( K ) + 1

C.....................................  Shift stored TMMIs down one for PSl (K)
                                    DO 121 JJ = NTIDX( K ), L+1, -1
                                      TIDXPSI( JJ,K )=TIDXPSI( JJ-1,K )
  121                               CONTINUE

C.....................................  Insert new TMMIT(J) for PSI(K)
                                    TIDXPSI ( L,K ) = TMMIT( J )

                                ENDIF  ! End overflow check

                            ENDIF  ! End check if TMMI in list already
  141                   CONTINUE   ! End loop through source's TMMIs
  151               CONTINUE       ! End loop through source's PSIs

  199           CONTINUE

C.................  Give warning if end of day is not hour 23

                IF( JTIME .NE. 230000 ) THEN 
                    WRITE( MESG, 94010 )
     &                     'Day not ending on HHMMSS 230000'

                    CALL M3WARN( 'PREDIUR', JDATE, JTIME, MESG )

                ENDIF

C.................  Write out min/max information for the current day
                IF( .NOT. WRITE3( 
C    &              ONAME, 'ALL', JDATE, 0, TMMI ) ) THEN
     &              ONAME, 'ALL', JDATE, FIRSTIME, TMMI ) ) THEN 

                    CALL M3EXIT( 'PREDIUR', JDATE, 0,
     &                           'Could not write timestep to "' //
     &                           ONAME( 1 : TRIMLEN( ONAME ) ) // 
     &                           '".', 2 )
                ENDIF

C.................................................................
C.............  Process first hour of a new day
C.................................................................
C           ELSEIF( JTIME .EQ. 0 .OR. H .EQ. 1 ) THEN
            ELSEIF( MOD( H, 24 ) .EQ. 1 ) THEN

C.................  Try to read current temperature file
                IF ( .NOT. READ3( TNAME( IFIL ),
     &               'TA', 1, MDATE, MTIME, TA ) ) THEN

C.....................  If last read failed, check if can go to next file
                    IF( IFIL .LT. NTFIL ) THEN
                        IFIL = IFIL + 1

C.........................  Read next file if possible
                        IF ( .NOT. READ3( TNAME( IFIL ),
     &                       'TA', 1, MDATE, MTIME, TA ) ) THEN
                            MESG = 'Could not read TA from ' // 
     &                             TNAME( IFIL )
                            CALL M3EXIT 
     &                           ( 'PREDIUR', MDATE, MTIME, MESG, 2 )
                        ENDIF

                    ELSE
                        MESG = 'Could not read TA from '//TNAME( IFIL )
                        CALL M3EXIT ( 'PREDIUR', MDATE, MTIME, MESG, 2 )
                    ENDIF
                END IF

                KU = 0
                DO 299 S = 1, NMSRC

C.....................  Apply ungridding matrix 
                    TASRC = 0.0
                    DO  201  J = 1, NU( S )
                        KU = KU + 1
                        TASRC = TASRC + TA( IU( KU ) ) * CU( KU )
  201               CONTINUE

C.....................  Initialize minimum and maximum per source
                    IF( NU( S ) .GT. 0 ) THEN

                        TKMIN( S ) = TASRC
                        TKMAX( S ) = TASRC

                    ENDIF

  299           CONTINUE

C.................................................................
C.............  Process hours in middle of day
C.................................................................
            ELSE  

C.................  Try to read current temperature file
                IF ( .NOT. READ3( TNAME( IFIL ),
     &               'TA', 1, MDATE, MTIME, TA ) ) THEN

C.....................  If last read failed, check if can go to next file
                    IF( IFIL .LT. NTFIL ) THEN
                        IFIL = IFIL + 1

C.........................  Read next file if possible
                        IF ( .NOT. READ3( TNAME( IFIL ),
     &                       'TA', 1, MDATE, MTIME, TA ) ) THEN
                            MESG = 'Could not read TA from ' // 
     &                             TNAME( IFIL )
                            CALL M3EXIT 
     &                           ( 'PREDIUR', MDATE, MTIME, MESG, 2 )
                        ENDIF

                    ELSE
                        MESG = 'Could not read TA from '//TNAME( IFIL )
                        CALL M3EXIT ( 'PREDIUR', MDATE, MTIME, MESG, 2 )
                    ENDIF
                END IF

                KU = 0
                DO 399 S = 1, NMSRC

C.....................  Apply ungridding matrix 
                    TASRC = 0.0
                    DO 301 J = 1, NU( S )
                        KU = KU + 1
                        TASRC = TASRC + TA( IU( KU ) ) * CU( KU )
  301               CONTINUE

C.....................  Update minimum and maximum per source
                    IF( NU( S ) .EQ. 0 ) THEN    !  source outside grid
 
                        TKMIN( S ) = 0.0
                        TKMAX( S ) = 0.0

                    ELSEIF( TASRC .LT. TKMIN( S ) ) THEN
                        TKMIN( S ) = TASRC

                    ELSEIF( TASRC .GT. TKMAX( S ) ) THEN
                        TKMAX( S ) = TASRC

                    ENDIF

  399           CONTINUE

            ENDIF  !  If end of day, last hour, first hour of day, or other

            CALL NEXTIME( JDATE, JTIME, TSTEP )
            CALL NEXTIME( MDATE, MTIME, TSTEP )

  499   CONTINUE   !  End loop on hours of temperature files
 
C......... Write temp combos indices for each PSI into ASCII file

        WRITE( *,92000 ) 
     &         ' ', 'Writing out EF-REF/TEMPERATURE file ...', ' '

        DO 521 K = 1, REFCNT

            WRITE( EDEV, 93010 ) REFLIST( K ), NTIDX( K ), 
     &                           ( TIDXPSI( J,K ), J = 1, NTIDX( K ) )

  521   CONTINUE

        IF( OFLAG ) THEN

            WRITE( MESG, 94010 )
     &         'During ungridding, ', OSRC,' sources excluded from grid'

            CALL M3WARN( 'PREDIUR', 0, 0, MESG )

        ENDIF


C......... End program sucessfully

        MESG = 'Normal completion of PREDIUR'
        CALL M3EXIT( 'PREDIUR', 0, 0, MESG, 0 )
 
C******************  FORMAT  STATEMENTS   ******************************
 
C...........   Informational (LOG) message formats... 92xxx
 
92000   FORMAT( 5X, A )
92010   FORMAT( 5X, A, 2( /, 5X, A, 1X, I7.7, 1X, I6.6, 1X, I2 ) )
 
C...........   Formatted file I/O formats............ 93xxx
 
93010   FORMAT( I5, ',', I8, ',', 20( I4, : ',' ), /, 
     &                   49( 15X, 20( I4, : ',' ), / ) )
 
C...........   Internal buffering formats............ 94xxx
 
94010   FORMAT( 10 ( A, :, 1X, I5, :, 1X ) )
94020   FORMAT( 10 ( A, 1X, I5, 2( 1X, A, 1X, F8.2 ) ) )
94030   FORMAT( A, F10.3, 1X, A )
94040   FORMAT( A5, I2.2 )
94050   FORMAT( 10 ( A, 1X, F8.2, 2( 1X, A, 1X, I5 ) ) )

        END



