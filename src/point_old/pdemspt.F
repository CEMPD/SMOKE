
         PROGRAM PDEMSPT

C***********************************************************************
C  program body starts at line 214
C
C  DESCRIPTION:
C       Process data contained in EMS-95 hour-specific source files,
C       match to inventory file, and write hour-specific inventory file 
C
C  PRECONDITIONS REQUIRED:
C       EMS-95 Input data for point sources.  Accepts multiple 
C       process.pt and pdayspec.pt files
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       Models-3 I/O
C       FINDC, GETYN, PROMPTFFILE, PROMPTMFILE, SORTIC, TRIMLEN
C       FIXSTK
C
C  REVISION  HISTORY:
C       Prototype  5/97 by MR Houyoux
C
C****************************************************************************/
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1998, MCNC--North Carolina Supercomputing Center
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

      IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'PTDIMS3.EXT'   !  point-source dimensioning parameters
        INCLUDE 'CHDIMS3.EXT'   !  emis chem info (both inventory and model)
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS and their descriptions:
        
        CHARACTER*2     CRLF
        LOGICAL         DSCGRID
        INTEGER         ENVINT
        LOGICAL         ENVYN
        INTEGER         FINDC   !  returns -1 for failure
        INTEGER         GETDATE
        INTEGER         GETNUM 
        LOGICAL         GETYN
        INTEGER         INDEX1
        INTEGER         JULIAN
        INTEGER         JUNIT
        INTEGER         LBLANK
        INTEGER         PROMPTFFILE
        CHARACTER*16    PROMPTMFILE
        INTEGER         SECSDIFF
        INTEGER         STR2INT
        REAL            STR2REAL
        INTEGER         TRIMLEN
        INTEGER         YEAR4

        EXTERNAL CRLF, DSCGRID, ENVINT, ENVYN, FINDC, GETNUM, GETDATE,
     &           GETYN, INDEX1, JULIAN, JUNIT, LAMBERT, LBLANK, 
     &           LL2LAM, UTM2LAM, PROMPTFFILE, PROMPTMFILE, SECSDIFF, 
     &           STR2INT, STR2REAL, TRIMLEN, YEAR4


C...........   PARAMETERS and their descriptions

        CHARACTER*5     BLANK5
        INTEGER         NPHRIP  ! Maximum number of hourly sources by inv pol
        INTEGER         NZONE   ! Maximum number of time zones

        PARAMETER     ( BLANK5  = ' ',
     &                  NPHRIP  = NPHR * NIPOL, 
     &                  NZONE   = 9  )

C...........   LOCAL VARIABLES and their descriptions:
C...............   Time zone tables:  FIP-independent; state-only; state-county

C...........   Point sources table for process.pt files  (unsorted; sorted)

        INTEGER        PS                ! Counter & number of process records
        INTEGER        INDXPA( NPEMAX )  ! subscript table
        CHARACTER*40   PSKEYA( NPEMAX )  ! FIP // SCC  // FCID // SKID
        CHARACTER*56   PSLSTA( NPEMAX )  ! FIP // FCID // SKID // DVID // PRID

        CHARACTER*40   PSKEY ( NPEMAX )  ! FIP // SCC  // FCID // SKID
        CHARACTER*56   PSLST ( NPEMAX )  ! FIP // FCID // SKID // DVID // PRID

C.......   Source list for matching 

        CHARACTER*40   SRCLST( NPSRC )

C...........   Unsorted record-days characteristics

        INTEGER        NPDPT ( MXHR )
        INTEGER        CODEA ( NPHRIP, MXHR )
        INTEGER        IDXSRC( NPHRIP, MXHR )
        INTEGER        SRCIDA( NPHRIP, MXHR )
        INTEGER        PRCIDA( NPHRIP, MXHR )
        REAL           EMISVA( NPHRIP, MXHR )

C...........   Hour-specific point Sources Table input 
C...........   sorted final version (image EMISREC of data record)
 
        INTEGER       NHRPT               ! current source-count
        INTEGER       INDXH( NPHR )       ! subscript into netCDF inventory file
        REAL          EMISH( NPHR, MXIPOL )! emissions values (tons/hr)
 
        COMMON / EMISREC / NHRPT, INDXH, EMISH

C............  Helper variables for hourly summing across multiple sources

        INTEGER       NPRHR              ! number of processes encounter per hour
        CHARACTER*62  PRHRLST( NPEMAX*NIPOL )! FIP//FCID//SKID//DVID//PRID//EINAM

C...........   For flexible output of only species for which there are inputs
        INTEGER       INDXO ( MXIPOL ) ! For each inv pol, how map to OUTVAR
        LOGICAL       OVFLAG( MXIPOL ) ! For each inv pol, will OUTVAR exist?
        CHARACTER*16  OUTVAR( MXIPOL ) ! Output variable names

C...........   File units and logical/physical names

        INTEGER         CDEV    ! Unit number for EMS95/SMOKE facility/stack ID
        INTEGER         DDEV    ! Unit number for process.pt or pdayspec.pt files
        INTEGER         IDEV    ! List of physical raw input file names
        INTEGER         LDEV    ! log-device

        CHARACTER*16    ENAME   ! point source inventory file
        CHARACTER*16    HNAME   !  hour-specific output emissions file

C...........   Other local variables

        REAL            EMIS( 24 )    !  Temporary emission value
                                
        INTEGER         I, H, J, K, L, P, S, T, V !  counters and indices
        INTEGER         HS
        INTEGER         L1, L2           !  counters and indices
        INTEGER         FIP     ! Temporary fip

        INTEGER         COD     !  Temporary inventory pollutant code
        INTEGER         CSS     !  Start of non-blank character string
        INTEGER         DAY     !  temporary day-of-week code
        INTEGER         IFIL    !  File counter
        INTEGER         IOS     !  I/O status
        INTEGER         IREC    !  input line (record) number
        INTEGER         JDATE   !  Current Julian date
        INTEGER         JTIME   !  Current HHMMSS time
        INTEGER         LS      !  Previous source ID
        INTEGER         MAXHS   !  maximum number of day-specific records needed
        INTEGER         MAXPTR  !  maximum time step pointer
        INTEGER         MINPTR  !  minimum time step pointer
        INTEGER         MONTH   !  temporary month code
        INTEGER         NHFIL   !  number of hour-specific data files
        INTEGER         NOUTV   !  number of output variables
        INTEGER         NPFIL   !  number of process data files
        INTEGER         NSTEPS  !  number of output time steps
        INTEGER         PTR     !  time step pointer
        INTEGER         RDATE   !  Reference date
        INTEGER         RTIME   !  Reference time
        INTEGER         SDATE   !  Starting output date 
        INTEGER         STIME   !  Starting output time
        INTEGER         TSTEP   !  time step (set to 10000)
        INTEGER         TZONE   !  Output time zone
        INTEGER         YEAR    !  temporary year code
        INTEGER         ZONE    !  temporary time zone code

        LOGICAL         EFLAG   !  TRUE iff processing error
        LOGICAL         PROMPTF !  TRUE iff PROMPTFLAG E.V. is Y
        LOGICAL         SFLAG   !  input verification:  report missing species
        LOGICAL         VFLAG   !  input verification:  TRUE iff ERROR

        CHARACTER*5     CFIP    !  Character FIP code
        CHARACTER*5     CPOL    !  Temporary pollutant code
        CHARACTER*8     CSCC    !  Character SCC code
        CHARACTER*12    DVID    !  Temporary device ID
        CHARACTER*12    PRID    !  Temporary process ID
        CHARACTER*12    SKID    !  Temporary stack ID
        CHARACTER*15    FCID    !  Temporary facility ID
        CHARACTER*56    PKEY
        CHARACTER*62    PVKEY
        CHARACTER*256   LINE    !  input line from POINT file
        CHARACTER*256   MESG    !  text for M3EXIT()
        CHARACTER*256   NAMTMP  !  Temporary buffer for input file names

        CHARACTER*3     ZONAM( NZONE ) 
        DATA   ZONAM / 'PST', 'PDT', 'MST', 'MDT', 'CST', 
     &                 'CDT', 'EST', 'EDT', 'GMT'         /

        INTEGER         ZONUM( NZONE )
        DATA   ZONUM / 8, 7, 7, 6, 6, 
     &                 5, 5, 4, 0     /

C***********************************************************************
C   begin body of program PDEMSPT

        LDEV = INIT3()

        CALL INITEM( LDEV )

        WRITE( *,92000 ) 
     &  ' ',
     &  'Program PDEMSPT to take sequence of EMS-95 point source',
     &  'foundation files and produce the POINT SOURCE EMISSIONS',
     &  'HOUR-SPECIFIC file. ',
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
            CALL M3EXIT( 'PDEMSPT', 0, 0, 'Ending program.', 2 )
        END IF

        SFLAG = ENVYN( 'RAW_SRC_CHECK',
     &                 'PDEMSPT check for missing species-records',
     &                 .FALSE., IOS )

        PROMPTF = ENVYN ( 'PROMPTFLAG', 'Prompt for inputs or not',
     &                    .TRUE., IOS )

        TZONE  = ENVINT( 'OUTZONE', 'Output time zone', 5, IOS )
        IF( IOS .NE. 0 .AND. .NOT. PROMPTF ) THEN
            EFLAG = .TRUE.
            MESG = 'Env var PROMPTFLAG is N, but OUTZONE is bad'
            CALL M3MSG2( MESG )
        ENDIF

        IF( EFLAG ) THEN
            CALL M3EXIT( 'PDEMSPT', 0, 0,
     &                    'Bad environment variable setting(s)', 2 )
        ENDIF

        IDEV = PROMPTFFILE( 
     &         'Enter the name of the RAW DATA FILENAME LISTING',
     &          .TRUE., .TRUE., 'PHLST', 'PDEMSPT' )

C.........  Get inventory file needed for setting up output file

        ENAME = PROMPTMFILE(
     &          'Enter logical name for POINT source inventory file',
     &          FSREAD3, 'PNTS', 'PDEMSPT' )

        CDEV = PROMPTFFILE(
     &          'Enter the name of the FACILITY/STACK input file',
     &          .TRUE., .TRUE., 'PSRC', 'PDEMSPT' )


C.........  Get output time zone
        TZONE = GETNUM( -12, 12, TZONE,
     &                  'Enter output time zone (0 for GMT, 5 for ' //
     &                  'EST, 4 for EDT)' )

C.............  Preprocess PHLST file here to be able to get good defaults
C.............  for output dates and times.

C.............  Loop through the PHLST file to make sure files are in the
C.............  right order, and to count how many files of each type.

        IFIL  =  0
        NHFIL =  0
        NPFIL =  0
        RDATE =  1980001  ! Jan 1, 1980 is the reference date
        RTIME =  0
        EFLAG = .FALSE.
44      CONTINUE    ! Head of file names loop

            READ( IDEV, 93000, END=55  ) NAMTMP  

C.............  Open files, and report status.
            DDEV   = JUNIT()
            OPEN( DDEV, ERR=1006, FILE=NAMTMP, STATUS='OLD' )

            IFIL = IFIL + 1

C.............  Read first line of file and determine its type
            READ( DDEV,93000 ) LINE

            REWIND( DDEV )

            L1 = TRIMLEN( LINE )

C.............  File type is process.pt
            IF( L1 .LE. 92 ) THEN 
                NPFIL = NPFIL + 1 

                IF( NHFIL .GT. 0 ) THEN
                    CALL M3EXIT( 'PDEMSPT', 0, 0, 'All process.pt ' //
     &                   'files must be at head of PHLST file.', 2 )
                ENDIF

C.............  File type is not process.pt but it is first file
            ELSEIF( IFIL .EQ. 1 ) THEN
                CALL M3EXIT( 'PDEMSPT', 0, 0, 'First file in PHLST ' //
     &                       'must be a process.pt file.', 2 )

C.............  File type is pdayspec.pt
            ELSE

                NHFIL = NHFIL + 1

                IREC = 0
49              CONTINUE   ! head of hourly preprocessing read loop

C.....................  Read line of file
                    READ( DDEV, 93000, END=51 ) LINE

                    IREC = IREC + 1

C.....................  Set Julian day from MMDDYY8 SAS format
                    MONTH = STR2INT( LINE( 62:63 ) )
                    DAY   = STR2INT( LINE( 65:66 ) )
                    YEAR  = YEAR4( STR2INT( LINE( 68:69 ) ) )

                    JDATE = 1000 * YEAR + JULIAN( YEAR, MONTH, DAY )
                    JTIME = 0
 
C.....................  Get time zone and convert date and time to output zone

                    CALL UPCASE( LINE( 70:72 ) )
                    ZONE= ZONUM( INDEX1( LINE( 70:72 ), NZONE, ZONAM ) )

                    IF( ZONE .LT. 0 ) THEN
                        EFLAG = .TRUE.
                        WRITE( MESG,94010 ) 
     &                    'Bad time zone "' // LINE(70:72) // 
     &                    '" at line', IREC, 'in file' // CRLF() //
     &                    BLANK5 // NAMTMP( 1:TRIMLEN( NAMTMP) )
                        CALL M3MESG( MESG )

                    ELSE
                        CALL NEXTIME( JDATE, JTIME, 
     &                                ( ZONE - TZONE ) * 10000 )

                    ENDIF

                    PTR = SECSDIFF( RDATE, RTIME, 
     &                              JDATE, JTIME ) / 3600 + 1

                    IF( NHFIL .EQ. 1 .AND. IREC .EQ. 1 ) THEN
                        MINPTR = PTR
                        MAXPTR = PTR
                    ELSEIF( PTR .LT. MINPTR ) THEN
                        MINPTR = PTR
                    ELSEIF( PTR + 23 .GT. MAXPTR ) THEN
                        MAXPTR = PTR + 23
                    ENDIF 

                GO TO 49

51              CONTINUE

            ENDIF 

            CLOSE ( DDEV )

        GO TO 44       ! to head of files preprocessing loop  

55      CONTINUE

        IF ( EFLAG ) THEN
            CALL M3EXIT( 'PDEMSPT', 0, 0,
     &           'Error in date formatting of day-specific data', 2 )
        END IF

        REWIND( IDEV )

C.........  Set default starting date and time, and default duration
C.........  Use reference time and the 1st hour in all files to
C.........  determine the starting time.
        SDATE = RDATE
        STIME = RTIME
        CALL NEXTIME( SDATE, STIME, ( MINPTR - 1 ) * 10000 )

        NSTEPS = MAXPTR - MINPTR + 1

C.........  Get output starting date
        SDATE = GETDATE( SDATE,
     &                   'Enter start date (YYYYDDD|YYYYMMDD) for ' //
     &                   'HOURLY-POINT file' )

C.........  Get output starting time
        STIME = GETNUM( 0, 240000, STIME,
     &          'Enter start time (HHMMSS) for HOURLY-POINT file' )

        TSTEP = 10000

        NSTEPS = GETNUM( 0, 9999999, NSTEPS, 
     &                   'Enter output duration (hours)' )

        IF( NSTEPS .GT. MXHR ) THEN
            WRITE( MESG, 94010 )
     &        'Maximum number of hours (MXHR) =', MXHR, 
     &         CRLF() // BLANK5 //
     &        'But duration requested was', NSTEPS, '. Must recompile.'

            CALL M3EXIT( 'PDEMSPT', 0, 0, MESG, 2 )
        ENDIF

C.........  Check dimensions of the input inventory file
        IF ( .NOT. DESC3( ENAME ) ) THEN
            CALL M3EXIT( 'PDEMSPT', 0, 0,
     &                   'Could not read file description for ' //
     &                   ENAME, 2 )
        ELSE IF ( NROWS3D .NE. NPSRC ) THEN
            WRITE( MESG, 94010 )
     &      'Dimension mismatch.  EMISSIONS file:', NROWS3D,
     &      'program:', NPSRC
            CALL M3EXIT( 'PDEMSPT', 0, 0, MESG, 2 )
        END IF

C.........  Read in source list FIPS code, SCC, Plant, Stack, EMS Plt, EMS Stk

        IREC = 0
66      CONTINUE

            READ( CDEV, 93010, END=69, IOSTAT=IOS ) 
     &            CFIP, CSCC, FCID, SKID
  
            IREC = IREC + 1 

            IF( IOS .GT. 0 ) THEN
                WRITE( MESG, 94010 ) 
     &              'Error ', IOS,  'reading PSRC file at line', IREC
                CALL M3EXIT( 'PDEMSPT', 0, 0, MESG, 2 )

            ELSEIF( IREC .LE. NPSRC ) THEN

                L1 = LBLANK( FCID )
                L2 = LBLANK( SKID )

                SRCLST( IREC ) = CFIP // CSCC // 
     &                           FCID( L1+1:TRIMLEN( FCID ) ) //
     &                           SKID( L2+1:TRIMLEN( SKID ) ) 

            ENDIF 

        GO TO 66  !  to head of read loop 

69      CONTINUE

        IF( IREC .GT. NPSRC ) THEN
            WRITE( MESG, 94010 )
     &      'Dimension mismatch. PSRC file:', IREC, 'program:', NPSRC
            CALL M3EXIT( 'PDEMSPT', 0, 0, MESG, 2 )
        ENDIF

C.........  Loop through all of the process files in the file list

        EFLAG = .FALSE.
        PS    =   0  !  For process.pt's
        DO 201 IFIL = 1, NPFIL

            READ( IDEV, 93000 ) NAMTMP  

C.............  Open files, and report status.
            DDEV   = JUNIT()
            OPEN( DDEV, ERR=1006, FILE=NAMTMP, STATUS='OLD' )

            WRITE( MESG,94010 ) 'Successful OPEN ' //
     &             'for PROCESS file:' // CRLF() // BLANK5 //
     &             NAMTMP( 1:TRIMLEN( NAMTMP ) )
            CALL M3MSG2( MESG )

C.............  Process each process.pt file and store unsorted arrays
C.............  for using in searching stage
            IREC = 0
133         CONTINUE

C.................  Read a line of process.pt file, and check input status

                READ( DDEV, 93000, END=134, IOSTAT=IOS ) LINE

                IREC = IREC + 1

                IF ( IOS .NE. 0 ) THEN
                    WRITE( MESG, 94010 ) 
     &                  'Error ', IOS,  'reading "' // 
     &                  NAMTMP( 1:TRIMLEN( NAMTMP ) ) //
     &                  '" at line', IREC
                    CALL M3EXIT( 'PDEMSPT', 0, 0, MESG, 2 )
                END IF

C.............  Convert and check SCC value

                PS = PS + 1

                IF ( PS .LE. NPEMAX ) THEN

                    WRITE ( CFIP,93020 ) 
     &                   1000 * STR2INT( LINE( 1:2 ) ) +
     &                          STR2INT( LINE( 3:5 ) )

                    CSCC = LINE( 57:64 )
                    IF( CSCC .EQ. ' ' ) CSCC = CMISS3

                    CSS  = LBLANK( LINE( 6:20 ) )
                    FCID = LINE( MIN(CSS+6,20):20 )
                    IF( FCID .EQ. ' ' ) FCID = CMISS3

                    CSS  = LBLANK( LINE( 21:32 ) )
                    SKID = LINE( MIN(CSS+21,32):32 )
                    IF( SKID .EQ. ' ' ) SKID = CMISS3

                    CSS  = LBLANK( LINE( 33:44 ) )
                    DVID = LINE( MIN(CSS+33,44):44 )
                    IF( DVID .EQ. ' ' ) DVID = CMISS3

                    CSS  = LBLANK( LINE( 45:56 ) )
                    PRID = LINE( MIN(CSS+45,56):56 )
                    IF( PRID .EQ. ' ' ) PRID = CMISS3

                    INDXPA( PS ) = PS

                    PSKEYA( PS ) = CFIP // 
     &                             CSCC //
     &                             FCID( 1:TRIMLEN( FCID ) ) //
     &                             SKID( 1:TRIMLEN( SKID ) )

                    PSLSTA( PS ) = CFIP // 
     &                             FCID( 1:TRIMLEN( FCID ) ) //
     &                             SKID( 1:TRIMLEN( SKID ) ) //
     &                             DVID( 1:TRIMLEN( DVID ) ) //
     &                             PRID( 1:TRIMLEN( PRID ) )

                END IF      !  if ps in bounds

            GO TO  133      !  to head of DDEV-read loop

134         CONTINUE        !  end of the DDEV-read loop

            CLOSE( DDEV )

201     CONTINUE            !  exit from loop on input process.pt files

C.........  Check for overflow of process records
        IF ( PS .GT. NPEMAX ) THEN		!  overflow
            EFLAG = .TRUE.

            WRITE( MESG,94010 ) 
     &          'Actual HOUR-SPECIFIC SOURCE  record-count', PS, 
     &           CRLF() // BLANK5 //
     &          'Max dimensioned record-count (NPEMAX)=', NPEMAX
            CALL M3MSG2( MESG )

            CALL M3MSG2( 'Max record-count exceeded in PROCESS files.' )

        END IF	!  if input process parameter table overflow or not

        IF ( EFLAG ) THEN
            CALL M3EXIT( 'PDEMSPT', 0, 0,
     &                   'Error processing of process.pt data', 2 )
        END IF


C.........  Sort process records (that are not overflowed ) so we can use 
C.........  quick search

        CALL SORTIC( PS, INDXPA,  PSLSTA )
        DO  222  I = 1, PS
            J = INDXPA( I )
            PSKEY( I ) = PSKEYA( J )
            PSLST( I ) = PSLSTA( J )
222     CONTINUE


C.........  Initialize main emissions storage array

        DO 235 T = 1, NSTEPS

            DO  233  S = 1, NPHR 
                EMISVA( S, T ) = BADVAL3
233         CONTINUE

235     CONTINUE

C.........  Process all of the day specific data.  There is no readily
C.........  apparent way of dealing with multiple input time zones, and 
C.........  multiple files, because each file can have sources in 
C.........  different time zones, and therefore covering different output 
C.........  hours.  So, the current design will store ALL emissions.
 
        EFLAG  =  .FALSE.
        MAXHS  = 0
        MAXPTR = 0

        DO 301 IFIL = 1, NHFIL

            READ( IDEV, 93000 ) NAMTMP  

C.............  Open files, and report status.
            DDEV   = JUNIT()
            OPEN( DDEV, ERR=1006, FILE=NAMTMP, STATUS='OLD' )

            WRITE( MESG,94010 ) 'Successful OPEN ' //
     &             'for PDAYSPEC file:' // CRLF() // BLANK5 //
     &              NAMTMP( 1:TRIMLEN( NAMTMP ) )
            CALL M3MSG2( MESG )

            IREC = 0
266         CONTINUE  !  Head of pdayspec file read loop

C.................  Read first line of file
                READ( DDEV, 93000, END=299 ) LINE

                IREC = IREC + 1

C.................  Even though converting the dates again is repetetive,
C.................  it is better than having to store the dates the first
C.................  time for _every_ hour in the input files.  Such an
C.................  approach would eliminate the advantage of being able to
C.................  reduce NPHR when only a limited part of the input file
C.................  is actually used. M Houyoux.

C.................  Set Julian day from MMDDYY8 SAS format
                MONTH = STR2INT( LINE( 62:63 ) )
                DAY   = STR2INT( LINE( 65:66 ) )
                YEAR  = YEAR4( STR2INT( LINE( 68:69 ) ) )

                JDATE = 1000 * YEAR + JULIAN( YEAR, MONTH, DAY )
                JTIME = 0

C.................  Get time zone and convert date and time to output zone

                CALL UPCASE( LINE( 70:72 ) )
                ZONE= ZONUM( INDEX1( LINE( 70:72 ), NZONE, ZONAM ) )

                IF( ZONE .LT. 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 
     &                'Bad time zone "' // LINE(70:72) // 
     &                '" at line', IREC, 'in file' // CRLF() //
     &                BLANK5 // NAMTMP( 1:TRIMLEN( NAMTMP) )
                    CALL M3MESG( MESG )

                ELSE
                    CALL NEXTIME( JDATE, JTIME, 
     &                            ( ZONE - TZONE ) * 10000 )

                ENDIF

                PTR = SECSDIFF( SDATE, STIME, 
     &                          JDATE, JTIME ) / 3600 + 1

C.................  Skip record if it is out of range of output file

                IF( PTR. LT. 1 .OR. PTR .GT. NSTEPS ) GO TO 266

C.................  Store maximum time step number of data

                IF( PTR + 23 .GT. MAXPTR ) MAXPTR = PTR + 23

C.................  Check pollutant code and set index I

                CSS  = LBLANK( LINE( 57:61 ) )
                CPOL = LINE( MIN(CSS+57,61):61 )
                COD  = INDEX1( CPOL( 1:TRIMLEN( CPOL ) ), 
     &                             NIPOL, EINAM )

                IF ( COD .LE. 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 
     &                'Bad line', IREC, 
     &                'Pollutant code "' // LINE(57:61) // '" in "'
     &                // NAMTMP( 1:TRIMLEN( NAMTMP ) ) // '"'
                    CALL M3MESG( MESG )
                    GO TO  266      !  to head of loop
                END IF

C.................  Check and set emissions values

                L1 = 66 
                L2 = 72 
                DO 277 J = 1, 24

                    L1 = L1 + 7
                    L2 = L2 + 7

                    EMIS( J )  = STR2REAL( LINE( L1:L2 ) )
                    IF ( EMIS( J ) .LT. 0.0 )  THEN
                        EFLAG = .TRUE.
                        WRITE( MESG,94010 ) 
     &                    'Bad line', IREC, '" in "' //
     &                    NAMTMP( 1:TRIMLEN( NAMTMP ) ) // 
     &                    '"--emis value "' // LINE( L1:L2 ) // '"'
                        CALL M3MESG( MESG )
                        GO TO  266  ! to head of read loop
                    END IF

277             CONTINUE

C.................  Set key for searching process array

                FIP  = 1000 * STR2INT( LINE( 1:2 ) ) +
     &                        STR2INT( LINE( 3:5 ) )
                WRITE( CFIP,93020 ) FIP

                CSS  = LBLANK( LINE( 6:20 ) )
                FCID = LINE( MIN(CSS+6,20):20 )
                IF( FCID .EQ. ' ' ) FCID = CMISS3 

                CSS  = LBLANK( LINE( 21:32 ) )
                SKID = LINE( MIN(CSS+21,32):32 )
                IF( SKID .EQ. ' ' ) SKID = CMISS3 

                CSS  = LBLANK( LINE( 33:44 ) )
                DVID = LINE( MIN(CSS+33,44):44 )
                IF( DVID .EQ. ' ' ) DVID = CMISS3

                CSS  = LBLANK( LINE( 45:56 ) )
                PRID = LINE( MIN(CSS+45,56):56 )
                IF( PRID .EQ. ' ' ) PRID = CMISS3

C.................  Search for this record in PROCESS records

                PKEY = CFIP // 
     &                 FCID( 1:TRIMLEN( FCID ) ) //
     &                 SKID( 1:TRIMLEN( SKID ) ) //
     &                 DVID( 1:TRIMLEN( DVID ) ) //
     &                 PRID( 1:TRIMLEN( PRID ) )

                P = FINDC( PKEY, PS, PSLST )

C.................  Search for source in source list
                IF( P .GT. 0 ) THEN

                    S = FINDC( PSKEY( P ), NPSRC, SRCLST )

                    IF( S .LE. 0 ) THEN
                        EFLAG = .TRUE.
                        WRITE( MESG,94010 ) 
     &                    'Bad line', IREC, ' in file' //
     &                    '--record not found in source list'
                        CALL M3MESG( MESG )
                        GO TO 266   !  to head of read loop
                    ENDIF

                ELSE

                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 
     &                'Bad line', IREC, 'in file' //
     &                '--record not found in process file records'
                    CALL M3MESG( MESG )
                    GO TO 266       !  to head of read loop

                ENDIF

C.................  Record needed data for this source and hour
                    
                H = 0 
                DO 288 T = PTR, MIN( PTR + 23, NSTEPS )

                    H = H + 1
                    NPDPT( T ) = NPDPT( T ) + 1

                    HS = NPDPT( T )

                    IF( HS .LE. NPHRIP ) THEN

                        IDXSRC( HS,T ) = HS
                        SRCIDA( HS,T ) = S
                        CODEA ( HS,T ) = COD
                        PRCIDA( HS,T ) = P
                        EMISVA( HS,T ) = EMIS( H )

                    ELSEIF( HS .GT. MAXHS ) THEN
                        MAXHS = HS

                    ENDIF

288             CONTINUE

            GO TO 266  ! To head of pdayspec read loop

299         CONTINUE

            CLOSE( DDEV )

301     CONTINUE 

        IF( MAXHS .GT. NPHRIP ) THEN
            EFLAG = .TRUE.

            WRITE( MESG, 94010 )
     &        'Maximum number of records (NPHR*NIPOL) =', NPHRIP, 
     &        CRLF() // BLANK5 //
     &        'But number needed was', MAXHS, '. Must recompile.'

            CALL M3MSG2( MESG )

        ENDIF

        IF ( EFLAG ) THEN
            CALL M3EXIT( 'PDEMSPT', 0, 0,
     &                   'Error processing day-specific data', 2 )
        END IF

C.........  Reset number of time steps, because not enough data
        IF( MAXPTR .LT. NSTEPS ) THEN

            WRITE( MESG,94010 ) 'Not enough data to support ',
     &                          NSTEPS, 'time steps.' //
     &                          CRLF() // BLANK5, MAXPTR,
     &                          'time steps used instead.'

            CALL M3WARN( 'PDEMSPT', 0, 0, MESG )

            NSTEPS = MAXPTR
            
        ENDIF

C.........  Loop through pollutant codes and flag all codes that are present
        DO 307 V = 1, NIPOL
           OVFLAG( V ) = .FALSE.
307     CONTINUE

        DO 311 T = 1, NSTEPS
            DO 309 J = 1, NPHR

                V = CODEA( J,T ) 
                IF( V .GT. 0 ) OVFLAG( V ) = .TRUE.

309         CONTINUE
311     CONTINUE

C.........  Set up output variable names in appropriate order
        J = 0
        DO 315 V = 1, NIPOL
            IF( OVFLAG( V ) ) THEN
                J = J + 1
                OUTVAR( J ) = EINAM( V )
                INDXO ( V ) = J
            ELSE
                INDXO ( V ) = 0
            ENDIF
315     CONTINUE 
        NOUTV = J  

C.........  Set up and open output hour-specific point source file
        IF ( .NOT. DESC3( ENAME ) ) THEN
            CALL M3EXIT( 'PDEMSPT', 0, 0,
     &                   'Could not read file description for ' //
     &                   ENAME, 2 )
        ENDIF

        FTYPE3D = IDDATA3
        SDATE3D = SDATE
        STIME3D = STIME
        TSTEP3D = TSTEP
        NVARS3D = NOUTV
        NROWS3D = NPHR      !  number of rows = # of point sources.

        DO 322 V = 1, NOUTV
            VNAME3D( V ) = OUTVAR( V )
            VTYPE3D( V ) = M3REAL
            UNITS3D( V ) = 'tons/hour'
            VDESC3D( V ) = 'specific period emissions totals'
322     CONTINUE

        FDESC3D( 1 ) = 'Point source hour-specific from EMS-95 data'

        DO 325 I = 2, MXDESC3 
            FDESC3D( I ) = ' '
325     CONTINUE

        HNAME = PROMPTMFILE(
     &          'Enter logical name for HOUR-POINT output file',
     &          FSUNKN3, 'PHOUR', 'PDEMSPT' )

C........................................................................
C.............  Head of the hourly output loop  .........................
C........................................................................

        JDATE   = SDATE
        JTIME   = STIME
        VFLAG   = .FALSE.

        DO 401 T = 1, NSTEPS

            DO  333  V = 1, NOUTV
            DO  332  S = 1, NPHR 
                EMISH( S,V ) = BADVAL3
332         CONTINUE
333         CONTINUE

            NPRHR = 0  ! Initialize the number of process per hour

            CALL SORTI1( NPDPT( T ), IDXSRC( 1,T ), SRCIDA( 1,T ) )

C.............  Store sorted records for this hour
            LS = 0  ! previous source
            K  = 0
            DO 344 I = 1, NPDPT( T )

                J = IDXSRC( I,T )
                S = SRCIDA( J,T )
                P = PRCIDA( J,T )
                V = INDXO( CODEA ( J,T ) )

                IF( V .LE. 0 ) THEN
                    MESG = 'Internal error indexing input pollutants '//
     &                     'to output pollutants'
                    CALL M3EXIT( 'PDEMSPT', JDATE, JTIME, MESG, 2 )
                ENDIF

                IF( S .NE. LS ) THEN

                    K = K + 1
                    INDXH( K ) = S
                    LS         = S

                ENDIF

                L1 = TRIMLEN( PSLST( P ) )
                PVKEY = PSLST( P )( 1:L1 ) // EINAM( V )
                L = INDEX1( PVKEY, NPRHR, PRHRLST )

                IF( EMISH( K,V ) .LT. AMISS3 ) THEN  ! EMISH is blank

                    EMISH( K,V ) = EMISVA( J,T )
                    NPRHR = NPRHR + 1
                    PRHRLST( NPRHR ) = PVKEY

C.................  Check if this process/device has been used this hour
                ELSEIF( L .LE. 0 ) THEN

                    EMISH( K,V ) = EMISH( K,V ) + EMISVA( J,T )
                    NPRHR = NPRHR + 1
                    PRHRLST( NPRHR ) = PVKEY

                ELSE
                    L2 = TRIMLEN( SRCLST( S ) )
                    VFLAG = .TRUE.
                    WRITE( MESG,94010 )
     &                  'Duplicate emissions record:  FIP:' //
     &                  SRCLST( S )(  1:5  ) // ' SCC:' //
     &                  SRCLST( S )(  6:13 ) // ' Plant&Stack:' //
     &                  SRCLST( S )( 14:L2 ) // ' Species:' //
     &                  EINAM( V )
                    CALL M3MESG( MESG )
                END IF

344         CONTINUE

            NHRPT = K

            IF ( VFLAG ) THEN
 
                CALL M3EXIT( 'PDEMSPT', 0, 0,
     &                 'Duplicate records found in input files', 2 )
 
            ENDIF

 
            IF ( SFLAG ) THEN       !  if missing records are fatal:
 
                VFLAG = .FALSE.
 
                DO  366  V = 1, NOUTV
                DO  365  K = 1, NHRPT
 
                   IF ( EMISH( K,V ) .LT. AMISS3 ) THEN
 
                        S = INDXH( K )

                        VFLAG = .TRUE.
                        WRITE( MESG,94010 )
     &                      'Missing emissions record:  FIP:' //
     &                      SRCLST( S )(  1:5  ) // ' SCC:' //
     &                      SRCLST( S )(  6:13 ) // ' Plant&Stack:' //
     &                      SRCLST( S )( 14:L2 ) // ' Species:',
     &                      OUTVAR( V )
                        CALL M3MESG( MESG )
 
                    END IF
 
365             CONTINUE
366             CONTINUE
 
                IF ( VFLAG ) THEN
                    CALL M3EXIT( 'PDEMSPT', 0, 0,
     &                       'Missing records in input file', 2 )
                END IF
 
            ELSE            !  "missing" ~~~> 0
 
                DO  368  V = 1, NOUTV
                DO  367  S = 1, NPHR 
                    IF ( EMISH( S,V ) .LT. AMISS3 ) EMISH( S,V ) = 0.0
367             CONTINUE
368             CONTINUE
 
            END IF          !  sflag or not:  check for non-set EMIS values

C.............  Write emissions for this hour

            IF ( .NOT. WRITE3( HNAME, ALLVAR3, JDATE, JTIME,
     &                         NHRPT) ) THEN
                CALL M3EXIT( 'PDEMSPT', JDATE, JTIME,
     &                       'Error writing output file "' //
     &                       HNAME( 1:TRIMLEN( HNAME ) ) // '"',  2 )
            END IF

            CALL NEXTIME( JDATE, JTIME, TSTEP )

401     CONTINUE     !  End of loop over hours


999     CONTINUE          !  exit program:  normal completion

        MESG = 'Successful completion of Program PDEMSPT'
        CALL M3EXIT( 'PDEMSPT', 0, 0, MESG, 0 )


C.........  Error opening raw input file
1006    MESG = 'Error opening file ' // NAMTMP( 1:TRIMLEN( NAMTMP ) )
        CALL M3EXIT( 'PDEMSPT', 0, 0, MESG, 2 )


C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )


C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

93010   FORMAT( A5, 1X, A8, 15X, A15, 1X, A12 ) 

93020   FORMAT( I5.5 )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I7, :, 1X ) )

        END

