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

         PROGRAM PDPOINT

C***********************************************************************
C  program body starts at line 133
C
C  DESCRIPTION:
C       Construct Models3/EDSS day- or hour-specific point source file 
C       sorted by source (on FIP,SCC,PLANT,STACK) from data contained in
C       EPS/AFS-style day specific point source file.
C
C  PRECONDITIONS REQUIRED:
C       Day- or hour-specific input data for point sources concatenated in
C       chronological order.
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       Models-3 I/O
C       FIND1, PROMPTFFILE, PROMPTMFILE, GETYN
C
C  REVISION  HISTORY:
C       Prototype  10/95 by CJC.
C
C***********************************************************************

      IMPLICIT NONE

C...........   INCLUDES:
        
        INCLUDE 'PTDIMS3.EXT'   !  point-source dimensioning parameters
        INCLUDE 'CHDIMS3.EXT'   !  emis chem info (both inventory and model)
        INCLUDE 'GRDIMS3.EXT'   !  grid parameters
        INCLUDE 'TMDIMS3.EXT'   !  time related parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.


C...........   PARAMETER

        INTEGER         NSMAX
        PARAMETER     ( NSMAX  = NPDAY * NIPOL )


C...........   EXTERNAL FUNCTIONS and their descriptions:
        
        CHARACTER*24    DT2STR
        LOGICAL         ENVYN
        INTEGER         FIND1   !  returns -1 for failure
        INTEGER         GETNUM
        LOGICAL         GETYN
        INTEGER         JULIAN
        INTEGER         PROMPTFFILE
        CHARACTER*16    PROMPTMFILE
        INTEGER         SECSDIFF
        INTEGER         SEC2TIME
        INTEGER         TIME2SEC
        INTEGER         STR2INT
        REAL            STR2REAL
        
        EXTERNAL    DT2STR, ENVYN, FIND1, GETNUM, GETYN, JULIAN, 
     &              PROMPTFFILE, PROMPTMFILE, SECSDIFF, SEC2TIME,
     &              STR2INT, STR2REAL, TIME2SEC


C...........   LOCAL VARIABLES and their descriptions:
C...........   NOTE that SSC (point-source-category) ID's are 10-digit
C...........   unsigned integers which may be treated as a leading 7-digit
C...........   field, and a trailing 3-digit field.  *7 and *3 arrays below
C...........   follow this scheme with parallel arrays

        CHARACTER*16  ENAME	!  logical name for PNTS file from RAWPOINT
        CHARACTER*16  DNAME	!  " DAY-POINTS output file
        CHARACTER*16  UNITS	!  " DAY-POINTS units for variables

        INTEGER     FIP
        INTEGER     SCC
        INTEGER     SIC
        INTEGER     PLANT
        INTEGER     STACK
        INTEGER     SDT, EDT, STM, ETM
        INTEGER     MM, DD, YY, HH, SECS, STEP
        REAL        EMIS
        INTEGER     IPOL
        INTEGER     JDATE, JTIME, TSTEP, EDATE, ETIME
                                 
C...........   Day-specific point Sources Table input unsorted copy (*A); 

        INTEGER       INDEXA( NSMAX )  !  subscript table for SORTI4()
        INTEGER       IFIPA ( NSMAX )  !  source FIPS (county) ID
        INTEGER       ISCCA ( NSMAX )  !  source SCC
        INTEGER       IPLANA( NSMAX )  !  Plant ID
        INTEGER       ISTACA( NSMAX )  !  Stack ID
        INTEGER       ICODEA( NSMAX )  !  spcecies-subscript
        REAL          EMISDA( NSMAX )  !  emissions values (tons/yr)

C...........   Other local variables
        
        INTEGER         S, I, V !  loop counters.
        INTEGER         PDEV    !  for point-source file
        INTEGER         LDEV    !  log-device
        INTEGER         IOS     !  I/O status
        INTEGER         IREC    !  input line (record) number
        CHARACTER*256   MESG    !  text for M3EXIT()
        CHARACTER*240   LINE    !  input line from POINT file
        LOGICAL         VFLAG   !  verification flag


C***********************************************************************
C   begin body of program PDPOINT
        
        LDEV = INIT3()
        
        WRITE( *,92000 ) 
     &' ',
     &'Program PDPOINT to take a day- or hour-specific EPS/AFS point',
     &'source file and the sorted PNTS file from RAWPOINT, and produce',
     &'DAY-SPECIFIC or HOUR-SPECIFIC POINT SOURCE EMISSIONS files for',
     &'use in TMPPOINT.',
     &'The time-specific input file should consist of records of the',
     &'appropriate type, sorted by date and time.',
     &' ',
     &'You will need to enter the logical names for the input and',
     &'output files (and to have set them prior to program launch,',
     &'using "setenv <logicalname> <pathname>").',
     &' ',
     &'You may use END_OF-FILE (control-D) to quit the program',
     &'during logical-name entry. Default responses are given in',
     &'brackets [LIKE THIS] and can be accepted by hitting the',
     &'<RETURN> key.',
     &' '
        
        IF ( .NOT. GETYN( 'Continue with program?', .TRUE. ) ) THEN
            WRITE( *,92000 ) ' ', 'Ending program PDPOINT', ' '
            CALL EXIT( 0 )
        END IF


C.......   Get file name; open input raw point source file

        ENAME = PROMPTMFILE( 
     &          'Enter logical name for POINTS inventory file',
     &          FSREAD3, 'PNTS', 'PDPOINT' )

        PDEV = PROMPTFFILE( 
     &         'Enter logical name for RAW DAY-POINT file',
     &         .TRUE., .TRUE., 'PRDAY', 'PDPOINT' )


C.......   Read the point source files
C.......   Get starting date (default borrowed from input file:
        
        WRITE( *,92000 ) ' ', 'Reading POINT SOURCE input files...', ' '

        READ( PDEV, 93000, END=122, IOSTAT=IOS )LINE
        IF ( IOS .NE. 0 ) THEN
            WRITE( MESG, 94010 ) 
     &          'Error ', IOS, 
     &          'reading POINT INPUT file at line', IREC
            CALL M3EXIT( 'PDPOINT', 0, 0, 2 )
        ELSE
            REWIND( PDEV )
        END IF

        SDT   = STR2INT( LINE( 60:67 ) )
        JTIME = 10000 * MOD( SDT,100 ) 
        SDT   = SDT / 100
        DD    = MOD( SDT , 100 )
        MM    = MOD( SDT / 100, 100 )
        YY    = 1900000 + 1000 * ( SDT / 10000 )
        SDT   = YY + JULIAN( YY, MM, DD )   !  as Julian date  YYYYDDD
        JDATE = SDT

        EDT   = STR2INT( LINE( 69:76 ) )
        ETIME = 10000 * + MOD( EDT , 100 )
        EDT   = EDT / 100
        DD    = MOD( EDT , 100 )
        MM    = MOD( EDT / 100, 100 )
        YY    = 1900000 + 1000 * ( EDT / 10000 )
        EDATE = YY + JULIAN( YY, MM, DD )   !  as Julian date  YYYYDDD
        CALL NEXTIME( EDATE, ETIME, 10000 ) !  ETIME:  *through* hr HH

        TSTEP = SEC2TIME( SECSDIFF( JDATE, JTIME, EDATE, ETIME ) )
        
        JDATE = GETNUM( SDT, 9999366, SDT, 
     &          'Enter starting date YYYYDDD  for DAY-POINT file' )
        JTIME = GETNUM( 0, 240000, JTIME, 
     &          'Enter starting time  HHMMSS  for DAY-POINT file' )

        TSTEP = GETNUM( 0, 240000, TSTEP, 
     &          'Enter time step      HHMMSS  for DAY-POINT file' )
        STEP  = TIME2SEC( TSTEP )


C.......   Set up file description (borrowing much of it from ENAME);
C.......   get file name; open output day-point file

        IF ( STEP .EQ. SEC2HR ) THEN
            UNITS = 'tons/hour'
            DNAME = 'PHOUR'
        ELSE IF ( STEP .EQ. SEC2DAY ) THEN
            UNITS = 'tons/day'
            DNAME = 'PDAY'
        ELSE
            DNAME = 'PSPEC'
            HH = STEP / SEC2HR
            WRITE( UNITS,94050 ) 'tons/(hr*', HH, ')'
            WRITE( MESG, 94010 ) 
     &          'Nonstandard period length:', HH, 'hours'
            CALL M3WARN( 'PDPOINT', 0,0, MESG )
        END IF

        IF ( .NOT. DESC3( ENAME ) ) THEN
            CALL M3EXIT( 'PDPOINT', 0, 0, 
     &                   'Could not read file description for ' // 
     &                   ENAME, 2 )
        END IF

        FTYPE3D = IDDATA3
        SDATE3D = JDATE
        STIME3D = JTIME
        TSTEP3D = TSTEP
        NVARS3D = NIPOL
        NROWS3D = NPDAY     !  number of rows = # of point sources.

        FDESC3D( 1 ) = 'NC Time-Specific Point Source emissions values'
        FDESC3D( 2 ) = 'Inventory  type "' // LINE(  1: 1 ) // '"'
        FDESC3D( 3 ) = 'Projection year "' // LINE(  2: 3 ) // '"'
        FDESC3D( 4 ) = 'Base year       "' // LINE(  5: 6 ) // '"'
        FDESC3D( 5 ) = 'Emission   type "' // LINE(  7: 8 ) // '"'

        DO  211  I = 6, MXDESC3
            FDESC3D( I ) = ' '
211     CONTINUE

        DO  222  V = 1, NIPOL
            VNAME3D( V ) = 'DAY' // EINAM( V )
            VTYPE3D( V ) = M3REAL
            UNITS3D( V ) = UNITS
            VDESC3D( V ) = 'specific period emissions totals'
222     CONTINUE

        DNAME = PROMPTMFILE( 
     &          'Enter logical name for DAY-POINT output file',
     &          FSUNKN3, DNAME, 'PDPOINT' )


C...........   Read the input file:

        VFLAG = .FALSE.
        IREC  = 0
        S     = 0

111     CONTINUE        !  head of the PDEV-read loop

            READ( PDEV, 93000, END=122, IOSTAT=IOS )LINE
            IREC = IREC + 1

            IF ( IOS .NE. 0 ) THEN

                WRITE( MESG, 94010 ) 
     &              'Error ', IOS, 
     &              'reading DAY-POINT file at line', IREC
                CALL M3EXIT( 'PDPOINT', JDATE, JTIME, 2 )

            END IF

            IF ( LINE( 57:57 ) .NE. 'S' ) THEN   !  special interval

                WRITE( MESG, 94010 ) 
     &              'Bad record type "' // LINE( 57:58 ) //
     &              '" in  DAY-POINT file at line', IREC
                CALL M3MESG( MESG )
                VFLAG = .TRUE.
                GO TO 111

            END IF

            SDT = STR2INT( LINE( 60:67 ) )
            EDT = STR2INT( LINE( 69:76 ) )
            HH  = MOD( SDT        , 100 )
            STM = 10000 * HH
            DD  = MOD( SDT / 100  , 100 )
            MM  = MOD( SDT / 10000, 100 )
            YY  = 1900000 + 1000 * ( SDT / 1000000 )
            SDT = YY + JULIAN( YY, MM, DD )   !  as Julian date  YYYYDDD

            IF ( SDT .LT. JDATE ) THEN

                WRITE( MESG, 94010 ) 
     &              'Records out of order:  from', JDATE, JTIME,
     &              'to', SDT, STM,
     &              'at line', IREC
                CALL M3MESG( MESG )
                VFLAG = .TRUE.
                GO TO 111

            ELSE IF ( SDT .GT. JDATE ) THEN       ! construct and write time step

                IF ( S .GT. NSMAX ) THEN
                    WRITE( MESG, 94020 ) 
     &                  'Number of sources', S, 
     &                  'Exceeds MAX', NSMAX,
     &                  'for', JDATE, JTIME,
     &                  'at line', IREC
                    CALL M3MESG( MESG )
                    S     = NSMAX
                    VFLAG = .TRUE.
                END IF

                CALL DAYREC( S, INDEXA, 
     &                       IFIPA, ISCCA, IPLANA, ISTACA, 
     &                       ICODEA, EMISDA, JDATE, JTIME,
     &                       ENAME, DNAME )

                SECS = SECSDIFF( JDATE, JTIME, SDT, STM )
                IF ( SECS .NE. STEP ) THEN
                    WRITE( MESG, 94040 ) 
     &                  'Bad time step:  from', JDATE, JTIME,
     &                  'to', SDT, STM,
     &                  'at line', IREC
                    CALL M3MESG( MESG )
                    VFLAG = .TRUE.
                END IF

                JDATE = SDT
                JTIME = STM
                S     = 0
                       
            END IF

            IPOL  = STR2INT ( LINE( 156:160 ) )
            IF ( IPOL .EQ. 42602 ) IPOL = 42603         !  fix NO2 ~~> NOx
            I = FIND1( IPOL, NIPOL, PCODE )
            IF ( I  .LT. 0 ) THEN

                WRITE( MESG,94010 ) 
     &              'Bad line', IREC, 
     &              'Pollutant code "' // LINE( 156:160 ) // '"'
                CALL  M3MESG( MESG )
                VFLAG = .TRUE.
                GO TO  111      !  to head of loop

            END IF

            EMIS  = STR2REAL( LINE( 162:171 ) )
            IF ( EMIS .LT. 0.0 )  THEN

                WRITE( MESG,94010 ) 
     &              'Bad line', IREC, 
     &              'Emis value "' // LINE( 162:171 ) // '"'
                CALL M3MESG( MESG )
                VFLAG = .TRUE.
                GO TO  111

            END IF

            HH  = MOD( EDT        , 100 )
            ETM = 10000 * HH
            DD  = MOD( EDT / 100  , 100 )
            MM  = MOD( EDT / 10000, 100 )
            YY  = 1900000 + 1000 * ( EDT / 1000000 )
            EDT = YY + JULIAN( YY, MM, DD )      !  as Julian date YYYYDDD
            CALL NEXTIME( EDT, ETM, 10000 )      ! ETM means *through* hr HH

            SECS  = SECSDIFF( SDT, STM, EDT, ETM )
            IF ( SECS .NE. STEP ) THEN

                WRITE( MESG, 94010 ) 
     &              'Bad duration ', SECS,
     &              'in  DAY-POINT file at line', IREC
                CALL M3MESG( MESG )
                VFLAG = .TRUE.
                GO TO 111

            END IF

            FIP   = STR2INT ( LINE(  12: 16 ) )
            SIC   = STR2INT ( LINE(  24: 27 ) )
            SCC   = STR2INT ( LINE(  29: 38 ) )
            PLANT = STR2INT ( LINE(  40: 44 ) )
            STACK = STR2INT ( LINE(  46: 48 ) )

            S = S + 1

            IF ( S .LE. NSMAX ) THEN

                INDEXA( S ) = S         !  index-table for use in SORTI4()
                IFIPA ( S ) = FIP
                ISCCA ( S ) = SCC
                IPLANA( S ) = PLANT
                ISTACA( S ) = STACK
                ICODEA( S ) = I         !  subscript into ACODE(*)
                EMISDA( S ) = EMIS

            ELSE

                WRITE( MESG, 94020 ) 
     &              'Number of sources', S, 
     &              'Exceeds MAX', NSMAX,
     &              'for', JDATE, JTIME, 
     &              'at line', IREC
                CALL M3MESG( MESG )
                VFLAG = .TRUE.

            END IF

            GO TO  111          !  to head of PDEV-read loop


122     CONTINUE        !  end of the PDEV-read loop
        
        IF ( S .GT. 0 ) THEN    !  if cleanup necessary
        
            IF ( S .GT. NSMAX ) THEN

                WRITE( MESG, 94020 )
     &              'Number of sources', S, 
     &              'Exceeds MAX', NSMAX,
     &              'for', JDATE, JTIME, 
     &              'at line', IREC
                CALL M3MESG( MESG )
                S     = NSMAX
                VFLAG = .TRUE.

            END IF

            CALL DAYREC( S, INDEXA, 
     &                   IFIPA, ISCCA, IPLANA, ISTACA, 
     &                   ICODEA, EMISDA, JDATE, JTIME,
     &                   ENAME, DNAME )

        END IF                  !  if s>0 (cleanup necessary)


999   CONTINUE          !  exit program:  normal completion

      IF ( VFLAG ) THEN
          MESG = 'Errors encountered in Program PDPOINT'
          CALL M3EXIT( 'PDPOINT', 0, 0, MESG, 2 )
      ELSE
          MESG = 'Successful completion of Program PDPOINT'
          CALL M3EXIT( 'PDPOINT', 0, 0, MESG, 0 )
      END IF


C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )

92010   FORMAT( 5X, A, :, I10 )


C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

93010   FORMAT( A16 )


C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I7, :, 1X ) )

94020   FORMAT( A, I5, X, A, I5, X, A, I9, ':', I6.6, :, X, A, I5  )

94030   FORMAT( 'Emis set to zero for source ', 
     &           I5, ':',I7, ':', I4, ':' I2, 2X, A, '--', A )

94040   FORMAT( 2( A, I9, ':', I6.6, 2X ) A, I5 )

94050   FORMAT( A, I5.5, A )


        END

