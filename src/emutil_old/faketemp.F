
        PROGRAM FAKETEMP

C***********************************************************************
C  subroutine body starts at line 138
C
C  FUNCTION:
C       Create an altered or fabricated gridded temperature data file using
C       one of a variety of methods, including: (1) uniform, (2) by state, 
C       (3) cell to value transform of existing data, (4) cell to cell 
C       transform of existing data. 
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C       Prototype 1/97 by MRH
C
C***************************************************************************
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

        INCLUDE 'GRDIMS3.EXT'    
        INCLUDE 'PARMS3.EXT'    
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations

C...........   EXTERNAL FUNCTIONS and their descriptions:

        LOGICAL         DSCGRID
        INTEGER         GETDATE
        INTEGER         GETMENU
        INTEGER         GETNUM
        REAL*4          GETREAL
        LOGICAL         GETYN
        INTEGER         INDEX1
        INTEGER         LEN2
        INTEGER         PROMPTFFILE
        CHARACTER*16    PROMPTMFILE
        INTEGER         TRIMLEN

        EXTERNAL        DSCGRID, GETDATE, GETMENU, GETNUM, GETREAL, 
     &                  GETYN, INDEX1, LEN2, PROMPTFFILE, PROMPTMFILE, 
     &                  TRIMLEN

C...........   Temperature arrays
        REAL       TA( NGRID ) ! storage for temperatures

C...........   Methods of faking the temperature file
        INTEGER         NMETHOD         ! Number of fake temp methods
        PARAMETER     ( NMETHOD = 4 ) 
        CHARACTER*40    MDESC( NMETHOD )

        DATA   MDESC / 'Uniform', 'Cell to value creation of new data',
     &                 'Cell to value transform of existing data', 
     &                 'Cell to cell  transform of existing data'  /

C...........   Variables for uniform method
        INTEGER    TPROMPT        ! Indicator for temporal resolution 

        CHARACTER*72    TRESOLUT( 3 )

        DATA  TRESOLUT / 'Per hour', 'Per Day', 'Constant' /

C...........   Variables for cell to value methods
        INTEGER     MXTVAL
        PARAMETER ( MXTVAL = 26 )  ! Number of letter in alphabet

        INTEGER     NTVAL
        REAL        TVAL ( MXTVAL )
        CHARACTER*1 TCLAS( MXTVAL )
        CHARACTER*1 TCGRD( NCOLS, NROWS )

C...........   Variables for cell to cell conversion method
        INTEGER    NCEL           ! Number of cell to cell conversions
        INTEGER    CINDEX1( NGRID )  ! Destination cell numbers
        INTEGER    CINDEX2( NGRID )  ! "From"      cell numbers

C...........   Unit numbers and file names
        INTEGER       CDEV        ! Cell to cell conversion file number
        INTEGER       LDEV        ! Log file unit number
        INTEGER       VDEV        ! Cell to value file number

        CHARACTER*16  ONAME       ! Output temperature file
        CHARACTER*16  TNAME       ! Input  temperature file

C...........   LOCAL VARIABLES and their descriptions:
        REAL       TMPR           ! Temorary temperature [K]

        INTEGER    I, J, N, T     ! Counters and Pointers
        INTEGER    IOS            ! I/O status
        INTEGER    IREC           ! Record number count
        INTEGER    JDATE          ! Current date
        INTEGER    JTIME          ! Current time
        INTEGER    METHOD         ! Temperature faking method
        INTEGER    NDAYS          ! Number of fake days
        INTEGER    NSTEP          ! Number of time steps
        INTEGER    SDATE          ! Starting date 
        INTEGER    STIME          ! Starting time
        INTEGER    TSTEP          ! Time step
        INTEGER    X1, X2, Y1, Y2 ! Temporary cell locations

        LOGICAL    EFLAG          ! Error flag
        LOGICAL    INFILE
        LOGICAL    CEL2CEL
        LOGICAL    CEL2VAL

        CHARACTER*16    SCRBUF          ! Scratch buffer
        CHARACTER*256   MESG 

C***********************************************************************
C   begin body of program FAKETEMP

        LDEV = INIT3()

        CALL INITEM( LDEV )
 
        WRITE( *,92000 )
     &  ' ',
     &  'Program FAKETEMP alters or fabricates a surface level NetCDF',
     &  'temperature file.  Several methods can be used to create this',
     &  'file:',
     &  '   (1) Uniform >> set a spatially uniform temperature for',
     &  '           each hour of a specific time period',
     &  '   (2) By state >> set a temporally uniform temperature for',
     &  '           each state in the domain',
     &  '   (3) Cell to value transform >> set a temporally uniform',
     &  '           temperature to groups of cells',
     &  '   (4) Cell to cell transform >> for each hour, set the',
     &  '           temperature of specific cells to the temperature',
     &  '           of other cells',
     &  ' ',
     &  'The program prompts for different inputs depending on the',
     &  'method selected.  The following data may be required: ',
     &  'existing lower-level temperature data, gridding surrogates',
     &  'and cross-reference files, cell to value transform data file,',
     &  'and/or a cell to cell transform data file.',
     &  ' ',
     &  'You may use END_OF-FILE (control-D) to quit the program',
     &  'during logical-name entry. Default responses are given in',
     &  'brackets [LIKE THIS] and can be accepted by hitting the',
     &  '<RETURN> key.',
     &  ' '

        IF ( .NOT. GETYN( 'Continue with program?', .TRUE. ) ) THEN
            WRITE( *,92000 ) ' ', 'Ending program FAKETEMP', ' '
            CALL EXIT( 0 )
        END IF

C.........  Set temperature fabrication method

   11   CONTINUE   
        METHOD = GETMENU( NMETHOD, 1, 'Which method to use?', MDESC )

C.........  Initialize program control variables
        TPROMPT = 0
        INFILE  = .FALSE.
        CEL2CEL = .FALSE.
        CEL2VAL = .FALSE.

C.........  Set program control variables based on method selected 

        IF( METHOD .EQ. 1 ) THEN 
            TPROMPT = GETMENU( 3, 1,
     &                'Temporal resolution for uniform temperatures?',
     &                TRESOLUT )

        ELSEIF( METHOD .EQ. 2 ) THEN
            CEL2VAL = .TRUE.

        ELSEIF( METHOD .EQ. 3 ) THEN
            INFILE  = .TRUE.
            CEL2VAL = .TRUE.

        ELSEIF( METHOD .EQ. 4 ) THEN
            INFILE  = .TRUE.
            CEL2CEL = .TRUE.

        ELSE
            WRITE( *,92010 ) 'Method', METHOD, 
     &                       'not currently supported.'
            GO TO 11

        ENDIF

C.........  Prompt for additional files, if needed...

C.........  Prompt for input temperature file
        IF( INFILE ) THEN
            TNAME = PROMPTMFILE(
     &              'Enter logical name for INPUT TEMPERATURE file',
     &              FSREAD3, 'TKDAT01', 'FAKETEMP' )

        ENDIF

C.........  Set output file header information

        IF( INFILE ) THEN

            IF ( .NOT. DESC3( TNAME ) ) THEN
                CALL M3EXIT( 'FAKETEMP', 0, 0,
     &          'Error reading HEADER from TEMPERATURE file', 2 )
            ENDIF

            NVARS3D = 1

            FDESC3D( 1 ) = 'Fake 1-layer temperature from "' //
     &                     MDESC( METHOD ) // '" method'

            VNAME3D( 1 ) = 'TA'
            VTYPE3D( 1 ) = M3REAL
            UNITS3D( 1 ) = 'deg K' 
            VDESC3D( 1 ) = 'Temperature'

            SDATE = SDATE3D
            STIME = STIME3D
            TSTEP = TSTEP3D
            NSTEP = MXREC3D

        ELSE
            NDAYS = 1
            NDAYS = GETNUM( 1, 9999, 1,
     &                      'Enter the number of days to fake' )

            SDATE = 1993201
            SDATE = GETDATE( SDATE,
     &                       'Enter start date (YYYYDDD)|(YYYYMMDD)' )
            STIME = 0 
            STIME = GETNUM( 0, 235959, STIME, 
     &                       'Enter start time  (HHMMSS)')
            TSTEP = 10000
            NSTEP = NDAYS * 24

C.............  Get grid header variables
            IF ( .NOT. DSCGRID( GRDNM, SCRBUF, GDTYP3D,
     &                  P_ALP3D, P_BET3D, P_GAM3D, XCENT3D, YCENT3D,
     &                  XORIG3D, YORIG3D, XCELL3D, YCELL3D,
     &                  NCOLS3D, NROWS3D, NTHIK3D ) ) THEN
 
                SCRBUF = GRDNM
                MESG   = '"' // SCRBUF( 1:TRIMLEN( SCRBUF ) ) //
     &                   '" not found in GRIDDESC file'
                CALL M3EXIT( 'FAKETEMP', 0, 0, MESG, 2 )
 
            END IF                  ! if dscgrid() failed

            IF( NROWS3D .NE. NROWS .OR. 
     &          NCOLS3D .NE. NCOLS      ) THEN

               WRITE( MESG,94010 ) 'In GRDIMS3.EXT, GRDNM ' // GRDNM //
     &                ' inconsistent w/ NCOLS=', NCOLS, 
     &                '& NROWS=', NROWS

            ENDIF

C.............  Set up remaining header variables
            FTYPE3D = GRDDED3
            SDATE3D = SDATE 
            STIME3D = STIME
            TSTEP3D = TSTEP
            NVARS3D = 1
            NLAYS3D = 1
            VGTYP3D = VGZVAL3
            VGTOP3D = AMISS3
            GDNAM3D = GRDNM

            FDESC3D( 1 ) = 'Fake 1-layer temperature from "' //
     &                     MDESC( METHOD ) // '" method'

            VNAME3D( 1 ) = 'TA'
            VTYPE3D( 1 ) = M3REAL
            UNITS3D( 1 ) = 'deg K' 
            VDESC3D( 1 ) = 'Temperature'

        ENDIF

C.........  Prompt for, read, and process cell to value file, which
C           can be used for either creating or modifying, depending
C           on status in INFILE

        IF( CEL2VAL ) THEN

            VDEV = PROMPTFFILE(
     &           'Enter name for CELL TO VALUE input file',
     &           .TRUE., .TRUE., 'TCTOV', 'FAKETEMP' )

            IF ( VDEV .LE. 0 ) THEN
                CALL M3EXIT( 'FAKETEMP', 0, 0,
     &             'Could not open CELL TO VALUE input file', 2 )
            END IF

            EFLAG = .FALSE.

C.............  Read past VDEV header 
            READ( VDEV, *, END=999 )
            READ( VDEV, *, END=999 )
            READ( VDEV, *, END=999 )

C.............  Read number of different temperatures used
            READ( VDEV, *, END=999 ) NTVAL

C.............  Ensure number of temperatures is not > max
            IF ( NTVAL .GT. MXTVAL ) THEN
                WRITE( MESG,94010 ) 
     &                'Number of temperature classes > MXTVAL=', MXTVAL
                CALL M3EXIT( 'FAKETEMP', 0, 0, MESG, 2 )
            ENDIF

C.............  Read classes (single chars) and temperature values (real)
C.............  A TVAL >= 9999. means do not set or adjust temperature
            IREC  = 4
            DO 15 I = 1, NTVAL

                READ( VDEV,*,END=999,IOSTAT=IOS ) TCLAS( I ), TVAL( I )

                IREC = IREC + 1
                IF( IOS .NE. 0 ) THEN
                    WRITE( MESG,94010 )
     &                  'I/O error', IOS,
     &                  'reading CELL TO VALUE file at line', IREC
                    CALL M3EXIT( 'FAKETEMP', 0, 0, MESG, 2 )
                ENDIF

   15       CONTINUE

C.............  Read grid of temperature classes
            DO  17 J = NROWS, 1, -1

                READ( VDEV, '(1000A1)', END=999, IOSTAT=IOS )
     &                        ( TCGRD( I,J ), I = 1,NCOLS )

                IREC = IREC + 1

                IF( IOS .NE. 0 ) THEN
                    WRITE( MESG,94010 )
     &                  'I/O error', IOS,
     &                  'reading CELL TO VALUE file at line', IREC
                    CALL M3EXIT( 'FAKETEMP', 0, 0, MESG, 2 )
                ENDIF

   17       CONTINUE  !  End of VDEV read 

        ENDIF

C.........  Prompt for, read, and process cell to cell conversion file
        IF( CEL2CEL ) THEN

            CDEV = PROMPTFFILE(
     &           'Enter name for CELL TO CELL CONVERSION file',
     &           .TRUE., .TRUE., 'TCTOC', 'FAKETEMP' )

            IF ( CDEV .LE. 0 ) THEN
                CALL M3EXIT( 'FAKETEMP', 0, 0,
     &             'Could not open CELL TO CELL CONVERSION file', 2 )
            END IF

            EFLAG = .FALSE.

C.............  Read past CDEV header 
            READ( CDEV, *, END=24 )
            READ( CDEV, *, END=24 )
            READ( CDEV, *, END=24 )

            I     = 0
            IREC  = 3
   22       CONTINUE  !  Head of CDEV read loop

                READ( CDEV, *, END=24, IOSTAT=IOS ) X1, Y1, X2, Y2
            
                IREC = IREC + 1

                IF( IOS .NE. 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 )
     &                  'I/O error', IOS,
     &                  'reading CELL TO CELL file at line', IREC
                    CALL M3MESG( MESG )
                    GO TO 22

                ELSEIF( X1 .GT. NCOLS .OR. 
     &                  Y1 .GT. NROWS      ) THEN
                    EFLAG = .TRUE.

                    SCRBUF = GRDNM
                    WRITE( MESG,94010 )
     &               'Destination cell location (', X1, ',', Y1,
     &               ') outside grid '// SCRBUF( 1:TRIMLEN( SCRBUF ) )//
     &               ' at line', IREC
                    CALL M3MESG( MESG )
                    GO TO 22

                ELSEIF( X2 .GT. NCOLS .OR.
     &                  Y2 .GT. NROWS      ) THEN
                    EFLAG = .TRUE.

                    SCRBUF = GRDNM
                    WRITE( MESG,94010 )
     &               '"From" cell location (', X2, ',', Y2,
     &               ') outside grid '// SCRBUF( 1:TRIMLEN( SCRBUF ) )//
     &               ' at line', IREC
                    CALL M3MESG( MESG )
                    GO TO 22

                ELSE
                    I = I + 1

                    CINDEX1( I ) = ( Y1-1 ) * NCOLS + X1
                    CINDEX2( I ) = ( Y2-1 ) * NCOLS + X2

                ENDIF

            GO TO 22  !  End of CDEV read loop

   24       CONTINUE

            IF( I .EQ. 0 ) THEN
                EFLAG = .TRUE.
                CALL M3MESG( 'No valid cell to cell entries in TCTOC' )

            ELSE
                NCEL = I

            ENDIF

            IF( EFLAG ) THEN
                MESG = 'Error reading CELL TO CELL CONVERSION file'
                CALL M3EXIT( 'FAKETEMP', 0, 0, MESG, 2 )
            ENDIF

        ENDIF

C.........  Get output file name using environment variable

        ONAME = PROMPTMFILE(
     &          'Enter logical name for output file',
     &          FSUNKN3, 'OUTPUT', 'FAKETEMP' )

C.........  Prompt for uniform, constant temperature
        IF( TPROMPT .EQ. 3 ) THEN
            WRITE( MESG,94010 ) 'Enter temperature value for output'
            TMPR = GETREAL( 0., 1000., 293., MESG )

        ENDIF

C.........  Read and write, but dimension on the fly

        JDATE = SDATE
        JTIME = STIME
        DO 111 T = 1, NSTEP

C.............  Intialize temperatures for current time step
            IF( INFILE ) THEN

                IF( .NOT. READ3( TNAME, 'TA', 1,
     &                           JDATE, JTIME, TA ) ) THEN

                    MESG = 'Error reading "' // 
     &                     SCRBUF( 1:TRIMLEN( SCRBUF ) ) //
     &                     '" from input TEMPERATURE file'

                    CALL M3EXIT( 'FAKETEMP', JDATE, JTIME, MESG, 2 )

                END IF
            ENDIF

C.............  Get day-specific temperature
            IF( TPROMPT .EQ. 2 .AND. MOD( T+23,24 ) .EQ. 0 ) THEN
                WRITE( MESG,94010 ) 
     &                 'Enter temperature value for date:', JDATE
                TMPR = GETREAL( 0., 1000., 293., MESG )

C.............  Get hour-specific temperature
            ELSEIF( TPROMPT .EQ. 1 ) THEN
                WRITE( MESG,94010 ) 
     &                 'Enter temperature value for date:', JDATE,
     &                 'and time:', JTIME
                TMPR = GETREAL( 0., 1000., 293., MESG )

            ENDIF

C.............  Populate temperature array
            IF( TPROMPT .GT. 0 ) THEN
                DO 44 I = 1, NGRID
                    TA( I ) = TMPR
   44           CONTINUE

            ELSEIF( CEL2VAL ) THEN
                
                DO 55 J = 1, NROWS
                    DO 54 I = 1, NCOLS
                        N = ( J-1 )*NCOLS + I

                        TMPR = TVAL( INDEX1( 
     &                               TCGRD( I,J ), NTVAL, TCLAS ) )

                        IF( TMPR .LT. 9999. ) THEN
                            TA( N ) = TMPR
                        ELSEIF( .NOT. INFILE ) THEN
                            TA( N ) = 0.
                        ENDIF 

   54               CONTINUE
   55           CONTINUE

            ELSEIF( CEL2CEL ) THEN

                DO 66 I = 1, NCEL
                     TA( CINDEX1( I ) ) = TA( CINDEX2( I ) )
   66           CONTINUE

            ENDIF

C.............  Write temperature array
            IF( .NOT. 
     &          WRITE3( ONAME, 'TA', JDATE, JTIME, TA(1) ) ) THEN

                WRITE( MESG, 94010 ) 'Could not write timestep to "'
     &                 // ONAME( 1:TRIMLEN( ONAME ) ) // '".' 
                CALL M3EXIT( 'FAKETEMP', JDATE, JTIME, MESG, 2 )

            ENDIF
                
            CALL NEXTIME( JDATE, JTIME, TSTEP )

  111   CONTINUE

        CALL M3EXIT( 'FAKETEMP', 0, 0,
     &               'Normal completion of PROGRAM FAKETEMP', 0 )

C.........  Abort if too early an end to TCTOV file
  999   WRITE( MESG,94010 ) 
     &        'Unexpected end to TCTOV at line', IREC

        CALL M3EXIT( 'FAKETEMP', 0, 0, MESG, 2 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx
 
92000   FORMAT( 5X, A )
92010   FORMAT( / 10( 5X, A, :, I8, :, 2X) )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 2X) )

        END


