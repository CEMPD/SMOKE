
C Version "@(#)$Id$ $Source$ $Date$ 

        PROGRAM ELEVPOINT

C***********************************************************************
C  program body starts at line  75
C
C  DESCRIPTION:
C       Uses PLUMRIS() and user-selected CUTOFF to determine which point
C       sources may possibly generate elevated plumes.  Outputs elevated
C       plume source FIP,SCC,PLANT,STACK,RISE to file.
C
C  PRECONDITIONS REQUIRED:
C       RAWPOINT point sources output file
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       PLUMRIS()
C
C  REVISION  HISTORY:
C       Prototype 10/95 by CJC
C       Version   12/95 by CJC uses stack parameters file to fix up missing
C       stack parameter values.
C***********************************************************************

      IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'PTDIMS3.EXT'   !  point-source dimensioning parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API constants
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations


C...........   EXTERNAL FUNCTIONS and their descriptions:

        REAL            GETREAL
        LOGICAL         GETYN
        REAL            PLUMRIS
        INTEGER         PROMPTFFILE
        CHARACTER*16    PROMPTMFILE

        EXTERNAL        GETREAL, GETYN, PLUMRIS, 
     &                  PROMPTFFILE, PROMPTMFILE


C...........   LOCAL VARIABLES and their descriptions:
C...........   Point Sources Table input sorted by FIP,SCC,PLANT,STACK

        INTEGER     IFIP  ( NPSRC )  !  source FIPS (county) ID
        INTEGER     ISCC  ( NPSRC )  !  source SCC
        INTEGER     IPLANT( NPSRC )  !  Plant ID
        INTEGER     ISTACK( NPSRC )  !  Stack ID
        REAL        STKHT ( NPSRC )  !  stack height   (m)
        REAL        STKDM ( NPSRC )  !  stack diameter (m)
        REAL        STKTK ( NPSRC )  !  exhaust temperature (deg K)
        REAL        STKVE ( NPSRC )  !  exhaust velocity    (m/s)

C...........   Other local variables

        INTEGER         S       !  loop counter.
        INTEGER         N       !  elevated sources counter.
        CHARACTER*16    ENAME   !  logical name for point-source input file
        INTEGER         FDEV    !  unit number for stack parameters file
        INTEGER         LDEV    !  unit number for log file
        INTEGER         PDEV    !  ... output elevated point source file
        REAL            CUTOFF  !  plume rise cutoff for elev pts
        REAL            RISE    !  calculated plume rise
        CHARACTER*60    MESG


C***********************************************************************
C   begin body of program ELEVPOINT

        LDEV = INIT3()

        WRITE( *,92000 ) 
     &  ' ',
     &  'Program ELEVPOINT to take the sorted POINT SOURCE EMISSIONS',
     &  'file, and produce the ASCII ELEVATED POINT SOURCE file. ',
     &  ' ',
     &  'You will need to enter the logical names for the input and',
     &  'output files (and to have set them prior to program launch,',
     &  'using "setenv <logicalname> <pathname>").  You will be',
     &  'prompted for the cutoff value to be used to determine which ',
     &  'sources may give rise to elevated plumes.',
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


C.......   Get file name; open input point source and output
C.......   elevated points files; get plume-rise cutoff for
C.......   elevated points file

        ENAME = PROMPTMFILE( 
     &          'Enter name for POINT SOURCE input file',
     &          FSREAD3, 'PNTS', 'TMPPOINT' )

        IF ( ENAME( 1:5 ) .EQ. 'NNONE ' ) THEN
            CALL M3EXIT( 'ELEVPOINT', 0, 0, 
     &           'Could not open POINT SOURCE input file', 2 )
        END IF

        FDEV = PROMPTFFILE( 'Enter name for STACK PARAMETER input file',
     &                      .TRUE., .TRUE., 'PSTK', 'ELEVPOINT' )
        IF ( FDEV .LE. 0 ) THEN
            CALL M3EXIT( 'ELEVPOINT', 0, 0, 
     &           'Could not open input STACK PARAMETER FILE', 2 )
        END IF

        PDEV = PROMPTFFILE(  
     &        'Enter name for ELEVATED POINT SOURCE output file',
     &        .FALSE., .TRUE., 'PELV', 'ELEVPOINT' )
        IF ( PDEV .LE. 0 ) THEN
            CALL M3EXIT( 'ELEVPOINT', 0, 0, 
     &           'Could not open ELEVATED POINT SOURCE file', 2 )
        END IF

        CUTOFF = GETREAL( 0.0, -BADVAL3, 75.0, 
     &           'Enter CUTOFF MIN for elevated plume rise (meters)' )


C...........   Read in emissions indexing data:  IFIP, ASC7, ASC3, 
C...........   IPLANT, ISTACK, TZONES

        WRITE( *,92000 ) 
     &  ' ',
     &  'Reading POINT SOURCE input file ' // ENAME,
     &  ' '

        IF ( .NOT. READ3( ENAME, 'IFIP', ALLAYS3,0,0, IFIP ) ) THEN
            CALL M3EXIT( 'ELEVPOINT', 0, 0, 
     &                   'Error reading FIP from file' // ENAME, 2)
        END IF

        IF ( .NOT. READ3( ENAME, 'ISCC', ALLAYS3,0,0, ISCC ) ) THEN
            CALL M3EXIT( 'ELEVPOINT', 0, 0, 
     &                   'Error reading ISCC from file' // ENAME, 2 )
        END IF

        IF ( .NOT. READ3( ENAME, 'IPLANT', ALLAYS3,0,0, IPLANT ) ) THEN
            CALL M3EXIT( 'ELEVPOINT', 0, 0, 
     &                   'Error reading IPLANT from file' // ENAME, 2 )
        END IF

        IF ( .NOT. READ3( ENAME, 'ISTACK', ALLAYS3,0,0, ISTACK ) ) THEN
            CALL M3EXIT( 'ELEVPOINT', 0, 0, 
     &                   'Error reading ISTACK from file' // ENAME, 2 )
        END IF

        IF ( .NOT. READ3( ENAME, 'STKHT', ALLAYS3,0,0, STKHT ) ) THEN
            CALL M3EXIT( 'ELEVPOINT', 0, 0, 
     &                   'Error reading STKHT from file' // ENAME, 2 )
        END IF

        IF ( .NOT. READ3( ENAME, 'STKDM', ALLAYS3,0,0, STKDM ) ) THEN
            CALL M3EXIT( 'ELEVPOINT', 0, 0, 
     &                   'Error reading STKDM from file' // ENAME, 2 )
        END IF

        IF ( .NOT. READ3( ENAME, 'STKTK', ALLAYS3,0,0, STKTK ) ) THEN
            CALL M3EXIT( 'ELEVPOINT', 0, 0, 
     &                   'Error reading STKTK from file' // ENAME, 2 )
        END IF

        IF ( .NOT. READ3( ENAME, 'STKVE', ALLAYS3,0,0, STKVE ) ) THEN
            CALL M3EXIT( 'ELEVPOINT', 0, 0, 
     &                   'Error reading STKVE from file' // ENAME, 2 )
        END IF


C...........   Fix up missing stack parameters:

        CALL FIXSTK( FDEV, IFIP, ISCC, STKHT, STKDM, STKTK, STKVE )


C...........   Process the stacks to determine elevated sources:

        WRITE( *,92000 ) 
     &  ' ',
     &  'Writing ELEVATED POINT SOURCE output file',
     &  ' '

        N = 0
        DO  11  S = 1, NPSRC

            RISE = PLUMRIS( STKHT( S ), 
     &                      STKTK( S ), 
     &                      STKVE( S ), 
     &                      STKDM( S ) )

            IF ( RISE .GT. CUTOFF ) THEN

                N = N + 1
                WRITE( PDEV,93010 )
     &                  IFIP  ( S ),
     &                  ISCC  ( S ),
     &                  IPLANT( S ),
     &                  ISTACK( S ),
     &                  RISE

            END IF      !  if rise > cutoff

11      CONTINUE        !  end loop on sources S

        WRITE( MESG,94010 ) 'Number of elevated sources:', N
        CALL M3MESG( MESG )

        MESG = 'Successful completion of program ELEVPOINT'
        CALL M3EXIT( 'ELEVPOINT', 0, 0, MESG, 0 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )


C...........   Formatted file I/O formats............ 93xxx

93010   FORMAT( 4 I10, F10.2 )


C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10 ( A, :, I5, :, 2X  ) )


        END

