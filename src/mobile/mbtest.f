
         PROGRAM MBTEST

C***********************************************************************
C  program body starts at line
C
C
C****************************************************************************
C
C...........   MODULES for public variables
C...........   This module is the inventory arrays
        USE MODSOURC

        IMPLICIT NONE

C...........   INCLUDES:
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS and their descriptions:
        INTEGER   GETFLINE
        INTEGER   PROMPTFFILE

        EXTERNAL  GETFLINE, PROMPTFFILE
        
C...........   File names and unit numbers:

        INTEGER         IDEV    !  Inventory file (various formats)
        INTEGER         LDEV    !  log-device

C...........   Other local variables
                                
        INTEGER         IOS     ! tmp i/o status
        INTEGER         MXMSRC  ! maximum number of mobile sources
        INTEGER         NDROP   ! number of dropped records
        INTEGER         NMSRC   ! actual number of mobile sources

        REAL            VDROP   ! dropped VMT from inventory file          

        LOGICAL      :: EFLAG = .FALSE.  !  TRUE iff ERROR

        CHARACTER*300   MESG    !  text for M3EXIT()

        CHARACTER*16  :: PROGNAME = 'MBTEST'   !  program name

C***********************************************************************
C   begin body of program MBTEST

        LDEV = INIT3()

C.........  Write out copywrite, version, web address, header info, and prompt
C           to continue running the program.
C        CALL INITEM( LDEV, SCCSW, PROGNAME )

C.........  Get file name for opening input raw point source file
        IDEV = PROMPTFFILE( 
     &         'Enter the name of the RAW MOBILE INVENTORY file',
     &          .TRUE., .TRUE., 'MBINV', PROGNAME )

C.........  Get number of lines in mobile file
        MXMSRC = GETFLINE( IDEV, 'Mobile input inventory')

        CALL M3MSG2( 'Setting up to read inventory data...' )

C.........  Allocate memory for inventory 
        ALLOCATE( CSCCA( MXMSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSCCA', PROGNAME )
        ALLOCATE( VMTA( MXMSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'VMTA', PROGNAME )
        ALLOCATE( CVTYPEA( MXMSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CVTYPEA', PROGNAME )
        ALLOCATE( ISPEEDA( MXMSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ISPEEDA', PROGNAME )
        ALLOCATE( INVYRA( MXMSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INVYRA', PROGNAME )
        ALLOCATE( IFIPA( MXMSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IFIPA', PROGNAME )

        CALL M3MSG2( 'Reading raw inventory data...' )

C.........  Read the raw inventory data to test the reader

        CALL RDIDAMB( IDEV, MXMSRC, NMSRC, EFLAG, NDROP, VDROP )

        IF( NDROP .GT. 0 ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94020 ) 'VMT dropped from inventory because ' //
     &             'of invalid records: ', VDROP
            CALL M3WARN( PROGNAME, 0, 0, MESG )

        END IF

        IF( EFLAG ) THEN
           MESG = 'Error reading raw inventory file(s)'         
           CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ENDIF

C.........  End program successfully
        MESG = ' '
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 0 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94020   FORMAT( A, F12.1 )
 
        END PROGRAM MBTEST
