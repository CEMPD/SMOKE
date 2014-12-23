
        SUBROUTINE RDORSDSC( FDEV )

C**************************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      This subroutine allocates memory for and reads the ORIS FIPS
C      codes and plant names
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Created 11/2001 by M. Houyoux
C
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
C**************************************************************************

C...........   Modules for public variables
C.........  This module contains the arrays for state and county summaries
        USE MODSTCY, ONLY: ORISFIP, ORISLST, ORISDSC, NORIS

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constat parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:

        LOGICAL        BLKORCMT
        LOGICAL        CHKINT
        INTEGER        GETFLINE
        INTEGER        STR2INT
        LOGICAL        USEEXPGEO

        EXTERNAL       BLKORCMT, CHKINT, GETFLINE, STR2INT, USEEXPGEO

C...........   Subroutine arguments
        INTEGER, INTENT (IN) :: FDEV          ! file unit number

C...........   Local paramaters
        INTEGER, PARAMETER :: NFIELD = 11     ! no. fields in input file

C...........   Local arrays
        CHARACTER(60) SEGMENT( NFIELD )        ! line parsing array

C...........   Local variables
        INTEGER         I                     ! indices and counters

        INTEGER         ENDLEN                ! end length for reading descriptn
        INTEGER         IOS                   ! i/o status
        INTEGER      :: IREC = 0              ! record number
        INTEGER      :: NLINES = 0            ! number of lines in input file

        LOGICAL      :: EFLAG = .FALSE.       ! true: error found

        CHARACTER(256)  LINE                  ! Read buffer for a line
        CHARACTER(256)  MESG                  ! Message buffer

        CHARACTER(ORSLEN3 ) CORS          ! tmp ORIS ID
        CHARACTER(DSCLEN3 ) PDSC          ! tmp plant description

        CHARACTER(16) :: PROGNAME = 'RDORSDSC'    !  program name

C***********************************************************************
C   Begin body of subroutine RDORSDSC

C.........  Write status message
        MESG = 'Reading ORIS descriptions file...'
        CALL M3MSG2( MESG )

C.........  Get the number of lines in the holidays file
        NLINES = GETFLINE( FDEV, 'ORIS Descriptions' )

C.........  Allocate memory for the SCC descriptions and initialize
        ALLOCATE( ORISFIP( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ORISFIP', PROGNAME )
        ALLOCATE( ORISLST( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ORISLST', PROGNAME )
        ALLOCATE( ORISDSC( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ORISDSC', PROGNAME )

C.........  Read the SCC descriptions, and store with SCC
        I = 0
        DO IREC = 1, NLINES

            READ ( FDEV, 93000, END=998, IOSTAT=IOS ) LINE

            IF ( IOS .GT. 0 ) THEN
                WRITE( MESG, 94010)
     &                'I/O error', IOS, 'reading ORIS '//
     &                'description file at line', IREC
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  Left adjust line
            LINE = ADJUSTL( LINE )

C.............  Skip blank and comment lines
            IF( BLKORCMT( LINE ) ) CYCLE

C.............  Get SCC line
            CALL PARSLINE( LINE, NFIELD, SEGMENT )

C.............  Set tmp ASCII fields into arrays of the correct length
            CORS = ADJUSTL( SEGMENT( 1 ) )
            CORS = ADJUSTR( CORS )
            
            PDSC = ADJUSTL( SEGMENT( 5 ) )

C.............  Check for integer field for FIPS code
            IF ( .NOT. USEEXPGEO() .AND.
     &           .NOT. CHKINT( SEGMENT( 2 ) ) ) THEN
                EFLAG = .TRUE.
                WRITE( MESG, 94010 ) 'ERROR: FIPS code not an ' //
     &                 'integer at line', IREC
                CALL M3MSG2( MESG )
                CYCLE
            END IF

C.............  Store entry
            I = I + 1
            ORISFIP( I ) = SEGMENT( 2 )
            ORISLST( I ) = CORS
            ORISDSC( I ) = PDSC

        END DO

        NORIS = I

        IF ( EFLAG ) THEN

            MESG = 'Problem reading ORIS descriptions.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF

C.........  Successful completion
        RETURN

C.........  Unexpected end of file
998     MESG = 'INTERNAL ERROR: Unexpected end of SCC description file'
        CALL M3MSG2( MESG )

        CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )


C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94020   FORMAT( A, 1X, I8, 1X, A, 1X, F10.6, 1X, A )

94030   FORMAT( A, 1X, I6.6, A, 100( ' SSC(', I2.2, '):', F10.6, : ) )

        END SUBROUTINE RDORSDSC
