
        SUBROUTINE RDCEMSUM

C***************************************************************************
C  subroutine body starts at line 86
C
C  DESCRIPTION:
C      This subroutine reads the CEM summary file produced by the utility
C      CEMScan.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C      Subroutines: I/O API subroutine
C
C  REVISION  HISTORY:
C      Created 06/05 by C. Seppanen
C
C***************************************************************************
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

C.........  MODULES for public variables
C.........  This module contains data for day- and hour-specific data
        USE MODDAYHR, ONLY: NOBRLIST, OBRLIST, ANNNOX, ANNSO2, ANNGLOAD,
     &                      ANNSLOAD, ANNHEAT

        IMPLICIT NONE

C.........  INCLUDES
        INCLUDE 'EMCNST3.EXT'   ! emissions constant parameters

C.........  EXTERNAL FUNCTIONS
        LOGICAL      BLKORCMT
        CHARACTER(2) CRLF
        INTEGER      FINDC
        INTEGER      PROMPTFFILE
        
        EXTERNAL     BLKORCMT, CRLF, FINDC, PROMPTFFILE

C.........  SUBROUTINE ARGUMENTS

C.........  File names and unit numbers
        INTEGER IDEV    ! input file unit
        
C.........  Other local variables
        INTEGER            I, N                   ! counters
        INTEGER            IOS                    ! I/O status
        
        REAL               NOXVAL                 ! tmp NOx value
        REAL               SO2VAL                 ! tmp SO2 value
        REAL               OPTIME                 ! tmp operating time (unused)
        REAL               GLOAD                  ! tmp gross load
        REAL               SLOAD                  ! tmp steam load
        REAL               HTINPUT                ! tmp heat input
        
        LOGICAL ::         EFLAG = .FALSE.        ! true: an error occurred
        
        CHARACTER(ORSLEN3) CORS                   ! tmp ORIS ID
        CHARACTER(BLRLEN3) BLID                   ! tmp boiler ID
        CHARACTER(256)     MESG                   ! message buffer

        CHARACTER(16) ::   PROGNAME = 'RDCEMSUM'  ! program name
        
C***********************************************************************
C   begin body of program RDCEMSUM

C.........  Prompt for input file
        MESG = 'Enter logical name of the CEM SUMMARY file'
        IDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'CEMSUM', PROGNAME )

C.........  Loop through file and count number of valid combinations
        NOBRLIST = 0
        DO
            READ( IDEV, *, IOSTAT=IOS ) MESG

C.............  Check for I/O errors
            IF( IOS > 0 ) THEN
                MESG = 'Problem reading CEM summary file'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  Check for end of file
            IF( IOS < 0 ) EXIT

C.............  Skip blank lines and comments
            IF( BLKORCMT( MESG ) ) THEN
                CYCLE
            ELSE
                BACKSPACE( IDEV )
            END IF
            
            READ( IDEV, 93010 ) CORS, BLID, NOXVAL, SO2VAL,
     &          OPTIME, GLOAD, SLOAD, HTINPUT

            NOBRLIST = NOBRLIST + 1
        END DO

        REWIND( IDEV )

C.........  Exit if error
        IF( EFLAG ) THEN
            MESG = 'Problem reading CEM summary file'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Allocate memory to store data
        ALLOCATE( OBRLIST( NOBRLIST ), STAT=IOS )
        CALL CHECKMEM( IOS, 'OBRLIST', PROGNAME )
        ALLOCATE( ANNNOX( NOBRLIST ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ANNNOX', PROGNAME )
        ALLOCATE( ANNSO2( NOBRLIST ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ANNSO2', PROGNAME )
        ALLOCATE( ANNGLOAD( NOBRLIST ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ANNGLOAD', PROGNAME )
        ALLOCATE( ANNSLOAD( NOBRLIST ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ANNSLOAD', PROGNAME )
        ALLOCATE( ANNHEAT( NOBRLIST ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ANNHEAT', PROGNAME )

C.........  Loop through file and store data        
        N = 0
        DO
            READ( IDEV, *, IOSTAT=IOS ) MESG

C.............  Check for end of file
            IF( IOS < 0 ) EXIT
        
            IF( BLKORCMT( MESG ) ) THEN
                CYCLE
            ELSE
                BACKSPACE( IDEV )
            END IF
            
            READ( IDEV, 93010 ) CORS, BLID, NOXVAL, SO2VAL,
     &          OPTIME, GLOAD, SLOAD, HTINPUT
            
            N = N + 1
            
            OBRLIST ( N ) = ADJUSTR( CORS ) // ADJUSTR( BLID )
            ANNNOX  ( N ) = NOXVAL
            ANNSO2  ( N ) = SO2VAL
            ANNGLOAD( N ) = GLOAD
            ANNSLOAD( N ) = SLOAD
            ANNHEAT ( N ) = HTINPUT
        
        END DO
        
        NOBRLIST = N

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

93010   FORMAT( A6, 1X, A6, 6( 1X, E17.10 ) )
        
        END SUBROUTINE RDCEMSUM
