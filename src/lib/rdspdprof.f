
       SUBROUTINE RDSPDPROF( FDEV )
       
C***********************************************************************
C  subroutine body starts at line 117
C
C  DESCRIPTION:
C      This subroutine reads the speed profiles file, then creates
C      sorted arrays by profile number.
C
C  PRECONDITIONS REQUIRED:
C      FDEV is opened
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
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
C***************************************************************************

C...........   Modules for public variables
C...........  This module is for mobile-specific data
        USE MODMOBIL, ONLY: SPDPROFS, SPDNUMS, NSPDPROF
        
        IMPLICIT NONE
        
C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        
C...........   EXTERNAL FUNCTIONS and their descriptions:
        LOGICAL         CHKINT
        LOGICAL         CHKREAL
        LOGICAL         BLKORCMT
        CHARACTER(2)    CRLF
        INTEGER         GETFLINE
        INTEGER         STR2INT
        REAL            STR2REAL
        
        EXTERNAL        BLKORCMT, CHKINT, CHKREAL, CRLF, GETFLINE, 
     &                  STR2INT, STR2REAL

C...........   Subroutine arguments
        INTEGER      , INTENT  (IN) :: FDEV       ! File unit number

C...........   Local parameters
        INTEGER    , PARAMETER :: MXCOL = 25

C...........   Array of input fields
        CHARACTER(SPDLEN3)  SEGMENT( MXCOL )

C...........   Local, allocatable arrays
        REAL,    ALLOCATABLE :: SPDPROFA( :,: )   ! unsorted hourly speed profiles
        INTEGER, ALLOCATABLE :: SPDNUMA ( : )     ! unsorted profile numbers
        INTEGER, ALLOCATABLE :: SPDIDXA ( : )     ! sorting index
                
C...........   Local variables
        INTEGER         I,J,K                     ! counters and indices
        
        INTEGER         IOS                       ! I/O status
        INTEGER         NLINES                    ! no. lines in file
        INTEGER         PREVNUM                   ! previous profile number
        
        LOGICAL      :: EFLAG = .FALSE.           ! true: error found
        
        CHARACTER(150)  LINE                      ! Read buffer for a line
        CHARACTER(256)  MESG                      ! message buffer

        CHARACTER(16) :: PROGNAME = 'RDSPDPROF'    ! program name

C***********************************************************************
C   Begin body of subroutine RDSPDPROF

C........  Count number of lines in file
       NLINES = GETFLINE( FDEV, 'Speed profiles file')

C........  Allocate unsorted arrays       
       ALLOCATE( SPDPROFA( NLINES, 24 ), STAT=IOS )
       CALL CHECKMEM( IOS, 'SPDPROFA', PROGNAME )
       ALLOCATE( SPDNUMA ( NLINES ), STAT=IOS )
       CALL CHECKMEM( IOS, 'SPDNUMA', PROGNAME )
       ALLOCATE( SPDIDXA ( NLINES ), STAT=IOS )
       CALL CHECKMEM( IOS, 'SPDIDXA', PROGNAME )

       K = 0

C........  Loop through file
       DO I = 1, NLINES
       
           READ( FDEV, 93000, IOSTAT=IOS ) LINE

C............  Check for I/O errors       
           IF( IOS > 0 ) THEN
               EFLAG = .TRUE.
               WRITE( MESG, 94010 ) 
     &            'I/O error', IOS, 
     &            'reading speed profiles file at line', I
               CALL M3MESG( MESG )
               CYCLE
           END IF
       
C............  Check for end of file
           IF( IOS < 0 ) THEN
               MESG = 'End of file reached unexpectedly. ' //
     &                'Check format of speed ' // CRLF() // BLANK5 //
     &                'profiles file.'
               CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
           END IF

C............  Skip blank lines and comment lines
           IF( BLKORCMT( LINE ) ) CYCLE

C............  Parse line into segments
           CALL PARSLINE( LINE, MXCOL, SEGMENT ) 

C............  Check that the profile number is a positive integer
           IF( .NOT. CHKINT( SEGMENT( 1 ) ) .OR.
     &         STR2INT( SEGMENT( 1 ) ) < 1 ) THEN
               EFLAG = .TRUE.
               WRITE( MESG,94010 ) 'ERROR: Profile number is not ' //
     &                'an integer or is less than 1 at line', I
               CALL M3MESG( MESG )
           END IF

C............  Check that hourly speeds are valid
           DO J = 2, 25
               IF( .NOT. CHKREAL( SEGMENT( J ) ) .OR.
     &             STR2REAL( SEGMENT( J ) ) < 0. ) THEN
                   EFLAG = .TRUE.
                   WRITE( MESG,94010 ) 'ERROR: Invalid hourly ' //
     &                    'speed for hour', J-1, 'at line', I
                   CALL M3MESG( MESG )
               END IF
           END DO
            
C............  Skip rest of loop if there has been an error
           IF( EFLAG ) CYCLE
           
           K = K + 1
           IF( K .GT. NLINES ) CYCLE  ! Ensure no overflow
           
C............  Store profile information
           SPDIDXA( K ) = K
           SPDNUMA( K ) = STR2INT( SEGMENT( 1 ) )
            
           DO J = 2, 25
               SPDPROFA( K,J-1 ) = STR2REAL( SEGMENT( J ) )
           END DO           
            
       END DO  ! End loop reading SPDPROF file

       NSPDPROF = K
       
       IF( NSPDPROF == 0 ) THEN
           EFLAG = .TRUE.
           MESG = 'ERROR: No valid entries in SPDPROF file!'
           CALL M3MSG2( MESG )
       
       ELSEIF( NSPDPROF > NLINES ) THEN
           EFLAG = .TRUE.
           WRITE( MESG,94010 ) 'INTERNAL ERROR: dimension for ' //
     &            'storing SPDPROF file was', NLINES,
     &            CRLF() // BLANK10 // 'but actually needed', NSPDPROF
           CALL M3MSG2( MESG )
           
       END IF

C........  Check for errors reading SPDPROF file and abort
       IF( EFLAG ) THEN
           MESG = 'Problem reading speed profiles file.'
           CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
       END IF

C........  Sort speed profiles by profile number
       CALL SORTI1( NSPDPROF, SPDIDXA, SPDNUMA )

C........  Allocate memory for sorted arrays
       ALLOCATE( SPDPROFS( NSPDPROF, 24 ), STAT=IOS )
       CALL CHECKMEM( IOS, 'SPDPROFS', PROGNAME )
       ALLOCATE( SPDNUMS ( NSPDPROF ), STAT=IOS )
       CALL CHECKMEM( IOS, 'SPDNUMS', PROGNAME )
       
C........  Create sorted profile arrays       
       PREVNUM = 0

       DO I = 1, NSPDPROF
           J = SPDIDXA( I )

C............  Check that duplicate profile numbers have not been used           
           IF( SPDNUMA( J ) == PREVNUM ) THEN
               EFLAG = .TRUE.
               WRITE( MESG,94010 ) 'ERROR: At least two speed ' //
     &                'profiles are assigned the profile number',
     &                PREVNUM
               CALL M3MESG( MESG )
               CYCLE
           END IF
           
           SPDNUMS( I ) = SPDNUMA( J )
           PREVNUM = SPDNUMA( J )
           
           SPDPROFS( I,: ) = SPDPROFA( J,: )
       END DO

       IF( EFLAG ) THEN
           MESG = 'Problem with speed profiles file.'
           CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
       END IF

C........  Deallocate local arrays
       DEALLOCATE( SPDIDXA, SPDNUMA, SPDPROFA )

       RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000  FORMAT( A )
       
C...........   Internal buffering formats............ 94xxx

94010  FORMAT( 10( A, :, I8, :, 1X ) )

       END SUBROUTINE RDSPDPROF
       
