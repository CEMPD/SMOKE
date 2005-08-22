
        SUBROUTINE RDMVREF( MDEV )

C***********************************************************************
C  subroutine body starts at line 99
C
C  DESCRIPTION:
C       Reads MVREF file, checks that settings are given for each 
C       reference county, ignores counties not specified in the MCREF
C       file, and sorts the data
C
C  PRECONDITIONS REQUIRED:
C       MDEV has been opened
C
C  SUBROUTINES AND FUNCTIONS CALLED:  none
C
C  REVISION  HISTORY:
C     10/01: Created by C. Seppanen
C
C***********************************************************************
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
C***********************************************************************

C.........  MODULES for public variables
C.........  This module is used for MOBILE6 setup information        
        USE MODMBSET, ONLY: NREFC, MVREFSORT, MCREFIDX, NREFFLAGS
                
        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        LOGICAL        CHKINT
        LOGICAL        BLKORCMT
        INTEGER        GETFLINE
        INTEGER        STR2INT
        INTEGER        FIND1
        CHARACTER(2)   CRLF    
        
        EXTERNAL  BLKORCMT, CHKINT, GETFLINE, STR2INT, FIND1, CRLF

C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT (IN) :: MDEV     ! MVREF file unit no.

C...........   Local allocatable arrays
        INTEGER, ALLOCATABLE :: MVREFRAW ( :,: )  ! raw MVREF data
          
        INTEGER, ALLOCATABLE :: IDX ( : )    ! index into MVREF data

C...........   Local arrays
        INTEGER       SETTINGS( NREFFLAGS )     ! converted settings flags

        CHARACTER(3)  SEGMENT( NREFFLAGS + 2 )  ! parsed input line
        
C...........   Other local variables
        INTEGER I, J, K, N                ! counters and indices                     
        
        INTEGER IOS                       ! I/O status
        INTEGER :: IREC = 0               ! record counter
        INTEGER NLINES                    ! number of lines

        INTEGER REFCOUNTY                 ! ref. county FIPS code
        INTEGER PRCOUNTY                  ! previous ref. county
                        
        LOGICAL      :: DUPFLAG = .FALSE.   ! true: duplicate entries found
        LOGICAL      :: EFLAG   = .FALSE.   ! true: error found    
        LOGICAL      :: SETFLAG = .FALSE.   ! true: error in settings values 
        LOGICAL      :: RFLAG   = .FALSE.   ! true: no. lines < no. ref. counties
        
        CHARACTER(100)     LINE     !  line buffer
        CHARACTER(300)     MESG     !  message buffer

        CHARACTER(16) :: PROGNAME = 'RDMVREF'   ! program name
        
C***********************************************************************
C   begin body of subroutine RDMVREF

C.........  Get the number of lines in the file     
        NLINES = GETFLINE( MDEV, 'Reference county settings file' )

C.........  Allocate memory to store settings information        
        ALLOCATE( MVREFRAW ( NLINES,NREFFLAGS + 1 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MVREFRAW', PROGNAME )
        ALLOCATE( IDX ( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IDX', PROGNAME )
        
C.........  If the no. of lines in the file is less than the no. of ref.
C              counties, then something is wrong, but we'll go through
C              the function to print out error messages
        IF( NLINES < NREFC ) THEN
            RFLAG = .TRUE.
            
            ALLOCATE( MVREFSORT ( NLINES,NREFFLAGS + 1 ), STAT=IOS )
            CALL CHECKMEM( IOS, 'MVREFSORT', PROGNAME )
        ELSE
            ALLOCATE( MVREFSORT ( NREFC,NREFFLAGS + 1 ), STAT=IOS ) 
            CALL CHECKMEM( IOS, 'MVREFSORT', PROGNAME )
        END IF
        
C.........  Initialize arrays
        MVREFRAW = 0.
        MVREFSORT = 0.
        
        DO I = 1, NLINES       
                
C.........  Read line
            READ( MDEV, 93000, END=999, IOSTAT=IOS ) LINE
            
            IREC = IREC + 1
            
            IF ( IOS /= 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG, 94010 )
     &              'I/O error', IOS, 'reading reference county ' //
     &              'settings file at line', IREC
                CALL M3MESG( MESG )
                CYCLE
            END IF
            
C.............  Skip blank lines
            IF( BLKORCMT( LINE ) ) CYCLE

C.............  Parse the line into segments
            CALL PARSLINE( LINE, NREFFLAGS + 2, SEGMENT )

C.............  Convert reference county to integer             
            IF( CHKINT(SEGMENT(1)) .AND. CHKINT(SEGMENT(2)) ) THEN
                CALL PADZERO( SEGMENT( 2 ) )
                REFCOUNTY = STR2INT( ADJUSTR(SEGMENT(1)) // 
     &                               SEGMENT(2) )
            ELSE
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: Bad reference county ' //
     &                     'FIPS code at line', IREC, 'of ' //
     &                     'reference county settings file.'
                CALL M3MESG( MESG )
                CYCLE
            END IF

C.............  Check settings flags and convert to integers
            DO J = 1, NREFFLAGS

                IF( CHKINT( SEGMENT( J + 2 ) ) ) THEN
                    SETTINGS( J ) = STR2INT( SEGMENT( J + 2 ) )
                ELSE
                    SETFLAG = .TRUE.
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 'ERROR: Bad settings value ' //
     &                          'at line', IREC, 'of reference ' //
     &                          CRLF() // BLANK10 // 'county ' //
     &                          'settings file.'
                    CALL M3MESG( MESG )
                    EXIT
                END IF
            END DO

C.............  Skip rest of loop if there was an error with the settings            
            IF( SETFLAG ) THEN
                SETFLAG = .FALSE.
                CYCLE
            END IF
                        
C.............  Store values in unsorted array
            MVREFRAW( I,1 ) = REFCOUNTY
            MVREFRAW( I,2:NREFFLAGS + 1 ) = SETTINGS( 1:NREFFLAGS )

        END DO  ! done reading MVREF file

C.........  Close MVREF file
        CLOSE( MDEV )
        
C.........  Abort if error found while reading settings file
        IF( EFLAG ) THEN
            MESG = 'Problem reading settings file'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.............  Sort MVREF index array
        DO I = 1, NLINES
            IDX( I ) = I
        END DO
                
        CALL SORTI1( NLINES, IDX, MVREFRAW( :,1 ) )

C.............  Check for duplicate references, then store sorted MVREF array
        PRCOUNTY = 0
        N = 0

        DO I = 1, NLINES
            J = IDX( I )
            
            REFCOUNTY = MVREFRAW( J,1 )

C.............  Skip any entries equal to zero due to blank lines
            IF( REFCOUNTY == 0 ) CYCLE
            
            IF( REFCOUNTY == PRCOUNTY ) THEN

                DUPFLAG = .TRUE.
                EFLAG   = .TRUE.
                
                WRITE( MESG,94010 ) 'ERROR: Duplicate entries in ' //
     &                 'reference county settings file for ' // 
     &                 CRLF() // BLANK10 // 'reference county', 
     &                 REFCOUNTY, '.'
                CALL M3MESG( MESG )
            ELSE

C.............  Check that current reference county is in the county cross-reference list
                K = FIND1( REFCOUNTY, NREFC, MCREFIDX( :,1 ) )
                
                IF( K < 0 ) THEN
                    WRITE( MESG,94010 ) 'WARNING: Ignoring ' //
     &                     'settings for reference county',
     &                     REFCOUNTY, ',' // CRLF() // BLANK10 // 
     &                     ' since it is not in the county ' //
     &                     ' cross-reference file.'
                    CALL M3MESG( MESG )
                    
                ELSE    
                    N = N + 1
        
                    MVREFSORT( N,1 ) = REFCOUNTY
                    MVREFSORT( N,2:NREFFLAGS + 1 ) = 
     &                              MVREFRAW( J,2:NREFFLAGS + 1 )
                END IF
                
            END IF
            
            PRCOUNTY = REFCOUNTY
            
        END DO

        IF( DUPFLAG ) THEN
            MESG = 'ERROR: Duplicate reference county settings ' //
     &             'found. ' // CRLF() // BLANK10 // 
     &             'Remove duplicate entries and try again.'
            CALL M3MSG2( MESG )
        END IF

C.............  Check that all reference counties have settings
        DO I = 1, NREFC
        
C.............  Get county from MCREF index array
            REFCOUNTY = MCREFIDX( I,1 )
           
            IF( RFLAG ) THEN
                J = FIND1( REFCOUNTY, NLINES, MVREFSORT( :,1 ) )
            ELSE
                J = FIND1( REFCOUNTY, NREFC, MVREFSORT( :,1 ) ) 
            END IF

            IF( J < 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG, 94010 ) 'ERROR: Missing settings ' //
     &                 ' for reference county', REFCOUNTY,
     &                 CRLF() // BLANK10 // 'in reference ' //
     &                 'county settings file.'
                CALL M3MESG( MESG )
            END IF
            
        END DO

        IF( EFLAG ) THEN
            MESG = 'Problem(s) found in reference county settings file.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF
        
        RETURN

999     MESG = 'End of file'
        MESG = 'End of file reached unexpectedly. ' //
     &         'Check format of MVREF' // CRLF() // BLANK5 //
     &         'input file.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )   
             
C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )  
      
C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )
94020   FORMAT( 3( A, 1X ), I8, 1X, A, 1X )
        
        END SUBROUTINE RDMVREF
