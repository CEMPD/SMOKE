
        SUBROUTINE RDGROUPS( GDEV, SRCARRAY, AVETYPE, NCNTY, LASTTIME )

C***********************************************************************
C  subroutine body starts at line 84
C
C  DESCRIPTION:
C       Reads the specified temporal averaging group file and stores
C       averaging type by source.
C
C  PRECONDITIONS REQUIRED:
C       GDEV must be open.
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
        
        USE MODINFO, ONLY: NSRC
        
        IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        INTEGER, EXTERNAL :: STR2INT, GETFLINE, FIND1FIRST
        
C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT (IN)  :: GDEV                ! GROUP file unit no.
        INTEGER, INTENT (OUT) :: SRCARRAY( NSRC,2 )  ! array to hold county codes
        INTEGER, INTENT (IN)  :: AVETYPE             ! temporal averaging type
        INTEGER, INTENT (OUT) :: NCNTY               ! no. counties for this aver. type
        LOGICAL, INTENT (IN)  :: LASTTIME            ! true: last time through routine

C...........   Local allocatable arrays
        INTEGER, ALLOCATABLE, SAVE :: SORTSRC( :,: ) ! sorted array of county codes
        INTEGER, ALLOCATABLE, SAVE :: IDX( : )       ! index to sort SRCARRAY
        
C...........   Local arrays
        CHARACTER(FIPLEN3) SEGMENT( 3 )          ! parsed input line
        
C...........   Other local variables
        INTEGER I, J, K                   ! counters and indices                     
        
        INTEGER IOS                       ! I/O status
        INTEGER :: IREC = 0               ! record counter
        INTEGER NLINES                    ! number of lines in GROUP file
        INTEGER CURRCOUNTY                ! current county

        LOGICAL :: EFLAG      = .FALSE.   ! true: error found
        LOGICAL :: INITIAL    = .TRUE.    ! true: first time through routine
        
        CHARACTER(100)     LINE     !  line buffer
        CHARACTER(300)     MESG     !  message buffer

        CHARACTER(16) :: PROGNAME = 'RDGROUPS'   ! program name
        
C***********************************************************************
C   begin body of subroutine RDGROUPS

C.........  Create sorted array first time through routine
        IF( INITIAL ) THEN
            ALLOCATE( SORTSRC( NSRC,2 ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SORTSRC', PROGNAME )
            ALLOCATE( IDX( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'IDX', PROGNAME )
            
            SORTSRC = 0
            
            DO I = 1, NSRC
                IDX( I ) = I
            END DO
            
            CALL SORTI1( NSRC, IDX, SRCARRAY( :,1 ) )
            
            DO I = 1, NSRC
               SORTSRC( I,1 ) = SRCARRAY( IDX( I ),1 )
            END DO
            
            INITIAL = .FALSE.
        END IF

C.........  Get number of lines in group file
        NLINES = GETFLINE( GDEV, 'GROUP file' )

        NCNTY = 0

C.........  Read through GROUP file for list of counties        
        DO I = 1, NLINES
        
C.........  Read line
            READ( GDEV, 93000, IOSTAT=IOS ) LINE
            
            IREC = IREC + 1
            
            IF( IOS /= 0 ) THEN
                EFLAG = .TRUE.
                
                IF( IOS == -1 ) THEN
                    MESG = 'End of file reached unexpectedly. ' //
     &              'Check format of GROUP file.'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )   
                END IF
                
                WRITE( MESG, 94010 )
     &              'I/O error', IOS,
     &              'reading GROUP file at line', IREC
                CALL M3MESG( MESG )
                CYCLE
            END IF
            
C.............  Skip blank lines
            IF( LINE == ' ' ) CYCLE        

C.............  Parse the line into 3 segments
            CALL PARSLINE( LINE, 3, SEGMENT )

C.............  Store current county            
            CURRCOUNTY = STR2INT( SEGMENT( 1 ) )

C.............  Increment county counter
            NCNTY = NCNTY + 1

C.............  Find start of current county in sorted source array            
            K = FIND1FIRST( CURRCOUNTY, NSRC, SORTSRC( :,1 ) )

C.............  For sources that match, store temporal averaging type            
            DO
            	
C.................  If county doesn't match, we're done
                IF( SORTSRC( K,1 ) /= CURRCOUNTY ) EXIT
                
                SORTSRC( K,2 ) = AVETYPE
                K = K + 1
                
C.................  Make sure we don't go outside the array
                IF( K > NSRC ) EXIT
                
            END DO

        END DO

C.........  Restore sorted array to original order
        IF( LASTTIME ) THEN
            DO I = 1, NSRC
                SRCARRAY( IDX( I ),2 ) = SORTSRC( I,2 )
            END DO
        END IF

C.........  Abort if error found while reading cross-reference file
        IF( EFLAG ) THEN
            MESG = 'Problem reading GROUP county list file'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF
        
        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )
93010   FORMAT( I6, 1X, I6, 1X, I1 ) 

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )     
        
        END SUBROUTINE RDGROUPS
        
