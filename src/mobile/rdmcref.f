
        SUBROUTINE RDMCREF( MDEV )

C.........  MODULES for public variables
C.........  This module contains the inventory arrays
        USE MODSOURC

C.........  This module contains the information about the source category
        USE MODINFO

C.........  This module contains the lists of unique source characteristics
        USE MODLISTS
                
        USE MODMBSET
                
        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        LOGICAL        CHKINT
        INTEGER        GETFLINE
        INTEGER        STR2INT
        INTEGER        FIND1
        CHARACTER*2    CRLF    
        
        EXTERNAL  CHKINT, GETFLINE, STR2INT, FIND1, CRLF

C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT (IN) :: MDEV     ! MCREF file unit no.

C...........   Local allocatable arrays
        INTEGER, ALLOCATABLE :: MCREFRAW ( :,: )  ! raw MCREF data
        INTEGER, ALLOCATABLE :: MCREFINV ( : )    ! sorted list of inventory counties in MCREF
          
        INTEGER, ALLOCATABLE :: IDX  ( : )    ! index into MCREF data
        INTEGER, ALLOCATABLE :: IDX2 ( : )    ! index into inv. county array

C...........   Local arrays
        CHARACTER*3   SEGMENT( 4 )          ! parsed input line
        
C...........   Other local variables
        INTEGER I, J, N, K                ! counters and indices                     
        
        INTEGER IOS                       ! I/O status
        INTEGER :: IREC = 0               ! record counter
        INTEGER NLINES                    ! number of lines
        
        INTEGER REFCOUNTY                 ! ref. county FIPS code
        INTEGER INVCOUNTY                 ! inv. county FIPS code
        INTEGER PRCOUNTY                  ! previous ref. county
        INTEGER PICOUNTY                  ! previous inv. county
        INTEGER :: NREF = 0               ! number of ref. counties
                
        LOGICAL      :: DUPFLAG = .FALSE.   ! true: duplicate entries found
        LOGICAL      :: EFLAG   = .FALSE.   ! true: error found    
        
        CHARACTER(LEN=100)     LINE     !  line buffer
        CHARACTER(LEN=300)     MESG     !  message buffer

        CHARACTER*16 :: PROGNAME = 'RDMCREF'   ! program name
        
C***********************************************************************
C   begin body of subroutine RDMCREF

C.........  Store no. inv. counties        	
        NINVC = NINVIFIP

C.........  Allocate arrays that use NINVC as dimension (sorted arrays and index)
        ALLOCATE( MCREFINV ( NINVC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MCREFINV', PROGNAME )
        ALLOCATE( IDX2 ( NINVC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IDX2', PROGNAME )
        
        ALLOCATE( MCREFSORT ( NINVC,2 ), STAT=IOS ) 
        CALL CHECKMEM( IOS, 'MCREFSORT', PROGNAME )
        
C.........  Get the number of lines in the file     
        NLINES = GETFLINE( MDEV, 'County cross-reference file' )

C.........  Allocate arrays that use NLINES as dimension (unsorted arrays and index)        
        ALLOCATE( MCREFRAW ( NLINES,2 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MCREFRAW', PROGNAME )
        ALLOCATE( IDX ( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IDX', PROGNAME )
        
C.........  Initialize arrays
        MCREFRAW = 0.
        MCREFINV = 0.
        MCREFSORT = 0.
        
        DO I = 1, NLINES       
            
            IDX( I ) = I
                
C.........  Read line
            READ( MDEV, 93000, END=999, IOSTAT=IOS ) LINE
            
            IREC = IREC + 1
            
            IF ( IOS /= 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG, 94010 )
     &              'I/O error', IOS,
     &              'reading county cross-reference file at line', IREC
                CALL M3MESG( MESG )
                CYCLE
            END IF
            
C.............  Skip blank lines
            IF( LINE == ' ' ) CYCLE

C.............  Parse the line into 4 segments
            CALL PARSLINE( LINE, 4, SEGMENT )

C.............  Convert inventory county to integer             
            IF( CHKINT(SEGMENT(1)) .AND. CHKINT(SEGMENT(2)) ) THEN
            	CALL PADZERO( SEGMENT( 2 ) )
                INVCOUNTY = STR2INT( ADJUSTR(SEGMENT(1)) // 
     &                               SEGMENT(2) )
            ELSE
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: Bad inventory county ' //
     &                     'FIPS code at line', IREC, 'of county ' //
     &                     'cross-reference file.'
                CALL M3MESG( MESG )
                CYCLE
            END IF
            
C.............  Convert reference county to integer 
            IF( CHKINT(SEGMENT(3)) .AND. CHKINT(SEGMENT(4)) ) THEN
                CALL PADZERO( SEGMENT( 4 ) )
                REFCOUNTY = STR2INT( ADJUSTR(SEGMENT(3)) // 
     &                               SEGMENT(4) )
            ELSE
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: Bad reference county ' //
     &                     'FIPS code at line', IREC, 'of county ' //
     &                     'cross-reference file.'
                CALL M3MESG( MESG )
                CYCLE
            END IF
                        
C.............  Store values in unsorted array
            MCREFRAW( I,1 ) = INVCOUNTY
            MCREFRAW( I,2 ) = REFCOUNTY

        END DO  ! done reading MCREF file

C.........  Close MCREF file
        CLOSE( MDEV )
        
C.........  Abort if error found while reading cross-reference file
        IF( EFLAG ) THEN
            MESG = 'Problem reading cross-reference file'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Sort MCREF index array by reference county
        CALL SORTI2( NLINES, IDX, MCREFRAW(:,2), MCREFRAW(:,1) )

C.........  Check for duplicate references, then store sorted MCREF array
        PRCOUNTY = 0
        PICOUNTY = 0
        N = 0

        DO I = 1, NLINES
            J = IDX( I )
            
            REFCOUNTY = MCREFRAW( J,2 )
            INVCOUNTY = MCREFRAW( J,1 )

C.............  Skip any entries equal to zero due to blank lines
            IF( REFCOUNTY == 0 .OR. INVCOUNTY == 0 ) CYCLE

            IF( REFCOUNTY == PRCOUNTY .AND. 
     &          INVCOUNTY == PICOUNTY ) THEN
     	
     	        DUPFLAG = .TRUE.
     	        EFLAG   = .TRUE.
                
                WRITE( MESG,94010 ) 'ERROR: Duplicate entries in ' //
     &                 'county cross-reference file for ' // CRLF() //
     &                 BLANK10 // 'inventory county', INVCOUNTY,
     &                 ', reference county', REFCOUNTY, '.'
                CALL M3MESG( MESG )
            ELSE
            	
C.............  Check that current county is in the inventory
                K = FIND1( INVCOUNTY, NINVC, INVIFIP )
                
                IF( K < 0 ) THEN
                    WRITE( MESG,94010 ) 'WARNING: Ignoring county ' //
     &                     'cross-reference for inventory county',
     &                     INVCOUNTY, ',' // CRLF() // BLANK10 // 
     &                     ' since it is not in the mobile inventory.'
                    CALL M3MESG( MESG )
                    
                ELSE    
                    N = N + 1
        
                    MCREFSORT( N,1 ) = INVCOUNTY
                    MCREFSORT( N,2 ) = REFCOUNTY
                END IF
                
            END IF
            
            PRCOUNTY = REFCOUNTY
            PICOUNTY = INVCOUNTY
            
        END DO

        IF( DUPFLAG ) THEN
            MESG = 'ERROR: Duplicate county cross-reference ' //
     &             'entries found. ' //CRLF()// BLANK10 // 
     &             'Remove duplicate entries and try again.'
            CALL M3MSG2( MESG )
        END IF

C.........  Initialize index for sorting
        DO I = 1, NINVC
            IDX2( I ) = I
        END DO
        	
C.........  Sort list of inv. counties in MCREF        	
        CALL SORTI1( NINVC, IDX2, MCREFSORT( :,1 ) )

C.........  Created sorted array
        DO I = 1, NINVC
            MCREFINV( I ) = MCREFSORT( IDX2( I ),1 )
        END DO

C.........  Check that all counties in the inventory have a reference county
        DO I = 1, NINVIFIP
        
C.............  Get county from inventory        
            INVCOUNTY = INVIFIP( I )
            
            J = FIND1( INVCOUNTY, NINVC, MCREFINV ) 

            IF( J < 0 ) THEN
                EFLAG = .TRUE.	
                WRITE( MESG, 94010 ) 'ERROR: Missing reference ' //
     &                 'county for inventory county', INVCOUNTY,
     &                CRLF() // BLANK10 // 'in county ' //
     &                'cross-reference file.'
                CALL M3MESG( MESG )
            END IF

        END DO

        IF( EFLAG ) THEN
            MESG = 'Problem(s) found in county cross-reference file.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Count total number of reference counties
        PRCOUNTY = 0

        DO I = 1, NINVC
            REFCOUNTY = MCREFSORT( I,2 )
            
            IF( REFCOUNTY /= PRCOUNTY ) THEN
               NREF = NREF + 1
            END IF
            
            PRCOUNTY = REFCOUNTY
               
        END DO

        NREFC = NREF

C.........  Create reference county index array
        ALLOCATE( MCREFIDX ( NREFC,2 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MCREFIDX', PROGNAME )
        
        PRCOUNTY = 0
        N = 0
        
        DO I = 1, NINVC
            REFCOUNTY = MCREFSORT( I,2 )
            
            IF( REFCOUNTY /= PRCOUNTY ) THEN
            	N = N + 1
            	
                MCREFIDX( N,1 ) = REFCOUNTY
                MCREFIDX( N,2 ) = I
            END IF
            
            PRCOUNTY = REFCOUNTY
        END DO

        RETURN

999     MESG = 'End of file'
        MESG = 'End of file reached unexpectedly. ' //
     &         'Check format of MCREF' // CRLF() // BLANK5 //
     &         'input file.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )   
             
C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )  
      
C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )
94020   FORMAT( 3( A, 1X ), I8, 1X, A, 1X )
        
        END SUBROUTINE RDMCREF