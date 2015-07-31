
        SUBROUTINE RDMXREF( MDEV, NCTY, GRIDCTY )

C***********************************************************************
C  subroutine body starts at line 107
C
C  DESCRIPTION:
C       Reads the MCREF file, checks that each inventory county is assigned
C       a reference county, ignores counties outside of the grid, and sorts the data
C
C  PRECONDITIONS REQUIRED:
C       MDEV must be opened
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
        USE MODMBSET, ONLY: NINVC, NREFC, MCREFSORT, MCREFIDX
                
        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        LOGICAL        BLKORCMT
        LOGICAL        CHKINT
        INTEGER        GETFLINE
        INTEGER        STR2INT
        INTEGER        FIND1
        CHARACTER(2)   CRLF    
        
        EXTERNAL  BLKORCMT, CHKINT, GETFLINE, STR2INT, FIND1, CRLF

C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT (IN) :: MDEV             ! MCREF file unit no.
        INTEGER, INTENT (IN) :: NCTY             ! no. counties
        INTEGER, INTENT (IN) :: GRIDCTY( NCTY )  ! counties in grid

C...........   Local allocatable arrays
        INTEGER, ALLOCATABLE :: MCREFRAW ( :,: )  ! raw MCREF data
        INTEGER, ALLOCATABLE :: IDX  ( : )    ! index into MCREF data

C...........   Local arrays
        CHARACTER(3)  SEGMENT( 6 )          ! parsed input line
        
C...........   Other local variables
        INTEGER I, J, N, K                ! counters and indices                     
        

        INTEGER :: CO = 0                 ! tmp country
        INTEGER :: ST = 0                 ! tmp state
        INTEGER :: CT = 0                 ! tmp county
        INTEGER :: IOS  = 0               ! I/O status
        INTEGER :: IREC = 0               ! record counter
        INTEGER :: NREF = 0               ! number of ref. counties
        INTEGER :: NLINES = 0             ! number of lines
        
        INTEGER REFCOUNTY                 ! ref. county FIPS code
        INTEGER INVCOUNTY                 ! inv. county FIPS code
        INTEGER PRCOUNTY                  ! previous ref. county
        INTEGER PICOUNTY                  ! previous inv. county
                
        LOGICAL      :: DUPFLAG = .FALSE.   ! true: duplicate entries found
        LOGICAL      :: EFLAG   = .FALSE.   ! true: error found    
        
        CHARACTER(100)     LINE     !  line buffer
        CHARACTER(300)     MESG     !  message buffer

        CHARACTER(16) :: PROGNAME = 'RDMCREF'   ! program name
        
C***********************************************************************
C   begin body of subroutine RDMCREF
        
C.........  Get the number of lines in the file     
        NLINES = GETFLINE( MDEV, 'County cross-reference file' )

C.........  Allocate arrays that use NLINES as dimension (unsorted arrays and index)        
        ALLOCATE( MCREFRAW ( NLINES,2 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MCREFRAW', PROGNAME )
        ALLOCATE( IDX ( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IDX', PROGNAME )
        
C.........  Initialize arrays
        MCREFRAW = 0.
        PICOUNTY = 0
        N = 0
                
        DO I = 1, NLINES       
            
            IDX( I ) = I
                
C.............  Read line
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
            
C.............  Skip blank or comment lines
            IF( BLKORCMT( LINE ) ) CYCLE

            CALL PARSLINE( LINE, 6, SEGMENT )

C.............  Check the format of input
            IF( .NOT. CHKINT( SEGMENT(1) ) ) EFLAG = .TRUE.
            IF( .NOT. CHKINT( SEGMENT(2) ) ) EFLAG = .TRUE.
            IF( .NOT. CHKINT( SEGMENT(3) ) ) EFLAG = .TRUE.
            IF( .NOT. CHKINT( SEGMENT(4) ) ) EFLAG = .TRUE.
            IF( .NOT. CHKINT( SEGMENT(5) ) ) EFLAG = .TRUE.
            IF( .NOT. CHKINT( SEGMENT(6) ) ) EFLAG = .TRUE.
            IF( SEGMENT(5)==' ' .OR. SEGMENT(6)==' ' ) EFLAG=.TRUE.

            IF( EFLAG ) THEN 
                WRITE( MESG,94010 ) 'ERROR: Bad inventory county ' //
     &                     'FIPS code at line', IREC, 'of county ' //
     &                     'cross-reference file.'
                CALL M3MESG( MESG )
                CYCLE
            END IF

C.............  Convert inventory county to integer             
            CO = STR2INT( SEGMENT( 1 ) )
            ST = STR2INT( SEGMENT( 2 ) )
            CT = STR2INT( SEGMENT( 3 ) )
            INVCOUNTY = CO*100000 + ST*1000 + CT

            CO = STR2INT( SEGMENT( 4 ) )
            ST = STR2INT( SEGMENT( 5 ) )
            CT = STR2INT( SEGMENT( 6 ) )
            REFCOUNTY = CO*100000 + ST*1000 + CT
           
C.............  Store values in unsorted array
            MCREFRAW( I,1 ) = INVCOUNTY
            MCREFRAW( I,2 ) = REFCOUNTY

C.............  Skip any entries equal to zero due to blank lines
            IF( REFCOUNTY == 0 .OR. INVCOUNTY == 0 ) CYCLE

C.............  Check if current inventory county is duplicate (match previous)
            IF( INVCOUNTY /= PICOUNTY ) THEN
 
C.............  Check that current county is inside the grid (and in the inventory)
                K = FIND1( INVCOUNTY, NCTY, GRIDCTY )

                IF( K > 0 ) N = N + 1

            END IF

            PICOUNTY = INVCOUNTY

        END DO  ! done reading MCREF file

        NINVC = N

C.........  Allocate arrays that use NINVC as dimension (sorted arrays and index)       
        ALLOCATE( MCREFSORT ( NINVC,2 ), STAT=IOS ) 
        CALL CHECKMEM( IOS, 'MCREFSORT', PROGNAME )
        MCREFSORT = 0.

C.........  Close MCREF file
        CLOSE( MDEV )
        
C.........  Abort if error found while reading cross-reference file
        IF( EFLAG ) THEN
            MESG = 'Problem reading cross-reference file'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Sort MCREF index array by reference county
        CALL SORTI2( NLINES, IDX, MCREFRAW(:,2), MCREFRAW(:,1) )

C.........  Check for duplicate entries and counties outside grid,
C           then store sorted MCREF array
        PICOUNTY = 0
        N = 0

        DO I = 1, NLINES
            J = IDX( I )
            
            REFCOUNTY = MCREFRAW( J,2 )
            INVCOUNTY = MCREFRAW( J,1 )

C.............  Skip any entries equal to zero due to blank lines
            IF( REFCOUNTY == 0 .OR. INVCOUNTY == 0 ) CYCLE

C.............  Check if current inventory county is duplicate (match previous)
            IF( INVCOUNTY == PICOUNTY ) THEN
 
                DUPFLAG = .TRUE.
                EFLAG   = .TRUE.
                
                WRITE( MESG,94010 ) 'ERROR: Duplicate entries in ' //
     &                 'county cross-reference file for ' // CRLF() //
     &                 BLANK10 // 'inventory county', INVCOUNTY,
     &                 ', reference county', REFCOUNTY, '.'
                CALL M3MESG( MESG )
            ELSE
 
C.............  Check that current county is inside the grid (and in the inventory)
                K = FIND1( INVCOUNTY, NCTY, GRIDCTY )
                
                IF( K < 0 ) THEN
                    WRITE( MESG,94010 ) 'WARNING: Ignoring county ' //
     &                     'cross-reference for inventory county',
     &                     INVCOUNTY, ',' // CRLF() // BLANK10 // 
     &                     ' since it is not inside the grid'
                    CALL M3MESG( MESG )
                    
                ELSE    
                    N = N + 1
        
                    MCREFSORT( N,1 ) = INVCOUNTY
                    MCREFSORT( N,2 ) = REFCOUNTY
                END IF
                
            END IF
            
            PICOUNTY = INVCOUNTY
            
        END DO

        IF( DUPFLAG ) THEN
            MESG = 'ERROR: Duplicate county cross-reference ' //
     &             'entries found. ' //CRLF()// BLANK10 // 
     &             'Remove duplicate entries and try again.'
            CALL M3MSG2( MESG )
        END IF

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

C.........  Deallocate local memory
        DEALLOCATE( MCREFRAW, IDX )
 
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
        
        END SUBROUTINE RDMXREF
