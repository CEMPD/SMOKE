
        SUBROUTINE RDMRCLIST( FDEV )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C       Reads the MRCLIST file (list of emission factor files for each 
C       reference county). Checks that each reference county has a
C       file and that all files can be opened. Sorts files by reference
C       county.
C
C  PRECONDITIONS REQUIRED:
C       FDEV must be opened
C
C  SUBROUTINES AND FUNCTIONS CALLED:  none
C
C  REVISION  HISTORY:
C     04/10: Created by C. Seppanen
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
C.........  This module is used for reference county information
        USE MODMBSET, ONLY: NREFC, MCREFIDX

C.........  This module contains data structures and flags specific to Movesmrg
        USE MODMVSMRG, ONLY: MRCLIST, MVFILDIR

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        LOGICAL       BLKORCMT
        LOGICAL       CHKINT
        INTEGER       GETFLINE
        INTEGER       STR2INT
        INTEGER       INDEXINT1
        CHARACTER(2)  CRLF
        
        EXTERNAL BLKORCMT, CHKINT, GETFLINE, STR2INT, INDEXINT1, CRLF

C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT (IN) :: FDEV             ! MRCLIST file unit no.

C...........   Local allocatable arrays
        INTEGER, ALLOCATABLE :: REFFIPA( : )     ! unsorted ref. county FIPs
        INTEGER, ALLOCATABLE :: REFFIP( : )      ! sorted ref. county FIPs
        
        INTEGER, ALLOCATABLE :: IDX( : )         ! sorting index

        CHARACTER(100), ALLOCATABLE :: FILESA( : )  ! unsorted files names
        CHARACTER(100), ALLOCATABLE :: FILES( : )   ! sorted file names

C...........   Local arrays
        CHARACTER(100)  SEGMENT( 3 )          ! parsed input line

C...........   Other local variables
        INTEGER         I, J, N     ! counters and indexes
        INTEGER         IOS         ! error status
        INTEGER      :: IREC = 0    ! record counter
        INTEGER         NLINES      ! number of lines
        INTEGER         PFIP        ! previous ref. county FIP
        INTEGER         TDEV        ! tmp. file unit
        
        LOGICAL      :: EFLAG = .FALSE.   ! true: error found

        CHARACTER(150)     LINE     ! line buffer
        CHARACTER(200)     FILENAME ! tmp. filename
        CHARACTER(300)     MESG     ! message buffer

        CHARACTER(16) :: PROGNAME = 'RDMRCLIST'   ! program name

C***********************************************************************
C   begin body of subroutine RDMRCLIST

C.........  Get the number of lines in the file
        NLINES = GETFLINE( FDEV, 'Reference county factors file list' )
        
        ALLOCATE( REFFIPA( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'REFFIPA', PROGNAME )
        ALLOCATE( FILESA( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'FILESA', PROGNAME )
        ALLOCATE( IDX( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IDX', PROGNAME )
        
        REFFIPA = 0
        FILESA = ' '
        IDX = 0
        
        DO I = 1, NLINES
        
            IDX( I ) = I
        
C.............  Read line
            READ( FDEV, 93000, END=999, IOSTAT=IOS ) LINE
            
            IREC = IREC + 1
            
            IF( IOS .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG, 94010 ) 'I/O error', IOS,
     &            'reading reference county factors file list at line',
     &            IREC
                CALL M3MESG( MESG )
                CYCLE
            END IF

C.............  Skip blank or comment lines
            IF( BLKORCMT( LINE ) ) CYCLE

C.............  Parse the line into 2 segments
            CALL PARSLINE( LINE, 2, SEGMENT )

C.............  Convert reference county to integer
            IF( .NOT. CHKINT( SEGMENT( 1 ) ) ) THEN
                EFLAG = .TRUE.
                WRITE( MESG, 94010 ) 'ERROR: Bad reference county ' //
     &            'FIPS code at line', IREC, 'of county factors ' //
     &            'file list.'
                CALL M3MESG( MESG )
                CYCLE
            END IF
            
            REFFIPA( I ) = STR2INT( ADJUSTR( SEGMENT( 1 ) ) )
            FILESA( I ) = SEGMENT( 2 )

C.............  Check that file can be opened
            FILENAME = TRIM( MVFILDIR ) // TRIM( FILESA( I ) )
            OPEN( TDEV, FILE=FILENAME, STATUS='OLD', IOSTAT=IOS )
            IF( IOS .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG, 94010 ) 'ERROR: Could not open file ' //
     &            FILESA( I ) // ' at line', IREC, 'of county ' //
     &            'factors file list.'
                CALL M3MESG( MESG )
            END IF
            CLOSE( TDEV )

        END DO
        
        CLOSE( FDEV )

C.........  Check for errors while reading file
        IF( EFLAG ) THEN
            MESG = 'Problem reading reference county factors list'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Sort list by reference county
        CALL SORTI1( NLINES, IDX, REFFIPA )

C.........  Check for duplicate entries and store sorted filenames
        ALLOCATE( REFFIP( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'REFFIPA', PROGNAME )
        ALLOCATE( FILES( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'FILESA', PROGNAME )
        
        REFFIP = 0
        FILES = ' '

        PFIP = -9
        N = 0
        DO I = 1, NLINES
            J = IDX( I )
            
            IF( REFFIPA( J ) == PFIP ) THEN
                EFLAG = .TRUE.
                WRITE( MESG, 94010 ) 'ERROR: Duplicate entries in ' //
     &            'reference county factors list for ' // CRLF() //
     &            BLANK10 // 'reference county', PFIP
                CALL M3MESG( MESG )
                CYCLE
            END IF

            N = N + 1
            REFFIP( N ) = REFFIPA( J )
            FILES( N ) = FILESA( J )
        
            PFIP = REFFIPA( J )
            
        END DO
        
        IF( EFLAG ) THEN
            MESG = 'Problem found in reference county factors list.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        ALLOCATE( MRCLIST( NREFC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MRCLIST', PROGNAME )
        
        MRCLIST = ' '

C.........  Check that each reference county has a factors file
        DO I = 1, NREFC
        
            J = INDEXINT1( MCREFIDX( I,1 ), NLINES, REFFIP )
            IF( J <= 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG, 94010 ) 'ERROR: No factor file found ' //
     &            'for reference county', MCREFIDX( I,1 ), 'in ' //
     &            'reference county factors list.'
                CALL M3MESG( MESG )
                CYCLE
            END IF
            
            MRCLIST( I ) = FILES( J )
        
        END DO
        
        IF( EFLAG ) THEN
            MESG = 'Problem found in reference county factors list.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF
        
        DEALLOCATE( REFFIPA, FILESA, IDX, REFFIP, FILES )

        RETURN

999     MESG = 'End of file'
        MESG = 'End of file reached unexpectedly. ' //
     &         'Check format of MRCLIST' // CRLF() // BLANK5 //
     &         'input file.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )   

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )  
      
C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )
        
        END SUBROUTINE RDMRCLIST
