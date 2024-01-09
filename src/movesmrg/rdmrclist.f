
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
C       MCREF file must be read already
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
        CHARACTER(2)  CRLF
        
        EXTERNAL BLKORCMT, CHKINT, GETFLINE, STR2INT, CRLF

C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT (IN) :: FDEV             ! MRCLIST file unit no.

C...........   Local allocatable arrays
        INTEGER           , ALLOCATABLE :: IDX( : )     ! sorting index
        INTEGER           , ALLOCATABLE :: REFFIPA( : ) ! unsorted ref. county FIPs
        CHARACTER(FIPLEN3), ALLOCATABLE :: REFFIP( : )  ! sorted ref. county FIPs
        CHARACTER(100)    , ALLOCATABLE :: FILESA( : )  ! unsorted files names
        CHARACTER(100)    , ALLOCATABLE :: FILES( : )   ! sorted file names
        
        INTEGER, ALLOCATABLE :: MONTHA( : )      ! unsorted month numbers
        INTEGER, ALLOCATABLE :: MONTH( : )       ! sorted month numbers

C...........   Local arrays
        CHARACTER(100)  SEGMENT( 3 )          ! parsed input line

C...........   Other local variables
        INTEGER         I, J, N     ! counters and indexes
        INTEGER         IOS         ! error status
        INTEGER      :: IREC = 0    ! record counter
        INTEGER         NLINES      ! number of lines
        INTEGER         PFIP        ! previous ref. county FIP
        INTEGER         PMONTH      ! previous fuel month
        INTEGER         TMONTH      ! tmp. fuel month
        INTEGER      :: TDEV = 0    ! tmp. file unit
        INTEGER         TIDX        ! tmp. index
        
        LOGICAL      :: EFLAG = .FALSE.   ! true: error found
        LOGICAL      :: FOUND = .FALSE.   ! true: found data for reference county

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
        ALLOCATE( MONTHA( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MONTHA', PROGNAME )
        ALLOCATE( FILESA( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'FILESA', PROGNAME )
        ALLOCATE( IDX( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IDX', PROGNAME )
        
        REFFIPA = 0
        MONTHA = 0
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

C.............  Parse the line into 3 segments
            CALL PARSLINE( LINE, 3, SEGMENT )

C.............  Convert reference county to integer
            IF( .NOT. CHKINT( SEGMENT( 1 ) ) ) THEN
                EFLAG = .TRUE.
                WRITE( MESG, 94010 ) 'ERROR: Bad reference county ' //
     &            'FIPS code at line', IREC, 'of county factors ' //
     &            'file list.'
                CALL M3MESG( MESG )
                CYCLE
            END IF
            
            REFFIPA( I ) = STR2INT( SEGMENT( 1 ) )

C.............  Convert month to integer
            IF( .NOT. CHKINT( SEGMENT( 2 ) ) ) THEN
                EFLAG = .TRUE.
                WRITE( MESG, 94010 ) 'ERROR: Bad month number ' //
     &            'at line', IREC, 'of county factors file list.'
                CALL M3MESG( MESG )
                CYCLE
            END IF
            
            TMONTH = STR2INT( SEGMENT( 2 ) )
            
            IF( TMONTH .LT. 1 .OR. TMONTH .GT. 12 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG, 94010 ) 'ERROR: Invalid fuel month ' //
     &            'at line', IREC, 'of county factors file list.'
                CALL M3MESG( MESG )
                CYCLE
            END IF
            
            MONTHA( I ) = TMONTH
            FILESA( I ) = SEGMENT( 3 )

C.............  Check that file can be opened
            FILENAME = TRIM( MVFILDIR ) // TRIM( FILESA( I ) )
            OPEN( TDEV, FILE=FILENAME, STATUS='OLD', IOSTAT=IOS )
            IF( IOS .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG, 94010 ) 'ERROR: Could not open file ' //
     &            TRIM( FILESA( I ) ) // ' at line', IREC, 'of county ' //
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

C.........  Sort list by reference county and fuel month
        CALL SORTI2( NLINES, IDX, REFFIPA, MONTHA )

C.........  Check for duplicate entries and store sorted filenames
        ALLOCATE( REFFIP( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'REFFIP', PROGNAME )
        ALLOCATE( MONTH( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MONTH', PROGNAME )
        ALLOCATE( FILES( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'FILES', PROGNAME )
        
        REFFIP = ' '
        MONTH = 0
        FILES = ' '

        PFIP = -9
        PMONTH = -9
        N = 0
        DO I = 1, NLINES
            J = IDX( I )
            
            IF( REFFIPA( J ) == PFIP .AND. MONTHA( J ) == PMONTH ) THEN
                EFLAG = .TRUE.
                WRITE( MESG, 94010 ) 'ERROR: Duplicate entries in ' //
     &            'reference county factors list for ' // CRLF() //
     &            BLANK10 // 'reference county', PFIP, 'and fuel ' //
     &            'month', PMONTH
                CALL M3MESG( MESG )
                CYCLE
            END IF

            N = N + 1
            WRITE( REFFIP( N ),'(I12.12)' ) REFFIPA( J )
            MONTH( N ) = MONTHA( J )
            FILES( N ) = FILESA( J )
        
            PFIP = REFFIPA( J )
            PMONTH = MONTHA( J )
            
        END DO
        
        IF( EFLAG ) THEN
            MESG = 'Problem found in reference county factors list.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        ALLOCATE( MRCLIST( NREFC,12 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MRCLIST', PROGNAME )
        MRCLIST = ' '

C.........  Check that each reference county has a factors file
        TIDX = 1
        DO I = 1, NREFC

C.............  Loop through sorted lines in MRCLIST file - the
C               FIPS codes are sorted in the same order as the 
C               reference counties in MCREFIDX
            FOUND = .FALSE.
            DO
                IF( REFFIP( TIDX ) == MCREFIDX( I,1 ) ) THEN
                    FOUND = .TRUE.
                    MRCLIST( I, MONTH( TIDX ) ) = FILES( TIDX )
                ELSE
                    IF( FOUND ) EXIT
                END IF

                TIDX = TIDX + 1
                IF( TIDX .GT. NLINES ) EXIT
            END DO

            IF( .NOT. FOUND ) THEN
                EFLAG = .TRUE.
c               WRITE( MESG, 94010 ) 'ERROR: No factor file found ' //
c    &            'for reference county', MCREFIDX( I,1 ), 'in ' //
c    &            'reference county factors list.'
                WRITE( MESG, '(A)' ) 'ERROR: No factor file found ' // ! UNC-IE: Jan 2024; correcting format follow Carlie's sugestion
     &            'for reference county ' // MCREFIDX( I,1 ) //
     &            ' in reference county factors list.'
                CALL M3MESG( MESG )
                CYCLE
            END IF
        
        END DO
        
        IF( EFLAG ) THEN
            MESG = 'Problem found in reference county factors list.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF
        
        DEALLOCATE( REFFIPA, MONTHA, FILESA, IDX, REFFIP, MONTH, FILES )

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
