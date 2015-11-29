
        SUBROUTINE RDSCCMAP( ADEV )

C***********************************************************************
C  subroutine body starts at line 99
C
C  DESCRIPTION:
C       Reads SCC_MAP file that maps aggregated SCC to full SCCs.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:  none
C
C  REVISION  HISTORY:
C     5/14: Created by B.H. Baek
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
        USE MODMOBIL, ONLY: NSCCMAP, SCCMAPLIST
                
        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        LOGICAL        BLKORCMT
        CHARACTER(2)   CRLF    
        INTEGER        GETFLINE
        
        EXTERNAL  BLKORCMT, CRLF, GETFLINE

C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT (IN) :: ADEV     ! COUNTY_FUELMONTH file unit no.

C...........   Local allocatable arrays
        INTEGER, ALLOCATABLE :: IDX ( : )          ! index of SCCs
        CHARACTER( SCCLEN3 ), ALLOCATABLE :: SCCMAPRAW ( :,: )  ! raw mapped SCCs

C...........   Local arrays
        CHARACTER( SCCLEN3 )  SEGMENT( 3 )  ! parsed input line
        
C...........   Other local variables
        INTEGER I, J, N, NS               ! counters and indices                     
        
        INTEGER    IOS                    ! I/O status
        INTEGER :: IREC = 0               ! record counter
        INTEGER :: NLINES = 0             ! number of lines
        INTEGER :: NSCC, NFSCC            ! no of SCCs

        LOGICAL      :: EFLAG   = .FALSE.   ! true: error found    

        CHARACTER(100)       LINE     !  line buffer
        CHARACTER(300)       MESG     !  message buffer
        CHARACTER( 8 )       CNFSCC   !  tmp buffer
        CHARACTER( SCCLEN3 ) CURSCC, PRVSCC, FULLSCC   ! current, previous, full SCCs

        CHARACTER(16) :: PROGNAME = 'RDSCCMAP'   ! program name

C***********************************************************************
C   begin body of subroutine RDSCCMAP

C.........  Get the number of lines in the file     
        MESG = 'Aggregated SCC to full SCC map input'
        NLINES = GETFLINE( ADEV,MESG )

C.........  Allocate memory to store settings information        
        ALLOCATE( IDX ( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IDX', PROGNAME )
        ALLOCATE( SCCMAPRAW ( NLINES,2 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SCCMAPRAW', PROGNAME )

C.........  Initialize arrays
        IDX = 0
        NSCC = 0
        SCCMAPRAW = ' '
        DO I = 1, NLINES       

            IDX( I ) = I

C.........  Read line
            READ( ADEV, 93000, END=999, IOSTAT=IOS ) LINE
            
            IREC = IREC + 1
            
            IF ( IOS /= 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG, 94010 )
     &              'I/O error', IOS, 'reading reference county ' //
     &              'fuel month settings file at line', IREC
                CALL M3MESG( MESG )
                CYCLE
            END IF
            
C.............  Skip blank/comment lines
            IF( BLKORCMT( LINE ) ) CYCLE

C.............  Parse the line into segments
            CALL PARSLINE( LINE, 2, SEGMENT )            

C.............  Store values in unsorted array
            NSCC = NSCC + 1
            SCCMAPRAW( I,1 ) = SEGMENT( 2 )    ! referenced SCCs
            SCCMAPRAW( I,2 ) = SEGMENT( 1 )    ! full SCCs

            CALL PADZERO( SCCMAPRAW( I,1 ) )
            CALL PADZERO( SCCMAPRAW( I,2 ) )

        END DO  ! done reading SCC map input file

        CLOSE( ADEV )

C.........  Abort if error found while reading settings file
        IF( EFLAG ) THEN
            MESG = 'Problem reading SCC mapping input file'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  If the no. of lines in the file is less than the no. of ref.
C              counties, then something is wrong, but we'll go through
C              the function to print out error messages
        ALLOCATE( SCCMAPLIST ( NSCC,3 ), STAT=IOS ) 
        CALL CHECKMEM( IOS, 'FMREFLIST', PROGNAME )
        SCCMAPLIST = ''

        CALL SORTIC( NLINES, IDX, SCCMAPRAW( :,1 ) )    ! sort ref SCCs

C.........  Count no of aggregated SCCs
        PRVSCC = ''
        NSCC   = 0
        DO I = 1, NLINES

            J = IDX( I )
            CURSCC  = SCCMAPRAW( J,1 )
            FULLSCC = SCCMAPRAW( J,2 )

            IF( CURSCC == ' ' ) CYCLE

            NSCC = NSCC + 1

            IF( PRVSCC /= CURSCC ) THEN
                NFSCC = 0
                NS = NSCC
            ELSE
                NFSCC = NFSCC + 1
            END IF
            
            WRITE( CNFSCC,'(I8)' ) NFSCC
            SCCMAPLIST(    NSCC,1 ) = CURSCC
            SCCMAPLIST(    NSCC,2 ) = FULLSCC
            SCCMAPLIST( NS:NSCC,3 ) = CNFSCC

            PRVSCC = CURSCC

        END DO

        NSCCMAP = NSCC

C.........  Deallocate local memory
        DEALLOCATE( IDX, SCCMAPRAW )

        RETURN

999     MESG = 'End of file reached unexpectedly. ' //
     &         'Check format of COUNTY_FUELMONTH' // CRLF() // BLANK5 //
     &         'input file.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )   
             
C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx
93000   FORMAT( A )  
      
C...........   Internal buffering formats............ 94xxx
94010   FORMAT( 10( A, :, I8, :, 1X ) )
        
        END SUBROUTINE RDSCCMAP
