
        SUBROUTINE RDFMREF( ADEV )

C***********************************************************************
C  subroutine body starts at line 99
C
C  DESCRIPTION:
C       Reads FUELMONTH_COUNTY file, checks that fuel month settings are given for each 
C       reference county, ignores counties not specified in the MCREF
C       file, and sorts the data
C
C  PRECONDITIONS REQUIRED:
C       ADEV has been opened
C
C  SUBROUTINES AND FUNCTIONS CALLED:  none
C
C  REVISION  HISTORY:
C     3/10: Created by B.H. Baek
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
        USE MODMBSET, ONLY: NREFC, MCREFIDX,  NREFFLAGS,
     &                      NREFF, FMREFSORT, NFUELC, FMREFLIST
                
        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        LOGICAL        CHKINT
        LOGICAL        BLKORCMT
        INTEGER        GETFLINE
        INTEGER        STR2INT
        INTEGER        FINDC
        CHARACTER(2)   CRLF    
        
        EXTERNAL  BLKORCMT, CHKINT, GETFLINE, STR2INT, FINDC, CRLF

C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT (IN) :: ADEV     ! COUNTY_FUELMONTH file unit no.

C...........   Local allocatable arrays
        INTEGER, ALLOCATABLE :: FMREFRAW ( :,: )  ! raw COUNTY_FULEMONTH data
        INTEGER, ALLOCATABLE :: IDX ( : )    ! index into COUNTY_FUELMONTH data

C...........   Local arrays
        CHARACTER(16)  SEGMENT( NREFFLAGS )  ! parsed input line
        
C...........   Other local variables
        INTEGER I, J, K, N, NF, NR        ! counters and indices                     
        
        INTEGER IOS                       ! I/O status
        INTEGER :: IREC = 0               ! record counter
        INTEGER NLINES                    ! number of lines

        INTEGER REFCOUNTY                 ! ref. county FIPS code
        INTEGER PRMONTH                   ! prev ref. county month
        INTEGER FMONTH                    ! current / no. of fuelmonth per ref. county
        INTEGER CMONTH                    ! current / no. of month(s) per ref. county

        LOGICAL      :: DUPFLAG = .FALSE.   ! true: duplicate entries found
        LOGICAL      :: EFLAG   = .FALSE.   ! true: error found    
        LOGICAL      :: SETFLAG = .FALSE.   ! true: error in settings values 
        LOGICAL      :: RFLAG   = .FALSE.   ! true: no. lines < no. ref. counties
      
        CHARACTER(FIPLEN3) PRVCNTY, REFCNTY, INVCNTY
        CHARACTER(100)     LINE     !  line buffer
        CHARACTER(300)     MESG     !  message buffer

        CHARACTER(16) :: PROGNAME = 'RDFMREF'   ! program name
        
C***********************************************************************
C   begin body of subroutine RDFMREF

C.........  Get the number of lines in the file     
        MESG = 'Reference county fuel month settings file'
        NLINES = GETFLINE( ADEV,MESG )

C.........  Allocate memory to store settings information        
        ALLOCATE( FMREFRAW ( NLINES,NREFFLAGS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'FMREFRAW', PROGNAME )
        ALLOCATE( IDX ( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IDX', PROGNAME )

C.........  Initialize arrays
        N = 0
        IDX = 0
        FMREFRAW = 0
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
            
C.............  Skip blank lines
            IF( BLKORCMT( LINE ) ) CYCLE

C.............  Parse the line into segments
            CALL PARSLINE( LINE, NREFFLAGS, SEGMENT )

C.............  Convert reference county to integer             
            IF( CHKINT( SEGMENT( 1 ) ) ) THEN
                REFCOUNTY = STR2INT( SEGMENT( 1 ) )
            ELSE
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: Bad reference county ' //
     &                     'FIPS code at line', IREC, 'of ' //
     &                     'reference county fuel month settings file.'
                CALL M3MESG( MESG )
                CYCLE
            END IF
            
            WRITE( REFCNTY,'(I12.12)' ) REFCOUNTY
            K = FINDC( REFCNTY, NREFC, MCREFIDX( :,1 ) )

            IF( K > 0 ) N = N + 1

            FMONTH = STR2INT( SEGMENT( 2 ) )
            CMONTH = STR2INT( SEGMENT( 3 ) )

C.............  Store values in unsorted array
            FMREFRAW( I,1 ) = REFCOUNTY
            FMREFRAW( I,2 ) = FMONTH
            FMREFRAW( I,3 ) = CMONTH

        END DO  ! done reading COUNTY_FULEMONTH file
        
        NREFF = N     ! no of fuel month entries
            
C.........  Close COUNTY_FUELMONTH file
        CLOSE( ADEV )
        
C.........  Abort if error found while reading settings file
        IF( EFLAG ) THEN
            MESG = 'Problem reading reference county fuel month ' //
     &             'setting file'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        CALL SORTI2( NLINES,IDX,FMREFRAW(:,1),FMREFRAW(:,2) )

C.............  Count no of ref counties
        NR = 0
        PRVCNTY = ' '
        DO I = 1, NLINES
            J = IDX( I )
            
            WRITE( REFCNTY,'(I12.12)' ) FMREFRAW( J,1 )

C.............  Skip any entries equal to zero due to blank lines
            IF( FMREFRAW( J,1 ) == 0 ) CYCLE
            
            K = FINDC( REFCNTY, NREFC, MCREFIDX( :,1 ) )

            IF( K < 1 ) CYCLE
            IF( PRVCNTY /= REFCNTY ) NR = NR + 1

            PRVCNTY = REFCNTY

        END DO

C.........  If the no. of lines in the file is less than the no. of ref.
C              counties, then something is wrong, but we'll go through
C              the function to print out error messages
        ALLOCATE( FMREFSORT ( NREFF,NREFFLAGS ), STAT=IOS ) 
        CALL CHECKMEM( IOS, 'FMREFSORT', PROGNAME )
        ALLOCATE( FMREFLIST ( NR,2 ), STAT=IOS ) 
        CALL CHECKMEM( IOS, 'FMREFLIST', PROGNAME )
        FMREFSORT = ' '
        FMREFLIST = ' ' 

C.........  Store sorted reference county fuel month setting array
        N = 0
        NR = 0
        PRVCNTY  = ' ' 
        DO I = 1, NLINES

            J = IDX( I )
            
            WRITE( REFCNTY,'(I12.12)' ) FMREFRAW( J,1 )
            CALL PADZERO( REFCNTY )
            FMONTH = FMREFRAW( J,2 )
            CMONTH = FMREFRAW( J,3 )

C.............  Skip any entries equal to zero due to blank lines
            IF( FMREFRAW( J,1 ) == 0 ) CYCLE
            
C.............  Check that current reference county is in the county cross-reference list
            K = FINDC( REFCNTY, NREFC, MCREFIDX( :,1 ) )

            IF( K < 0 ) CYCLE 

            N = N + 1
            FMREFSORT( N,1 ) = REFCNTY
            WRITE( FMREFSORT( N,2 ),'(I8)' ) FMONTH
            WRITE( FMREFSORT( N,3 ),'(I8)' ) CMONTH

            IF( PRVCNTY == REFCNTY ) THEN
                IF( PRMONTH /= CMONTH ) THEN
                    NF = NF + 1
                    WRITE( FMREFLIST( NR,2 ),'(I8)' ) NF
                END IF
            ELSE
                NR = NR + 1
                NF = 1
                FMREFLIST( NR,1 ) = REFCNTY
            END IF
                
            PRMONTH = CMONTH
            PRVCNTY = REFCNTY
            
        END DO
        
        NFUELC = NR

        IF( EFLAG ) THEN
            MESG = 'Problem(s) found in reference county fuel month ' //
     &             'settings file.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Deallocate local memory
        DEALLOCATE( FMREFRAW, IDX )

        RETURN

999     MESG = 'End of file'
        MESG = 'End of file reached unexpectedly. ' //
     &         'Check format of COUNTY_FUELMONTH' // CRLF() // BLANK5 //
     &         'input file.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )   
             
C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )  
      
C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )
94020   FORMAT( 3( A, 1X ), I8, 1X, A, 1X )
        
        END SUBROUTINE RDFMREF
