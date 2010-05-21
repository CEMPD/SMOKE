
        SUBROUTINE RDMETMOVES( FDEV )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C       Reads the Met4moves output file (minimum and maximum temperatures
C       for each county). Checks that each inventory county has an entry.
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
C.........  This module contains data structures and flags specific to Movesmrg
        USE MODMVSMRG, ONLY: AVGMIN, AVGMAX

C.........  This module contains the lists of unique source characteristics
        USE MODLISTS, ONLY: NINVIFIP, INVIFIP

C.........  This module is used for reference county information
        USE MODMBSET, ONLY: NINVC, MCREFSORT

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        LOGICAL       BLKORCMT
        LOGICAL       CHKINT
        LOGICAL       CHKREAL
        INTEGER       GETFLINE
        INTEGER       FIND1
        INTEGER       STR2INT
        REAL          STR2REAL
        CHARACTER(2)  CRLF
        
        EXTERNAL BLKORCMT, CHKINT, CHKREAL, FIND1, GETFLINE, 
     &           STR2INT, STR2REAL, CRLF

C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT (IN) :: FDEV             ! Met4moves file unit no.

C...........   Local allocatable arrays

C...........   Local arrays
        CHARACTER(100)  SEGMENT( 7 )          ! parsed input line

C...........   Other local variables
        INTEGER         I, K        ! counters and indexes
        INTEGER         IOS         ! error status
        INTEGER      :: IREC = 0    ! record counter
        INTEGER         NLINES      ! number of lines
        INTEGER         CNTY        ! current FIPS code
        
        REAL            MINVAL      ! minimum temperature value
        REAL            MAXVAL      ! maximum temperature value

        LOGICAL      :: EFLAG = .FALSE.   ! true: error found

        CHARACTER(150)     LINE     ! line buffer
        CHARACTER(300)     MESG     ! message buffer

        CHARACTER(16) :: PROGNAME = 'RDMETMOVES'   ! program name

C***********************************************************************
C   begin body of subroutine RDMETMOVES

C.........  Allocate memory to store min and max temperatures
        ALLOCATE( AVGMIN( NINVIFIP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'AVGMIN', PROGNAME )
        ALLOCATE( AVGMAX( NINVIFIP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'AVGMAX', PROGNAME )
        AVGMIN = BADVAL3  ! array
        AVGMAX = BADVAL3  ! array

C.........  Get the number of lines in the file
        NLINES = GETFLINE( FDEV, 'Met4moves output file' )

C.........  Read through file and match to FIPS from inventory
        DO I = 1, NLINES
        
            READ( FDEV, 93000, END=999, IOSTAT=IOS ) LINE
            
            IREC = IREC + 1
            
            IF( IOS .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG, 94010 ) 'I/O error', IOS,
     &            'reading Met4moves output file at line', IREC
                CALL M3MESG( MESG )
                CYCLE
            END IF

C.............  Skip blank or comment lines
            IF( BLKORCMT( LINE ) ) CYCLE

C.............  Parse the line into 7 segments
            CALL PARSLINE( LINE, 7, SEGMENT )

C.............  Convert county to integer
            IF( .NOT. CHKINT( SEGMENT( 1 ) ) ) THEN
                EFLAG = .TRUE.
                WRITE( MESG, 94010 ) 'ERROR: Bad FIPS code ' //
     &            'at line', IREC, 'of Met4moves output file.'
                CALL M3MESG( MESG )
                CYCLE
            END IF

C.............  Find county in inventory list
            CNTY = STR2INT( ADJUSTR( SEGMENT( 1 ) ) )
            K = FIND1( CNTY, NINVIFIP, INVIFIP )
            
            IF( K .LE. 0 ) THEN
                WRITE( MESG, 94010 ) 'NOTE: Skipping line',
     &            IREC, 'of Met4moves output file because FIPS code',
     &            CNTY, 'is not in the inventory.'
                CALL M3MESG( MESG )
                CYCLE
            END IF
            
C.............  Check min and max temperature values
            IF( .NOT. CHKREAL( SEGMENT( 6 ) ) ) THEN
                EFLAG = .TRUE.
                WRITE( MESG, 94010 ) 'ERROR: Bad minimum ' //
     &            'temperature value at line', IREC, 
     &            'of Met4moves output file.'
                CALL M3MESG( MESG )
                CYCLE
            END IF
            
            IF( .NOT. CHKREAL( SEGMENT( 7 ) ) ) THEN
                EFLAG = .TRUE.
                WRITE( MESG, 94010 ) 'ERROR: Bad maximum ' //
     &            'temperature value at line', IREC,
     &            'of Met4moves output file.'
                CALL M3MESG( MESG )
                CYCLE
            END IF
            
            MINVAL = STR2REAL( ADJUSTR( SEGMENT( 6 ) ) )
            MAXVAL = STR2REAL( ADJUSTR( SEGMENT( 7 ) ) )
            
            IF( MINVAL .GT. MAXVAL ) THEN
                EFLAG = .TRUE.
                WRITE( MESG, 94010 ) 'ERROR: Minimum temperature ' //
     &            'is greater than maximum temperature at line',
     &            IREC, 'of Met4moves output file.'
                CALL M3MESG( MESG )
                CYCLE
            END IF
            
C.............  Check for duplicate entries
            IF( AVGMIN( K ) .GT. AMISS3 .OR.
     &          AVGMAX( K ) .GT. AMISS3 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG, 94010 ) 'ERROR: Duplicate county',
     &            CNTY, 'at line', IREC, 'of Met4moves output file.'
                CALL M3MESG( MESG )
                CYCLE
            END IF
            
            AVGMIN( K ) = MINVAL
            AVGMAX( K ) = MAXVAL

        END DO

        CLOSE( FDEV )
        
        IF( EFLAG ) THEN
            MESG = 'Problem found in Met4moves output file.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF
        
C.........  Check that every county in the domain has a min and max temp
        DO I = 1, NINVC

            K = FIND1( MCREFSORT( I,1 ), NINVIFIP, INVIFIP )
        
            IF( AVGMIN( K ) .LT. AMISS3 .OR.
     &          AVGMAX( K ) .LT. AMISS3 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG, 94010 ) 'ERROR: Missing temperature ' //
     &            'data for county', INVIFIP( K ), 'in ' //
     &            'Met4moves output file.'
                CALL M3MESG( MESG )
                CYCLE
            END IF
        
        END DO
        
        IF( EFLAG ) THEN
            MESG = 'Missing data in Met4moves output file.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        RETURN

999     MESG = 'End of file'
        MESG = 'End of file reached unexpectedly. ' //
     &         'Check format of METMOVES' // CRLF() // BLANK5 //
     &         'input file.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )   

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )  
      
C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )
        
        END SUBROUTINE RDMETMOVES
