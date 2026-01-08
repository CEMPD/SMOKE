
        SUBROUTINE RDMEDSINFO

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C       Reads MEDS column/row & GAI lookup tables for assigning 
C       Lat/Lon and 12-digit CO-ABS-DIST code 
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:  none
C
C  REVISION  HISTORY:
C     10/13: Created by BH Baek
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
C       Updated with USE M3UTILIO by Huy Tran UNC-IE on 2026-01
C***********************************************************************

C.........  MODULES for public variables
        USE M3UTILIO

        USE MODSOURC, ONLY: NMEDGRD, CMEDGRD, NMEDGAI, COABDST

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
C        INCLUDE 'PARMS3.EXT'    !  I/O API parameters


C...........   EXTERNAL FUNCTIONS and their descriptions:
        LOGICAL       BLKORCMT
        INTEGER       GETFLINE
C       INTEGER       PROMPTFFILE
C       INTEGER       STR2INT
C       REAL          STR2REAL
        
C        EXTERNAL BLKORCMT, PROMPTFFILE, GETFLINE, STR2INT, STR2REAL
        EXTERNAL     BLKORCMT, GETFLINE

C...........   Local allocatable arrays

C...........   Local arrays
        CHARACTER( 20 )  SEGMENT( 6 )          ! parsed input line

C...........   Other local variables
        INTEGER         I, J, K, L, N    ! counters and indexes
        INTEGER         IOS         ! error status
        INTEGER      :: IREC = 0    ! record counter
        INTEGER         IDEV, JDEV  
        INTEGER         NLINES      ! number of lines
        INTEGER         ROW, COL, FIPS
        
        CHARACTER( 3 )     ARBN, DSTR
        CHARACTER(500)     LINE     ! line buffer
        CHARACTER(300)     MESG     ! message buffer

        CHARACTER(16) :: PROGNAME = 'RDMEDSINFO'   ! program name

C***********************************************************************
C   begin body of subroutine RDMEDSINFO

C.........  Open Row/col and GAI lookup tables
        IDEV = PROMPTFFILE(
     &      'Enter logical name for input ROWCOL_LATLON file',
     &      .TRUE., .TRUE., 'ROWCOL_LATLON', PROGNAME )

        JDEV = PROMPTFFILE(
     &      'Enter logical name for input GAI lookup input file',
     &      .TRUE., .TRUE., 'GAI_LOOKUP_TABLE', PROGNAME )

C.........  Count no of ROWCOL_LATLON entries
        NLINES = GETFLINE( IDEV, 'Reading ROWCOL_LATLON input'
     &      // ' file for MEDS gridded inventory data processing' )
        IREC = 0
        NMEDGRD = 0
        DO I = 1, NLINES

            READ( IDEV, 93000, IOSTAT=IOS ) LINE
            IREC = IREC + 1

            IF( IOS .NE. 0 ) THEN
                WRITE( MESG, 94010 ) 'I/O error', IOS,
     &            'reading control factor file at line', IREC
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  Skip blank or comment lines
            IF( BLKORCMT( LINE ) ) CYCLE

C.............  Parse line into fields
            NMEDGRD = NMEDGRD + 1

        END DO

        REWIND( IDEV )

C.........  Allocate ARRAYS
        ALLOCATE( CMEDGRD( NMEDGRD,3 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CMEDGRD', PROGNAME )
        CMEDGRD = ''

        IREC = 0
        N = 0
        DO I = 1, NLINES

            READ( IDEV, 93000, IOSTAT=IOS ) LINE
            IREC = IREC + 1

            IF( IOS .NE. 0 ) THEN
                WRITE( MESG, 94010 ) 'I/O error', IOS,
     &            'reading control factor file at line', IREC
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  Skip blank or comment lines
            IF( BLKORCMT( LINE ) ) CYCLE

C.............  Parse line into fields
            CALL PARSLINE( LINE, 4, SEGMENT )

C.............  Store data
            N = N + 1

            COL = STR2INT( SEGMENT( 1 ) )
            ROW = STR2INT( SEGMENT( 2 ) )

            WRITE( CMEDGRD( N,1 ),'( 2I3.3 )' ) COL, ROW   ! store col_row
            CMEDGRD( N,2 ) = ADJUSTL( SEGMENT( 3 ) )       ! store longitude
            CMEDGRD( N,3 ) = ADJUSTL( SEGMENT( 4 ) )       ! store latitude

        END DO
        
        CLOSE( IDEV )

C.........  Count no of GAI entries
        NLINES = GETFLINE( JDEV, 'Reading GAI_LOOKUP_TABLE input'
     &      // ' file for MEDS gridded inventory data processing' )
        IREC = 0
        NMEDGAI = 0
        DO I = 1, NLINES

            READ( JDEV, 93000, IOSTAT=IOS ) LINE
            IREC = IREC + 1

            IF( IOS .NE. 0 ) THEN
                WRITE( MESG, 94010 ) 'I/O error', IOS,
     &            'reading control factor file at line', IREC
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  Skip blank or comment lines
            IF( BLKORCMT( LINE ) ) CYCLE

C.............  Parse line into fields
            NMEDGAI = NMEDGAI + 1

        END DO

        REWIND( JDEV )

C.........  Allocate ARRAYS
        ALLOCATE( COABDST( NMEDGAI,2 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'COABDST', PROGNAME )
        COABDST = ''

        IREC = 0
        N = 0
        DO I = 1, NLINES

            READ( JDEV, 93000, IOSTAT=IOS ) LINE
            IREC = IREC + 1

            IF( IOS .NE. 0 ) THEN
                WRITE( MESG, 94010 ) 'I/O error', IOS,
     &            'reading control factor file at line', IREC
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  Skip blank or comment lines
            IF( BLKORCMT( LINE ) ) CYCLE

C.............  Parse line into fields
            CALL PARSLINE( LINE, 4, SEGMENT )

C.............  Store data
            N = N + 1

            COABDST( N,1 ) = ADJUSTL( SEGMENT( 1 ) )

            ARBN = ADJUSTR( TRIM( SEGMENT(2) ) )
            FIPS = STR2INT( SEGMENT(3) )
            DSTR = ADJUSTR( TRIM( SEGMENT(4) ) )
            CALL PADZERO( ARBN )
            CALL PADZERO( DSTR )
            
            WRITE( COABDST( N,2 ),'(A3,I6.6,A3)' ) ARBN, FIPS, DSTR 

        END DO

        CLOSE( JDEV )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )  
94010   FORMAT( 10( A, :, I8, :, 1X ) )
        
        END SUBROUTINE RDMEDSINFO

