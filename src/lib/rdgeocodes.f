
        SUBROUTINE RDGEOCODES( NCODES, INCODES )

C***********************************************************************
C
C  DESCRIPTION:
C       This subroutine opens and reads the four files that contain
C       the codes and names of the geographic regions.
C
C  REVISION HISTORY:
C       Created 12/13 by C. Seppanen
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
C.........  This module contains the arrays for state and county summaries
        USE MODSTCY, ONLY: NGEOLEV1, NCOUNTRY, NSTATE, NCOUNTY,
     &                     GEOLEV1COD, CTRYCOD, STATCOD, CNTYCOD,
     &                     GEOLEV1NAM, CTRYNAM, STATNAM, CNTYNAM,
     &                     CNTYTZON, CNTYTZNM, USEDAYLT

        IMPLICIT NONE

C...........   INCLUDES:
        INCLUDE 'EMCNST3.EXT'  ! emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        LOGICAL      BLKORCMT
        INTEGER      ENVINT
        INTEGER      GETFLINE
        INTEGER      INDEX1
        INTEGER      PROMPTFFILE
        
        EXTERNAL     BLKORCMT, ENVINT, GETFLINE, INDEX1, PROMPTFFILE

C...........   Subroutine arguments
        INTEGER,            INTENT(IN) :: NCODES            ! num. of input codes
        CHARACTER(FIPLEN3), INTENT(IN) :: INCODES( NCODES ) ! input codes array or empty

C...........   Local parameters
        INTEGER, PARAMETER :: MXCOL = 4  ! max. number of columns in file

C...........   Local arrays
        CHARACTER(20)      :: SEGMENT( MXCOL )  ! segments of parsed line

C...........   Local allocatable arrays
        INTEGER,            ALLOCATABLE :: INDEXA( : )   ! sorting index
        CHARACTER(FIPLEN3), ALLOCATABLE :: CODEA( : )    ! unsorted codes
        CHARACTER(20),      ALLOCATABLE :: NAMEA( : )    ! unsorted names
        INTEGER,            ALLOCATABLE :: TZONEA( : )   ! unsorted time zones
        CHARACTER(3),       ALLOCATABLE :: TZNAMA( : )   ! unsorted time zone names
        LOGICAL,            ALLOCATABLE :: USEDYA( : )   ! unsorted daylight flag

C...........   Other local variables
        INTEGER          I, J, K                 ! counters and indices
        INTEGER          CODELEN                 ! length of current code
        INTEGER          FILENUM                 ! current file number
        INTEGER          IOS                     ! i/o status
        INTEGER          LDEV                    ! file unit number
        INTEGER          NLINES                  ! number of lines in current file
        INTEGER          TZONE0                  ! default time zone
        
        LOGICAL ::       EFLAG = .FALSE.         ! error flag

        CHARACTER(16)    ENVVAR                  ! env. variable name for current file
        CHARACTER(100)   FILEDESC                ! description of current file
        CHARACTER(300)   LINE                    ! line buffer
        CHARACTER(300)   MESG                    ! message buffer
        CHARACTER(FIPLEN3) PREVCODE              ! previous code

        CHARACTER(16) :: PROGNAME = 'RDGEOCODES' ! program name

C***********************************************************************
C   begin body of subroutine RDGEOCODES

C.........  Get default time zone
        MESG = 'Default time zone for sources'
        TZONE0 = ENVINT( 'SMK_DEFAULT_TZONE', MESG, 5, IOS )

        MESG = 'Reading geographic region names and time zones...'
        CALL M3MSG2( MESG )

        DO FILENUM = 1, 4
        
            WRITE( FILEDESC, '(A,1X,I1,1X,A)' ) 'geographic level', FILENUM, 'codes'

C.............  Open code file
            MESG = 'Enter logical name for ' // TRIM( FILEDESC ) // ' file'
            WRITE( ENVVAR, '(A,I1)' ) 'GEOCODE_LEVEL', FILENUM
            LDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., ENVVAR, PROGNAME )

            NLINES = GETFLINE( LDEV, FILEDESC )
    
            ALLOCATE( INDEXA( NLINES ), STAT=IOS )
            CALL CHECKMEM( IOS, 'INDEXA', PROGNAME )
            ALLOCATE( CODEA( NLINES ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CODEA', PROGNAME )
            ALLOCATE( NAMEA( NLINES ), STAT=IOS )
            CALL CHECKMEM( IOS, 'NAMEA', PROGNAME )
            ALLOCATE( TZONEA( NLINES ), STAT=IOS )
            CALL CHECKMEM( IOS, 'TZONEA', PROGNAME )
            ALLOCATE( TZNAMA( NLINES ), STAT=IOS )
            CALL CHECKMEM( IOS, 'TZNAMA', PROGNAME )
            ALLOCATE( USEDYA( NLINES ), STAT=IOS )
            CALL CHECKMEM( IOS, 'USEDYA', PROGNAME )
        
            K = 0
            DO I = 1, NLINES
            
                READ( LDEV, '(A)', IOSTAT=IOS ) LINE
    
C.................  Check for I/O errors
                IF( IOS > 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG, 94010 ) 
     &                'I/O error', IOS, 'reading', TRIM( FILEDESC ), 'at line', I
                    CALL M3MESG( MESG )
                    CYCLE
                END IF

C.................  Check for end of file
                IF( IOS < 0 ) THEN
                    MESG = 'End of file reached unexpectedly. ' //
     &                     'Check format of ' // TRIM( FILEDESC ) // ' file.'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF

C.................  Skip blank and comment lines
                IF( BLKORCMT( LINE ) ) CYCLE

C.................  Parse line into segments
                CALL PARSLINE( LINE, MXCOL, SEGMENT )

C.................  Store fields in unsorted arrays
                K = K + 1
                INDEXA( K ) = K
                
                CODEA( K ) = ADJUSTL( SEGMENT( 1 )( 1:FIPLEN3 ) )

C.................  Fill remaining places with zeroes
                CODELEN = FILENUM * 3
                IF( CODELEN < FIPLEN3 ) THEN
                    CODEA( K )( CODELEN+1:FIPLEN3 ) = REPEAT( '0', FIPLEN3-CODELEN )
                END IF

                NAMEA( K ) = ADJUSTL( SEGMENT( 2 )( 1:20 ) )

C.................  Read time zone from level 4 file
                IF( FILENUM == 4 ) THEN
                    TZNAMA( K ) = ADJUSTL( SEGMENT( 3 )( 1:3 ) )

C.....................  Find time zone name in master list                    
                    J = INDEX1( TZNAMA( K ), MXTZONE, TZONNAM )

                    IF( J > 0 ) THEN
                        TZONEA( K ) = TZONNUM( J )

C.....................  Use default time zone if name is not found
                    ELSE
                        WRITE( MESG,94010 ) 
     &                    'WARNING: Applying default time zone', TZONE0,
     &                    'to geographic code: ' // CODEA( K )
                        CALL M3MESG( MESG )

                        TZONEA( K ) = TZONE0
                    END IF
                    
                    USEDYA( K ) = .TRUE.
                    IF( SEGMENT( 4 ) /= '' ) USEDYA( K ) = .FALSE.
                END IF

            END DO
        
            IF( EFLAG ) THEN
                MESG = 'Problem reading ' // TRIM( FILEDESC ) // ' file.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  Store number of codes and allocate memory for sorted data
            SELECT CASE( FILENUM )
            CASE( 1 )
                NGEOLEV1 = K
                ALLOCATE( GEOLEV1COD( NGEOLEV1 ), STAT=IOS )
                CALL CHECKMEM( IOS, 'GEOLEV1COD', PROGNAME )
                ALLOCATE( GEOLEV1NAM( NGEOLEV1 ), STAT=IOS )
                CALL CHECKMEM( IOS, 'GEOLEV1NAM', PROGNAME )
            CASE( 2 )
                NCOUNTRY = K
                ALLOCATE( CTRYCOD( NCOUNTRY ), STAT=IOS )
                CALL CHECKMEM( IOS, 'CTRYCOD', PROGNAME )
                ALLOCATE( CTRYNAM( NCOUNTRY ), STAT=IOS )
                CALL CHECKMEM( IOS, 'CTRYNAM', PROGNAME )
            CASE( 3 )
                NSTATE = K
                ALLOCATE( STATCOD( NSTATE ), STAT=IOS )
                CALL CHECKMEM( IOS, 'STATCOD', PROGNAME )
                ALLOCATE( STATNAM( NSTATE ), STAT=IOS )
                CALL CHECKMEM( IOS, 'STATNAM', PROGNAME )
            CASE( 4 )
                NCOUNTY = K
                ALLOCATE( CNTYCOD( NCOUNTY ), STAT=IOS )
                CALL CHECKMEM( IOS, 'CNTYCOD', PROGNAME )
                ALLOCATE( CNTYNAM( NCOUNTY ), STAT=IOS )
                CALL CHECKMEM( IOS, 'CNTYNAM', PROGNAME )
                ALLOCATE( CNTYTZON( NCOUNTY ), STAT=IOS )
                CALL CHECKMEM( IOS, 'CNTYTZON', PROGNAME )
                ALLOCATE( CNTYTZNM( NCOUNTY ), STAT=IOS )
                CALL CHECKMEM( IOS, 'CNTYTZNM', PROGNAME )
                ALLOCATE( USEDAYLT( NCOUNTY ), STAT=IOS )
                CALL CHECKMEM( IOS, 'USEDAYLT', PROGNAME )
            END SELECT

C.............  Sort codes
            CALL SORTIC( K, INDEXA, CODEA )

C.............  Check for duplicates and store sorted codes
            PREVCODE = ' '
            DO I = 1, K
            
                J = INDEXA( I )
                
                IF( CODEA( J ) == PREVCODE ) THEN
                    EFLAG = .TRUE.
                    MESG = 'Duplicate code ' // PREVCODE // ' found in ' //
     &                     TRIM( FILEDESC ) // ' file.'
                    CALL M3MESG( MESG )
                    CYCLE
                END IF
                
                SELECT CASE( FILENUM )

                CASE( 1 )
                    GEOLEV1COD( I ) = CODEA( J )
                    GEOLEV1NAM( I ) = NAMEA( J )
                CASE( 2 )
                    CTRYCOD( I ) = CODEA( J )
                    CTRYNAM( I ) = NAMEA( J )
                CASE( 3 )
                    STATCOD( I ) = CODEA( J )
                    STATNAM( I ) = NAMEA( J )
                CASE( 4 )
                    CNTYCOD( I ) = CODEA( J )
                    CNTYNAM( I ) = NAMEA( J )
                    CNTYTZON( I ) = TZONEA( J )
                    CNTYTZNM( I ) = TZNAMA( J )
                    USEDAYLT( I ) = USEDYA( J )
                END SELECT
                
                PREVCODE = CODEA( J )
                
            END DO
        
            IF( EFLAG ) THEN
                MESG = 'Problem with ' // TRIM( FILEDESC ) // ' file.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  Deallocate local arrays
            DEALLOCATE( INDEXA, CODEA, NAMEA, TZONEA, TZNAMA, USEDYA )
        
            CLOSE( LDEV )
            
        END DO  ! end loop over levels

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )
        
        END SUBROUTINE RDGEOCODES
