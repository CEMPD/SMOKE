
        SUBROUTINE RDGRPLIST( GDEV, NLINES, GRPLIST )
        
        IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        
C...........   SUBROUTINE ARGUMENTS
        INTEGER,                INTENT (IN)  :: GDEV                ! GROUP file unit no.
        INTEGER,                INTENT (IN)  :: NLINES              ! no. lines in file
        CHARACTER(LEN=FIPLEN3), INTENT (OUT) :: GRPLIST( NLINES,2 ) ! contents of GROUP file

C...........   Local arrays
        CHARACTER(LEN=FIPLEN3) SEGMENT( 2 )          ! parsed input line
        
C...........   Other local variables
        INTEGER I, J, K                   ! counters and indices                     
        
        INTEGER IOS                       ! I/O status
        INTEGER :: IREC = 0               ! record counter

        LOGICAL :: EFLAG      = .FALSE.   ! true: error found
        
        CHARACTER(LEN=100)     LINE     !  line buffer
        CHARACTER(LEN=300)     MESG     !  message buffer

        CHARACTER*16 :: PROGNAME = 'RDGRPLIST'   ! program name
        
C***********************************************************************
C   begin body of subroutine RDGRPLIST

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
     &              'reading county list file at line', IREC
                CALL M3MESG( MESG )
                CYCLE
            END IF
            
C.............  Skip blank lines
            IF( LINE == ' ' ) CYCLE        

C.............  Parse the line into 2 segments
            CALL PARSLINE( LINE, 2, SEGMENT )
            
            GRPLIST( I,1 ) = SEGMENT( 1 )
            GRPLIST( I,2 ) = SEGMENT( 2 )

        END DO

C.........  Abort if error found while reading cross-reference file
        IF( EFLAG ) THEN
            MESG = 'Problem reading GROUP county list file'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF
        
        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )     
        
        END SUBROUTINE RDGRPLIST
        