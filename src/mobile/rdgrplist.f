
        SUBROUTINE RDGRPLIST( GDEV, NLINES, GRPLIST )

C***********************************************************************
C  subroutine body starts at line 72
C
C  DESCRIPTION:
C       Reads list of counties in current temporal averaging group
C
C  PRECONDITIONS REQUIRED:
C       GDEV must be open
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
C COPYRIGHT (C) 2002, MCNC Environmental Modeling Center
C All Rights Reserved
C
C See file COPYRIGHT for conditions of use.
C
C Environmental Modeling Center
C MCNC
C P.O. Box 12889
C Research Triangle Park, NC  27709-2889
C
C smoke@emc.mcnc.org
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C***********************************************************************

        IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        INTEGER, EXTERNAL :: STR2INT
        
C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT (IN)  :: GDEV                ! GROUP file unit no.
        INTEGER, INTENT (IN)  :: NLINES              ! no. lines in file
        INTEGER, INTENT (OUT) :: GRPLIST( NLINES,3 ) ! contents of GROUP file

C...........   Local arrays
        CHARACTER(LEN=FIPLEN3) SEGMENT( 3 )          ! parsed input line
        
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
            CALL PARSLINE( LINE, 3, SEGMENT )
            
            GRPLIST( I,1 ) = STR2INT( SEGMENT( 1 ) )
            GRPLIST( I,2 ) = STR2INT( SEGMENT( 2 ) )
            GRPLIST( I,3 ) = STR2INT( SEGMENT( 3 ) )

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
93010   FORMAT( I6, 1X, I6, 1X, I1 ) 

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )     
        
        END SUBROUTINE RDGRPLIST
        