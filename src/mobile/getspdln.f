
        SUBROUTINE GETSPDLN( SDEV, COUNTY, NLINES, CURRLINE )
        
C***********************************************************************
C  subroutine body starts at line 68
C
C  DESCRIPTION:
C       Finds line number in SPDSUM file where current county information
C       begins.
C
C  PRECONDITIONS REQUIRED:
C       SDEV must be open.
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

C...........   FUNCTION ARGUMENTS
        INTEGER,                INTENT (IN)   :: SDEV     ! SPDSUM file unit no.
        CHARACTER(LEN=FIPLEN3), INTENT (IN)   :: COUNTY   ! county to search for
        INTEGER,                INTENT (IN)   :: NLINES   ! number of lines in SPDSUM file
        INTEGER,                INTENT (INOUT):: CURRLINE ! current line number in SPDSUM file

C...........   Other local variables
        INTEGER I                         ! counters and indices                     

        INTEGER IOS                       ! I/O status
        
        LOGICAL :: FNDLINE = .FALSE.      ! true: found starting line for county
        
        CHARACTER(LEN=FIPLEN3) CURRCOUNTY  ! current county in SPDSUM file
        CHARACTER(LEN=100)     LINE     !  line buffer
        CHARACTER(LEN=300)     MESG     !  message buffer

        CHARACTER*16 :: PROGNAME = 'GETSPDLN'   ! program name

C***********************************************************************
C   begin body of subroutine GETSPDLN
        
        FNDLINE = .FALSE.

C.........  Loop through remaining lines in SPDSUM file
        DO I = CURRLINE, NLINES
        
C.............  Read line
            READ( SDEV, 93000, IOSTAT=IOS ) LINE
            
            IF( IOS /= 0 ) THEN
                IF( IOS == -1 ) THEN
                    MESG = 'End of file reached unexpectedly. ' //
     &              'Check format of SPDSUM file.'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )   
                END IF
                
                WRITE( MESG, 94010 )
     &              'I/O error', IOS,
     &              'reading speed summary file at line', I
                CALL M3MESG( MESG )
                CYCLE
            END IF
                
            CURRCOUNTY = LINE( 1:6 )
            CALL PADZERO( CURRCOUNTY )
            
            IF( CURRCOUNTY == COUNTY ) THEN
                FNDLINE = .TRUE.
                BACKSPACE( SDEV )
                CURRLINE = I
                EXIT
            END IF

        END DO
        	
        IF( .NOT. FNDLINE ) THEN
            REWIND( SDEV )
            CURRLINE = 0
        END IF
        	
        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )
       
        END SUBROUTINE GETSPDLN
        