
        SUBROUTINE CHKM6SCN( FILENAME, SCENARIO, NLINES, EFYEAR,
     &                       FLATFLAG )
     
C***********************************************************************
C  subroutine body starts at line 77
C
C  DESCRIPTION:
C       Checks the MOBILE6 input scenario for various commands. Comments
C       out any commands that are not allowed or replaced by SMOKE. Checks
C       that calendar year matches EF_YEAR.
C
C  PRECONDITIONS REQUIRED:
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
        
        IMPLICIT NONE        

C...........   INCLUDES:

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        
C...........   EXTERNAL FUNCTIONS and their descriptions:
        LOGICAL        CHKINT
        INTEGER        STR2INT
        CHARACTER*2    CRLF    
        
        EXTERNAL  CHKINT, STR2INT, CRLF

C...........   SUBROUTINE ARGUMENTS
        CHARACTER*300, INTENT (IN)    :: FILENAME             ! scenario filename
        CHARACTER*150, INTENT (INOUT) :: SCENARIO( NLINES )   ! M6 scenario
        INTEGER,       INTENT (IN)    :: NLINES               ! no. lines in scenario (size of array)
        INTEGER,       INTENT (IN)    :: EFYEAR               ! emission factor year
        LOGICAL,       INTENT (IN)    :: FLATFLAG             ! true: use flat hourly VMT profile

C...........   Other local variables
        INTEGER I, J, K                       ! counters and indices                     
        INTEGER YEAR                          ! calendar year

        CHARACTER(LEN=150) CURRLINE           ! current line from scenario
        CHARACTER(LEN=150) RPLCLINE           ! replacement line
        CHARACTER(LEN=19)  COMMAND            ! Mobile 6 command
        CHARACTER(LEN=130) YEARSTR            ! calendar year from scenario
        CHARACTER(LEN=300) MESG               ! message buffer

        CHARACTER*16 :: PROGNAME = 'CHKM6SCN'   ! program name
        
C***********************************************************************
C   begin body of subroutine CHKM6SCN

C.........  Loop through lines of scenario
        DO I = 1, NLINES
            CURRLINE = SCENARIO( I )
            
C.............  Skip comment and blank lines
            IF( CURRLINE( 1:1 ) == '*' ) CYCLE
            IF( CURRLINE( 1:1 ) == '>' ) CYCLE
            IF( CURRLINE == ' ' ) CYCLE

C.............  Get Mobile6 command                   
            COMMAND = CURRLINE( 1:19 )

C.............  Comment out unused commands (this can be added to as needed)
            IF( INDEX( COMMAND, 'SCENARIO RECORD' ) > 0 .OR.
     &          INDEX( COMMAND, 'PARTICLE SIZE' ) > 0 .OR.
     &          INDEX( COMMAND, 'AVERAGE SPEED' ) > 0 .OR.
     &          INDEX( COMMAND, 'SPEED VMT' ) > 0 ) THEN
                RPLCLINE( 1:1 ) = '*'
                RPLCLINE( 2:150 ) = CURRLINE( 1:149 )
                SCENARIO( I ) = RPLCLINE
            END IF

C.............  Check for VMT BY HOUR command
            IF( FLATFLAG ) THEN
                IF( INDEX( COMMAND, 'VMT BY HOUR' ) > 0 ) THEN
                    RPLCLINE( 1:1 ) = '*'
                    RPLCLINE( 2:150 ) = CURRLINE( 1:149 )
                    SCENARIO( I ) = RPLCLINE
                END IF
            END IF

C.............  Search command for year command
            IF( INDEX( COMMAND,'CALENDAR YEAR' ) > 0 ) THEN
                YEARSTR = ADJUSTL( CURRLINE( 21:150 ) )
                IF( CHKINT( YEARSTR ) ) THEN
                    YEAR = STR2INT( YEARSTR )
                    
                    IF( YEAR /= EFYEAR ) THEN
                       WRITE( MESG,94010 ) 'WARNING: Changing CALENDAR'
     &                  // ' YEAR in MOBILE6 scenario ' // CRLF()
     &                  // BLANK5 // FILENAME( 1:LEN_TRIM( FILENAME ) )
     &                  // CRLF() // BLANK5 // 'to match EF_YEAR',
     &                     EFYEAR
                       CALL M3MESG( MESG )
                       
                       WRITE( RPLCLINE,94020 ) 
     &                        'CALENDAR YEAR      : ', EFYEAR
                       SCENARIO( I ) = RPLCLINE
                    END IF
                END IF
            END IF
 
        END DO

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )
94020   FORMAT( A21, I4 )
        
        END SUBROUTINE CHKM6SCN
        
