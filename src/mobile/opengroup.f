
        SUBROUTINE OPENGROUP( DDEV, WDEV, MDEV, EDEV )

C***********************************************************************
C  subroutine body starts at line 67
C
C  DESCRIPTION:
C       Opens any available temperature averaging group files. Gets
C       file name from environment variable, then checks to see if file
C       exists before opening.
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

C...........   EXTERNAL FUNCTIONS and their descriptions:   
        INTEGER, EXTERNAL :: GETEFILE

C...........   FUNCTION ARGUMENTS       
        INTEGER, INTENT (OUT) :: DDEV  ! unit no. for daily group file
        INTEGER, INTENT (OUT) :: WDEV  ! unit no. for weekly group
        INTEGER, INTENT (OUT) :: MDEV  ! unit no. for monthly group
        INTEGER, INTENT (OUT) :: EDEV  ! unit no. for episode group
        
C...........   LOCAL VARIABLES and their descriptions:
        INTEGER    I              ! counter
        INTEGER    IOS            ! temporary I/O status
        INTEGER    FILEDEV        ! unit no. for current file

        LOGICAL :: FEXIST   = .FALSE.  !  true: file exists

        CHARACTER(16)  :: LOGICNAM  !  logical file name of current file
        CHARACTER(300) :: FILENAME  !  name of current group file
        CHARACTER(300) :: MESG      !  message buffer 

        CHARACTER(16)  :: PROGNAME = 'OPENGROUP' ! program name
        
C***********************************************************************
C   begin body of program OPENGROUP

        DDEV = 0
        WDEV = 0
        MDEV = 0
        EDEV = 0

C.........  Set message and logical file name         
        MESG = 'Daily group file name'
        LOGICNAM = 'DAILYGROUP'

C.........  Loop through the four time periods        
        DO I = 1, 4
            FILEDEV = 0

C.............  Get file name from environment variable        
            CALL ENVSTR( LOGICNAM, MESG, 'NONE', FILENAME, IOS )

C.............  If file exists, open it
            IF( FILENAME /= 'NONE' ) THEN
                INQUIRE( FILE=FILENAME, EXIST=FEXIST )
                IF( FEXIST ) THEN
                    FILEDEV = GETEFILE( LOGICNAM, .TRUE., .TRUE., 
     &                                    PROGNAME )
                    IF( FILEDEV == -1 ) THEN
                        MESG = 'Could not open ' // 
     &                   LOGICNAM( 1:LEN_TRIM( LOGICNAM ) ) // 'file'
                        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                    END IF 
                END IF
            END IF

C.............  Store unit number and set message for next time period            
            SELECT CASE( I )
            CASE ( 1 )
                DDEV = FILEDEV
                MESG = 'Weekly group file name'
                LOGICNAM = 'WEEKLYGROUP'
            CASE ( 2 )
                WDEV = FILEDEV
                MESG = 'Monthly group file name'
                LOGICNAM = 'MONTHLYGROUP'
            CASE ( 3 )
                MDEV = FILEDEV
                MESG = 'Episode group file name'
                LOGICNAM = 'EPISODEGROUP'
            CASE ( 4 )
                EDEV = FILEDEV
            END SELECT            
            
        END DO

C.........  Make sure at least one file was opened
        IF( DDEV < 1 .AND. WDEV < 1 .AND. MDEV < 1 .AND. EDEV < 1 ) THEN
            MESG = 'No group files available.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        END SUBROUTINE OPENGROUP
        
