
        SUBROUTINE RDSPDSRC( SDEV, NSRC, SRCARRAY )

C***********************************************************************
C  subroutine body starts at line 80
C
C  DESCRIPTION:
C       Reads the SPDSUM file
C
C  PRECONDITIONS REQUIRED:
C       SDEV must be opened
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
        INTEGER       STR2INT
        CHARACTER(2)  CRLF    
        
        EXTERNAL  STR2INT, CRLF

C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT (IN)  :: SDEV               ! SPDSUM file unit number
        INTEGER, INTENT (IN)  :: NSRC               ! no. of sources
        INTEGER, INTENT (OUT) :: SRCARRAY( NSRC,2 ) ! array to hold county codes

C...........   Local arrays
        INTEGER       SOURCES( 7 )          ! line of sources from SPDSUM file

C...........   Other local variables
        INTEGER I, J                      ! counters and indices                     
        
        INTEGER IOS                       ! I/O status
        INTEGER :: IREC = 0               ! record counter
        INTEGER COUNTY                    ! county from SPDSUM file
        INTEGER ROADTYPE                  ! roadtype from SPDSUM file
         
        REAL SPEED                        ! speed from SPDSUM file

        LOGICAL   :: EFLAG     = .FALSE.  !  true: error found
                
        CHARACTER   CONTCHAR          ! continuation character from SPDSUM file
        
        CHARACTER(100)     INTFMT     ! SPDSUM format string with integer
        CHARACTER(100)     REALFMT    ! SPDSUM format string with real
        CHARACTER(300)     MESG       !  message buffer

        CHARACTER(16) :: PROGNAME = 'RDSPDSRC'   ! program name

C***********************************************************************
C   begin body of subroutine RDSPDSRC

        CALL GETSPDFMT( INTFMT, REALFMT )

        DO
        
            SOURCES = 0   ! array
        
C.............  Read line from SPDSUM file
            READ( SDEV, REALFMT, IOSTAT=IOS, END=10 ) COUNTY, ROADTYPE, 
     &            SPEED, SOURCES, CONTCHAR

C.............  Exit if we've reached the end of the file
            IF( IOS == -1 ) EXIT

            IREC = IREC + 1

C.............  Check for other I/O errors        
            IF( IOS /= 0 ) THEN
                EFLAG = .TRUE.

                WRITE( MESG, 94010 )
     &                 'I/O error', IOS,
     &                 'reading speed summary file at line', IREC
                CALL M3MESG( MESG )
                CYCLE
            END IF

C.............  Store county code for each source
            DO J = 1, 7
                IF( SOURCES( J ) /= 0 ) THEN
                    SRCARRAY( SOURCES( J ),1 ) = COUNTY
                END IF
            END DO
                        
        END DO

C.........  Abort if error found while reading SPDSUM file
        IF( EFLAG ) THEN
            MESG = 'Problem reading SPDSUM file'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

10      RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )
94020   FORMAT( 10( A, :, F6.2, :, 1X ) )
        
        END SUBROUTINE RDSPDSRC
        
