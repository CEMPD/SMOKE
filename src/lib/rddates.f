
        SUBROUTINE RDDATES( PFLAG, FDEV, NTPERIOD )

C**************************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      This subroutine allocates memory for and reads the time periods.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Created 8/2005 by B. Baek 
C
C**************************************************************************
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
C**************************************************************************

C...........   Modules for public variables
C...........   This module contains the lists of unique source characteristics

        USE MODTMPRL, ONLY: STDATE, STTIME, RUNLEN, ITDATE

        IMPLICIT NONE

C...........   EXTERNAL FUNCTIONS and their descriptions:

        INTEGER        GETFLINE
        LOGICAL        BLKORCMT
        INTEGER        STR2INT
        
        EXTERNAL       BLKORCMT, GETFLINE, STR2INT

C...........   Subroutine arguments

        LOGICAL, INTENT(IN   ) :: PFLAG
        INTEGER, INTENT(IN   ) :: FDEV       ! file unit number
        INTEGER, INTENT(  OUT) :: NTPERIOD   ! No of time periods

C...........   Local variables
        INTEGER         I, N                  ! indices and counters
        INTEGER         TSTEP

        INTEGER         ENDLEN                ! end length for reading descriptn
        INTEGER         IOS                   ! i/o status
        INTEGER      :: IREC = 0              ! record number
        INTEGER      :: NLINES = 0            ! number of lines in input file

        CHARACTER(256)  LINE                  ! Read buffer for a line
        CHARACTER(300)  MESG                  ! Message buffer
        CHARACTER(60)   SEGMENT( 3 )          ! line parsing array

        CHARACTER(16) :: PROGNAME = 'RDDATES'    !  program name

C***********************************************************************
C   Begin body of subroutine RDDATES
        
        IF ( .NOT. PFLAG ) THEN
            NTPERIOD = 1
            ALLOCATE( ITDATE( 1 ),
     &                STDATE( 1 ),
     &                STTIME( 1 ),
     &                RUNLEN( 1 ), STAT=IOS )
            CALL CHECKMEM( IOS, 'ITDATE,STDATE,STTIME,RUNLEN', PROGNAME )
            RETURN
        END IF
        
        
        REWIND( FDEV )  ! In case of multiple calls

C.........  Get the number of lines in the PROCDATES file

        NLINES = GETFLINE( FDEV, 'PROCDATES Descriptions' )

C.........  Allocate memory for the PROCDATES descriptions and initialize
        ALLOCATE( ITDATE( NLINES ),
     &            STDATE( NLINES ),
     &            STTIME( NLINES ),
     &            RUNLEN( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ITDATE,STDATE,STTIME,RUNLEN', PROGNAME )

C.........  Read the PROCDATE file and store STDATE, ETDATE, and RUNLEN

        N = 0
        IREC = 0       
        DO I = 1, NLINES

            READ ( FDEV, 93000, END=998, IOSTAT=IOS ) LINE
            IREC = IREC + 1

            IF ( IOS .GT. 0 ) THEN
                WRITE( MESG, 94010)
     &                'I/O error', IOS, 'reading PROCDATES '//
     &                'description file at line', IREC
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  Left adjust line
            LINE = ADJUSTL( LINE )

C.............  Skip blank and comment lines
            IF( BLKORCMT( LINE ) ) CYCLE

C.............  Get line
            CALL PARSLINE( LINE, 3, SEGMENT )

            N = N + 1      ! actual line #s after skipping blank and comments

            STDATE( N ) = STR2INT( SEGMENT( 1 ) )
            STTIME( N ) = STR2INT( SEGMENT( 2 ) )
            RUNLEN( N ) = STR2INT( SEGMENT( 3 ) )
            
        END DO  

        NTPERIOD = N

C.........  Successful completion
        RETURN

C.........  Unexpected end of file
998     MESG = 'INTERNAL ERROR : Unexpected end of PROCDATES desc. file'
        CALL M3MSG2( MESG )

        CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

C******************  FORMAT  STATEMENTS   ******************************

93000   FORMAT( A )

94010   FORMAT( 10( A, :, I8, :, 1X ) )


        END SUBROUTINE RDDATES
