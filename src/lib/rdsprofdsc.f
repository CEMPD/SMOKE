
        SUBROUTINE RDSPROFDSC( FDEV )

C**************************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      This subroutine allocates memory for and reads the GSPRO descriptions.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Created 4/2004 by M. Houyoux
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
C.........  This module contains the speciation-profiles matrix, among other things.
        USE MODSPRO,  ONLY : NSPROF, SPROFN, SPCDESC

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constat parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:

        INTEGER        FINDC 
        INTEGER        GETFLINE
        LOGICAL        BLKORCMT

        EXTERNAL       FINDC, GETFLINE, BLKORCMT

C...........   Subroutine arguments
        INTEGER, INTENT (IN) :: FDEV          ! file unit number

C.........  Local parameters
        INTEGER, PARAMETER :: MXSEG = 5       ! max # of potential line segments

C.........  Other arrays
        CHARACTER( SDSLEN3 ) SEGMENT( MXSEG ) ! Segment of parsed lines

C...........   Local variables
        INTEGER         J, L, N, NS          ! indices and counters

        INTEGER         IOS                   ! i/o status
        INTEGER      :: IREC = 0              ! record number
        INTEGER      :: NLINES = 0            ! number of lines in input file

        LOGICAL      :: EFLAG = .FALSE.       ! true: error found

        CHARACTER(256)  LINE                  ! Read buffer for a line
        CHARACTER(300)  MESG                  ! Message buffer
        
        CHARACTER(SPNLEN3) CSPF               ! tmp GSPRO 

        CHARACTER(16) :: PROGNAME = 'RDSPROFDSC'    !  program name

C***********************************************************************
C   Begin body of subroutine RDSPROFDSC

        REWIND( FDEV )  ! In case of multiple calls

C.........  Get the number of lines in the file
        NSPROF = GETFLINE( FDEV, 'GSPRO Speciation Profile Descriptions' )

C.........  Allocate memory for the GSPRO descriptions and initialize
        ALLOCATE( SPROFN( NSPROF ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SPROFN', PROGNAME )
        ALLOCATE( SPCDESC( NSPROF ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SPCDESC', PROGNAME )
        SPROFN  = ''
        SPCDESC = 'Description unavailable'          ! array

C.........  Read the GSPRO descriptions, and store with GSPRO
        NS = 0
        DO N = 1, NSPROF

            READ ( FDEV, 93000, END=998, IOSTAT=IOS ) LINE
            IREC = IREC + 1

            IF ( IOS .GT. 0 ) THEN
                WRITE( MESG, 94010)
     &                'I/O error', IOS, 'reading GSPRO '//
     &                'description file at line', IREC
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  If it is not header, Skip blank and comment lines
            IF( BLKORCMT( LINE ) ) CYCLE

C.............  Left adjust line
            LINE = ADJUSTL( LINE )

C.............  Get GSPRODESC line
            CALL PARSLINE( LINE, 2, SEGMENT )
            CSPF = ADJUSTR( SEGMENT( 1 )( 1:SPNLEN3 ) )

C.............  Look for matched GSPRO 
            J = FINDC( CSPF, NSPROF, SPROFN )
            IF ( J < 1 ) THEN
                NS = NS + 1
                SPROFN ( NS ) = CSPF
                SPCDESC( NS ) = ADJUSTL( SEGMENT( 2 ) )
            END IF

        END DO

C.........  Successful completion
        RETURN

C.........  Unexpected end of file
998     MESG = 'INTERNAL ERROR: Unexpected end of GSPRO description file'
        CALL M3MSG2( MESG )

        CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )


C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE RDSPROFDSC 
