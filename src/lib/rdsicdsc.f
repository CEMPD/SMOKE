
        SUBROUTINE RDSICDSC( FDEV )

C**************************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      This subroutine allocates memory for and reads the SIC descriptions.
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
C...........   This module contains the lists of unique source characteristics
        USE MODLISTS, ONLY: SICDESC, NINVSIC, INVSIC

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constat parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:

        INTEGER        FIND1 
        INTEGER        GETFLINE
        INTEGER        STR2INT

        EXTERNAL       FIND1, GETFLINE, STR2INT

C...........   Subroutine arguments
        INTEGER, INTENT (IN) :: FDEV           ! file unit number

C...........   Local variables
        INTEGER         J, N               ! indices and counters

        INTEGER         ENDLEN                ! end length for reading descriptn
        INTEGER         IOS                   ! i/o status
        INTEGER      :: IREC = 0              ! record number
        INTEGER      :: NLINES = 0            ! number of lines in input file
        INTEGER         SIC                   ! tmp SIC

        LOGICAL      :: EFLAG = .FALSE.       ! true: error found


        CHARACTER*256   LINE                  ! Read buffer for a line
        CHARACTER*300   MESG                  ! Message buffer

        CHARACTER*16 :: PROGNAME = 'RDSICDSC'    !  program name

C***********************************************************************
C   Begin body of subroutine RDSICDSC

        REWIND( FDEV )  ! In case of multiple calls

C.........  Get the number of lines in the holidays file
        NLINES = GETFLINE( FDEV, 'SIC Descriptions' )

C.........  Allocate memory for the SCC descriptions and initialize
        ALLOCATE( SICDESC( NINVSIC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SICDESC', PROGNAME )
        SICDESC = 'Description unavailable'          ! array

C.........  Read the SCC descriptions, and store with SCC
        ENDLEN = SICLEN3 + SDSLEN3
        DO N = 1, NLINES

            READ ( FDEV, 93000, END=998, IOSTAT=IOS ) LINE
            IREC = IREC + 1

            IF ( IOS .GT. 0 ) THEN
                WRITE( MESG, 94010)
     &                'I/O error', IOS, 'reading SIC '//
     &                'description file at line', IREC
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  Skip blank lines
            IF( LINE .EQ. ' ' ) CYCLE

C.............  Left adjust line
            LINE = ADJUSTL( LINE )

C.............  Get SIC line
            SIC = STR2INT( LINE( 1:SICLEN3 ) )

C.............  Find SCC in inventory list, and if it's in the inventory, 
C               store the description.
            J = FIND1( SIC, NINVSIC, INVSIC )

            IF ( J .GT. 0 ) THEN
                SICDESC( J ) = ADJUSTL( LINE( SICLEN3+1:ENDLEN ) )
            END IF

        END DO

C.........  Successful completion
        RETURN

C.........  Unexpected end of file
998     MESG = 'INTERNAL ERROR: Unexpected end of SIC description file'
        CALL M3MSG2( MESG )

        CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )


C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94020   FORMAT( A, 1X, I8, 1X, A, 1X, F10.6, 1X, A )

94030   FORMAT( A, 1X, I6.6, A, 100( ' SSC(', I2.2, '):', F10.6, : ) )

        END SUBROUTINE RDSICDSC
