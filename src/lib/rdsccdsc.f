
        SUBROUTINE RDSCCDSC( FDEV )

C**************************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      This subroutine allocates memory for and reads the SCC descriptions.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Created 7/2001 by M. Houyoux
C
C**************************************************************************
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
C**************************************************************************

C...........   Modules for public variables
C...........   This module contains the lists of unique source characteristics
        USE MODLISTS

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constat parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:

        INTEGER        FINDC 
        INTEGER        GETFLINE

        EXTERNAL       FINDC, GETFLINE

C...........   Subroutine arguments
        INTEGER, INTENT (IN) :: FDEV           ! file unit number

C...........   Local variables
        INTEGER         J, N               ! indices and counters

        INTEGER         ENDLEN                ! end length for reading descriptn
        INTEGER         IOS                   ! i/o status
        INTEGER      :: IREC = 0              ! record number
        INTEGER      :: NLINES = 0            ! number of lines in input file

        LOGICAL      :: EFLAG = .FALSE.       ! true: error found


        CHARACTER*256   LINE                  ! Read buffer for a line
        CHARACTER*300   MESG                  ! Message buffer

        CHARACTER(LEN=SCCLEN3 ) TSCC          ! tmp SCC

        CHARACTER*16 :: PROGNAME = 'RDSCCDSC'    !  program name

C***********************************************************************
C   Begin body of subroutine RDSCCDSC

C.........  Get the number of lines in the holidays file
        NLINES = GETFLINE( FDEV, 'SCC Descriptions' )

C.........  Allocate memory for the SCC descriptions and initialize
        ALLOCATE( SCCDESC( NINVSCC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SCCDESC', PROGNAME )
        SCCDESC = 'Description unavailable'          ! array

C.........  Read the SCC descriptions, and store with SCC
        ENDLEN = SCCLEN3 + SDSLEN3
        DO N = 1, NLINES

            READ ( FDEV, 93000, END=998, IOSTAT=IOS ) LINE
            IREC = IREC + 1

            IF ( IOS .GT. 0 ) THEN
                WRITE( MESG, 94010)
     &                'I/O error', IOS, 'reading SCC '//
     &                'description file at line', IREC
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  Skip blank lines
            IF( LINE .EQ. ' ' ) CYCLE

C.............  Left adjust line
            LINE = ADJUSTL( LINE )

C.............  Get SCC line
            TSCC = LINE( 1:SCCLEN3 )
            CALL PADZERO( TSCC )

C.............  Find SCC in inventory list, and if it's in the inventory, 
C               store the description.
            J = FINDC( TSCC, NINVSCC, INVSCC )

            IF ( J .GT. 0 ) THEN
                SCCDESC( J ) = ADJUSTL( LINE( SCCLEN3+1:ENDLEN ) )
            END IF

        END DO

C.........  Successful completion
        RETURN

C.........  Unexpected end of file
998     MESG = 'INTERNAL ERROR: Unexpected end of SCC description file'
        CALL M3MSG2( MESG )

        CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )


C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94020   FORMAT( A, 1X, I8, 1X, A, 1X, F10.6, 1X, A )

94030   FORMAT( A, 1X, I6.6, A, 100( ' SSC(', I2.2, '):', F10.6, : ) )

        END SUBROUTINE RDSCCDSC
