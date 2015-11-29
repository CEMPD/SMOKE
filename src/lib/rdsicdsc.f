
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
        INTEGER         J, L, N               ! indices and counters

        INTEGER         ENDLEN                ! end length for reading descriptn
        INTEGER         IOS                   ! i/o status
        INTEGER      :: IREC = 0              ! record number
        INTEGER      :: NLINES = 0            ! number of lines in input file

        LOGICAL      :: DFLAG = .FALSE.       ! true: processing delimited format
        LOGICAL      :: FFLAG = .FALSE.       ! true: processing fixed format
        LOGICAL      :: EFLAG = .FALSE.       ! true: error found

        CHARACTER(256)  LINE                  ! Read buffer for a line
        CHARACTER(300)  MESG                  ! Message buffer
        
        CHARACTER(SICLEN3) CSIC               ! tmp SIC

        CHARACTER(16) :: PROGNAME = 'RDSICDSC'    !  program name

C***********************************************************************
C   Begin body of subroutine RDSICDSC

        REWIND( FDEV )  ! In case of multiple calls

C.........  Get the number of lines in the file
        NLINES = GETFLINE( FDEV, 'SIC Descriptions' )

C.........  Allocate memory for the SIC descriptions and initialize
        ALLOCATE( SICDESC( NINVSIC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SICDESC', PROGNAME )
        SICDESC = 'Description unavailable'          ! array

C.........  Read the SIC descriptions, and store with SIC
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

C.............  Check header to define the format (comma-delimited or fixed)
C               If it is not header, Skip blank and comment lines
            IF( INDEX( LINE, '#FIXED' ) > 0 ) THEN
                FFLAG = .TRUE.    ! processing fixed format
                CYCLE
            ELSE IF( INDEX( LINE, '#DELIMITED' ) > 0 ) THEN
                DFLAG = .TRUE.    ! processing delimited format
                CYCLE
            ELSE IF( BLKORCMT( LINE ) ) THEN
                CYCLE
            END IF

C.............  Left adjust line
            LINE = ADJUSTL( LINE )

C.............  Get SIC line
            IF( DFLAG ) THEN
                CALL PARSLINE( LINE, 2, SEGMENT )
                CSIC = SEGMENT( 1 )( 1:SICLEN3 )
                CALL PADZERO( CSIC )

C.................  Find SIC in inventory list, and if it's in the
C                   inventory, store the description.
                J = FINDC( CSIC, NINVSIC, INVSIC )

                IF ( J .GT. 0 ) THEN
                    SICDESC( J ) = ADJUSTL( SEGMENT( 2 ) )
                END IF
   
            ELSE IF( FFLAG ) THEN
                CSIC = LINE( 1:SICLEN3-SICEXPLEN3 )
                CALL PADZERO( CSIC )

C.................  Find SIC in inventory list, and if it's in the
C                   inventory, store the description.
                J = FINDC( CSIC, NINVSIC, INVSIC )

                IF ( J .GT. 0 ) THEN
                    SICDESC( J ) = ADJUSTL( LINE( SICLEN3-SICEXPLEN3+1:ENDLEN ) )
                END IF

            ELSE
                MESG = 'ERROR: Missing file format header in the SIC '//
     &              'description file. Refer to Chapter 8 of the '//
     &              'SMOKE manual for information on the SICDESC format'
                CALL M3MSG2( MESG )
                CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

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

        END SUBROUTINE RDSICDSC
