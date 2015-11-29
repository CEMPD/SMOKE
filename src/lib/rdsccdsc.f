
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
        USE MODLISTS, ONLY: SCCDESC, NINVSCC, INVSCC, SCCDLEV

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constat parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:

        INTEGER        FINDC 
        INTEGER        GETFLINE
        LOGICAL        BLKORCMT

        EXTERNAL       BLKORCMT, FINDC, GETFLINE

C...........   Subroutine arguments
        INTEGER, INTENT (IN) :: FDEV           ! file unit number

C.........  Local parameters
        INTEGER, PARAMETER :: MXSEG = 5       ! max # of potential line segments

C.........  Other arrays
        CHARACTER( SDSLEN3 ) SEGMENT( MXSEG ) ! Segment of parsed lines

C...........   Local variables
        INTEGER         J, K, L, N, P         ! indices and counters

        INTEGER         ENDLEN                ! end length for reading descriptn
        INTEGER         IOS                   ! i/o status
        INTEGER      :: IREC = 0              ! record number
        INTEGER      :: NLINES = 0            ! number of lines in input file

        LOGICAL      :: DFLAG = .FALSE.       ! true: processing delimited format
        LOGICAL      :: FFLAG = .FALSE.       ! true: processing fixed format
        LOGICAL      :: EFLAG = .FALSE.       ! true: error found

        CHARACTER(1)    CHAR                  ! single character buffer
        CHARACTER(256)  LINE                  ! Read buffer for a line
        CHARACTER(300)  MESG                  ! Message buffer

        CHARACTER( SCCLEN3 ) TSCC          ! tmp SCC

        CHARACTER(16) :: PROGNAME = 'RDSCCDSC'    !  program name

C***********************************************************************
C   Begin body of subroutine RDSCCDSC

        REWIND( FDEV )  ! In case of multiple calls

C.........  Get the number of lines in the holidays file
        NLINES = GETFLINE( FDEV, 'SCC Descriptions' )

C.........  Deallocate SCCDLEV if it was allocated 
        IF (ALLOCATED(SCCDLEV)) DEALLOCATE(SCCDLEV)

C.........  Allocate memory for the SCC descriptions and initialize
        ALLOCATE( SCCDESC( NINVSCC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SCCDESC', PROGNAME )
        ALLOCATE( SCCDLEV( NINVSCC, NSCCLV3-1 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SCCDLEV', PROGNAME )

        SCCDESC = 'Description unavailable'          ! array
        SCCDLEV = 23                                 ! array

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

C.............  Get SCC line
            IF( DFLAG ) THEN
                CALL PARSLINE( LINE, 2, SEGMENT )
                TSCC = SEGMENT( 1 )( 1:SCCLEN3 )
                CALL PADZERO( TSCC )

C.................  Find SCC in inventory list, and if it's in the
C                   inventory, store the description.
                J = FINDC( TSCC, NINVSCC, INVSCC )
                IF ( J .GT. 0 ) THEN
                    SCCDESC( J ) = ADJUSTL( SEGMENT( 2 ) )
                END IF
   
            ELSE IF( FFLAG ) THEN
                TSCC = LINE( 1:SCCLEN3-SCCEXPLEN3 )
                CALL PADZERO( TSCC )

C.................  Find SCC in inventory list, and if it's in the
C                   inventory, store the description.
                J = FINDC( TSCC, NINVSCC, INVSCC )

                IF ( J .GT. 0 ) THEN
                    SCCDESC( J ) = ADJUSTL( LINE( SCCLEN3-SCCEXPLEN3+1:ENDLEN ) )
                END IF

            ELSE
                MESG = 'ERROR: Missing file format header in the SCC '//
     &              'description file. Refer to Chapter 8 of the '//
     &              'SMOKE manual for information on the SCCDESC format'
                CALL M3MSG2( MESG )
                CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

            END IF
            
C.................  Determine and store the position of semi-colons that indicate
C                   the parts of the SCCs for the different levels
            IF( J > 0 ) THEN
                L = LEN_TRIM( SCCDESC( J ) )
                P = 0
                CHAR = ' '
                DO K = 1, L
                    CHAR = SCCDESC( J )( K:K )
                    IF( CHAR == ';' ) THEN
                        P = P + 1
                        SCCDLEV( J,P ) = K - 1
                    END IF

C.....................  Get out of loop after 2nd-to-last - all other semi-colons
C                       are included in final section of name
                    IF( P == NSCCLV3-2 ) EXIT
                END DO
                SCCDLEV( J,NSCCLV3-1 ) = L
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
