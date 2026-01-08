
        SUBROUTINE RDPELV( FDEV, NSRC, ASCIFLAG, NMAJOR, NPING )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C       Allocates memory for and reads in the PELV file output from ELEVPOINT.
C       This file contains a column for Major sources, plume-in-grid, sources, 
C       and for stack group numbers.  Only the source codes that are either 
C       identified as major source or PinG source are listed. 
C
C  PRECONDITIONS REQUIRED:
C
C  REVISION  HISTORY:
C       Written  1/99 by M. Houyoux
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
C       Updated with USE M3UTILIO by Huy Tran UNC-IE on 2026-01
C****************************************************************************

C.........  MODULES for public variables
C.........  This module contains arrays for plume-in-grid and major sources
        USE M3UTILIO

        USE MODELEV, ONLY: LMAJOR, LPING, GROUPID

        IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'EMCNST3.EXT'     ! emissions constant parameters
C        INCLUDE 'PARMS3.EXT'      ! I/O API constants
C        INCLUDE 'FDESC3.EXT'      ! I/O API file description data structure
C        INCLUDE 'IODECL3.EXT'     ! I/O API function declarations

C...........   EXTERNAL FUNCTIONS:
C       CHARACTER(2)  CRLF
C       LOGICAL       ENVYN

C        EXTERNAL      CRLF, ENVYN

C...........   ARGUMENTS and their descriptions: actually-occurring ASC table

        INTEGER, INTENT (IN)  :: FDEV    !  unit number for elev srcs file 
        INTEGER, INTENT (IN)  :: NSRC    !  no. sources
        LOGICAL, INTENT (IN)  :: ASCIFLAG!  true: outputing an ASCII elevated
        INTEGER, INTENT (OUT) :: NMAJOR  !  number of major sources
        INTEGER, INTENT (OUT) :: NPING   !  number of PinG sources

C...........   Arrays dimensioned by subroutine arguments
        INTEGER      SINDX( NSRC ) ! sorting index for group ID

C...........   OTHER LOCAL VARIABLES and their descriptions:

        INTEGER         I, J, K, L, L2, S   !  counters and indices
        INTEGER         IOS              !  I/O Status
        INTEGER         IGRP             !  group number for PinG source, or 0
        INTEGER         IMAJR            !  src ID for major source, or 0
        INTEGER         IPING            !  src ID for PinG source, or 0
        INTEGER         IREC             !  input line counter
        INTEGER         PGRP             !  group from previous iteration

        LOGICAL      :: EFLAG = .FALSE.  !  error flag

        CHARACTER(32 )  SRGFMT           !  buffer for grid format (MODEL-3) 
        CHARACTER(300)  BUFFER           !  buffer for formatted source chars
        CHARACTER(300)  MESG             !  message buffer

        CHARACTER(16) :: PROGNAME = 'RDPELV' ! program name

C***********************************************************************
C   begin body of subroutine RDPELV

C.........   Get settings from environment variables
C.........   This variable enables have no "major" sources identified in the
C            PELV file, and actually not doing plume rise on any sources instead
C            of the default behaviour, which is to do plume rise on all sources

C.........   Check gridded information is consistent with current grid info
        IF( FDEV > 0 ) CALL RDSRGHDR( .FALSE., FDEV, SRGFMT )   ! CHKGRID may be initialized

C.........  Allocate the MODELEV arrays for identifying major/PinG sources
        ALLOCATE( LMAJOR( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LMAJOR', PROGNAME )
        ALLOCATE( LPING( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LPING', PROGNAME )
        ALLOCATE( GROUPID( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'GROUPID', PROGNAME )

C.........  Initialize arrays
        LMAJOR   = .FALSE.   ! arrays 
        LPING    = .FALSE.
        GROUPID  = 0 

        MESG = 'Determining elevated/plume-in-grid sources...'
        CALL M3MSG2( MESG )

C.........  Read in lines knowing that they are formatted in the CSOURC
C           spacing
        NMAJOR = 0
        NPING  = 0
        IREC   = 0
        DO           !  head of the FDEV-read loop

C.............  If no input file, then end read loop
            IF( FDEV .LE. 0 ) EXIT

            READ( FDEV, *, END=23, IOSTAT=IOS ) IMAJR, IPING, IGRP
            IREC = IREC + 1

            IF( IOS .GT. 0 ) THEN

                EFLAG = .TRUE.
                WRITE( MESG,94010 )
     &              'ERROR: I/O error', IOS,
     &              'reading elevated sources file at line', IREC
                CALL M3MESG( MESG )

            END IF

C.............  Set sources that are major sources
            IF( IMAJR .GT. 0 .AND. IMAJR .LE. NSRC ) THEN
                NMAJOR = NMAJOR + 1
                LMAJOR  ( IMAJR ) = .TRUE.
                GROUPID ( IMAJR ) = IGRP
            END IF 

C.............  Set sources that are PinG sources
            IF( IPING .GT. 0 .AND. IPING .LE. NSRC ) THEN
                NPING = NPING + 1                
                LPING   ( IPING ) = .TRUE.
                GROUPID ( IPING ) = IGRP

C.................  When outputting ASCII output file, flag PinG sources
C                   as elevated sources.
                IF ( ASCIFLAG ) THEN
                    NMAJOR = NMAJOR + 1
                    LMAJOR  ( IPING ) = .TRUE.
                END IF

            END IF 

C.............  If index is out of range, ELEVPOINT needs rerunning
            IF( IMAJR .GT. NSRC .OR. IPING .GT. NSRC ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: Source ID', NSRC, 
     &                 'in elevated sources file is greater than' //
     &                 CRLF()// BLANK10 // 
     &                 'the maximum number of sourced!' //
     &                 CRLF()// BLANK10 // 'Elevpoint should be rerun.'
                CALL M3MESG( MESG )

            END IF 

        END DO       !  to head of elevated sources loop

23      CONTINUE    !  end of read loop

C.........  If there are no major sources specifically identified...
        IF( NMAJOR .EQ. 0 ) THEN

C.............  Error if ASCII elevated file being created and no PELV file
C.............  Note that when called from Laypoint, ASCIFLAG is never true,
C               even if Laypoint is being run for UAM-style explicit plumerise
            IF( ASCIFLAG .AND. FDEV .LE. 0 ) THEN
                EFLAG = .TRUE.
                MESG = 'ERROR: No PELV file input, but it is ' //
     &                 'required for UAM-style processing.'

C.............  Error if ASCII elevated file being created and major sources
            ELSE IF( ASCIFLAG ) THEN
                EFLAG = .TRUE.
                MESG = 'ERROR: No elevated sources identified in ' //
     &                 'PELV input file.'

C.............  For CMAQ-style, change array to make all sources potentially
C               elevated when there is no PELV file used
            ELSE IF( FDEV .LE. 0 ) THEN
                NMAJOR = NSRC
                LMAJOR = .TRUE.  ! array
                MESG = 'NOTE: All non-PinG sources are potentially '//
     &                 'elevated'

C.............  For CMAQ-style, with a PELV file, change all to major, but
C               give a warning.
            ELSE
                NMAJOR = NSRC
                LMAJOR = .TRUE.  ! array
                MESG = 'WARNING: No major/PinG sources in PELV file, '//
     &                 CRLF() // BLANK10 // 'All non-PinG sources ' //
     &                 'potentially elevated'

            END IF

        ELSE
            WRITE( MESG,94010 ) 'NOTE:', NMAJOR,' sources ' //
     &             'will be elevated, and the rest are' //
     &             CRLF()// BLANK10 // 'low-level (layer 1)'
     &             
        END IF
        CALL M3MSG2( MESG )

        IF( NPING .EQ. 0 ) THEN
            MESG = 'NOTE: No PinG sources will be modeled'
        ELSE
            WRITE( MESG,94010 ) 'NOTE:', NPING,' PinG sources will ' //
     &             'be modeled'
        END IF
        CALL M3MSG2( MESG )

C.........  Abort if error occured on read
        IF( EFLAG ) THEN
            MESG = 'Problem with elevated sources file'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

93020   FORMAT( I7, I3 )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10 ( A, :, I10, :, 2X ) )


        END SUBROUTINE RDPELV

