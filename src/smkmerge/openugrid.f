
        SUBROUTINE OPENUGRID( RLZN )

C***********************************************************************
C  subroutine OPNMRGIN body starts at line
C
C  DESCRIPTION:
C      The purpose of this subroutine is to open all of the necessary
C      files for the merge routine and set the episode information 
C      for the calling program.
C
C  PRECONDITIONS REQUIRED:  
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C       Created 2/99 by M. Houyoux
C       Modified 6/02 by G. Cano
C
C***********************************************************************
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
C****************************************************************************

C.........  MODULES for public variables
C.........  This module contains the major data structure and control flags
        USE MODMERGE

C.........  This module contains the global variables for the 3-d grid
        USE MODGRID

C.........  This module contains the global variables for uncertainty
        USE MODUNCERT

        IMPLICIT NONE

C.........  INCLUDES:
        
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file desc. data structures

C.........  EXTERNAL FUNCTIONS and their descriptions:
        
        CHARACTER*16    PROMPTMFILE  
        LOGICAL         SETENVVAR

        EXTERNAL  PROMPTMFILE, SETENVVAR

C...........   Subroutine arguments
        INTEGER      , INTENT (IN)  :: RLZN      ! realization to merge

C.........  Other local variables

        INTEGER         I, J, N, V       ! counters and indices

        CHARACTER*16    GUNAME       ! set grid file name
        CHARACTER*300   MESG         ! message buffer
        CHARACTER*300   LNAME        ! physical file name for input grid
        CHARACTER(LEN=IOVLEN3) COORD3D    ! coordinate system name 
        CHARACTER(LEN=IOVLEN3) COORUN3D   ! coordinate system projection units
        CHARACTER(LEN=IOVLEN3) PROJTYPE   ! projection type

        CHARACTER*16 :: PROGNAME = 'OPENUGRID' ! program name

C***********************************************************************
C   begin body of subroutine OPENUGRID

C.........  Open gridding matrix for are or mobile sources
        IF( AFLAG .OR. MFLAG) THEN

C.............  initialize uncertainty grid settings if needed
            GUNAME = UCAT( 1:1 ) // 'GMATU'

            IF( RLZN .GT. 1 ) THEN
                IF( .NOT.( CLOSE3( GUNAME ) ) ) THEN
                    WRITE( MESG,94010 ) 
     &                     'WARNING: Unable to close file' 
     &                     // GUNAME, RLZN
                    CALL M3MSG2( MESG )
                END IF
            END IF

            CALL SETUGRID

            IF( .NOT.( SETENVVAR( GUNAME, LNAME ) ) ) THEN

                WRITE( MESG, 94010 ) 
     &                 'Unable to assign setenv for' //
     &                 LNAME, RLZN
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            END IF

C.............  Open gridding matrix, compare number of sources, and 
C               compare or initialize grid information.
            IF( AFLAG ) THEN

                AGNAME = PROMPTMFILE( 
     &              'Enter logical name for the AREA GRIDDING MATRIX',
     &              FSREAD3, GUNAME, PROGNAME )
            END IF

C.............  If we have mobile sources 
            IF( MFLAG ) THEN

                MGNAME = PROMPTMFILE( 
     &              'Enter logical name for the MOBILE GRIDDING MATRIX',
     &               FSREAD3, GUNAME, PROGNAME )

            END IF

        END IF  

C...........   Internal buffering formats.............94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )


C********************************************************************** 

        CONTAINS

            SUBROUTINE SETUGRID

            INTEGER        L1, L2, L3 

            CHARACTER*16  FNAM      ! GMATU file name buffer 
            CHARACTER*16  TMPBUF    ! temporary string buffer   
            CHARACTER*300 DIRBUF    ! directory buffer 

C**********************************************************************

C.............  Set up name for uncertainty output files
            CALL NAMEVAL( UCAT( 1:1 ) // 'GMATU_ODIR', DIRBUF )
            L1 = LEN_TRIM( DIRBUF )

            CALL NAMEVAL( UCAT( 1:1 ) // 'GMATU_FNAM', FNAM )
            L2 = LEN_TRIM( FNAM )

            WRITE( TMPBUF,94010 ) '', RLZN
            TMPBUF = ADJUSTL( TMPBUF )
            CALL PADNZERO( RMXLEN, TMPBUF )
            L3 = LEN_TRIM( TMPBUF )

            LNAME = DIRBUF( 1:L1 ) // '/' // 
     &              FNAM( 1: L2 )  // '_' //
     &              TMPBUF( 1: L3 )  // '.ncf' 

C...........   Internal buffering formats.............94xxx

94010       FORMAT( 10( A, :, I8, :, 1X ) )

            END SUBROUTINE SETUGRID

        END SUBROUTINE OPENUGRID
