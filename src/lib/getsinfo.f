
        SUBROUTINE GETSINFO

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION: 
C     This subroutine get the source-category-specific source information from
C     the inventory header, which should be read prior to calling this
C     routine.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:  M3EXIT
C
C  REVISION  HISTORY:
C       Created 3/99 by M Houyoux
C
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1998, MCNC--North Carolina Supercomputing Center
C All Rights Reserved
C
C See file COPYRIGHT for conditions of use.
C
C Environmental Programs Group
C MCNC--North Carolina Supercomputing Center
C P.O. Box 12889
C Research Triangle Park, NC  27709-2889
C
C env_progs@mcnc.org
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C****************************************************************************

C.........  MODULES for public variables

C.........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE

C...........   INCLUDES:
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS:
        CHARACTER*2     CRLF
        INTEGER         GETIFDSC

        EXTERNAL        CRLF, GETIFDSC

C...........   LOCAL VARIABLES their descriptions:

        INTEGER       I, J        ! counters and indices 
        INTEGER       IOS         ! memory allocation status
        INTEGER       NVAR        ! number of non-pollutant variables 

        LOGICAL      :: ALLOCPOL = .FALSE.   ! true: allocate & store pollutants

        CHARACTER*300 MESG         ! Message buffer

        CHARACTER*16 :: PROGNAME = 'GETSINFO'    ! Program name

C***********************************************************************
C   begin body of subroutine GETSINFO

C.........  Set the number of sources 
        NSRC = NROWS

C.........  Set source-category-specific information
        SELECT CASE( CATEGORY )

        CASE( 'AREA' )

            ALLOCPOL = .TRUE.
            LSCCEND  = 4
            RSCCBEG  = 5
            NPPOL    = NARPPOL3
            NCHARS   = 3  ! or 2?
            JSCC     = 0
            PLTIDX = 0  ! Change this to be area max (used in XREFTBL)

        CASE( 'MOBILE' )

            ALLOCPOL = .FALSE.
c            NPPOL    = NMBPPOL3


        CASE( 'POINT' )

            ALLOCPOL = .TRUE.
            LSCCEND  = 5
            RSCCBEG  = 6
            PLTIDX   = 2
            NPPOL    = NPTPPOL3
            NCHARS   = GETIFDSC( FDESC3D, '/NUMBER CHARS/', .TRUE. )
            JSCC     = GETIFDSC( FDESC3D, '/SCC POSITION/', .TRUE. )

            ALLOCATE( SC_BEGP( NCHARS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SC_BEGP', PROGNAME )
            ALLOCATE( SC_ENDP( NCHARS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SC_ENDP', PROGNAME )

            DO I = 1, NCHARS
                SC_BEGP( I ) = PTBEGL3( I )
                SC_ENDP( I ) = PTENDL3( I )
            END DO
            
        CASE DEFAULT
            MESG = 'Category ' // CATEGORY // ' not known in program.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END SELECT

C.............  Allocate memory for pollutant names
        IF( ALLOCPOL ) THEN

            NVAR     = GETIFDSC( FDESC3D, '/NON POLLUTANT/', .TRUE. )
            NPPOL    = GETIFDSC( FDESC3D, '/PER POLLUTANT/', .TRUE. )

            NIPOL = ( NVARS3D - NVAR ) / NPPOL

C.............  Allocate memory for and store pollutant names 
            ALLOCATE( EINAM( NIPOL ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EINAM', PROGNAME )

            J = NVAR + 1
            DO I = 1, NIPOL

                EINAM( I ) = VNAME3D( J )
                J = J + NPPOL   ! skip over other pollutant-spec variables

            ENDDO

        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats.............94xxx
 
94010   FORMAT( 10( A, :, I6, :, 2X ) )

        END SUBROUTINE GETSINFO
