
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
C     O3SEASON_YN e.v. to have been checked to set the value of INVPIDX
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

C.........  This module is required by the FileSetAPI
        USE MODFILESET

C.........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE

C...........   INCLUDES:
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
c        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
c        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS:
        CHARACTER*2     CRLF
        LOGICAL         ENVYN
        INTEGER         GETIFDSC

        EXTERNAL        CRLF, ENVYN, GETIFDSC

C...........   LOCAL VARIABLES their descriptions:

        INTEGER       I, J, K     ! counters and indices 
        INTEGER       IOS         ! memory allocation status
        INTEGER       NVAR        ! number of non-pollutant variables 

        LOGICAL       LO3SEAS     ! true: use ozone season emissions

        CHARACTER*300 MESG         ! Message buffer

        CHARACTER*16 :: PROGNAME = 'GETSINFO'    ! Program name

C***********************************************************************
C   begin body of subroutine GETSINFO

C.........  Set the number of sources 
        NSRC = NROWS3D

C.........  Set the number of source characteristics and allocate memory for
C           the field positions
        NCHARS   = GETIFDSC( FDESC3D, '/NUMBER CHARS/', .TRUE. )
        JSCC     = GETIFDSC( FDESC3D, '/SCC POSITION/', .TRUE. )
        JSTACK   = GETIFDSC( FDESC3D, '/STACK POSITION/', .FALSE. )

        IF( .NOT. ALLOCATED( SC_BEGP ) ) THEN
            ALLOCATE( SC_BEGP( NCHARS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SC_BEGP', PROGNAME )
            ALLOCATE( SC_ENDP( NCHARS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SC_ENDP', PROGNAME )
        END IF

C.........  Set source-category-specific information
        SELECT CASE( CATEGORY )

        CASE( 'AREA' )

            LSCCEND  = 7
            RSCCBEG  = 8
            NPPOL    = NARPPOL3
            PLTIDX   = MXARCHR3
            MXCHRS   = MXARCHR3

            DO I = 1, NCHARS
                SC_BEGP( I ) = ARBEGL3( I )
                SC_ENDP( I ) = ARENDL3( I )
            END DO

            SCCLEV1 = 2
            SCCLEV2 = 4
            SCCLEV3 = 7
            SCCLEV4 = 10

        CASE( 'MOBILE' )

            LSCCEND  = SCCLEN3 - VIDLEN3   ! For reset SCCs CRWT // VCID
            RSCCBEG  = LSCCEND + 1
            NPPOL    = NMBPPOL3
            PLTIDX   = 2
            MXCHRS   = MXMBCHR3

            DO I = 1, NCHARS
                SC_BEGP( I ) = MBBEGL3( I )
                SC_ENDP( I ) = MBENDL3( I )
            END DO

            SCCLEV1 = 2
            SCCLEV2 = 4
            SCCLEV3 = 7
            SCCLEV4 = 10
            
        CASE( 'POINT' )

            LSCCEND  = 5
            RSCCBEG  = 6
            PLTIDX   = 2
            NPPOL    = NPTPPOL3
            MXCHRS   = MXPTCHR3

            DO I = 1, NCHARS
                SC_BEGP( I ) = PTBEGL3( I )
                SC_ENDP( I ) = PTENDL3( I )
            END DO
            
            SCCLEV1 = 3   ! Assumes right-justified 10-digit w/ first 2 zero
            SCCLEV2 = 5
            SCCLEV3 = 8
            SCCLEV4 = 10
            
        CASE DEFAULT
            MESG = 'Category ' // CATEGORY // ' not known in program.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END SELECT

C.............  Get sizes of inventory from the header

        NVAR     = GETIFDSC( FDESC3D, '/NON POLLUTANT/', .TRUE. )
        NIPOL    = GETIFDSC( FDESC3D, '/POLLUTANTS/', .FALSE. )
        NPPOL    = GETIFDSC( FDESC3D, '/PER POLLUTANT/', .FALSE. )
        NIACT    = GETIFDSC( FDESC3D, '/ACTIVITIES/', .FALSE. )
        NPACT    = GETIFDSC( FDESC3D, '/PER ACTIVITY/', .FALSE. )

        NIPOL = MAX( 0, NIPOL )
        NPPOL = MAX( 0, NPPOL )
        NIACT = MAX( 0, NIACT )
        NPACT = MAX( 0, NPACT )
        NIPPA = NIPOL + NIACT

C............. Retrieve base year information
        BYEAR = GETIFDSC( FDESC3D, '/BASE YEAR/', .FALSE. )

C............. Allocate memory and store the source attributes units
        IF( .NOT. ALLOCATED( ATTRUNIT ) ) THEN

            ALLOCATE( ATTRUNIT( NVAR ), STAT=IOS )
            CALL CHECKMEM( IOS, 'ATTRUNIT', PROGNAME )

        END IF

        DO I = 1, NVAR
            ATTRUNIT( I ) = VUNITSET( I )
        END DO

        IF( .NOT. ALLOCATED( EANAM ) ) THEN

C............. Allocate memory for the array that stores pollutant
C............. names and activity names, units, and descriptions. Populate 
C              this array in the loops below that fill EINAM and NIACT
            ALLOCATE( EANAM( NIPPA ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EANAM', PROGNAME )
            ALLOCATE( EAREAD( NIPPA ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EAREAD', PROGNAME )
            ALLOCATE( EAUNIT( NIPPA ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EAUNIT', PROGNAME )
            ALLOCATE( EADESC( NIPPA ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EADESC', PROGNAME )

C.............  Allocate memory for and store pollutant names 
            ALLOCATE( EINAM( NIPOL ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EINAM', PROGNAME )

        END IF

        K = 0
        J = NVAR + 1
        DO I = 1, NIPOL

            K = K + 1
            EINAM ( I ) = VNAMESET( J )
	    EANAM ( K ) = VNAMESET( J )
            EAREAD( K ) = VNAMESET( J + INVPIDX )
            EAUNIT( K ) = VUNITSET( J + INVPIDX )
            EADESC( K ) = VDESCSET( J + INVPIDX )

            J = J + NPPOL   ! skip over other pollutant-spec variables

        END DO

C.........  Allocate memory for and store activity names
        IF( .NOT. ALLOCATED( ACTVTY ) ) THEN
            ALLOCATE( ACTVTY( NIACT ), STAT=IOS )
            CALL CHECKMEM( IOS, 'ACTVTY', PROGNAME )
        ENDIF

        DO I = 1, NIACT

            K = K + 1
            ACTVTY( I ) = VNAMESET( J )
            EANAM ( K ) = VNAMESET( J )
            EAREAD( K ) = VNAMESET( J )
            EAUNIT( K ) = VUNITSET( J )
            EADESC( K ) = VDESCSET( J )
            J = J + NPACT   ! skip over other activity-spec variables

        END DO


        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats.............94xxx
 
94010   FORMAT( 10( A, :, I6, :, 2X ) )

        END SUBROUTINE GETSINFO
