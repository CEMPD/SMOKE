
        SUBROUTINE INITINFO( FILFMT )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      This subroutine initializes the source-category-specific inventory
C      characteristics stored in MODINFO.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Created 4/99 by M. Houyoux
C      09/2025 by HT UNC-IE:  Use M3UTILIO; remove unused FMT statement
C
C****************************************************************************/
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
C***************************************************************************
        USE M3UTILIO

C...........   MODULES for public variables
C.........  This module contains the information about the source category
        USE MODINFO, ONLY: CATEGORY, NCHARS, JSCC, JSTACK,
     &                     NPACT, NPPOL, MXCHRS, LSCCEND, RSCCBEG,
     &                     PLTIDX, SCCLEV1, SCCLEV2, SCCLEV3, SCCLEV4,
     &                     NEM, NDY, NCE, NRE, NRP, NEF, NC1, NC2,
     &                     SC_BEGP, SC_ENDP

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT (IN) :: FILFMT            ! inv format code

C...........   Variables for reading dummy names of emission output

        INTEGER,            ALLOCATABLE :: IDUMARR( : ) !  int dummy array
        CHARACTER(IOVLEN3), ALLOCATABLE :: ENAMES ( : ) !  dummy names
        CHARACTER(IODLEN3), ALLOCATABLE :: CDUMARR( : ) !  char dummy array
    
C...........   Other local variables
        INTEGER         I, IOS               ! memory allocation status

        LOGICAL      :: EFLAG = .FALSE.   ! true: error detected
        LOGICAL,SAVE :: FIRST = .TRUE.    ! true: first time through subroutine

        CHARACTER(300)  MESG

        CHARACTER(16) :: PROGNAME =  'INITINFO' ! program name

C***********************************************************************
C   begin body of subroutine INITINFO

C.........  Skip the rest of the routine if it has been called before           
        IF( .NOT. FIRST ) RETURN
           
C.........  Set category-independent variables
        NPACT = 1   ! no. variables per activity

C.........  Choose the source category for which to set the source information
        SELECT CASE( CATEGORY )

C.........  For area sources ...
        CASE ( 'AREA' )

            NPPOL  = NARPPOL3
            MXCHRS = MXARCHR3
            NCHARS = MXCHRS
            JSCC   = 2
            LSCCEND  = SCCEXPLEN3 + 7         ! first 7
            RSCCBEG  = SCCEXPLEN3 + 8         ! last 3
            PLTIDX   = MXARCHR3  ! needed for ar-to-point
        
            SCCLEV1 = SCCEXPLEN3 + 2
            SCCLEV2 = SCCEXPLEN3 + 4
            SCCLEV3 = SCCEXPLEN3 + 7
            SCCLEV4 = SCCEXPLEN3 + 10

C.........  For mobile sources ...
        CASE ( 'MOBILE' ) 

            NPPOL  = NMBPPOL3
            MXCHRS = MXMBCHR3
            NCHARS = MXCHRS  ! FIPS / SCC (veh type &/or road class) / veh type / lnk
            JSCC   = 0
            LSCCEND  = SCCEXPLEN3 + 7
            RSCCBEG  = SCCEXPLEN3 + 8

            SCCLEV1 = SCCEXPLEN3 + 2
            SCCLEV2 = SCCEXPLEN3 + 4
            SCCLEV3 = SCCEXPLEN3 + 7
            SCCLEV4 = SCCEXPLEN3 + 10

C.........  For point sources ...
        CASE ( 'POINT' )

            NPPOL  = NPTPPOL3
            MXCHRS = MXPTCHR3
            NCHARS = 6
            JSCC   = 6
            JSTACK = 4
            LSCCEND  = SCCEXPLEN3 + 5
            RSCCBEG  = SCCEXPLEN3 + 6

            SCCLEV1 = SCCEXPLEN3 + 3   ! Assumes right-justified 10-digit w/ first 2 zero
            SCCLEV2 = SCCEXPLEN3 + 5
            SCCLEV3 = SCCEXPLEN3 + 8
            SCCLEV4 = SCCEXPLEN3 + 10
          
        END SELECT

C.........  Allocate memory for source characteristic positions (+1 is for
C           pollutant position)
        ALLOCATE( SC_BEGP( MXCHRS + 1 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SC_BEGP', PROGNAME )
        ALLOCATE( SC_ENDP( MXCHRS + 1 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SC_ENDP', PROGNAME )

C.........  Store source characteristic positions
        DO I = 1, MXCHRS + 1
            SELECT CASE( CATEGORY )
            CASE ( 'AREA' )
                SC_BEGP( I ) = ARBEGL3( I )
                SC_ENDP( I ) = ARENDL3( I )
            CASE ( 'MOBILE' ) 
                SC_BEGP( I ) = MBBEGL3( I )
                SC_ENDP( I ) = MBENDL3( I )
            CASE ( 'POINT' )
                SC_BEGP( I ) = PTBEGL3( I )
                SC_ENDP( I ) = PTENDL3( I )
            END SELECT
        END DO

C.........  Allocate memory for dummy variables for determining positions
C           of the fields in the pollutant-specific variables
        ALLOCATE( ENAMES( NPPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ENAMES', PROGNAME )
        ALLOCATE( CDUMARR( NPPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CDUMARR', PROGNAME )
        ALLOCATE( IDUMARR( NPPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IDUMARR', PROGNAME )

C.........  Get dummy output pollutant variable names to use for determining
C           the positions of the fields in the pollutant-specific variable(s)

        CALL BLDENAMS( CATEGORY, 1, NPPOL, 'DUM', ENAMES, CDUMARR,
     &                 IDUMARR, CDUMARR )

C.........  Based on the order of the output names, find the positions in the
C           second dimension of the pollutant-specific data (e.g., POLVLA
C           and POLVAL)
        DO I = 1, NPPOL
            IF( ENAMES( I )(1:IOVLEN3)  .EQ. 'DUM' )    NEM = I
            IF( ENAMES( I )(1:CPRTLEN3) .EQ. AVEDAYRT ) NDY = I
            IF( ENAMES( I )(1:CPRTLEN3) .EQ. CTLEFFRT ) NCE = I
            IF( ENAMES( I )(1:CPRTLEN3) .EQ. RULEFFRT ) NRE = I
            IF( ENAMES( I )(1:CPRTLEN3) .EQ. RULPENRT ) NRP = I
            IF( ENAMES( I )(1:CPRTLEN3) .EQ. EMISFCRT ) NEF = I
            IF( ENAMES( I )(1:CPRTLEN3) .EQ. CECOD1RT ) NC1 = I
            IF( ENAMES( I )(1:CPRTLEN3) .EQ. CECOD2RT ) NC2 = I
        END DO

C.........  Ensure that all of the expect output variable names are present
        SELECT CASE( CATEGORY )
        CASE( 'AREA' )

            IF( NEM .EQ. 0 .OR. NDY .EQ. 0 .OR. NCE .EQ. 0 .OR.
     &          NRE .EQ. 0 .OR. NRP .EQ. 0 .OR. NEF .EQ. 0      ) 
     &          EFLAG = .TRUE.

        CASE( 'MOBILE' )

             IF( NEM .EQ. 0 ) EFLAG = .TRUE.

        CASE( 'POINT' ) 

            IF( NEM .EQ. 0 .OR. NDY .EQ. 0 .OR. NCE .EQ. 0 .OR.
     &          NRE .EQ. 0 .OR. NEF .EQ. 0 .OR. NC1 .EQ. 0 .OR.
     &          NC2 .EQ. 0 ) EFLAG = .TRUE.
 
        END SELECT

        IF( EFLAG ) THEN
            MESG = 'INTERNAL ERROR: Could not find position of ' //
     &             'one or more of the pollutant-specific data.'
            CALL M3MSG2( MESG )
            CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )
        END IF

C.........  Deallocate memory used in this program only
        DEALLOCATE( ENAMES, CDUMARR, IDUMARR )

        FIRST = .FALSE.

        RETURN

        END SUBROUTINE INITINFO
