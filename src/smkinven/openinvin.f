
        SUBROUTINE OPENINVIN( CATEGORY, IDEV, RDEV, PDEV, VDEV,
     &                        ZDEV, INNAME )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      This subroutine opens the appropriate input files for inventory import
C      given the source category.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C      Subroutines: I/O API subroutines
C
C  REVISION  HISTORY:
C      Created 4/99 by M. Houyoux
C
C****************************************************************************/
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1999, MCNC--North Carolina Supercomputing Center
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
C***************************************************************************

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constat parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS and their descriptionsNRAWIN
        LOGICAL         ENVYN
        INTEGER         INDEX1
        INTEGER         PROMPTFFILE

        EXTERNAL        ENVYN, INDEX1, PROMPTFFILE

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN) :: CATEGORY  ! source category
        INTEGER     , INTENT(OUT) :: IDEV      ! unit number for inven file
        INTEGER     , INTENT(OUT) :: RDEV      ! unit no. for stack replacements
        INTEGER     , INTENT(OUT) :: PDEV      ! unit no. for pol codes & names
        INTEGER     , INTENT(OUT) :: VDEV      ! unit no. for activity names
        INTEGER     , INTENT(OUT) :: ZDEV      ! unit no. for time zones
        CHARACTER(*), INTENT(OUT) :: INNAME    ! tmp file name buffer 

C...........   LOCAL PARAMETERS
        INTEGER    , PARAMETER :: NCAT = 3          ! number src categories
        CHARACTER*6, PARAMETER :: CATLIST( NCAT ) = ! src categories
     &                         ( / 'AREA  ', 'MOBILE', 'POINT ' / )

        CHARACTER*5, PARAMETER :: NAMLIST( NCAT ) = ! input names
     &                         ( / 'ARINV' , 'MBINV' , 'PTINV'  / )

C...........   Other local variables

        INTEGER       IOS    ! i/o status
        INTEGER       J      ! counter and indices
        INTEGER       LCAT   ! length of CATEGORY string

        LOGICAL       CFLAG  ! true: use SIPOLS file
        LOGICAL       VFLAG  ! true: use ACTVNAMS file

        CHARACTER*300          MESG        ! message buffer 

        CHARACTER*16 :: PROGNAME = 'OPENINVIN' ! program name

C***********************************************************************
C   begin body of subroutine OPENINVIN

C.........  Set controls for reading the pollutants and activities files
C.........  Default is for mobile to read in activities and not pollutants
C           and for other source categories to read in pollutants and not
C           activities
        IF( CATEGORY .EQ. 'MOBILE' ) THEN
            CFLAG = .FALSE.
            VFLAG = .TRUE.
        ELSE
            CFLAG = .TRUE.
            VFLAG = .FALSE.
        END IF

C.........  Get value of these controls from the environment
        MESG = 'Indicator for using pollutants list'
        CFLAG = ENVYN( 'SMK_USE_SIPOLS', MESG, CFLAG, IOS )

        MESG = 'Indicator for using activities list'
        VFLAG = ENVYN( 'SMK_USE_ACTVNAMS', MESG, VFLAG, IOS )

C.........  Get the length of the category name
        LCAT = LEN_TRIM( CATEGORY )

C.........  Find name name of raw inventory file
        J = INDEX1( CATEGORY, NCAT, CATLIST )
        IF( J .LE. 0 ) THEN
            MESG = 'INTERNAL ERROR: Do not know about category ' //
     &             CATEGORY( 1:LCAT ) // ' in program ' // PROGNAME
            CALL M3MSG2( MESG )
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ELSE
 
            INNAME = NAMLIST( J )

        END IF

C.........  Set prompt for raw inventory file
        MESG = 'Enter logical name of the RAW ' // CATEGORY( 1:LCAT ) // 
     &         ' INVENTORY ' // 'file'

C.........  Get file name an open raw input inventory file
        IDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., INNAME, PROGNAME )

C.........  Open category-specific inputs
        SELECT CASE( CATEGORY )
        CASE( 'MOBILE' ) 

C.............  Get file name for converting road-class to road type & 
C               vehicle type name to vehicle type number.
            MESG = 'Enter logical name for MOBILE CODES file'
            RDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'MCODES', 
     &                          PROGNAME )

        CASE( 'POINT' )

C.............  Get file name for input replacement stack parameters file
            MESG = 'Enter logical name for REPLACEMENT STACK ' //
     &             'PARAMETERS file'
            RDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'PSTK', PROGNAME )

        END SELECT

C.........  Get file name for country, state, and county file, with time zones
        ZDEV = PROMPTFFILE(
     &             'Enter logical name for COUNTRY, STATE, AND ' //
     &             'COUNTY file', .TRUE., .TRUE., 'COSTCY', PROGNAME )

C.........  Get file name for inventory pollutants codes/names
        IF( CFLAG ) THEN 
            MESG = 'Enter logical name for POLLUTANT CODES & NAMES file'
            PDEV = PROMPTFFILE( MESG,.TRUE.,.TRUE.,'SIPOLS',PROGNAME )
        ELSE
            PDEV = 0
        END IF

C.........  Get file name for inventory pollutants codes/names
        IF( VFLAG ) THEN
            MESG = 'Enter logical name for ACTIVITY CODES & NAMES file'
            VDEV = PROMPTFFILE( MESG,.TRUE.,.TRUE.,'ACTVNAMS',PROGNAME )
        ELSE
            VDEV = 0
        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )
 
        END SUBROUTINE OPENINVIN

