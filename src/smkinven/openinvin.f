
        SUBROUTINE OPENINVIN( CATEGORY, IDEV, RDEV, PDEV, ZDEV, INNAME )

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
C***************************************************************************

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constat parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS and their descriptionsNRAWIN
        INTEGER         INDEX1
        INTEGER         PROMPTFFILE

        EXTERNAL        INDEX1, PROMPTFFILE

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN) :: CATEGORY  ! source category
        INTEGER     , INTENT(OUT) :: IDEV      ! unit number for inven file
        INTEGER     , INTENT(OUT) :: RDEV      ! unit no. for stack replacements
        INTEGER     , INTENT(OUT) :: PDEV      ! unit no. for pol codes & names
        INTEGER     , INTENT(OUT) :: ZDEV      ! unit no. for time zones
        CHARACTER(*), INTENT(OUT) :: INNAME    ! tmp file name buffer 

C...........   LOCAL PARAMETERS
        INTEGER    , PARAMETER :: NCAT = 3          ! number src categories
        CHARACTER*6, PARAMETER :: CATLIST( NCAT ) = ! src categories
     &                         ( / 'AREA  ', 'MOBILE', 'POINT ' / )

        CHARACTER*5, PARAMETER :: NAMLIST( NCAT ) = ! input names
     &                         ( / 'ARINV' , 'MBINV' , 'PTINV'  / )

C...........   Other local variables

        INTEGER       J    ! counter and indices
        INTEGER       LCAT   ! length of CATEGORY string

        CHARACTER*300          MESG        ! message buffer 

        CHARACTER*16 :: PROGNAME = 'OPENINVIN' ! program name

C***********************************************************************
C   begin body of subroutine OPENINVIN

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
        MESG = 'Enter name of the RAW ' // CATEGORY( 1:LCAT ) // 
     &         ' INVENTORY ' // 'file'

C.........  Get file name an open raw input inventory file
        IDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., INNAME, PROGNAME )

C.........  Open category-specific inputs
        SELECT CASE( CATEGORY )
        CASE( 'POINT' )

C.............  Get file name for input replacement stack parameters file
            RDEV = PROMPTFFILE(
     &                 'Enter REPLACEMENT STACK PARAMETERS file',
     &                 .TRUE., .TRUE., 'PSTK', PROGNAME )

        END SELECT

C.........  Get file name for time zones files
        ZDEV = PROMPTFFILE( 
     &             'Enter logical name for TIME ZONE file',
     &             .TRUE., .TRUE., 'ZONES', PROGNAME )

C.........  Get file name for inventory pollutants codes/names
        MESG = 'Enter logical name for POLLUTANT CODES & NAMES file'
        PDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'SIPOLS', PROGNAME )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )
 
        END SUBROUTINE OPENINVIN

