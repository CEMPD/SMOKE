
        SUBROUTINE OPENINVIN( CATEGORY, IDEV, DDEV, HDEV, RDEV, SDEV, 
     &                        PDEV, VDEV, ZDEV, ENAME, INNAME, 
     &                        IDNAME, IHNAME )

C***********************************************************************
C  subroutine body starts at line 114
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
        LOGICAL                ENVYN
        INTEGER                INDEX1
        INTEGER                PROMPTFFILE
        CHARACTER(LEN=NAMLEN3) PROMPTMFILE

        EXTERNAL        ENVYN, INDEX1, PROMPTFFILE, PROMPTMFILE

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN) :: CATEGORY  ! source category
        INTEGER     , INTENT(OUT) :: IDEV      ! unit number for inven file
        INTEGER     , INTENT(OUT) :: DDEV      ! unit no. for day-specific file
        INTEGER     , INTENT(OUT) :: HDEV      ! unit no. for hr-specific file
        INTEGER     , INTENT(OUT) :: RDEV      ! unit no. for stack replacements
        INTEGER     , INTENT(OUT) :: SDEV      ! unit no. for optional inven in
        INTEGER     , INTENT(OUT) :: PDEV      ! unit no. for pol codes & names
        INTEGER     , INTENT(OUT) :: VDEV      ! unit no. for activity names
        INTEGER     , INTENT(OUT) :: ZDEV      ! unit no. for time zones
        CHARACTER(*), INTENT(OUT) :: ENAME     ! optional netCDF inven input
        CHARACTER(*), INTENT(OUT) :: INNAME    ! average inventory name
        CHARACTER(*), INTENT(OUT) :: IDNAME    ! day-specific inventory 
        CHARACTER(*), INTENT(OUT) :: IHNAME    ! hour-specific inventory name 

C...........   LOCAL PARAMETERS
        INTEGER    , PARAMETER :: NCAT = 3          ! number src categories
        CHARACTER*6, PARAMETER :: CATLIST( NCAT ) = ! src categories
     &                         ( / 'AREA  ', 'MOBILE', 'POINT ' / )

        CHARACTER*5, PARAMETER :: ANAMLIST( NCAT ) = ! ave data input names
     &                         ( / 'ARINV' , 'MBINV' , 'PTINV'  / )

        CHARACTER*5, PARAMETER :: DNAMLIST( NCAT ) = ! daily data input names
     &                         ( / 'ARDAY' , 'MBDAY' , 'PTDAY'  / )

        CHARACTER*6, PARAMETER :: HNAMLIST( NCAT ) = ! hourly data input names
     &                         ( / 'ARHOUR' , 'MBHOUR' , 'PTHOUR'  / )

C...........   Other local variables

        INTEGER       IOS    ! i/o status
        INTEGER       J      ! counter and indices
        INTEGER       LCAT   ! length of CATEGORY string

        LOGICAL       CFLAG  ! true: use SIPOLS file
        LOGICAL       DFLAG  ! true: import day-specific file
        LOGICAL       HFLAG  ! true: import hour-specific file
        LOGICAL       IFLAG  ! true: import annual/average inventory
        LOGICAL       VFLAG  ! true: use ACTVNAMS file

        CHARACTER(LEN=NAMLEN3) ANAME
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

        MESG = 'Import average inventory data'
        IFLAG = ENVYN ( 'IMPORT_AVEINV_YN', MESG, .TRUE., IOS )

        MESG = 'Import day-specific data'
        DFLAG = ENVYN ( 'DAY_SPECIFIC_YN', MESG, .FALSE., IOS )

        MESG = 'Import hour-specific data'
        HFLAG = ENVYN ( 'HOUR_SPECIFIC_YN', MESG, .FALSE., IOS )

C.........  Abort if no settings set to read data
        IF( .NOT. IFLAG .AND. 
     &      .NOT. DFLAG .AND.
     &      .NOT. HFLAG       ) THEN

            MESG = 'ERROR: Environment settings indicate no files ' //
     &             'are to be read'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF

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
 
            INNAME = ANAMLIST( J )
            IDNAME = DNAMLIST( J )
            IHNAME = HNAMLIST( J )

        END IF

C.........  Initialize unit numbers and logical file names
        ENAME = ' '
        IDEV = 0
        DDEV = 0
        HDEV = 0
        RDEV = 0
        SDEV = 0
        PDEV = 0
        VDEV = 0
        ZDEV = 0

C.........  Get file name and open average input inventory file when inventory
C           is to be imported
        IF( IFLAG ) THEN
            MESG = 'Enter logical name of the RAW ' // 
     &             CATEGORY( 1:LCAT ) // ' AVERAGE INVENTORY ' // 'file'

            IDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., INNAME, PROGNAME )

C.........  If inventory already exists, open files for reading later
        ELSE

C.............  Get input inventory file names given source category
            CALL GETINAME( CATEGORY, ENAME, ANAME )

            ENAME = PROMPTMFILE( 
     &              'Enter logical name for the I/O API INVENTORY file',
     &              FSREAD3, ENAME, PROGNAME )

            SDEV = PROMPTFFILE( 
     &               'Enter logical name for the ASCII INVENTORY file',
     &               .TRUE., .TRUE., ANAME, PROGNAME )
        END IF

C.........  Get file name and open daily input inventory file
        IF( DFLAG ) THEN
            MESG = 'Enter logical name of the RAW ' //
     &             CATEGORY( 1:LCAT ) // ' DAILY INVENTORY ' // 'file'

            DDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., IDNAME, PROGNAME )
        END IF

C.........  Get file name and open daily input inventory file
        IF( HFLAG ) THEN
            MESG = 'Enter logical name of the RAW ' //
     &             CATEGORY( 1:LCAT ) // ' HOURLY INVENTORY ' // 'file'

            HDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., IHNAME, PROGNAME )
        END IF

        IF( IFLAG ) THEN

C.............  Open category-specific inputs
            SELECT CASE( CATEGORY )
            CASE( 'MOBILE' ) 

C.................  Get file name for converting road-class to road type & 
C                   vehicle type name to vehicle type number.
        	MESG = 'Enter logical name for MOBILE CODES file'
         	RDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'MCODES', 
     &                              PROGNAME )

            CASE( 'POINT' )

C.................  Get file name for input replacement stack parameters file
                MESG = 'Enter logical name for REPLACEMENT STACK ' //
     &                 'PARAMETERS file'
                RDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'PSTK', 
     &                              PROGNAME )

            END SELECT

C.............  Get file name for country, state, and county file, with time 
C               zones
            ZDEV = PROMPTFFILE(
     &             'Enter logical name for COUNTRY, STATE, AND ' //
     &             'COUNTY file', .TRUE., .TRUE., 'COSTCY', PROGNAME )

C.............  Get file name for inventory pollutants codes/names
            IF( CFLAG ) THEN 
        	MESG = 'Enter logical name for POLLUTANT CODES & ' //
     &                 'NAMES file'
        	PDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'SIPOLS',
     &                              PROGNAME )
            END IF

C.............  Get file name for inventory pollutants codes/names
            IF( VFLAG ) THEN
                MESG = 'Enter logical name for ACTIVITY CODES & ' //
     &                 'NAMES file'
                VDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'ACTVNAMS',
     &                              PROGNAME )
            END IF

        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )
 
        END SUBROUTINE OPENINVIN

