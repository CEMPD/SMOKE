
        SUBROUTINE OPENINVIN( CATEGORY, IDEV, DDEV, HDEV, RDEV, SDEV, 
     &                        XDEV, EDEV, PDEV, VDEV, ZDEV, CDEV, ODEV, 
     &                        ENAME, INNAME, IDNAME, IHNAME )

C***********************************************************************
C  subroutine body starts at line 119
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
C**************************************************************************
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
C***************************************************************************

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
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
        INTEGER     , INTENT(OUT) :: XDEV      ! unit no. for vmt mix file
        INTEGER     , INTENT(OUT) :: EDEV      ! unit no. for speeds file
        INTEGER     , INTENT(OUT) :: PDEV      ! unit no. for pol codes & names
        INTEGER     , INTENT(OUT) :: VDEV      ! unit no. for activity names
        INTEGER     , INTENT(OUT) :: ZDEV      ! unit no. for time zones
        INTEGER     , INTENT(OUT) :: CDEV      ! unit no. for SCCs description
        INTEGER     , INTENT(OUT) :: ODEV      ! unit no. for ORIS description
        CHARACTER(*), INTENT(OUT) :: ENAME     ! optional netCDF inven input
        CHARACTER(*), INTENT(OUT) :: INNAME    ! average inventory name
        CHARACTER(*), INTENT(OUT) :: IDNAME    ! day-specific inventory 
        CHARACTER(*), INTENT(OUT) :: IHNAME    ! hour-specific inventory name 

C...........   Other local variables

        INTEGER       IOS    ! i/o status
        INTEGER       J      ! counter and indices
        INTEGER       LCAT   ! length of CATEGORY string

        LOGICAL       CFLAG  ! true: use SIPOLS file
        LOGICAL       DFLAG  ! true: import day-specific file
        LOGICAL       GFLAG  ! true: import gridded I/O API inventory
        LOGICAL       HFLAG  ! true: import hour-specific file
        LOGICAL       IFLAG  ! true: import annual/average inventory
        LOGICAL       SFLAG  ! true: import speeds file
        LOGICAL       VFLAG  ! true: use ACTVNAMS file
        LOGICAL       XFLAG  ! true: import VMT mix file

        CHARACTER(LEN=NAMLEN3) ANAME
        CHARACTER(LEN=NAMLEN3) NAMBUF      ! file name buffer
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

        IF ( CATEGORY .EQ. 'POINT' ) THEN
            MESG = 'Import day-specific data'
            DFLAG = ENVYN ( 'DAY_SPECIFIC_YN', MESG, .FALSE., IOS )

            MESG = 'Import hour-specific data'
            HFLAG = ENVYN ( 'HOUR_SPECIFIC_YN', MESG, .FALSE., IOS )
        END IF

        IF ( CATEGORY .EQ. 'AREA' ) THEN
            MESG = 'Import gridded I/O API inventory data'
            GFLAG = ENVYN ( 'IMPORT_GRDIOAPI_YN', MESG, .FALSE., IOS )
        END IF

        IF ( CATEGORY .EQ. 'MOBILE' ) THEN
            MESG = 'Import VMT mix data'
            XFLAG = ENVYN ( 'IMPORT_VMTMIX_YN', MESG, .FALSE., IOS )

            MESG = 'Import mobile speeds data'
            SFLAG = ENVYN ( 'IMPORT_SPEEDS_YN', MESG, .FALSE., IOS )
        END IF

C.........  Make sure VMT mix and speeds will only be imported for mobile 
C           sources
        IF( CATEGORY .NE. 'MOBILE' ) THEN
            XFLAG = .FALSE.
            SFLAG = .FALSE.

C.........  Make sure gridded point source file is not attempted
        ELSE IF ( ( CATEGORY .EQ. 'POINT' .OR. 
     &              CATEGORY .EQ. 'MOBILE'    ) .AND. 
     &            GFLAG                               ) THEN
            MESG = 'Cannot import gridded mobile or point source data.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF

C.........  When gridded data are imported, override other settings
        IF( GFLAG ) THEN
            CFLAG = .FALSE.
            VFLAG = .FALSE.
            IFLAG = .FALSE.
            DFLAG = .FALSE.
            HFLAG = .FALSE.
            XFLAG = .FALSE.
            SFLAG = .FALSE.
        END IF

C.........  Abort if no settings set to read data
        IF( .NOT. IFLAG .AND. 
     &      .NOT. DFLAG .AND.
     &      .NOT. HFLAG .AND.
     &      .NOT. GFLAG      ) THEN

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
            ENAME  = GNAMLIST( J )
            IDNAME = DNAMLIST( J )
            IHNAME = HNAMLIST( J )

        END IF

C.........
        IF( GFLAG ) THEN

            MESG = 'Enter logical name of the GRIDDED ' // 
     &             CATEGORY( 1:LCAT ) // ' INVENTORY ' // 'file'

            ENAME = PROMPTMFILE( MESG, FSREAD3, ENAME, PROGNAME )

C.........  Get ASCII file name and open average input inventory file when 
C           inventory is to be imported
        ELSE IF( IFLAG ) THEN
            MESG = 'Enter logical name of the RAW ' // 
     &             CATEGORY( 1:LCAT ) // ' AVERAGE INVENTORY ' // 'file'

            IDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., INNAME, PROGNAME )

C.........  If SMOKE inventory already exists, open files for reading later
        ELSE

C.............  Get input inventory file names given source category
C.............  Use NAMBUF for HP safety
            CALL GETINAME( CATEGORY, ENAME, ANAME )

            NAMBUF = PROMPTMFILE( 
     &              'Enter logical name for the I/O API INVENTORY file',
     &              FSREAD3, ENAME, PROGNAME )
            ENAME = NAMBUF

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

C.........  Get VMT Mix file
        IF( XFLAG ) THEN

            MESG = 'Enter logical name for VMT MIX file'
            XDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'VMTMIX', 
     &                          PROGNAME )

       END IF

C.........  Get speeds file
        IF( SFLAG ) THEN

            MESG = 'Enter logical name for MOBILE SPEEDS file'
            EDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'MSPEEDS', 
     &                          PROGNAME )

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

        END IF

C.........  Get file name for country, state, and county file, with time 
C           zones
        IF( .NOT. GFLAG ) THEN
            ZDEV = PROMPTFFILE(
     &             'Enter logical name for COUNTRY, STATE, AND ' //
     &             'COUNTY file', .TRUE., .TRUE., 'COSTCY', PROGNAME )
        END IF

C.........  Get list of powerplant SCCs, in case needed for reporting CEM
C           matches with the inventory
        IF ( HFLAG ) THEN
            CDEV = PROMPTFFILE(
     &             'Enter logical name for SCC DESCRIPTION file ',
     &             .TRUE., .TRUE., 'SCCDESC', PROGNAME )

            ODEV = PROMPTFFILE(
     &             'Enter logical name for ORIS DESCRIPTION file ',
     &             .TRUE., .TRUE., 'ORISDESC', PROGNAME )
        END IF

C.........  Get file name for inventory pollutants codes/names
        IF( CFLAG ) THEN 
            MESG = 'Enter logical name for POLLUTANT CODES & ' //
     &             'NAMES file'
            PDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'SIPOLS',
     &                          PROGNAME )
        END IF

C.........  Get file name for inventory pollutants codes/names
        IF( VFLAG ) THEN
            MESG = 'Enter logical name for ACTIVITY CODES & ' //
     &             'NAMES file'
            VDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'ACTVNAMS',
     &                          PROGNAME )
        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )
 
        END SUBROUTINE OPENINVIN

