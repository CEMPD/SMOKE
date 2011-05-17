
        SUBROUTINE OPENINVIN( CATEGORY, ADEV, DDEV, HDEV, RDEV, SDEV, 
     &                        XDEV, EDEV, PDEV, ZDEV, CDEV, ODEV, UDEV,
     &                        YDEV, ENAME, INNAME, IDNAME, IHNAME )

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
C...........  This module contains the information about the source category
        USE MODINFO, ONLY: NMAP, MAPNAM, MAPFIL

C............ This module contains the cross-reference tables
        USE MODXREF, ONLY: PROC_HAPS 

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS and their descriptionsNRAWIN
        CHARACTER(2)       CRLF
        LOGICAL            ENVYN
        INTEGER            INDEX1
        INTEGER            PROMPTFFILE
        CHARACTER(NAMLEN3) PROMPTMFILE

        EXTERNAL        CRLF, ENVYN, INDEX1, PROMPTFFILE, PROMPTMFILE

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN) :: CATEGORY  ! source category
        INTEGER     , INTENT(OUT) :: ADEV      ! unit number for inven file
        INTEGER     , INTENT(OUT) :: DDEV      ! unit no. for day-specific file
        INTEGER     , INTENT(OUT) :: HDEV      ! unit no. for hr-specific file
        INTEGER     , INTENT(OUT) :: RDEV      ! unit no. for stack replacements
        INTEGER     , INTENT(OUT) :: SDEV      ! unit no. for optional inven in
        INTEGER     , INTENT(OUT) :: XDEV      ! unit no. for vmt mix file
        INTEGER     , INTENT(OUT) :: EDEV      ! unit no. for speeds file
        INTEGER     , INTENT(OUT) :: PDEV      ! unit no. for inven data table
        INTEGER     , INTENT(OUT) :: ZDEV      ! unit no. for time zones
        INTEGER     , INTENT(OUT) :: CDEV      ! unit no. for SCCs description
        INTEGER     , INTENT(OUT) :: ODEV      ! unit no. for ORIS description
        INTEGER     , INTENT(OUT) :: UDEV      ! unit no. for non-HAP inclusions/exclusions
        INTEGER     , INTENT(OUT) :: YDEV      ! unit no. for area-to-point
        CHARACTER(*), INTENT(OUT) :: ENAME     ! optional netCDF inven input
        CHARACTER(*), INTENT(OUT) :: INNAME    ! average inventory name
        CHARACTER(*), INTENT(OUT) :: IDNAME    ! day-specific inventory 
        CHARACTER(*), INTENT(OUT) :: IHNAME    ! hour-specific inventory name 

C...........   Other local variables
        INTEGER       IDEV   ! tmp unit number if ENAME is map file
        INTEGER       IOS    ! i/o status
        INTEGER       J      ! counter and indices
        INTEGER       LCAT   ! length of CATEGORY string

        LOGICAL    :: CFLAG = .FALSE.  ! true: open area-to-point file
        LOGICAL    :: DFLAG = .FALSE.  ! true: open day-specific file
        LOGICAL    :: GFLAG = .FALSE.  ! true: open gridded I/O API inventory
        LOGICAL    :: HFLAG = .FALSE.  ! true: open hour-specific file
        LOGICAL    :: IFLAG = .FALSE.  ! true: open annual/average inventory
        LOGICAL    :: NFLAG = .FALSE.  ! true: open non-HAP inclusions/exclusions
        LOGICAL    :: MFLAG = .FALSE.  ! true: treat all sources as treated
        LOGICAL    :: SFLAG = .FALSE.  ! true: open speeds file
        LOGICAL    :: XFLAG = .FALSE.  ! true: open VMT mix file

        CHARACTER(NAMLEN3) ANAME
        CHARACTER(NAMLEN3) NAMBUF      ! file name buffer
        CHARACTER(NAMLEN3) INAME       ! tmp name for inven file of unknown fmt
        CHARACTER(256)     MESG        ! message buffer 

        CHARACTER(16) :: PROGNAME = 'OPENINVIN' ! program name

C***********************************************************************
C   begin body of subroutine OPENINVIN

C.........  Get value of these controls from the environment
        MESG = 'Import average inventory data'
        IFLAG = ENVYN ( 'IMPORT_AVEINV_YN', MESG, .TRUE., IOS )

        IF ( CATEGORY .EQ. 'POINT' ) THEN
            MESG = 'Import day-specific data'
            DFLAG = ENVYN ( 'DAY_SPECIFIC_YN', MESG, .FALSE., IOS )

            MESG = 'Import hour-specific data'
            HFLAG = ENVYN ( 'HOUR_SPECIFIC_YN', MESG, .FALSE., IOS )
        END IF
        
        MESG = 'Define the processing method of combining haradous ' //
     &         'air pollutants with criteria VOC.'
        CALL ENVSTR( 'SMK_PROCESS_HAPS', MESG, ' ', PROC_HAPS, IOS )

        SELECT CASE( PROC_HAPS )
        CASE( 'ALL' )
            MESG = 'Treat all sources as integrate sources to compute'//
     &             ' NONHAP[VOC|TOG].'
            NFLAG = .FALSE.

        CASE( 'NONE' )
            MESG = 'Treat all sources as non-integrate sources.'
            NFLAG = .FALSE.

        CASE( 'PARTIAL' )
            MESG = 'Partially treat sources as either integrate or ' //
     &             'non-integrate sources.'
            NFLAG = .TRUE.

        CASE DEFAULT

C.............  Check older flag setting to support backward comparability
            NFLAG = ENVYN ( 'SMK_NHAPEXCLUDE_YN', ' ', .FALSE., IOS )
            IF( NFLAG ) THEN
                MESG = 'WARNING: SMK_NHAPEXCLUDE_YN is no longer '//
     &             'supported.'//CRLF()//BLANK10 //'Please use SMK_'//
     &             'PROCESS_HAPS [ALL|NONE|PARTIAL] instead'
                CALL M3MSG2( MESG )

                MESG = 'Partially treat sources as either integrate '//
     &                 'or non-integrate sources.'
                PROC_HAPS = 'PARTIAL'
                NFLAG = .TRUE.

            ELSE
                MESG = 'No processing of combining criteria VOC with '//
     &                 'hazardous air pollutants (HAP).'
                NFLAG = .FALSE.

            END IF

        END SELECT
        
        CALL M3MSG2( MESG )
        
        IF ( CATEGORY .EQ. 'AREA' ) THEN
            MESG = 'Read and use area-to-point factors file'
            CFLAG = ENVYN ( 'SMK_ARTOPNT_YN', MESG, .FALSE., IOS )
            
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

C.........   Get NetCDF gridded inventory file (without SCCs)
        IF( GFLAG ) THEN

            MESG = 'Enter logical name of the GRIDDED ' // 
     &             CATEGORY( 1:LCAT ) // ' INVENTORY ' // 'file'

            ENAME = PROMPTMFILE( MESG, FSREAD3, ENAME, PROGNAME )

C.........  Get ASCII file name and open average input inventory file when 
C           inventory is to be imported
        ELSE IF( IFLAG ) THEN
            MESG = 'Enter logical name of the RAW ' // 
     &             CATEGORY( 1:LCAT ) // ' AVERAGE INVENTORY ' // 'file'

            ADEV = PROMPTFFILE( MESG, .TRUE., .TRUE., INNAME, PROGNAME )

C.........  If SMOKE inventory already exists, open files for reading later
        ELSE

C.............  Get input inventory file names given source category
C.............  Use NAMBUF for HP safety
            CALL GETINAME( CATEGORY, ENAME, ANAME )

C.........  Prompt for and open inventory file 
            INAME = ENAME
            MESG = 'Enter logical name for the MAP INVENTORY file'
            IDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., INAME, PROGNAME )

C.............  Open and read map file
            CALL RDINVMAP( INAME, IDEV, ENAME, ANAME, SDEV )

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

C.........  Get SCC descriptions if needed for reporting purposes
        IF ( HFLAG .OR. CFLAG ) THEN
            CDEV = PROMPTFFILE(
     &             'Enter logical name for SCC DESCRIPTION file ',
     &             .TRUE., .TRUE., 'SCCDESC', PROGNAME )
        END IF

C.........  Get ORIS descriptions file
        IF( HFLAG ) THEN
            ODEV = PROMPTFFILE(
     &             'Enter logical name for ORIS DESCRIPTION file ',
     &             .TRUE., .TRUE., 'ORISDESC', PROGNAME )
        END IF

C.........  Get file name for non-HAP exclusion file
        IF( NFLAG ) THEN
            MESG = 'Enter logical name for NHAPEXCLUDE file'
            UDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'NHAPEXCLUDE',
     &                          PROGNAME )
        END IF

C.........  Get file name for area-to-point factors file
        IF( CFLAG ) THEN
            MESG = 'Enter logical name for AREA-TO-POINT FACTORS file'
            YDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'ARTOPNT',
     &                          PROGNAME )
        END IF

C.........  Get file name for inventory pollutants codes/names
        MESG = 'Enter logical name for INVENTORY DATA TABLE file'
        PDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'INVTABLE',
     &                      PROGNAME )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )
 
        END SUBROUTINE OPENINVIN

