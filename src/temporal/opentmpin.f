
        SUBROUTINE OPENTMPIN( MODELNAM, UFLAG, ENAME, ANAME, DNAME, 
     &                        HNAME, GNAME, SDEV, XDEV, RDEV, CDEV, 
     &                        HDEV, TDEV, MDEV, EDEV, PYEAR )

C***********************************************************************
C  subroutine body starts at line 123
C
C  DESCRIPTION:
C      This subroutine opens the file or files for output from the tmppoint
C      program. It also populates the MODINFO module with the source-category
C      specific information. It also determines the name of the temperature
C      variable.
C 
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Created 7/99 by M. Houyoux
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

C...........   MODULES for public variables   
C...........   This module is the derived meteorology data for emission factors
        USE MODMET, ONLY:

C...........  This module contains the information about the source category
        USE MODINFO, ONLY: CATEGORY, CRL, NSRC, NIACT, INVPIDX

C.........  This module contains the global variables for the 3-d grid
        USE MODGRID, ONLY:

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   ! emissions constat parameters
        INCLUDE 'PARMS3.EXT'    !  i/o api parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.
        INCLUDE 'SETDECL.EXT'   !  FileSetAPI variables and functions
        INCLUDE 'FLTERR.EXT'    ! error filter statement function

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER*2            CRLF
        LOGICAL                ENVYN
        CHARACTER(LEN=IOVLEN3) GETCFDSC
        INTEGER                GETIFDSC
        INTEGER                INDEX1
        INTEGER                PROMPTFFILE
        CHARACTER(LEN=NAMLEN3) PROMPTMFILE

        EXTERNAL        CRLF, ENVYN, GETIFDSC, GETCFDSC, INDEX1, 
     &                  PROMPTFFILE, PROMPTMFILE

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT    (IN) :: MODELNAM ! name for EF model
        LOGICAL     , INTENT    (IN) :: UFLAG    ! use uniform temporal profile
        CHARACTER(*), INTENT(IN OUT) :: ENAME ! name for I/O API inven input
        CHARACTER(*), INTENT(IN OUT) :: ANAME ! name for ASCII inven input 
        CHARACTER(*), INTENT   (OUT) :: DNAME ! day-spec file
        CHARACTER(*), INTENT   (OUT) :: HNAME ! hour-spec file
        CHARACTER(*), INTENT   (OUT) :: GNAME ! ungridding matrix
        INTEGER     , INTENT   (OUT) :: SDEV  ! unit no.: ASCII inven file
        INTEGER     , INTENT   (OUT) :: XDEV  ! unit no.: x-ref file
        INTEGER     , INTENT   (OUT) :: RDEV  ! unit no.: tmprl profile file
        INTEGER     , INTENT   (OUT) :: CDEV  ! unit no.: region codes file
        INTEGER     , INTENT   (OUT) :: HDEV  ! unit no.: holidays file
        INTEGER     , INTENT   (OUT) :: TDEV  ! unit no.: emissions process file
        INTEGER     , INTENT   (OUT) :: MDEV  ! unit no.: mobile codes file
        INTEGER     , INTENT   (OUT) :: EDEV  ! unit no.: emission factor file list
        INTEGER     , INTENT   (OUT) :: PYEAR ! projected year

C...........   Other local variables
        INTEGER         IDEV        ! tmp unit number if ENAME is map file
        INTEGER         IOS         ! status from environment variables
        INTEGER         J           ! index
        INTEGER         L           ! string length

        LOGICAL         DFLAG       ! true: day-specific  file available
        LOGICAL      :: EFLAG = .FALSE.  ! true: error found
        LOGICAL         HFLAG       ! true: hour-specific file available
        LOGICAL         OFLAG       ! true: average day emissions needed
        LOGICAL         XFLAG       ! true: use daylight time exemptions file

        CHARACTER*16    INAME       ! tmp name for inven file of unknown fmt
        CHARACTER*256   MESG        ! message buffer 

        CHARACTER(LEN=NAMLEN3)  NAMBUF ! file name buffer

        CHARACTER*16 :: PROGNAME = 'OPENTMPIN' ! program name

C***********************************************************************
C   begin body of subroutine OPENTMPIN

C.........  Get environment variables that control program behavior
        IF ( CATEGORY .EQ. 'POINT' ) THEN
            DFLAG = ENVYN( 'DAY_SPECIFIC_YN', 'Use day-specific data',
     &                      .FALSE., IOS )

            HFLAG = ENVYN( 'HOUR_SPECIFIC_YN', 'Use hour-specific data',
     &                     .FALSE., IOS )
        END IF

        OFLAG = ENVYN( 'SMK_AVEDAY_YN', MESG, .FALSE., IOS )

C.........  Prompt for and open inventory file
        INAME = ENAME 
        MESG = 'Enter logical name for the MAP INVENTORY file'
        IDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., INAME, PROGNAME )

C.........  Open and read map file
        CALL RDINVMAP( INAME, IDEV, ENAME, ANAME, SDEV )

        IF( DFLAG ) THEN
            NAMBUF = PROMPTMFILE( 
     &               'Enter logical name for DAY-SPECIFIC file',
     &               FSREAD3, CRL // 'DAY', PROGNAME )
            DNAME = NAMBUF
        END IF

        IF( HFLAG ) THEN
            NAMBUF = PROMPTMFILE( 
     &               'Enter logical name for HOUR-SPECIFIC file',
     &               FSREAD3, CRL // 'HOUR', PROGNAME )
            HNAME = NAMBUF
        END IF

        IF( .NOT. UFLAG ) THEN
            XDEV = PROMPTFFILE( 
     &           'Enter logical name for TEMPORAL CROSS-REFERENCE file',
     &           .TRUE., .TRUE., CRL // 'TREF', PROGNAME )

            RDEV = PROMPTFFILE( 
     &           'Enter logical name for TEMPORAL PROFILES file',
     &           .TRUE., .TRUE., CRL // 'TPRO', PROGNAME )
        END IF

C.........  Store source-category-specific header information, 
C           including the inventory pollutants in the file (if any).  Note that 
C           the I/O API head info is passed by include file and the
C           results are stored in module MODINFO.
C.........  Set average day emissions flag (INVPIDX)
        IF( OFLAG ) INVPIDX = 1
        CALL GETSINFO( ENAME )

        PYEAR   = GETIFDSC( FDESC3D, '/PROJECTED YEAR/', .FALSE. )

C.............  Store non-category-specific header information
        NSRC = NROWS3D

C.........  Open region codes file for determining daylight savings time status
        CDEV = PROMPTFFILE(
     &             'Enter logical name for COUNTRY, STATE, AND ' //
     &             'COUNTY file', .TRUE., .TRUE., 'COSTCY', PROGNAME )
        
C.........  Open holidays file for determining holidays by region
        HDEV = PROMPTFFILE(
     &             'Enter logical name for HOLIDAYS file',
     &             .TRUE., .TRUE., 'HOLIDAYS', PROGNAME )

C.........  Open additional files for when activity data are in the inventory.
C.........  NOTE - this structure currently assumes that all of the
C           files needed for using MOBILE5 emission factors would be needed
C           in all cases.  For driving other emission factor models, other
C           logic would need to be used that evaluates the emission factor
C           model assigned to each activity, and opens files depending on the
C           emission factor model.
C.........  Use NAMBUF for the HP
        IF( NIACT .GT. 0 ) THEN

            L = LEN_TRIM( MODELNAM )

            NAMBUF= PROMPTMFILE( 
     &              'Enter logical name for UNGRIDDING MATRIX file',
     &              FSREAD3, CRL // 'UMAT', PROGNAME )
            GNAME = NAMBUF
 
C.............  Get the header description from the ungridding matrix file
            CALL RETRIEVE_IOAPI_HEADER( GNAME )

C.............  Check the number of sources in the ungridding matrix
            CALL CHKSRCNO( 'mobile', 'MUMAT', NROWS3D, NSRC, EFLAG )
            
            TDEV = PROMPTFFILE( 
     &             'Enter logical name for EMISSION PROCESSES file',
     &             .TRUE., .TRUE., CRL // 'EPROC', PROGNAME )
     
            EDEV = PROMPTFFILE(
     &             'Enter logical name for EMISSION FACTORS LIST file',
     &             .TRUE., .TRUE., CRL // 'EFLIST', PROGNAME )
     
        END IF

C.........  Open files that are specific to mobile sources
        IF( CATEGORY .EQ. 'MOBILE' ) THEN

            MDEV = PROMPTFFILE( 
     &             'Enter logical name for MOBILE CODES file',
     &             .TRUE., .TRUE., 'MCODES', PROGNAME )

        END IF

C.........  Abort if error was found
        IF ( EFLAG ) THEN

            MESG = 'Problem with input files'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94000   FORMAT( I2.2 )
 
94010   FORMAT( 10( A, :, I8, :, 1X ) )

C******************  INTERNAL SUBPROGRAMS   ******************************

        CONTAINS

C.............  This internal subprogram tries to retrieve the I/O API header
C               and aborts if it was not successful
            SUBROUTINE RETRIEVE_IOAPI_HEADER( FILNAM )

C.............  Subprogram arguments
            CHARACTER(*) FILNAM

C----------------------------------------------------------------------

            IF ( .NOT. DESC3( FILNAM ) ) THEN

                MESG = 'Could not get description of file "' //
     &                 FILNAM( 1:LEN_TRIM( FILNAM ) ) // '"'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            END IF

            END SUBROUTINE RETRIEVE_IOAPI_HEADER

        END SUBROUTINE OPENTMPIN

