
        SUBROUTINE OPENTMPIN( MODELNAM, UFLAG, ENAME, ANAME, DNAME, 
     &                        HNAME, FNAME, NNAME, MNAME, GNAME, WNAME,
     &                        TVARNAME, SDEV, XDEV, RDEV, UDEV, FDEV,
     &                        TDEV, MDEV )

C***********************************************************************
C  subroutine body starts at line
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

C...........   MODULES for public variables   
C...........   This module is the derived meteorology data for emission factors
        USE MODMET

C...........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   ! emissions constat parameters
        INCLUDE 'PARMS3.EXT'    ! I/O API parameters
        INCLUDE 'IODECL3.EXT'   ! I/O API function declarations
        INCLUDE 'FDESC3.EXT'    ! I/O API file description data structures.
        INCLUDE 'FLTERR.EXT'    ! error filter statement function

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER*2            CRLF
        LOGICAL                ENVYN
        CHARACTER(LEN=IOVLEN3) GETCFDSC
        INTEGER                INDEX1
        INTEGER                PROMPTFFILE
        CHARACTER(LEN=NAMLEN3) PROMPTMFILE

        EXTERNAL        CRLF, ENVYN, GETCFDSC, INDEX1, PROMPTFFILE, 
     &                  PROMPTMFILE

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT    (IN) :: MODELNAM ! name for EF model
        LOGICAL     , INTENT    (IN) :: UFLAG    ! use uniform temporal profile
        CHARACTER(*), INTENT(IN OUT) :: ENAME ! name for I/O API inven input
        CHARACTER(*), INTENT(IN OUT) :: ANAME ! name for ASCII inven input 
        CHARACTER(*), INTENT   (OUT) :: DNAME ! day-spec file
        CHARACTER(*), INTENT   (OUT) :: HNAME ! hour-spec file
        CHARACTER(*), INTENT   (OUT) :: FNAME ! non-diurnal EF file
        CHARACTER(*), INTENT   (OUT) :: NNAME ! diurnal EF file
        CHARACTER(*), INTENT   (OUT) :: MNAME ! surface temperature file
        CHARACTER(*), INTENT   (OUT) :: GNAME ! ungridding matrix
        CHARACTER(*), INTENT   (OUT) :: WNAME ! ungridded min/max temperatures
        CHARACTER(*), INTENT   (OUT) :: TVARNAME ! tmpr variable name
        INTEGER     , INTENT   (OUT) :: SDEV  ! unit no.: ASCII inven file
        INTEGER     , INTENT   (OUT) :: XDEV  ! unit no.: x-ref file
        INTEGER     , INTENT   (OUT) :: RDEV  ! unit no.: tmprl profile file
        INTEGER     , INTENT   (OUT) :: UDEV  ! unit no.: optional UAM elev srcs
        INTEGER     , INTENT   (OUT) :: FDEV  ! unit no.: EF x-ref file
        INTEGER     , INTENT   (OUT) :: TDEV  ! unit no.: speciation list file
        INTEGER     , INTENT   (OUT) :: MDEV  ! unit no.: mobile codes file

C...........   LOCAL PARAMETERS
        CHARACTER*50  SCCSW          ! SCCS string with version number at end

        PARAMETER   ( SCCSW   = '@(#)$Id$'
     &              )

C...........   Other local variables

        INTEGER         IOS         ! status from environment variables
        INTEGER         J           ! index
        INTEGER         L           ! string length

        REAL            MAX1, MAX2  ! tmp maximum temperature values
        REAL            MIN1, MIN2  ! tmp minimum temperature values

        LOGICAL         DFLAG       ! day-specific  file available
        LOGICAL         EFLAG       ! true: error found
        LOGICAL         HFLAG       ! hour-specific file available

        CHARACTER*16    MNAME0      ! default gridded temperature file name
        CHARACTER*300   MESG        ! message buffer 

        CHARACTER*16 :: PROGNAME = 'OPENTMPIN' ! program name

C***********************************************************************
C   begin body of subroutine OPENTMPIN

C.........  Get environment variables that control program behavior
        DFLAG = ENVYN ( 'DAY_SPECIFIC_YN', 'Use day-specific data',
     &                   .FALSE., IOS )

        HFLAG = ENVYN ( 'HOUR_SPECIFIC_YN', 'Use hour-specific data',
     &                   .FALSE., IOS )

C.........  Prompt for and open input I/O API and ASCII files
        ENAME = PROMPTMFILE( 
     &          'Enter logical name for the I/O API INVENTORY file',
     &          FSREAD3, ENAME, PROGNAME )

        SDEV = PROMPTFFILE( 
     &           'Enter logical name for the ASCII INVENTORY file',
     &           .TRUE., .TRUE., ANAME, PROGNAME )

        IF( DFLAG ) DNAME = PROMPTMFILE( 
     &          'Enter logical name for DAY-SPECIFIC file',
     &          FSREAD3, CRL // 'DAY', PROGNAME )

        IF( HFLAG ) HNAME = PROMPTMFILE( 
     &          'Enter logical name for HOUR-SPECIFIC file',
     &          FSREAD3, CRL // 'HOUR', PROGNAME )

        IF( .NOT. UFLAG ) THEN
            XDEV = PROMPTFFILE( 
     &           'Enter logical name for TEMPORAL CROSS-REFERENCE file',
     &           .TRUE., .TRUE., CRL // 'TREF', PROGNAME )

            RDEV = PROMPTFFILE( 
     &           'Enter logical name for TEMPORAL PROFILES file',
     &           .TRUE., .TRUE., CRL // 'TPRO', PROGNAME )
        END IF

C.........  Get source category information from the inventory files
C.........  Get header description of inventory file
C.........  Exit if getting the description fails 
        IF( .NOT. DESC3( ENAME ) ) THEN
            L = LEN_TRIM( ENAME )
            MESG = 'Could not get description of file "' //
     &             ENAME( 1:L ) // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C.........  Otherwise, store source-category-specific header information, 
C           including the inventory pollutants in the file (if any).  Note that 
C           the I/O API head info is passed by include file and the
C           results are stored in module MODINFO.
        ELSE

            CALL GETSINFO

C.............  Store non-category-specific header information
            NSRC = NROWS3D

        ENDIF        

C.........  Open additional files for when activity data are in the inventory.
C.........  NOTE - this structure currently assumes that all of the
C           files needed for using MOBILE5 emission factors would be needed
C           in all cases.  For driving other emission factor models, other
C           logic would need to be used that evaluates the emission factor
C           model assigned to each activity, and opens files depending on the
C           emission factor model.

        IF( NIACT .GT. 0 ) THEN

            L = LEN_TRIM( MODELNAM )

            MESG = 'Enter logical name for ' // MODELNAM( 1:L ) // 
     &             ' NON-DIURNAL EMISSION FACTORS file'
            FNAME = PROMPTMFILE( MESG, FSREAD3, CRL//'EFSND', PROGNAME )
    
            MESG = 'Enter logical name for ' // MODELNAM( 1:L ) // 
     &             ' DIURNAL EMISSION FACTORS file'
            NNAME = PROMPTMFILE( MESG, FSREAD3, CRL//'EFSD', PROGNAME )

            GNAME = PROMPTMFILE( 
     &              'Enter logical name for UNGRIDDING MATRIX file',
     &              FSREAD3, CRL // 'UMAT', PROGNAME )

            WNAME = PROMPTMFILE( 
     &              'Enter logical name for UNGRIDDED MIN/MAX ' //
     &              'TEMPERATURE file', FSREAD3, 'MINMAXT', PROGNAME )

C.............  Get the header description from the min/max temperatures file
            IF( .NOT. DESC3( WNAME ) ) THEN
                L = LEN_TRIM( WNAME )
        	MESG = 'Could not get description of file "' //
     &                 WNAME( 1:L ) // '"'
        	CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  Determine the temperature variable that was used to create the
C               min/max temperature file.
            TVARNAME = GETCFDSC( FDESC3D, '/T_VNAME/', .TRUE. )

C.............  Based on the temperature variable name, set the default name 
C               for the gridded temperature file
            MNAME0 = 'MET_CRO_2D'
            IF ( TVARNAME .EQ. 'TA' ) MNAME0 = 'MET_CRO_3D'

            MNAME = PROMPTMFILE( 
     &              'Enter logical name for SURFACE TEMPERATURE file',
     &              FSREAD3, MNAME0, PROGNAME )

C.............  Get the header of the gridded temperature file
            IF( .NOT. DESC3( MNAME ) ) THEN
                L = LEN_TRIM( MNAME )
        	MESG = 'Could not get description of file "' //
     &                 MNAME( 1:L ) // '"'
        	CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  Check to make sure the temperature variable of interest is in
C               the file.
            J = INDEX1( TVARNAME, NVARS3D, VNAME3D )             

C.............  If not, write a warning and get the temperature variable 
C               name from the environment
            IF( J .LE. 0 ) THEN

                CALL TEMPERATURE_WARNING

                MESG = 'NOTE: Getting temperature variable name from '//
     &                 'the environment...'
                CALL M3MSG2( MESG )

        	MESG = 'Temperature variable name'
        	CALL ENVSTR( 'TVARNAME', MESG, 'TEMP1P5', TVARNAME, IOS )

C.................  Write message if TVARNAME environment variable is undefined
                IF( IOS .LT. 0 ) THEN
                    MESG = 'NOTE: Using default temperature '//
     &                     'variable name from the environment.'
                    CALL M3MSG2( MESG )
                END IF

C.................  Ensure that the new temperature variable name of interest 
C                   is in the gridded temperature file
                J = INDEX1( TVARNAME, NVARS3D, VNAME3D )             

                IF( J .LE. 0 ) THEN

                    CALL TEMPERATURE_WARNING

                    MESG = 'ERROR: Could not get a variable name '//
     &                     'to use for gridded temperature file.'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF

            END IF

C.............  Compare the min/max temperature information in the min/max 
C               temperature file and in the emission factors files...
C.............  Retrieve header of min/max temperature file
            IF( .NOT. DESC3( WNAME ) ) THEN
                MESG = 'Could not get description for file ' // WNAME
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  Retrieve temperature ranges from min/max file header
C.............  Populate table of valid min/max temperatures in MODMET
            CALL TMPRINFO( .FALSE., 'BOTH' )

C.............  Store min/max temperatures for comparison
            MIN1 = MINT_MIN
            MIN2 = MINT_MAX
            MAX1 = MAXT_MIN
            MAX2 = MAXT_MAX

C.............  Retrieve header of non-diurnal emission factors file
            IF( .NOT. DESC3( FNAME ) ) THEN
                MESG = 'Could not get description for file ' // FNAME
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  Retrieve temperature ranges from non-diurnal EFs file header
            CALL TMPRINFO( .FALSE., 'NOMINMAX' )

C.............  Compare mint_min and maxt_max
            CALL COMPARE_TMPRS( 'NOMINMAX' )

C.............  Retrieve header of diurnal emission factors file
            IF( .NOT. DESC3( NNAME ) ) THEN
                MESG = 'Could not get description for file ' // NNAME
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  Retrieve temperature ranges from diurnal EFs file header
            CALL TMPRINFO( .FALSE., 'BOTH' )

C.............  Compare all min/max temperatures
            CALL COMPARE_TMPRS( 'BOTH' )

            FDEV = PROMPTFFILE( 
     &             'Enter logical name for EMISSION FACTORS X-REF file',
     &             .TRUE., .TRUE., CRL // 'PLIST', PROGNAME )

            TDEV = PROMPTFFILE( 
     &             'Enter logical name for EMISSION PROCESSES file',
     &             .TRUE., .TRUE., CRL // 'EPROC', PROGNAME )

        END IF

C.........  Open files that are specific to mobile sources
        IF( CATEGORY .EQ. 'MOBILE' ) THEN

            MDEV = PROMPTFFILE( 
     &             'Enter logical name for MOBILE CODES file',
     &             .TRUE., .TRUE., 'MCODES', PROGNAME )

        END IF

C.........  Report the name of the temperature variable
        IF( TVARNAME .NE. ' ' ) THEN

            L = LEN_TRIM( TVARNAME )
            MESG = 'NOTE: Using temperature variable name "' //
     &             TVARNAME( 1:L ) // '".'
            CALL M3MSG2( MESG )

        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94000   FORMAT( I2.2 )
 
94010   FORMAT( 10( A, :, I8, :, 1X ) )

C******************  INTERNAL SUBPROGRAMS   ******************************

        CONTAINS

C.............  This subroutine writes a warning message that the temperature
C               variable name is not consistent with the file
            SUBROUTINE TEMPERATURE_WARNING

C.............  Local variables
            INTEGER L, L2

C..........................................................................
            L  = LEN_TRIM( TVARNAME )
            L2 = LEN_TRIM( MNAME )
            MESG = 'WARNING: temperature variable "' // 
     &             TVARNAME( 1:L ) // '" is not in gridded' // 
     &             CRLF()// BLANK10// 'temperature file "' //
     &             MNAME( 1:L2 )// '".'
            CALL M3MSG2( MESG )

            END SUBROUTINE TEMPERATURE_WARNING

C----------------------------------------------------------------------------
C----------------------------------------------------------------------------

C.............  This subroutine compares the minimum/maximum temperatures
C               and sets an error flag
            SUBROUTINE COMPARE_TMPRS( CHECKTYP )

            INCLUDE 'FLTERR.EXT'    ! error filter statement function

            CHARACTER(*), INTENT( IN ) :: CHECKTYP

C..........................................................................            

            IF( FLTERR( MIN1, MINT_MIN ) ) THEN
                EFLAG = .TRUE.
            END IF

            IF( FLTERR( MAX2, MAXT_MAX ) ) THEN
                EFLAG = .TRUE.
            END IF

            IF( CHECKTYP .NE. 'NOMINMAX' ) THEN
                
        	IF( FLTERR( MIN2, MINT_MAX ) ) THEN
                    EFLAG = .TRUE.
        	END IF

        	IF( FLTERR( MAX1, MAXT_MIN ) ) THEN
                    EFLAG = .TRUE.
        	END IF

            END IF

            END SUBROUTINE COMPARE_TMPRS

        END SUBROUTINE OPENTMPIN

