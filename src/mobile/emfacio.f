
        SUBROUTINE EMFACIO( MODELNAM, NNAME, DNAME, JYEAR, EDEV, IDEV, 
     &                      RDEV, SDEV, NDREUSE, DIREUSE )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      This subroutine opens the existing or new diurnal and non-diurnal
C      emission factor files. It uses environment variables to control
C      whether existing files are to be updated or new files are to
C      be created.
C 
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Created 7/99 by M. Houyoux
C
C************************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2000, MCNC--North Carolina Supercomputing Center
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
C.........  This module is for mobile-specific data
        USE MODMOBIL

C...........   This module contains emission factor tables and related
        USE MODEMFAC

C...........   This module contains the information about the source category
        USE MODINFO

C...........   This module is the derived meteorology data for emission factors
        USE MODMET

        IMPLICIT NONE

C...........   INCLUDES:
        INCLUDE 'EMCNST3.EXT'    !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'     !  I/O API parameters
        INCLUDE 'IODECL3.EXT'    !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'     !  I/O API file description data structures

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER*2            CRLF
        INTEGER                ENVINT
        LOGICAL                ENVYN
        CHARACTER(LEN=IODLEN3) GETCFDSC
        INTEGER                PROMPTFFILE
        CHARACTER*16           PROMPTMFILE
        CHARACTER*16           VERCHAR

        EXTERNAL     CRLF, GETCFDSC, PROMPTFFILE, PROMPTMFILE, VERCHAR
     &               

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT    (IN) :: MODELNAM ! name of EF model
        CHARACTER(*), INTENT   (OUT) :: NNAME    ! name of non-diurnal EF file
        CHARACTER(*), INTENT   (OUT) :: DNAME    ! name of     diurnal EF file
        INTEGER     , INTENT   (OUT) :: JYEAR    ! year of EFs
        INTEGER     , INTENT   (OUT) :: EDEV     ! unit no. of EF-code/temp file
        INTEGER     , INTENT   (OUT) :: IDEV     ! unit no. of EF x-ref file
        INTEGER     , INTENT   (OUT) :: RDEV     ! unit no. of EF model inputs
        INTEGER     , INTENT   (OUT) :: SDEV     ! unit no. of update file
        LOGICAL     , INTENT   (OUT) :: NDREUSE  ! true: reusing non-diur EFs
        LOGICAL     , INTENT   (OUT) :: DIREUSE  ! true: reusing diur EFs

C...........   LOCAL PARAMETERS
        CHARACTER*50, PARAMETER :: SCCSW = ! SCCS string with version no. at end
     &               '@(#)$Id$'

C...........   Local file names and unit numbers
        INTEGER        CDEV        ! mobile codes for vehicle types and rtypes
 
        CHARACTER*16   DNAME_IN    ! input diurnal EF file name
        CHARACTER*16   NNAME_IN    ! input non-diurnal EF file name

C...........   Other local variables
        INTEGER        IOS     ! E.V. return status

        LOGICAL        TFLAG   ! true: not reusing non-diurnal EF file
        LOGICAL        UFLAG   ! true: use EF forced update file

        CHARACTER*300  MESG    !  message buffer

        CHARACTER(LEN=IODLEN3)  IFDESC2, IFDESC3 ! fields 2 & 3 from INVEN FDESC

        CHARACTER*16 :: PROGNAME = 'EMFACIO' ! program name

C***********************************************************************
C   begin body of subroutine EMFACIO

C.........  Get environment variables that control the opening of files...

C.........  Set whether non-diurnal emission factors will be reused
        MESG = 'Reuse setting for nondiurnal emission factors'
        NDREUSE = ENVYN ( 'REUSE_NONDIURNAL', MESG, .FALSE., IOS )

C.........  Set whether diurnal emission factors will be reused
        MESG = 'Reuse setting for diurnal emission factors'
        DIREUSE = ENVYN ( 'REUSE_DIURNAL', MESG, .FALSE., IOS )

C.........  Set whether a targeted emission factor update file will be used
        MESG = 'Indicator for using emission factor update file'
        UFLAG = ENVYN ( 'EF_FORCE_UPDATE', MESG, .FALSE., IOS )

C.........  Get a year for computing the emission factors
        MESG = 'Emission factors year'
        JYEAR = ENVINT( 'EF_YEAR', MESG, 1988, IOS )

C.........  Prompt for inputs files...

C.........  Prompt for diurnal master list EFs input file for modeling year
        IF( NDREUSE ) 
     &      NNAME_IN = PROMPTMFILE(
     &              'Enter logical name for input NON-DIURNAL EF file',
     &              FSRDWR3, CRL // 'EFSND_IN', PROGNAME )

C.........  Prompt for diurnal master list EFs input file for modeling year
        IF( DIREUSE ) 
     &      DNAME_IN = PROMPTMFILE(
     &              'Enter logical name for input DIURNAL EF file',
     &              FSRDWR3, CRL // 'EFSD_IN', PROGNAME )

C.........  Prompt for emission factors file that force updating of certain EFs
        IF( UFLAG ) 
     &      SDEV = PROMPTFFILE(
     &                'Enter logical name for TARGETED EF UPDATE file',
     &                .TRUE., .TRUE., CRL // 'EFUPD', PROGNAME )

C.........  Prompt for EF-based file of min/max temperature combinations
        EDEV  = PROMPTFFILE(
     &          'Enter logical name for EF-REF/TEMPERATURE file',
     &          .TRUE., .TRUE., CRL // 'EFTEMP', PROGNAME )

C.........  Prompt for scenario specific emission factor cross-reference file
        IDEV = PROMPTFFILE(
     &         'Enter logical name for EMISSION FACTOR X-REF file',
     &         .TRUE., .TRUE., CRL // 'PLIST', PROGNAME )
        
C.........  Prompt foremission factor data file
        RDEV = PROMPTFFILE( 
     &         'Enter logical name for EMISSION FACTOR DATA file',
     &         .TRUE., .TRUE., CRL // 'PREF', PROGNAME )

C.........  Get file name for getting the number of vehicle types
        MESG = 'Enter logical name for MOBILE CODES file'
        CDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'MCODES', PROGNAME )

C.........  Get, among other things, the number of vehicle types from MCODES
        CALL RDMVINFO( CDEV )

C............................................................................
C           Non-diurnal file section
C............................................................................

C.........  Get description of non-diurnal emission factors input file
C.........  Store various things from the header
        IF( NDREUSE ) THEN

            IF( .NOT. DESC3( NNAME_IN ) ) THEN
                MESG = 'Could not get description of ' // NNAME_IN
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

        END IF

C.........  Get temperature information from non-diurnal emission factors
C           file or from environment variables, and store for comparison
C           with non-diurnal file
C.........  Build the allowable non-diurnal temperature array
        TFLAG = ( .NOT. NDREUSE )
        CALL TMPRINFO( TFLAG, 'NOMINMAX' )

C.........  If starting with existing non-diurnal emission factors file,
C           set parameters that would otherwise get from prompt
        IF( NDREUSE ) THEN

C.............  Set output file name to same as input
            NNAME = NNAME_IN

C.............  Set modeling year (override the environment variable setting)
            JYEAR = SDATE3D / 1000

C.............  Set the number of non diurnal emission factors
            NNDI = NVARS3D

C.........  If not starting with an existing non-diurnal emission factors file,
C           set up the header, prompt for the file, and open it
        ELSE

C.............  Initialize emissions header to missing
            CALL HDRMISS3

C.............  Get emission factor names, units, and descriptions for the  
C               emission factor model of interest
            CALL EFSETUP( MODELNAM, 'NONDIURNAL', MXVARS3, NNDI, 
     &                    VNAME3D, UNITS3D, VDESC3D )

C.............  Set I/O API variable information
            VTYPE3D( 1:NNDI ) = M3REAL           

C.............  Set up non-default header options for non-diurnal EF file 
            XORIG3D = -DBLE( MINT_MIN ) ! Negative starting temp (- for PAVE)
            YORIG3D = 26.D0           ! Set to prevent PAVE crash
            XCELL3D = DBLE( TMMINVL ) ! Temperature increment
            YCELL3D = 1.D0            ! Set to prevent PAVE crash
            SDATE3D = JYEAR * 1000 + 1
            STIME3D = 1
            TSTEP3D = 1
            NVARS3D = NNDI
            NCOLS3D = NTMPR           ! note: Rename this!
            NROWS3D = NVTYPE 
            GDTYP3D = LATGRD3         ! Set to prevent PAVE crash

            FDESC3D( 1 ) = CATEGORY( 1:LEN_TRIM( CATEGORY ) ) //
     &                     ' non-diurnal emission factors'
            FDESC3D( 2 ) =               '/FROM/ '    // PROGNAME
            FDESC3D( 3 ) =               '/VERSION/ '// VERCHAR( SCCSW )
            WRITE( FDESC3D( 4 ), 94030 ) '/MINT_MIN/', MINT_MIN
            WRITE( FDESC3D( 5 ), 94030 ) '/MAXT_MAX/', MAXT_MAX
            WRITE( FDESC3D( 6 ), 94030 ) '/T_INTERVAL/', TMMINVL
            FDESC3D( 7 ) =               '/T_UNITS/ "K"'
            FDESC3D( 8 ) =               '/FROM MODEL/' // MODELNAM

            MESG = 'Enter logical name for output NON-DIURNAL ' //
     &             'EMISSION FACTORS file'
            NNAME = CRL // 'EFSND'
            NNAME = PROMPTMFILE( MESG, FSNEW3, NNAME, PROGNAME ) 

        END IF  ! File exists already or not

C............................................................................
C           Diurnal file section
C............................................................................

C.........  Get description of diurnal emission factors input file
C.........  Store various things from the header
        IF( DIREUSE ) THEN

            IF( .NOT. DESC3( DNAME_IN ) ) THEN
                MESG = 'Could not get description of ' // DNAME_IN
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

        END IF

C.........  Get temperature information from diurnal emission factors
C           file or from environment variables
C.........  Build the allowable min/max diurnal temperature arrays
        TFLAG = ( .NOT. DIREUSE )
        CALL TMPRINFO( TFLAG, 'MINMAX' )

C.........  If starting with existing non-diurnal emission factors file,
C           set parameters that would otherwise get from prompt
        IF ( DIREUSE ) THEN 

C.............  Set output file name to same as input
            DNAME  = DNAME_IN

C.............  Make sure that the non-diurnal and diurnal files are consistent
            CALL CHKEFHDR( NNAME, DNAME )

C.............  Set the number of diurnal emission factors
            NDIU = NVARS3D

C.........  If not starting with an existing diurnal emission factors file,
C           set up the header, prompt for the file, and open it
        ELSE

C.............  Initialize emissions header to missing
            CALL HDRMISS3

C.............  Get emission factor names, units, and descriptions for the  
C               emission factor model of interest
            CALL EFSETUP( MODELNAM, 'DIURNAL', MXVARS3, NDIU, 
     &                   VNAME3D, UNITS3D, VDESC3D )

C.............  Set I/O API variable type
            VTYPE3D( 1:NDIU ) = M3REAL           
 
C.............  Set up diurnal master List file description

            FTYPE3D = GRDDED3
            XORIG3D = -DBLE( MINT_MIN ) ! Negative starting temp (- for PAVE)
            YORIG3D = 26.D0             ! Set for PAVE
            XCELL3D = DBLE( TMMINVL )     ! Temperature increment
            YCELL3D = 1.0               ! Set for PAVE
            SDATE3D = JYEAR * 1000 + 1
            STIME3D = 1
            TSTEP3D = 1
            NVARS3D = NDIU
            NCOLS3D = NVLDTMM
            NROWS3D = NVTYPE 
            GDTYP3D = LATGRD3           ! Set for PAVE

            FDESC3D( 1 ) = CATEGORY( 1:LEN_TRIM( CATEGORY ) ) //
     &                     ' diurnal emission factors'
            FDESC3D( 2 ) =               '/FROM/ '    // PROGNAME
            FDESC3D( 3 ) =               '/VERSION/ '// VERCHAR( SCCSW )
            WRITE( FDESC3D( 4 ), 94030 ) '/MINT_MIN/', MINT_MIN
            WRITE( FDESC3D( 5 ), 94030 ) '/MINT_MAX/', MINT_MAX
            WRITE( FDESC3D( 6 ), 94030 ) '/MAXT_MIN/', MAXT_MIN
            WRITE( FDESC3D( 7 ), 94030 ) '/MAXT_MAX/', MAXT_MAX
            WRITE( FDESC3D( 8 ), 94030 ) '/T_INTERVAL/', TMMINVL
            WRITE( FDESC3D( 9 ), 94030 ) '/T_MAXINTVL/', TMXINVL
            FDESC3D( 10 ) =              '/T_UNITS/ "K"'
            FDESC3D( 11 ) =              '/FROM MODEL/' // MODELNAM

            MESG = 'Enter logical name for output DIURNAL ' //
     &             'EMISSION FACTORS file'
            DNAME = CRL // 'EFSD'
            DNAME = PROMPTMFILE( MESG, FSNEW3, DNAME, PROGNAME )

        ENDIF

C.........  Write message about which year emission factors will be for
        WRITE( MESG,94010 ) 
     &         'NOTE: Emission factors are for year', JYEAR
        CALL M3MSG2( MESG )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94030   FORMAT( A, F15.9, 1X, A )

        END SUBROUTINE EMFACIO
