
        SUBROUTINE OPENTMPIN( ENAME, ANAME, DNAME, HNAME, FNAME,
     &                        NNAME, MNAME, GNAME, WNAME, SDEV, XDEV,
     &                        RDEV, UDEV, FDEV, TDEV )

c note: delete UDEV from here, inputs definition, and temporal.f

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      This subroutine opens the file or files for output from the tmppoint
C      program
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

C.........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constat parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS and their descriptions:
        INTEGER                PROMPTFFILE
        CHARACTER(LEN=NAMLEN3) PROMPTMFILE

        EXTERNAL        PROMPTFFILE, PROMPTMFILE

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT(IN OUT) :: ENAME !  name for I/O API inven input
        CHARACTER(*), INTENT(IN OUT) :: ANAME !  name for ASCII inven input 
        CHARACTER(*), INTENT   (OUT) :: DNAME !  day-spec file
        CHARACTER(*), INTENT   (OUT) :: HNAME !  hour-spec file
        CHARACTER(*), INTENT   (OUT) :: FNAME !  non-diurnal EF file
        CHARACTER(*), INTENT   (OUT) :: NNAME !  diurnal EF file
        CHARACTER(*), INTENT   (OUT) :: MNAME !  surface temperature file
        CHARACTER(*), INTENT   (OUT) :: GNAME !  ungridding matrix
        CHARACTER(*), INTENT   (OUT) :: WNAME !  ungridded min/max temperatures
        INTEGER     , INTENT   (OUT) :: SDEV  !  unit no.: ASCII inven file
        INTEGER     , INTENT   (OUT) :: XDEV  !  unit no.: x-ref file
        INTEGER     , INTENT   (OUT) :: RDEV  !  unit no.: tmprl profile file
        INTEGER     , INTENT   (OUT) :: UDEV  !  unit no.: optional elev srcs
        INTEGER     , INTENT   (OUT) :: FDEV  !  unit no.: emis factors x-ref
        INTEGER     , INTENT   (OUT) :: TDEV  !  unit no.: speciation list file
 
C...........   LOCAL PARAMETERS
        CHARACTER*50  SCCSW          ! SCCS string with version number at end

        PARAMETER   ( SCCSW   = '@(#)$Id$'
     &              )

C...........   Other local variables

        INTEGER         IOS         ! status from environment variables

        LOGICAL         DFLAG       ! day-specific  file available
        LOGICAL         HFLAG       ! hour-specific file available

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

        XDEV = PROMPTFFILE( 
     &           'Enter logical name for TEMPORAL XREF file',
     &           .TRUE., .TRUE., CRL // 'TREF', PROGNAME )

        RDEV = PROMPTFFILE( 
     &           'Enter logical name for TEMPORAL PROFILES file',
     &           .TRUE., .TRUE., CRL // 'TPRO', PROGNAME )

        IF( CATEGORY .EQ. 'MOBILE' ) THEN

            FNAME = PROMPTMFILE( 
     &              'Enter logical name for NON-DIURNAL EMISSIONS ' //
     &              'FACTOR file', FSREAD3, CRL // 'EFSND', PROGNAME )
    
            NNAME = PROMPTMFILE( 
     &              'Enter logical name for DIURNAL EMISSIONS ' //
     &              'FACTOR file', FSREAD3, CRL // 'EFSD', PROGNAME )

            MNAME = PROMPTMFILE( 
     &              'Enter logical name for SURFACE TEMPERATURE file',
     &              FSREAD3, 'MET_CRO_3D', PROGNAME )

            GNAME = PROMPTMFILE( 
     &              'Enter logical name for UNGRIDDING MATRIX file',
     &              FSREAD3, CRL // 'UMAT', PROGNAME )

            WNAME = PROMPTMFILE( 
     &              'Enter logical name for UNGRIDDED MIN/MAX ' //
     &              'TEMPERATURE file', FSREAD3, 'MINMAXT', PROGNAME )

            FDEV = PROMPTFFILE( 
     &             'Enter logical name for EMISSION FACTORS XREF ' // 
     &             'file', .TRUE., .TRUE., CRL // 'PLIST', PROGNAME )

            TDEV = PROMPTFFILE( 
     &             'Enter logical name for SPECIATION LIST file',
     &             .TRUE., .TRUE., 'SPCS', PROGNAME )

        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94000   FORMAT( I2.2 )
 
94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE OPENTMPIN

