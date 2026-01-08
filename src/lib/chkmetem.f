
        LOGICAL FUNCTION CHKMETEM( GNAM2D, MNAM2D, DNAM2D, 
     &                             GNAM3D, MNAM3D, DNAM3D )

C***********************************************************************
C  function body starts at line
C
C  DESCRIPTION:
C      This function compares the headers of the 6 or fewer meteorology 
C      files used in the emissions processing.  If the file names provided
C      as subroutine arguments are "NONE", then that file is not available
C      in the current call to the subroutine, and the associated part of the
C      checking for that file is skipped. The return value is TRUE for 
C      files that are consistent and FALSE otherwise
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Created 3/99 by M. Houyoux
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
C       Updated with USE M3UTILIO by Huy Tran UNC-IE on 2026-01
C***************************************************************************
        USE M3UTILIO

        USE MODGRDLIB       !! for dblerr() generic

        IMPLICIT NONE

C.........  INCLUDE FILES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
C        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
C        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
C        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C.........  EXTERNAL FUNCTIONS
C       CHARACTER(2)  CRLF
        CHARACTER(50) GETCFDSC

C        EXTERNAL      CRLF, GETCFDSC
        EXTERNAL     GETCFDSC

C.........  SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN) :: GNAM2D ! name of grid cross-point 2d file
        CHARACTER(*), INTENT (IN) :: MNAM2D ! name of met cross-point 2d file
        CHARACTER(*), INTENT (IN) :: DNAM2D ! name of grid dot-point 2d file
        CHARACTER(*), INTENT (IN) :: GNAM3D ! name of grid cross-point 3d file
        CHARACTER(*), INTENT (IN) :: MNAM3D ! name of met cross-point 3d file
        CHARACTER(*), INTENT (IN) :: DNAM3D ! name of met dot-point 3d file

C.........  Grid information
        INTEGER :: NCOLS = 0      ! number of columns in grid
        INTEGER :: NGRID = 0      ! number of cells in grid
        INTEGER :: NLAYS = 0      ! number of layers
        INTEGER :: NROWS = 0      ! number of rows in grid
        INTEGER :: GDTYP = -1     ! i/o api grid type code
        REAL(8) :: P_ALP = 0.D0   ! projection alpha
        REAL(8) :: P_BET = 0.D0   ! projection beta
        REAL(8) :: P_GAM = 0.D0   ! projection gamma
        REAL(8) :: XCENT = 0.D0   ! x-center of projection
        REAL(8) :: YCENT = 0.D0   ! y-center of projection
        REAL(8) :: XORIG = 0.D0   ! x-origin of grid
        REAL(8) :: YORIG = 0.D0   ! y-origin of grid
        REAL(8) :: XCELL = 0.D0   ! x-dim of cells
        REAL(8) :: YCELL = 0.D0   ! y-dim of cells
        CHARACTER(IOVLEN3) :: GRDNM = ' '  ! grid name

C.........  Vertical structure information
        INTEGER :: VGTYP  = -1     ! type of vertical coordinates
        REAL    :: VGTOP  = 0.0    ! model-top, for sigma coord types
        REAL    :: VGLVS( 0:MXLAYS3 ) ! vertical coordinate values

C.........  Other met file information for comparison
        CHARACTER(50) :: METSCEN  = ' '  ! Name of met scenario
        CHARACTER(50) :: CLOUDSHM = ' '  ! Name of cloud scheme

C.........  Other local variables
        INTEGER         I, J, J1, J2, K, L      ! indicies and counters

        LOGICAL      :: EFLAG     = .FALSE.  ! true: error in comparing files
        LOGICAL      :: DOT_BASIS = .FALSE.  ! true: the comparison file is dot
        LOGICAL      :: THREE_D   = .FALSE.  ! true: one of the inputs is 3d
        LOGICAL      :: VFLAG     = .FALSE.  ! true: error in vertical

        CHARACTER(IOVLEN3) FILNAM 
        CHARACTER(300)     MESG 

        CHARACTER(16) :: PROGNAME = 'CHKMETEM' ! program name

C***********************************************************************
C   begin body of function CHKMETEM

C.........  Find the file to use as a base-line file for checking the others...
C.........  Determine if any of the 3-d files are used. If so, find one to
C           use as the base-line
        IF( GNAM3D .NE. 'NONE' .OR. 
     &      MNAM3D .NE. 'NONE'      ) THEN
c     &      MNAM3D .NE. 'NONE' .OR. 
c     &      DNAM3D .NE. 'NONE'      ) THEN

            THREE_D = .TRUE.

            IF( MNAM3D .NE. 'NONE'  ) THEN
                FILNAM = MNAM3D

c            ELSE IF( DNAM3D .NE. 'NONE' ) THEN
c                DOT_BASIS = .TRUE.
c                FILNAM = DNAM3D

            ELSE IF( GNAM3D .NE. 'NONE' ) THEN
                FILNAM = GNAM3D

            END IF

C.........  If only 2-d files, then find one to use as the base-line
        ELSE IF( MNAM2D .NE. 'NONE' ) THEN
            FILNAM = MNAM2D
            
c        ELSE IF( DNAM2D .NE. 'NONE' ) THEN
c            DOT_BASIS = .TRUE.
c            FILNAM = DNAM2D

        ELSE IF( GNAM2D .NE. 'NONE' ) THEN
            FILNAM = GNAM2D

C.........  Internal error if none of the files are valid
        ELSE
            MESG = 'INTERNAL ERROR: No opened meteorology files ' //
     &             'provided in ' // PROGNAME // ' call!'
            CALL M3MSG2( MESG )
            CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

        END IF

C.........  Set the reference local grid and other meteorology file settings
        IF( .NOT. DESC3( FILNAM ) ) THEN

            MESG = 'Could not get description of file "' //
     &             FILNAM( 1:LEN_TRIM( FILNAM ) ) // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        ELSE
            GRDNM = GDNAM3D
            GDTYP = GDTYP3D
            P_ALP = P_ALP3D
            P_BET = P_BET3D
            P_GAM = P_GAM3D
            XCENT = XCENT3D
            YCENT = YCENT3D
            XORIG = XORIG3D
            YORIG = YORIG3D
            XCELL = XCELL3D
            YCELL = YCELL3D
            NCOLS = NCOLS3D
            NROWS = NROWS3D

C.............  Adjust rows and columns if file being used as basis is dot-point
            IF( DOT_BASIS ) THEN
                NCOLS = NCOLS - 1
                NROWS = NROWS - 1
c                XORIG = XORIG + 0.5 * XCELL
c                YORIG = YORIG + 0.5 * YCELL
            END IF

C.............  Set 3d file settigs
            IF( THREE_D ) THEN
                NLAYS = NLAYS3D
                VGTYP = VGTYP3D
                VGTOP = VGTOP3D

                J = LBOUND( VGLVS3D,1 )
                DO K = 0, NLAYS
                    VGLVS( K ) = VGLVS3D( J )
                    J = J + 1
                END DO
            END IF

C.............  Try to find scenario name a cloud scheme file description info
            METSCEN  = GETCFDSC( FDESC3D, '/MET SCENARIO/', .FALSE. )
            CLOUDSHM = GETCFDSC( FDESC3D, '/CLOUD SCHEME/', .FALSE. )

C.............  Write message about which file used for initializing...
            L = LEN_TRIM( FILNAM )
            MESG = 'NOTE: file "' // FILNAM( 1:L ) // '" was used ' //
     &             'to initialize meteorology file header checks.'
            CALL M3MSG2( MESG )

        END IF

C.........  Checking for grid cross-point 2d file
        IF ( GNAM2D .NE. 'NONE' ) CALL CHECK_MET_INFO( GNAM2D, .FALSE. )

C.........  Checking for meteorology cross-point 2d file
        IF ( MNAM2D .NE. 'NONE' ) CALL CHECK_MET_INFO( MNAM2D, .FALSE. ) 

C.........  Checking for grid dot-point 2d file
        IF ( DNAM2D .NE. 'NONE' ) CALL CHECK_MET_INFO( DNAM2D, .TRUE. )

C.........  Checking for grid cross-point 3d file
        IF ( GNAM3D .NE. 'NONE' ) CALL CHECK_MET_INFO( GNAM3D, .FALSE. )

C.........  Checking for meteorology cross-point 3d file
        IF ( MNAM3D .NE. 'NONE' ) CALL CHECK_MET_INFO( MNAM3D, .FALSE. )

C.........  Checking for meteorology dot-point 3d file
        IF ( DNAM3D .NE. 'NONE' ) CALL CHECK_MET_INFO( DNAM3D, .TRUE. )

        IF( EFLAG ) THEN

            CHKMETEM = .FALSE.

        ELSE

            CHKMETEM = .TRUE.

        END IF

        RETURN

 
C******************  INTERNAL SUBPROGRAMS  *****************************

        CONTAINS
 
C.............  This internal subprogram tries to retrieve the I/O API header
C               and aborts if it was not successful
            SUBROUTINE CHECK_MET_INFO( FILNAM, DOT_STATUS )

C.............  Subprogram arguments
            CHARACTER(*) FILNAM      ! Name of file to check
            LOGICAL      DOT_STATUS  ! Status of file as dot-point or cross-pt

C.............  Local variables
            INTEGER       I, J1, J2, K, L, L2

            CHARACTER(50)  CVAL        ! Current file's METSCEN or CLOUDSHM
            CHARACTER(300) MESG        ! Output message

C----------------------------------------------------------------------

C.............  Get file description from I/O API header
            IF ( .NOT. DESC3( FILNAM ) ) THEN

                MESG = 'Could not get description of file "' //
     &                 FILNAM( 1:LEN_TRIM( FILNAM ) ) // '"'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            END IF

C.............  Adjust columns and rows for comparison if file is a dot-point
            IF( DOT_STATUS ) THEN
                NCOLS3D = NCOLS3D - 1
                NROWS3D = NROWS3D - 1
                XORIG3D = XORIG
                YORIG3D = YORIG
c                XORIG3D = XORIG3D + 0.5D0 * XCELL3D
c                YORIG3D = YORIG3D + 0.5D0 * YCELL3D
            END IF

            L2 = LEN_TRIM( FILNAM )

C.............  Check horizontal parameters from header
            IF ( NCOLS .NE. NCOLS3D .OR.
     &           NROWS .NE. NROWS3D .OR.
     &           DBLERR( XCELL, XCELL3D ) .OR.
     &           DBLERR( YCELL, YCELL3D ) .OR.
     &           DBLERR( XORIG, XORIG3D ) .OR.
     &           DBLERR( YORIG, YORIG3D ) .OR.
     &           DBLERR( XCENT, XCENT3D ) .OR.
     &           DBLERR( YCENT, YCENT3D ) .OR.
     &           DBLERR( P_ALP, P_ALP3D ) .OR.
     &           DBLERR( P_BET, P_BET3D ) .OR.
     &           DBLERR( P_GAM, P_GAM3D )      ) THEN

                EFLAG = .TRUE.
                MESG = 'ERROR: Horizontal grid parameters in file "' // 
     &                 FILNAM( 1:L2 ) // '"'// CRLF() // BLANK10//
     &                 'are inconsistent with initialized values.'
                CALL M3MSG2( MESG )

            END IF

C.............  Check vertical parameters from header
            IF( NLAYS3D .GT. 1 ) THEN

            
                IF ( NLAYS .NE. NLAYS3D .OR.
     &               VGTYP .NE. VGTYP3D .OR.
     &               VGTOP .NE. VGTOP3D      ) THEN

                    EFLAG = .TRUE.
                    VFLAG = .TRUE.

                ELSE

                    J = LBOUND( VGLVS3D,1 )
                    DO K = 0, NLAYS

                        IF( FLTERR( VGLVS( K ), VGLVS3D( J ) ) ) THEN
                            EFLAG = .TRUE.
                            VFLAG = .TRUE.
                        END IF
                        J = J + 1

                    END DO

                END IF

                IF( VFLAG ) THEN
                    MESG = 'ERROR: Vertical grid parameters in file "'// 
     &                     FILNAM( 1:L2 ) // '"'//CRLF() // BLANK10//
     &                     'are inconsistent with initialized values.'
                    CALL M3MSG2( MESG )
                END IF

            END IF

C.............  Look for met scenario and cloud scheme indicators in I/O API
C               file header. If they are not there, do nothing. If they are 
C               there, compare to the original settings.
            CVAL = GETCFDSC( FDESC3D, '/MET SCENARIO/', .FALSE. )
            IF( CVAL .NE. ' ' .AND. CVAL .NE. METSCEN ) THEN
                EFLAG = .TRUE.
                MESG = 'ERROR: Meteorology scenario in file "' // 
     &                 FILNAM( 1:L2 ) // '" is inconsistent '//
     &                 'with initialized value.'
                CALL M3MSG2( MESG )
            END IF

            CVAL = GETCFDSC( FDESC3D, '/CLOUD SCHEME/', .FALSE. )
            IF( CVAL .NE. ' ' .AND. CVAL .NE. CLOUDSHM ) THEN
                EFLAG = .TRUE.
                MESG = 'ERROR: Cloud scheme in file "' // 
     &                 FILNAM( 1:L2 ) // '" is inconsistent '//
     &                 CRLF() // BLANK10// 'with initialized value.'
                CALL M3MSG2( MESG )
            END IF

            RETURN

            END SUBROUTINE CHECK_MET_INFO

        END FUNCTION CHKMETEM
