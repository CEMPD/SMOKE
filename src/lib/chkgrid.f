
        SUBROUTINE CHKGRID( CATDESC, FTYPE, EFLAG )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      This subroutine updates the grid information and compares to 
C      the existing information, if it has been previously set.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C
C**************************************************************************
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

C.........  MODULES for public variables
C.........  This module contains the major data structure and control flags
        USE MODMERGE

        IMPLICIT NONE

C.........  INCLUDES:
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file desc. data structures
        INCLUDE 'FLTERR.EXT'    !  error filter statement function

C.........  EXTERNAL FUNCTIONS and their descriptions:
        INTEGER      GETIFDSC  
        EXTERNAL     GETIFDSC

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT(IN) :: CATDESC  ! category descriptions
        CHARACTER(*), INTENT(IN) :: FTYPE    ! GMAT|GRID file type of interest
        LOGICAL     , INTENT(OUT):: EFLAG    ! true: comparison failed

C...........   Local variables
        INTEGER       L       ! length of file description
        INTEGER       NC      ! tmp number of columns
        INTEGER       NR      ! tmp number of rows

        LOGICAL, SAVE :: GFLAG = .FALSE.  ! true: grid settings have been init

        CHARACTER*20    FILDESC  ! description of input file
        CHARACTER*300   MESG     ! message buffer

        CHARACTER*16 :: PROGNAME = 'CHKGRID' ! program name

C***********************************************************************
C   begin body of function CHKGRID

C.............  Set tmp rows, columns, and total cells depending on file type
        IF( FTYPE .EQ. 'GMAT' ) THEN
            NC = GETIFDSC( FDESC3D, '/NCOLS3D/', .TRUE. )
            NR = GETIFDSC( FDESC3D, '/NROWS3D/', .TRUE. )
            FILDESC = 'gridding matrix'

        ELSEIF( FTYPE .EQ. 'GRID' ) THEN
            NC = NCOLS3D
            NR = NROWS3D
            FILDESC = 'gridded file'

        ELSE
            MESG = 'INTERNAL ERROR: File type "' // FTYPE // 
     &              '" not known in call to ' // PROGNAME
            CALL M3MSG2( MESG )
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        ENDIF

        L = LEN_TRIM( FILDESC )

C.............  If grid information has already been initialized, then compare
C               existing to this file.
        IF( GFLAG ) THEN

            IF ( NCOLS .NE. NC      .OR.
     &           NROWS .NE. NR      .OR.
     &           GDTYP .NE. GDTYP3D .OR.
     &           FLTERR( XCELL, SNGL( XCELL3D ) ) .OR.
     &           FLTERR( YCELL, SNGL( YCELL3D ) ) .OR.
     &           FLTERR( XORIG, SNGL( XORIG3D ) ) .OR.
     &           FLTERR( YORIG, SNGL( YORIG3D ) ) .OR.
     &           FLTERR( XCENT, SNGL( XCENT3D ) ) .OR.
     &           FLTERR( YCENT, SNGL( YCENT3D ) ) .OR.
     &           FLTERR( P_ALP, SNGL( P_ALP3D ) ) .OR.
     &           FLTERR( P_BET, SNGL( P_BET3D ) ) .OR.
     &           FLTERR( P_GAM, SNGL( P_GAM3D ) )      ) THEN

                EFLAG = .TRUE.
                MESG = 'Grid parameters in ' // CATDESC // ' ' //
     &                 FILDESC( 1:L ) // ' are not consistent ' //
     &                 'with initialized values.'
                CALL M3MSG2( MESG )

            END IF

C.........  Initialize grid information
        ELSE

            GFLAG = .TRUE.
            GRDNM = GDNAM3D
            GDTYP = GDTYP3D
            P_ALP = SNGL( P_ALP3D )
            P_BET = SNGL( P_BET3D )
            P_GAM = SNGL( P_GAM3D )
            XCENT = SNGL( XCENT3D )
            YCENT = SNGL( YCENT3D )
            XORIG = SNGL( XORIG3D )
            YORIG = SNGL( YORIG3D )
            XCELL = SNGL( XCELL3D )
            YCELL = SNGL( YCELL3D )
            NCOLS = NC
            NROWS = NR
            NGRID = NCOLS * NROWS

            MESG = 'NOTE: Grid settings initialized using ' // 
     &             CATDESC // ' ' // FILDESC( 1:L ) // '.'

            CALL M3MSG2( MESG )

        ENDIF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END

