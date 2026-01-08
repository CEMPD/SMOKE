
        SUBROUTINE CHKGRID( DATDESC, FTYPE, CHKLEVEL, EFLAG )

C***********************************************************************
C  subroutine body starts at line 91
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
C Project Title: EDSS Tools Library
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

C.........  MODULES for public variables
C.........  This module contains the global variables for the 3-d grid
        USE M3UTILIO

        USE MODGRID, ONLY: NCOLS, NROWS, XORIG, YORIG, XOFF, YOFF,
     &                     GDTYP, XCELL, YCELL, XCENT, YCENT,
     &                     P_ALP, P_BET, P_GAM, OFFLAG, GRDNM,
     &                     XOFF_A, YOFF_A, XDIFF, YDIFF, NGRID

        IMPLICIT NONE

C.........  INCLUDES:
        INCLUDE 'IOCNST3.EXT'   !  emissions constant parameters
C        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
C        INCLUDE 'FDESC3.EXT'    !  I/O API file desc. data structures

C.........  EXTERNAL FUNCTIONS and their descriptions:
C       CHARACTER(2) CRLF
        INTEGER      GETIFDSC  

C        EXTERNAL     CRLF, GETIFDSC
        EXTERNAL     GETIFDSC

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT(IN   ) :: DATDESC  ! data descriptions
        CHARACTER(*), INTENT(IN   ) :: FTYPE    ! GMAT|GRID file type of interest
        INTEGER     , INTENT(IN   ) :: CHKLEVEL ! strigency of check
        LOGICAL     , INTENT(INOUT):: EFLAG    ! true: comparison failed

C...........   Local parameters
        INTEGER, PARAMETER :: CHK_ALL = 0     ! check all of the grid settings
        INTEGER, PARAMETER :: CHK_SUBGRID = 1 ! check all but allow subgrids
        INTEGER, PARAMETER :: CHK_TMPSUBG = 2 ! check all, allow temporary subgrids

C...........   Local variables
        INTEGER       L       ! length of file description
        INTEGER       NC      ! tmp number of columns
        INTEGER       NR      ! tmp number of rows
        INTEGER       XO      ! tmp x-offset  
        INTEGER       YO      ! tmp y-offset  

        REAL(8)       CHK_X   ! tmp val for checking subgrid even with grid
        REAL(8)       CHK_Y   ! tmp val for checking subgrid even with grid
        REAL(8)       X0      ! tmp x origin
        REAL(8)       Y0      ! tmp y origin

        LOGICAL, SAVE :: GFLAG  = .FALSE. ! true: grid settings have been init
        LOGICAL       :: SFLAG  = .FALSE. ! true: local error

        CHARACTER(25)   FILDESC  ! description of input file
        CHARACTER(300)  MESG     ! message buffer

        CHARACTER(16) :: PROGNAME = 'CHKGRID' ! program name

C***********************************************************************
C   begin body of function CHKGRID

C.............  Initialize local error flag
        SFLAG = .FALSE.

C.............  Set tmp rows, columns, and total cells depending on file type
        X0 = XORIG3D
        Y0 = YORIG3D

        IF( FTYPE .EQ. 'GMAT' ) THEN
            NC = GETIFDSC( FDESC3D, '/NCOLS3D/', .TRUE. )
            NR = GETIFDSC( FDESC3D, '/NROWS3D/', .TRUE. )
            FILDESC = 'gridding matrix'

        ELSEIF( FTYPE .EQ. 'GROUPS' ) THEN
            NC = GETIFDSC( FDESC3D, '/NCOLS3D/', .TRUE. )
            NR = GETIFDSC( FDESC3D, '/NROWS3D/', .TRUE. )
            FILDESC = 'stack groups file'

        ELSEIF( FTYPE .EQ. 'GRID' ) THEN
            NC = NCOLS3D
            NR = NROWS3D
            FILDESC = 'gridded file'

        ELSEIF( FTYPE .EQ. 'GRIDDESC' ) THEN
            NC = NCOLS3D
            NR = NROWS3D
            FILDESC = 'grid description file'

        ELSEIF( FTYPE .EQ. 'SURROGATES' ) THEN
            NC = NCOLS3D
            NR = NROWS3D
            FILDESC = 'surrogates file'

        ELSEIF( FTYPE .EQ. 'LANDUSE' ) THEN
            NC = NCOLS3D
            NR = NROWS3D
            FILDESC = 'landuse file'

        ELSEIF( FTYPE .EQ. 'DOT' ) THEN
            NC = NCOLS3D - 1
            NR = NROWS3D - 1
            X0 = XORIG
            Y0 = YORIG
c            X0 = XORIG3D + 0.5 * XCELL3D
c            Y0 = YORIG3D + 0.5 * YCELL3D
            FILDESC = 'dot-gridded file'

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

C.............  Check settings that must be consistent for exact grid match
            IF( CHKLEVEL .EQ. CHK_ALL ) THEN

                IF( NCOLS .NE. NC       .OR.
     &              NROWS .NE. NR       .OR.
     &              DABS( XORIG-X0 ) > 0.001D0  .OR.
     &              DABS( YORIG-Y0 ) > 0.001D0  )  THEN
                    SFLAG = .TRUE.
                    MESG = 'ERROR: Columns, rows, x-origin, or ' //
     &                     'y-origin for ' // DATDESC // ' in ' //
     &                     CRLF() // BLANK10 // FILDESC( 1:L ) // 
     &                     ' are inconsistent with initialized values.'
                    CALL M3MSG2( MESG ) 

                END IF

                XOFF = 0
                YOFF = 0

            END IF

C.............  Check settings that must be consistent for grids and subgrids
            IF( CHKLEVEL .LE. CHK_TMPSUBG ) THEN

                IF( GDTYP .NE. GDTYP3D .OR.
     &              DABS( XCELL-XCELL3D  ) > 0.001D0  .OR.
     &              DABS( YCELL-YCELL3D  ) > 0.001D0  .OR.
     &              DABS( XCENT-XCENT3D  ) > 0.001D0  .OR.
     &              DABS( YCENT-YCENT3D  ) > 0.001D0  .OR.
     &              DABS( P_ALP-P_ALP3D  ) > 0.001D0  .OR.
     &              DABS( P_BET-P_BET3D  ) > 0.001D0  .OR.
     &              DABS( P_GAM-P_GAM3D  ) > 0.001D0  ) THEN
                    SFLAG = .TRUE.
                    MESG = 'ERROR: Grid type, cell sizes, or ' //
     &                     'grid projection for ' // DATDESC // ' in '//
     &                     CRLF() // BLANK10 // FILDESC( 1:L ) // 
     &                     ' are inconsistent with initialized values.'
                    CALL M3MSG2( MESG ) 

                END IF

C.................  Ensure that origins are compatible with each other by
C                   making sure they line up based on the cell sizes
                CHK_X  = ( X0 - XORIG ) / XCELL
                CHK_X  = CHK_X - INT( CHK_X )
                CHK_Y  = ( Y0 - YORIG ) / YCELL
                CHK_Y  = CHK_Y - INT( CHK_Y )
                IF( DABS( CHK_X ) > 0.001D0  .OR.
     &              DABS( CHK_Y ) > 0.001D0       ) THEN
                    SFLAG = .TRUE.
                    MESG = 'ERROR: Grid origins not compatible ' //
     &                     'between ' // DATDESC // ' in ' // 
     &                     CRLF() // BLANK10 // FILDESC( 1:L ) // 
     &                     ' and initialized values.'
                    CALL M3MSG2( MESG ) 

                END IF

C.................  If offset has been set, then check to ensure its the same
                IF( OFFLAG ) THEN

C.....................  If file has different origin from the subgrid...
                    IF( X0 .NE. XORIG .OR. 
     &                  Y0 .NE. YORIG       ) THEN

                        XO = INT( ( X0 - XORIG ) / XCELL )
                        YO = INT( ( Y0 - YORIG ) / YCELL )
                        IF( XOFF .NE. XO .OR.
     &                      YOFF .NE. YO      ) THEN

                            SFLAG = .TRUE.
                            MESG = 'ERROR: Subgrid offset for ' //
     &                          DATDESC // ' in ' // CRLF() // BLANK10// 
     &                          FILDESC( 1:L ) // 'is ' //
     &                         'inconsistent with initialized values.'
                            CALL M3MSG2( MESG ) 


                        END IF

C.....................  If file has same origin as subgrid
                    ELSE

C.........................  Check that current subgrid is the same as the 
C                           previous subgrid
                        IF ( NCOLS .NE. NC       .OR.
     &                       NROWS .NE. NR       .OR.
     &                       DABS( XORIG-X0 ) > 0.001D0 .OR.
     &                       DABS( YORIG-Y0 ) > 1.00D0     ) THEN

                             SFLAG = .TRUE.
                             MESG = 'ERROR: Columns, rows, x-origin, '//
     &                          'or y-origin for ' //DATDESC //' in ' 
     &                          //CRLF() //BLANK10 //FILDESC( 1:L ) // 
     &                          'are inconsistent with values from ' // 
     &                          GRDNM
                             CALL M3MSG2( MESG ) 

                        END IF

                    END IF

C.................  If offset for final subgrid hasn't been set yet...
                ELSE

C.....................  Compute possible offset from upper right hand corner,
C                       and if there is one, set flag
                    XOFF_A = INT( ( XORIG + NCOLS * XCELL   ) - 
     &                            ( X0    + NC    * XCELL3D ) ) / XCELL
                    YOFF_A = INT( ( YORIG + NROWS * YCELL   ) - 
     &                            ( Y0    + NR    * YCELL3D ) ) / YCELL

C.....................  Compute possible offset from origin, and if there is 
C                       one, set flag
                    XOFF_A = INT( ( X0 - XORIG ) / XCELL )
                    YOFF_A = INT( ( Y0 - YORIG ) / YCELL )
                    
C.....................  Reset origin and number of cells to latest grid
                    IF( FTYPE /= 'DOT' ) GRDNM = GDNAM3D

C.....................  Only store grid and offset parameters if the subgrid is
C                       not temporary
                    IF ( CHKLEVEL .LE. CHK_SUBGRID ) THEN
                        XOFF = XOFF_A
                        YOFF = YOFF_A
                        IF( XOFF .NE. 0 .OR. YOFF .NE. 0 ) 
     &                      OFFLAG = .TRUE.
                        XDIFF = NCOLS - NC
                        YDIFF = NROWS - NR
                        XORIG = X0
                        YORIG = Y0
                        NCOLS = NC
                        NROWS = NR
                        NGRID = NCOLS * NROWS
                    END IF

                END IF

            END IF

C.........  Initialize grid information
        ELSE

            GFLAG = .TRUE.
            GRDNM = GDNAM3D
            GDTYP = GDTYP3D
            P_ALP = P_ALP3D
            P_BET = P_BET3D
            P_GAM = P_GAM3D
            XCENT = XCENT3D
            YCENT = YCENT3D
            XORIG = X0
            YORIG = Y0
            XCELL = XCELL3D
            YCELL = YCELL3D
            NCOLS = NC
            NROWS = NR
            NGRID = NCOLS * NROWS

            MESG = 'NOTE: Grid settings initialized using ' // 
     &             DATDESC // ' in ' // CRLF()// BLANK10 // 
     &             FILDESC( 1:L ) // '.'

            CALL M3MSG2( MESG )

        ENDIF

        IF( SFLAG ) THEN

            MESG = 'ERROR: Grid parameters for ' // DATDESC // ' in ' //
     &             CRLF() // BLANK10 // FILDESC( 1:L ) //
     &             ' are inconsistent with initialized values.'
            CALL M3MSG2( MESG )

        END IF

        IF( SFLAG ) EFLAG = SFLAG

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END

