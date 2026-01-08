
        SUBROUTINE HDRMISS3

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION: 
C     This subroutine initializes the I/O API header to missing and default 
C     values appropriate for emission files.
C
C  PRECONDITIONS REQUIRED:
C     
C  SUBROUTINES AND FUNCTIONS CALLED: 
C
C  REVISION  HISTORY:
C       Created 3/99 by M Houyoux
C
C***********************************************************************
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
C****************************************************************************

        USE M3UTILIO

        IMPLICIT NONE

C...........   INCLUDES:
C        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
C        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

        CHARACTER(16) :: PROGNAME = 'HDRMISS3'    ! Program name

C***********************************************************************
C   begin body of subroutine HDRMISS3

        FTYPE3D = GRDDED3
        SDATE3D = 0
        STIME3D = 0
        TSTEP3D = 0
        NCOLS3D = 1
        NROWS3D = 1
        NLAYS3D = 1
        NTHIK3D = 1
        NVARS3D = 0
        GDTYP3D = IMISS3
        P_ALP3D = AMISS3
        P_BET3D = AMISS3
        P_GAM3D = AMISS3
        XCENT3D = AMISS3 
        YCENT3D = AMISS3
        XORIG3D = AMISS3
        YORIG3D = AMISS3
        XCELL3D = AMISS3
        YCELL3D = AMISS3
        VGTYP3D = IMISS3
        VGTOP3D = AMISS3
        VGLVS3D = 0.      ! array
        GDNAM3D = ' '

        VNAME3D = ' '     ! array
        VTYPE3D = 0       ! array
        UNITS3D = ' '     ! array
        VDESC3D = ' '     ! array
        FDESC3D = ' '     ! array


C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats.............94xxx
 
94010   FORMAT( 10( A, :, I6, :, 2X ) )

        END SUBROUTINE HDRMISS3
