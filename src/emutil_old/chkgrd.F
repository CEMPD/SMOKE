C ****************************************************************************
C *
C * Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C *                System
C * File: @(#)$Id$
C *
C * COPYRIGHT (C) 1996, MCNC--North Carolina Supercomputing Center
C * All Rights Reserved
C *
C * See file COPYRIGHT for conditions of use.
C *
C * Environmental Programs Group
C * MCNC--North Carolina Supercomputing Center
C * P.O. Box 12889
C * Research Triangle Park, NC  27709-2889
C *
C * env_progs@mcnc.org
C *
C * Pathname: $Source$
C * Last updated: $Date$ 
C *
C ****************************************************************************/
 
C.... SUBROUTINE CHKGRD compares input NetCDF file header grid info with 
C     output UAM info based on USERIN file

      SUBROUTINE CHKGRD( IERR )

      include 'PARMS3.EXT'
      include 'FDESC3.EXT'
      include 'param.inc'
      include 'uamcntl.inc'
      include 'uamregn.inc'

C.... Local variables
      INTEGER*4 IERR     ! Error flag
 
      IERR = 0

C.... Coordinate system
      IF    ( GDTYP3D .NE. UTMGRD3 ) THEN
        WRITE (*,*) 'Input grid type not UTM'
        GOTO 9999
      ELSEIF( IZONE .NE. P_ALP3D ) THEN
        WRITE (*,*) 'UTM Zones do not match'
        GOTO 9999
      ELSEIF( (XORIG+UTMX) .NE. XORIG3D ) THEN
        WRITE (*,*) 'XORIG+UTM is not equal to XORIG3D'
        GOTO 9999
      ELSEIF( (YORIG+UTMY) .NE. YORIG3D ) THEN
        WRITE (*,*) 'YORIG+UTM is not equal to YORIG3D'
        GOTO 9999
      ENDIF
C.... Grid parameters
      IF     (DXCELL .NE. XCELL3D ) THEN
        WRITE (*,*) 'DXCELL is not equal to XCELL3D'
        GOTO 9999
      ELSEIF( DYCELL .NE. YCELL3D ) THEN
        WRITE (*,*) 'DYCELL is not equal to YCELL3D'
        GOTO 9999
      ELSEIF( NXCELL .NE. NCOLS3D ) THEN
        WRITE (*,*) 'NXCELL is not equal to NCOLS3D'
        GOTO 9999
      ELSEIF( NYCELL .NE. NROWS3D ) THEN
        WRITE (*,*) 'NYCELL is not equal to NROWS3D'
        GOTO 9999
      ENDIF

      RETURN

9999  IERR = 1
      RETURN

      END
