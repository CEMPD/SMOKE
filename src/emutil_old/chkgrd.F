 
C Version @(#)$Id$ $Source$ $Date$ 
 
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
