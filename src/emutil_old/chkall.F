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

C.... SUBROUTINE CHKALL compares input NetCDF file header info with 
C     output UAM info based on USERIN file

      SUBROUTINE CHKALL( IERR )

      include 'PARMS3.EXT'
      include 'FDESC3.EXT'
      include 'param.inc'
      include 'uamcntl.inc'
      include 'uamregn.inc'

      COMMON /TIMES/ JDATE, EDDATE, JTIME, EDTIME, TSTEP

C.... Local variables
      INTEGER*4 EDATE3D  ! End date in NetCDF file (YYYYDDD)
      INTEGER*4 ETIME3D  ! End time in NetCDF file (HHMMSS)
      INTEGER*4 EDDATE   ! End date of output file (YYYYDDD)
      INTEGER*4 EDTIME   ! End time of output file (HHMMSS)
      INTEGER*4 IERR     ! Error flag
      INTEGER*4 JDATE    ! Date variable
      INTEGER*4 JTIME    ! Time variable
      INTEGER*4 TSTEP    ! Time step = 1 hour
 

C.... Functions and their descriptions
      INTEGER NOS_TO_S   ! Converts HHMM real to HHMMSS integer
      REAL    TIME
      REAL    S_TO_NOS   ! Converts HHMMSS integer to HHMM real
      INTEGER ITIME

      NOS_TO_S(TIME ) = (INT(TIME/100.))*10000 +
     &                  (INT(TIME)-INT(TIME/100.)*100)*1000

      S_TO_NOS(ITIME) = REAL(ITIME/100)

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

C.... Time information
C.... Initialize output times and dates for NCD2UAM in HHMMSS and YYYYDDD 
      JDATE  = 1900000 + IBGDAT
      EDDATE = 1900000 + IENDAT
      JTIME  = NOS_TO_S( BEGTIM )
      EDTIME = NOS_TO_S( ENDTIM ) 

C.... Adjust output times and dates for HR=24 in HHMMSS, HHMM, YYYYDDD, YYDDD
      CALL NEXTIME (JDATE , JTIME , 0)
      CALL NEXTIME (EDDATE, EDTIME, 0)
      IBGDAT = JDATE  - 1900000
      IENDAT = EDDATE - 1900000
      BEGTIM = S_TO_NOS( JTIME  )
      ENDTIM = S_TO_NOS( EDTIME )

C.... Calculate end date and time for input file
      EDATE3D = SDATE3D
      ETIME3D = STIME3D
      CALL NEXTIME (EDATE3D, ETIME3D, MXREC3D*TSTEP3D)

C.... Compare output dates and times with input dates and times
      IF    ( JDATE .LT. SDATE3D ) THEN
        WRITE (*,*) 'IBGDAT starts before first date in input file'
        GOTO 9999
      ELSEIF( JDATE .EQ. SDATE3D ) THEN
        IF  ( JTIME .LT. STIME3D ) THEN
          WRITE (*,*) 
     &          'BEGTIM starts before first time on day 1 of input file'
          GOTO 9999
        ENDIF
      ENDIF

      IF    ( EDDATE .GT. EDATE3D ) THEN
        WRITE (*,*) 'IENDAT ends after last date in input file'
        GOTO 9999
      ELSEIF( EDDATE .EQ. EDATE3D ) THEN
        IF  ( EDTIME .GT. ETIME3D ) THEN
          WRITE (*,*)
     &          'ENDTIM ends after last time on last day of input file'
          GOTO 9999
        ENDIF
      ENDIF

      IF ( TSTEP3D .NE. 10000 ) THEN
        WRITE (*,*) ' Input file time step is not equal to 1 hour'
        GOTO 9999
      ENDIF

C.... Initialize time step
      TSTEP  = TSTEP3D

      RETURN

9999  IERR = 1
      RETURN

      END
