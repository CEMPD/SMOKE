
      SUBROUTINE CHKTMP( IERR, JDATE, EDDATE, JTIME, EDTIME, TSTEP )
      
C***************************************************************************
C  subroutine body starts at line 66
C
C  DESCRIPTION:
C     SUBROUTINE CHKTMP compares input NetCDF file header time info with 
C     output UAM info based on USERIN file
C
C***************************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1998, MCNC--North Carolina Supercomputing Center
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

      IMPLICIT NONE

      INCLUDE 'PARMS3.EXT'
      INCLUDE 'FDESC3.EXT'
      INCLUDE 'param.inc'
      INCLUDE 'uamcntl.inc'

C.... Local variables
      INTEGER*4 EDATE3D  ! End date in NetCDF file (YYYYDDD)
      INTEGER*4 ETIME3D  ! End time in NetCDF file (HHMMSS)
      INTEGER*4 EDDATE   ! Last day of output file (YYYYDDD)
      INTEGER*4 EDTIME   ! Start of last hour of output file (HHMMSS)
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

C***********************************************************************
C   begin body of subroutine CHKTMP

      IERR = 0

C.... Time information
C.... Initialize output times and dates 
      JDATE  = 1900000 + IBGDAT    ! in YYYYDDD
      EDDATE = 1900000 + IENDAT    ! in YYYYDDD
      JTIME  = NOS_TO_S( BEGTIM )  ! in HHMMSS
      EDTIME = NOS_TO_S( ENDTIM )  ! in HHMMSS

C.... Adjust output times and dates in case HR=24
      CALL NEXTIME (JDATE , JTIME , 0)
      CALL NEXTIME (EDDATE, EDTIME, 0)
      IBGDAT = JDATE  - 1900000    ! in YYDDD
      IENDAT = EDDATE - 1900000    ! in YYDDD
      BEGTIM = S_TO_NOS( JTIME  )  ! in HHMM
      ENDTIM = S_TO_NOS( EDTIME )  ! in HHMM

      CALL NEXTIME( EDDATE, EDTIME, -10000 )

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
