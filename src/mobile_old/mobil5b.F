
      SUBROUTINE MOBILE(INERR , ITEMP, MODEV, NSAV, EFACT, DFACT )
 
cmh   SUBROUTINE MOBILE(INERR, *)
C
C  All "cmh" comments are modifications for the SMOKE modeling system
C  made by Marc Houyoux at MCNC.
C
C***************************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1996, MCNC--North Carolina Supercomputing Center
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
C****************************************************************************
C
C  Called by user's program (example: Air Quality Analysis
C  System).
C
C  Calls ADJUST, BDSAVE, CONSEC, EFCALX, ONESEC, OUTNEW, OUTPUT,
C        PARSEC, QUITER, REINIT, RESTOR, REGMOD and SAVER.
C
C  The following variables are passed in argument lists or common
C  blocks and in MOBILE:
C
C    argument lists: ICY,INERR
C    common blocks:
C    /EVAL/   MEVAL
C    /FLAGS3/ OUTFMT
C    /IOUCOM/ IOUERR
C    /REGION/ IREJN
C    /SAVE01/ JCALL
C    /YEARS4/ IYEND
C
C  Output on return:
C
C    parameter list: INERR
C
C  Local variable dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  ICYOLD   I    calendar year from previous pass of scenario loop
C  IROLD    I    region from previous pass of scenario loop
C                      0 = message has not yet been given
C                      1 = message has been given: not again!
C  IY       I    loop counting variable for interpolation runs
C  LEAST1   I    flag: 0 = OUTPUT has not yet been called.
C                      1 = OUTPUT has been called at least 1 time.
C  MEVOLD   I    calendar month from previous pass of scenario loop
C  NIY      I    number of year interations for by-month interpolation
C  NIYOLD   I    number of year interations from previous pass
C                of scenario loop
C
C
C  Notes:
C
C  For MOBILE4.1 the subroutine MOBILE was created from MAIN.
C  The /BYMYCd/ CBs were added.
C  May-18-1993 @ ARC-bsg Subtask 238 the common block IMPAR8
C  was removed, it contained the variable ICALC
C  June-15-1993 @ ARC-bk Subtask 244 the common block IMPAR5
C  was modified, the HDGV I/M credit was expanded to include
C  an initial value for NOx.
C  24-August-1994 @CSC-pme request 2-446 Include file BASE12.I was removed
C  (24-May-1996) @DynTel-yc Request 2-621 The calendar year was expanded
C                 from year 2020 to year 2051
C
      INCLUDE 'EVAL.I'
      INCLUDE 'FLAGS2.I'
      INCLUDE 'FLAGS3.I'
      INCLUDE 'FLAGS4.I'
      INCLUDE 'IOUCOM.I'
      INCLUDE 'REGION.I'
      INCLUDE 'SAVE01.I'
      INCLUDE 'YEARS4.I'
c
cmh add includes from Mobile5b
      INCLUDE 'RESUL1.I'
      INCLUDE 'RESUL2.I'
      INCLUDE 'EVAPGR.I'
      INCLUDE 'VDATA.I'
c
cmh add includes from SMOKE
      INCLUDE 'MBDIMS3.EXT'
      INCLUDE 'IODECL3.EXT'
c
cmh add external functions
      CHARACTER*2   CRLF
      EXTERNAL      CRLF

cmh add parameter statements
      CHARACTER*5   BLANK5
      PARAMETER   ( BLANK5 = ' ' )
c
cmh add LMOS definition of ndi_di
      CHARACTER*1 ndi_di
c
cmh add declarations
      REAL    EFACT ( MXNDET, NVTYPE )
      REAL    DFACT ( MXDET, NVTYPE )
      REAL    EFSAVE( NVTYPE, MXM5PSI, SAVTMP, MXEMIS )
      REAL    TMPRTR  ! Ambient temperature
      REAL    TMIN    ! Minimum daily temperature
      REAL    TMAX    ! Maximum daily temperature
      INTEGER MODEV                ! Unit number of MOBILE5a input file
      INTEGER ITEMP   ! When == 0, tmpr is only change from previous call
      INTEGER MCNT
      INTEGER NSAV    ! Number of PSIs to save during call to subroutine
      INTEGER TCNT
      LOGICAL WFLAG
      CHARACTER*8   M5VRSION
      CHARACTER*256 MESG
cmh
cmh add commons
      COMMON /WARNING/ WFLAG
      COMMON /TOVRIDE/ TMPRTR, TMIN, TMAX
      COMMON /MULTISC/ EFSAVE, M5VRSION
c
cmh add save
      SAVE MCNT, TCNT
C
cmh add lines
c
cmh  Add LMOS and SCENE3 common blocks
c
      COMMON /LMOS/ EFDNL(9),EFEVND(9),VDNL,DNLGPM
      COMMON /SCENE3/ ndi_di
 
      IOUGEN = MODEV
      IOUGET = 5
      IOUREP = INIT3()
      IOUASK = INIT3()
      IOUERR = INIT3()
      IOBAT = 19
 
      WFLAG = .FALSE.     ! Initialize warning message write to off
 
      IF( ITEMP .EQ. 1 ) THEN
          IREAD = 0 ! Reset MOBILE5's indicator for reading TECH12 file
          ITEMP = 0 ! Indicate that only temperature is changing
          WFLAG = .TRUE.  ! Turn warning messages on
          TCNT  = 1 ! Reset temperature counter for EFSAVE 
      ELSE
          TCNT  = TCNT + 1
      ENDIF
c
      MCNT  = 0 ! Reset scenario counter for EFSAVE 

cmh end add

cmh add comparison of MOBILE5 versions
      INCLUDE 'VNAME.I'

      IF( M5VRSION .NE. VERSON( 1:8 ) .AND. 
     &    M5VRSION .NE. ' '                 ) THEN
          MESG = 'In MPREF, MOBILE5A packet found, but program ' //
     &           'compiled using ' // VERSON( 1:8 )
          CALL M3EXIT( 'MOBILE', 0, 0, MESG, 2 )
      ENDIF
C
      IF(JCALL.EQ.0) GOTO 10
C
C  Subroutine call, setup run.
C
      CALL REINIT
      IF(JCALL.EQ.1) CALL BDSAVE
      IF(JCALL.EQ.2) CALL RESTOR
      JCALL=2
C
C  Initialize non-COMMON variables (INERR & local variables).
C
   10 INERR=0
      ICYOLD=-1
      IROLD=-1
      MEVOLD=-1
      NIYOLD=-1
      LEAST1=0
C
C  Control and One-Time Data Sections:
C
C  Get the execution control flags and parameters, and then use them to read
C  in data, if any, that is read in once and applied to all scenarios.  Note:
C  alternate RETURN1 is the fatal read error exit: execution stops.
C
      IF(IOUASK.EQ.0.OR.IOUASK.EQ.9) WRITE(IOUASK,110)
  110 FORMAT(' Reading information.')
C
      CALL CONSEC(*99,*999)
      CALL ONESEC(INERR,*99)
C
C  ljn hc=3 3/5/93 Warning MSG in QUITER
C
      IF(HCFLAG.EQ.3.AND.OUTFMT.NE.3.AND.OUTFMT.NE.5)
     *                    CALL QUITER(0.,0,153,INERR)
C
C  only get error diagnostics out of the run.
C
cmh20 CALL PARSEC(ICY,INERR,*98)

cmh add if statement which prevents PARSEC from having to
cmh use the end of file as the indicator to stop running

   20 IF( MCNT .EQ. NSAV ) GO TO 98  

      CALL PARSEC(ICY,INERR,*98)
C
C  50+ errors is (arbritrary) run limit => stop execution.  Otherwise,
C  one or more input processing errors => abort current scenario & then do
C  input processing for error diagnostics on the remaining scenarios.  INERR
C  is never reset to 0, but instead accumulates until the run ends or the
C  cutoff value (50) is reached.
C
cmh comment out this check, as per LMOS Mobile5a change
c
c     IF(INERR.GT.50) GOTO 30
      IF(INERR.GT.0) GOTO 20
C
      IF(IOUASK.EQ.0.OR.IOUASK.EQ.9) WRITE(IOUASK,115)
  115 FORMAT(' Performing calculations.')
C
C  For interpolated (non Jan 1st) emissions, loop through ICY and ICY+1,
C  save intermediate results, then interpolate to the evaluation month.
C
      NIY=2
      IF(MEVAL.EQ.1.OR.ICY.EQ.IYEND) NIY=1
C
      DO 35 IY=1,NIY
      ICY1=ICY+IY-1
      IF(NIY.GT.1 .OR.
     *   NIYOLD.NE.NIY .OR.
     *   ICYOLD.NE.ICY1 .OR.
     *   IROLD.NE.IREJN .OR.
     *   MEVOLD.NE.MEVAL) CALL REGMOD(ICY1,INERR)
c
cmh comment out this check, as per LMOS Mobile5a change
c
c     IF(INERR.GT.50) GOTO 30
      IF(INERR.GT.0) GOTO 20
C
      ICYOLD=ICY1
      IROLD=IREJN
      MEVOLD=MEVAL
      NIYOLD=NIY
C
      CALL EFCALX(ICY1,INERR,IY)
      IF(IDLFLG.EQ.2) CALL IDLCAL(ICY)  ! cmh: this line is new for Mobile5b
      IF(INERR.GT.50) GOTO 30
      IF(INERR.GT.0) GOTO 20
C
      IF(NIY.GT.1) CALL SAVER(IY)
C
  35  CONTINUE
C
      IF(NIY.GT.1) CALL ADJUST
C
      IF(IOUASK.EQ.0.OR.IOUASK.EQ.9) WRITE(IOUASK,120)
  120 FORMAT(' Preparing output.')
C
cmh   CALL OUTPUT(ICY)
c
cmh add lines
      MCNT = MCNT + 1

      IF( MCNT .EQ. 1 .AND. TCNT .EQ. 1 .AND. NSAV .GT. 1 ) THEN

          DO II = 1, MXEMIS
              DO T = 1, SAVTMP
                  DO M = 1, MXM5PSI 
                      DO JJ = 1, NVTYPE 
                          EFSAVE( JJ, M, T, II ) = 0.0
                      ENDDO
                  ENDDO
              ENDDO
          ENDDO

      ENDIF

      IF( MCNT .EQ. 1 ) THEN   ! NOTE: NOT and ELSEIF

          DO JJ = 1, NVTYPE
              DO II = 1, MXNDET
                  EFACT ( II,JJ ) = 0.0
              ENDDO
              DO II = 1, MXDET
                  DFACT( II,JJ ) = 0.0
              ENDDO
          ENDDO

          DO JJ = 1,NVTYPE
              EFACT( 1,JJ ) = EFFTP ( 2,JJ )   ! Exhaust CO
              EFACT( 2,JJ ) = EFFTP ( 3,JJ )   ! Exhaust NOX
              EFACT( 3,JJ ) = EFEXH ( JJ )     ! Exhaust VOC
              EFACT( 4,JJ ) = EFEVAP( JJ )     ! Evaporative VOC
              EFACT( 5,JJ ) = EFRUNL( JJ )     ! Running Loss VOC
              EFACT( 6,JJ ) = EFRSTL( JJ )     ! Resting Loss VOC
              EFACT( 7,JJ ) = EFREFL( JJ )     ! Refueling Loss VOC
cmh           EFACT( 8, JJ ) = RLGGAL( JJ )    ! Refueling loss VOC g/gal
cmh           EFACT( 11,JJ ) = RSTGPH( JJ )    ! Resting Loss VOC g/hr
          ENDDO

          DO JJ = 1,4
              DFACT( 1, JJ ) = GREVP ( 2,JJ )  ! Weighted Diurnal Evap VOC
              DFACT( 2, JJ ) = EFDNL (   JJ )  ! LMOS Diurnal Evap VOC
              DFACT( 3, JJ ) = GREVP ( 1,JJ )  ! Hot Soak Evap VOC
              DFACT( 4, JJ ) = GREVP ( 4,JJ )  ! Crankcase Evap VOC
          ENDDO

          DFACT( 1, 8 ) = GREVP ( 2,8 )        ! Note: 8 is MC, 9 is combined trucks
          DFACT( 2, 8 ) = EFDNL (   8 )  
          DFACT( 3, 8 ) = GREVP ( 1,8 )
          DFACT( 4, 8 ) = GREVP ( 4,8 )        ! May always be zero (?)

      ENDIF

      IF( MCNT .LE. MXM5PSI .AND. 
     &    TCNT .LE. SAVTMP  .AND. 
     &    NSAV .GT. 1             ) THEN
 
          DO JJ = 1,NVTYPE
              EFSAVE( JJ, MCNT, TCNT, 1 ) = EFFTP ( 2,JJ )
              EFSAVE( JJ, MCNT, TCNT, 2 ) = EFFTP ( 3,JJ )
              EFSAVE( JJ, MCNT, TCNT, 3 ) = EFEXH ( JJ ) 
              EFSAVE( JJ, MCNT, TCNT, 4 ) = EFEVAP( JJ )  
              EFSAVE( JJ, MCNT, TCNT, 5 ) = EFRUNL( JJ )
              EFSAVE( JJ, MCNT, TCNT, 6 ) = EFRSTL( JJ )
              EFSAVE( JJ, MCNT, TCNT, 7 ) = EFREFL( JJ )
          ENDDO
 
          DO JJ = 1,4
              EFSAVE( JJ, MCNT, TCNT, 8 ) = GREVP ( 2,JJ )
              EFSAVE( JJ, MCNT, TCNT, 9 ) = EFDNL (   JJ )  
              EFSAVE( JJ, MCNT, TCNT, 10) = GREVP ( 1,JJ )
              EFSAVE( JJ, MCNT, TCNT, 11 )= GREVP ( 4,JJ )
          ENDDO
 
          EFSAVE( 8, MCNT, TCNT, 8  ) = GREVP ( 2,8 )
          EFSAVE( 8, MCNT, TCNT, 9  ) = EFDNL (   8 )
          EFSAVE( 8, MCNT, TCNT, 10 ) = GREVP ( 1,8 )
          EFSAVE( 8, MCNT, TCNT, 11 ) = GREVP ( 4,8 )

      ENDIF   ! Overflow check for EFSAVE
cmh end add

      LEAST1=1
      GOTO 20
C
   30 CALL QUITER(0.0,0,28,INERR)
C
C  If user is supplying replacement basic FTP ef parameter records, always
C  print the records contents and disposition table, even when input errors
C  have prevented the normal output routines from being called.  The only
C  exception is when a READ error / end-of-file occurs in CONSEC or ONESEC.
C
C  This section can be expanded to print other optional input.
C
   98 IF(LEAST1.EQ.0.AND.
     *  (NEWFLG.EQ.2.OR.NEWFLG.EQ.4.OR.NEWFLG.EQ.6)) CALL OUTNEW
C
c
cmh add bounding exit errors
      IF( MCNT .GT. MXM5PSI ) THEN

            WRITE( MESG, 94010 ) 
     &            'Number of scenarios in single MOBILE5 run:', MCNT, 
     &             CRLF() // BLANK5 //
     &            'Max permitted MOBILE5 scenarios (MXM5PSI):', MXM5PSI,
     &             CRLF() // BLANK5 //
     &             'Maximum exceeded.'
            CALL M3MSG2( MESG )
            CALL M3EXIT( 'MOBILE', 0, 0, 'Error in MOBILE5 inputs', 2 )

      ENDIF

      IF( NSAV .GT. 1 .AND. TCNT .GT. SAVTMP ) THEN
            WRITE( MESG, 94010 ) 
     &            'Number of saved temperatures             :', TCNT, 
     &             CRLF() // BLANK5 //
     &            'Max permitted saved temperatures (SAVTMP):', SAVTMP,
     &             CRLF() // BLANK5 //
     &             'Maximum exceeded.'
            CALL M3MSG2( MESG )
            CALL M3EXIT( 'MOBILE', 0, 0, 'Redimension and compile', 2 )

      ENDIF

cmh end add
c
   99 RETURN
C
C  End-of-file on 1st read of call => input for call missing entirely.
C  Calling program can branch to attach a new input file or terminate
C  calls to MOBILE.
C
cmh  999 RETURN 1
  999 RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx
 
94010   FORMAT( 10 ( A, :, I5, :, 2X ) )

      END
