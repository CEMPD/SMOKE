
      SUBROUTINE MOBILE( INERR, TI, TMMI, MODEV, NSAV, RFLAG, 
     &                   EFACT, DFACT, EFSAVND, EFSAVDI )
 
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
C Portions COPYRIGHT (C) 1998, MCNC--North Carolina Supercomputing Center
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
C****************************************************************************

C...........   MODULES for public variables
C.........  This module contains emission factor tables and related
        USE MODEMFAC

C.........  This module contains the information about the source category
        USE MODINFO

C...........   This module is the derived meteorology data for emission factors
        USE MODMET

        IMPLICIT NONE

      INCLUDE 'EVAL.I'
      INCLUDE 'FLAGS2.I'
      INCLUDE 'FLAGS3.I'
      INCLUDE 'FLAGS4.I'
      INCLUDE 'IMPAR7.I'
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
      INCLUDE 'EMCNST3.EXT'   !  Emissions constants
      INCLUDE 'CONST3.EXT'    !  Physical constants
      INCLUDE 'M5CNST3.EXT'   !  Mobile5a/b constants
      INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
      REAL, PARAMETER :: CTOF = 9./5.        ! Celcius to Farhenheit 

cmh add functions
      CHARACTER*2   CRLF
      EXTERNAL      CRLF
c
cmh add LMOS definition of ndi_di
      CHARACTER*1 ndi_di
c
cmh add declarations

C.......  Subroutine arguments
      INTEGER, INTENT (IN OUT):: INERR
      INTEGER, INTENT (IN)    :: TI     ! ambient tmpr index
      INTEGER, INTENT (IN)    :: TMMI   ! min/max tmpr index
      INTEGER, INTENT (IN)    :: MODEV  ! Unit number of MOBILE5a input file
      INTEGER, INTENT (IN)    :: NSAV   ! No. PSIs to save during call
      LOGICAL, INTENT (IN)    :: RFLAG  ! true: reset saved EFs
      REAL   , INTENT (OUT)   :: EFACT( NVTYPE, NNDI )
      REAL   , INTENT (OUT)   :: DFACT( NVTYPE, NDIU )
      REAL   , INTENT (OUT)   :: EFSAVND(NVTYPE,MXPPGRP,NTMPR,NNDI)
      REAL   , INTENT (OUT)   :: EFSAVDI(NVTYPE,MXPPGRP,NVLDTMM,NDIU)

C.......  Other local variables because of SMOKE
      INTEGER       II, IOS, IY, JJ, M, T   ! indices and counters
      INTEGER       MCNT

      REAL          TMPRTR  ! Ambient temperature
      REAL          TMIN    ! Minimum daily temperature
      REAL          TMAX    ! Maximum daily temperature

      LOGICAL       WFLAG
      LOGICAL    :: FIRSTTEM = .TRUE.
      LOGICAL    :: FIRSTIME = .TRUE.

      CHARACTER*8   M5VRSION
      CHARACTER*300 MESG
cmh
cmh add commons
      COMMON /WARNING/ WFLAG
      COMMON /TOVRIDE/ TMPRTR, TMIN, TMAX
      COMMON /MULTISC/ M5VRSION
c
cmh add other declarations for MOBILE5
      INTEGER   ICY, ICY1, ICYOLD, IROLD, LEAST1, MEVOLD, NIY, NIYOLD

cmh add save
      SAVE MCNT, FIRSTTEM, FIRSTIME
C
cmh add declarations of LMOS common
      REAL EFDNL, EFEVND, VDNL, DNLGPM
c
cmh  Add LMOS and SCENE3 common blocks
c
      COMMON /LMOS/ EFDNL(9),EFEVND(9),VDNL,DNLGPM
      COMMON /SCENE3/ ndi_di

cmh add program name
      CHARACTER*16 :: PROGNAME = 'MOBILE'

cmh ***************************************************************************
cmh start of subroutine MOBIL5B

cmh add firstime check of EMS_LOC
      IF( FIRSTIME ) THEN

          FIRSTIME = .FALSE.
          CALL ENVSTR( 'EMS_LOC', 'Path for I/M tables', ' ',
     &                  MESG, IOS )
          IF( IOS .NE. 0 ) THEN

              MESG = 'ERROR: Environment variable EMS_LOC not set!' //
     &                CRLF() // BLANK10 // 'Must set it to a ' //
     &               'directory with the TECH12.D and IMDATA.D files.'
              CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
          END IF

      END IF
 
cmh add conversion of temperatures from Kelvin to Farhenheit
      TMPRTR = MAX( ( TMPRTR - CTOK ) * CTOF + 32., 0. )
      TMPRTR = MIN( TMPRTR, 110. )
      TMIN   = MAX( ( TMIN   - CTOK ) * CTOF + 32., 0. )
      TMAX   = MIN( ( TMAX   - CTOK ) * CTOF + 32., 110. )

      IOUGEN = MODEV
      IOUREP = INIT3()
      IOUASK = INIT3()
      IOUERR = INIT3()
cmh   IOUGET = 5
cmh   IOBAT = 19
 
      WFLAG = .FALSE.     ! Initialize warning message write to off
 
      IF( RFLAG ) THEN
          IREAD = 0 ! Reset MOBILE5's indicator for reading TECH12 file
          WFLAG = .TRUE.  ! Turn warning messages on
      END IF
c
      MCNT  = 0 ! Reset scenario counter for EFSAVE 

C.......  Initialize saved emission factor variables
      IF( RFLAG .AND. NSAV .GT. 1 ) THEN 
          EFSAVND = 0.0   ! 4-d array
          EFSAVDI = 0.0   ! 4-d array
      ENDIF

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
      CALL CONSEC(*901,*902)

cmh add reset of TEMFLG for use with SMOKE.  Always want it to be set to assume
c   mobile temperatures will be hour-based.  TEMFLG in FLAGS3.I.
      IF( TEMFLG .NE. 2 ) THEN

          IF( FIRSTTEM ) THEN
              MESG = 'Resetting TEMFLG in MOBILE5 inputs to 2 for ' //
     &               'use in SMOKE' // CRLF() // BLANK5
              CALL M3WARN( 'MOBILE', 0, 0, MESG )
              FIRSTTEM = .FALSE.
          ENDIF

          TEMFLG = 2
      ENDIF

      CALL ONESEC(INERR,*903)
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

   20 IF( MCNT .EQ. NSAV ) GO TO 88    ! To safe exit

      CALL PARSEC(ICY,INERR,*98)       ! *98 is error exit
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

      IF( MCNT .EQ. 1 ) THEN   ! NOTE: NOT an ELSEIF !

          EFACT = 0.0  ! array
          DFACT = 0.0  ! array

          DO JJ = 1,NVTYPE
              EFACT( JJ,1 ) = EFFTP ( 2,JJ )   ! Exhaust CO
              EFACT( JJ,2 ) = EFFTP ( 3,JJ )   ! Exhaust NOX
              EFACT( JJ,3 ) = EFEXH ( JJ )     ! Exhaust VOC
              EFACT( JJ,4 ) = EFEVAP( JJ )     ! Evaporative VOC
              EFACT( JJ,5 ) = EFRUNL( JJ )     ! Running Loss VOC
              EFACT( JJ,6 ) = EFRSTL( JJ )     ! Resting Loss VOC
              EFACT( JJ,7 ) = EFREFL( JJ )     ! Refueling Loss VOC
cmh           EFACT( JJ,8 ) = RLGGAL( JJ )     ! Refueling loss VOC g/gal
cmh           EFACT( JJ,11) = RSTGPH( JJ )     ! Resting Loss VOC g/hr
          END DO

          DO JJ = 1,4
              DFACT( JJ,1 ) = GREVP ( 2,JJ )  ! Weighted Diurnal Evap VOC
              DFACT( JJ,2 ) = EFDNL (   JJ )  ! LMOS Diurnal Evap VOC
              DFACT( JJ,3 ) = GREVP ( 1,JJ )  ! Hot Soak Evap VOC
              DFACT( JJ,4 ) = GREVP ( 4,JJ )  ! Crankcase Evap VOC
          END DO

          DFACT( 8, 1 ) = GREVP ( 2,8 )   ! Note: 8 is MC, 9 is combined trucks
          DFACT( 8, 2 ) = EFDNL (   8 )  
          DFACT( 8, 3 ) = GREVP ( 1,8 )
          DFACT( 8, 4 ) = GREVP ( 4,8 )   ! May always be zero (?)
          DFACT( 9, 1 ) = GREVP ( 2,9 )   ! Note: 9 is MC, 9 is combined trucks
          DFACT( 9, 2 ) = EFDNL (   9 )  
          DFACT( 9, 3 ) = GREVP ( 1,9 )
          DFACT( 9, 4 ) = GREVP ( 4,9 )   ! May always be zero (?)

      END IF

C.......  Make sure that the number of scenarios is within the compiled limits
C         and that we are needing to save the factors for this call.
      IF( MCNT .LE. MXM5SCEN .AND. 
     &    NSAV .GT. 1              ) THEN

C...........  Screen for cases where we are not running MOBILE for non-diurnal
C             factors
          IF( TI .GT. 0 ) THEN
 
              DO JJ = 1,NVTYPE
        	  EFSAVND( JJ, MCNT, TI, 1 ) = EFFTP ( 2,JJ )
        	  EFSAVND( JJ, MCNT, TI, 2 ) = EFFTP ( 3,JJ )
        	  EFSAVND( JJ, MCNT, TI, 3 ) = EFEXH ( JJ ) 
        	  EFSAVND( JJ, MCNT, TI, 4 ) = EFEVAP( JJ )  
        	  EFSAVND( JJ, MCNT, TI, 5 ) = EFRUNL( JJ )
        	  EFSAVND( JJ, MCNT, TI, 6 ) = EFRSTL( JJ )
        	  EFSAVND( JJ, MCNT, TI, 7 ) = EFREFL( JJ )
              END DO

          END IF

C...........  Screen for cases where we are not running MOBILE for diurnal
C             factors
          IF( TMMI .GT. 0 ) THEN
 
              DO JJ = 1,4
        	  EFSAVDI( JJ, MCNT, TMMI, 1 ) = GREVP ( 2,JJ )
        	  EFSAVDI( JJ, MCNT, TMMI, 2 ) = EFDNL (   JJ )  
        	  EFSAVDI( JJ, MCNT, TMMI, 3 ) = GREVP ( 1,JJ )
        	  EFSAVDI( JJ, MCNT, TMMI, 4 ) = GREVP ( 4,JJ )
              END DO

              EFSAVDI( 8, MCNT, TMMI, 1 ) = GREVP ( 2,8 )
              EFSAVDI( 8, MCNT, TMMI, 2 ) = EFDNL (   8 )
              EFSAVDI( 8, MCNT, TMMI, 3 ) = GREVP ( 1,8 )
              EFSAVDI( 8, MCNT, TMMI, 4 ) = GREVP ( 4,8 )
              EFSAVDI( 9, MCNT, TMMI, 1 ) = GREVP ( 2,9 )
              EFSAVDI( 9, MCNT, TMMI, 2 ) = EFDNL (   9 )
              EFSAVDI( 9, MCNT, TMMI, 3 ) = GREVP ( 1,9 )
              EFSAVDI( 9, MCNT, TMMI, 4 ) = GREVP ( 4,9 )

          END IF

      END IF   ! Overflow check for EFSAVND and EFSAVDI
cmh end add

      LEAST1=1
      GOTO 20    ! to head of scenario read loop
C
   30 CALL QUITER(0.0,0,28,INERR)
c
cmh add bounding exit errors
      IF( MCNT .GT. MXM5SCEN ) THEN

            WRITE( MESG, 94010 ) 
     &         'Number of scenarios in single MOBILE5 run:', MCNT, 
     &          CRLF() // BLANK5 //
     &         'Max permitted MOBILE5 scenarios (MXM5SCEN):', MXM5SCEN,
     &          CRLF() // BLANK5 //
     &          'Maximum exceeded.'
            CALL M3MSG2( MESG )
            CALL M3EXIT( 'MOBILE', 0, 0, 'Error in MOBILE5 inputs', 2 )

      ENDIF

88    CONTINUE

      RETURN

cmh end add
c
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
cmh   99 RETURN
C
C  End-of-file on 1st read of call => input for call missing entirely.
C  Calling program can branch to attach a new input file or terminate
C  calls to MOBILE.
C
cmh  999 RETURN 1
c
cmh add
      MESG = 'Problem reading scenario section of MOBILE5b inputs!'
      CALL M3MSG2( MESG )
      GOTO 999

901   MESG = 'Error exit status 1 while reading control section ' //
     &       'of MOBILE5b inputs!'
      CALL M3MSG2( MESG )
      GOTO 999

902   MESG = 'Error exit status 2 while reading control section ' //
     &       'of MOBILE5b inputs!'
      CALL M3MSG2( MESG )
      GOTO 999

903   MESG = 'Problem reading one-time section of MOBILE5b inputs!'
      CALL M3MSG2( MESG )
      GOTO 999

999   CONTINUE      ! error abort because of misformatted inputs

      MESG = 'Error or misformatted inputs encountered ' //
     &       'during MOBILE5b call.'
      CALL M3EXIT( 'MOBILE', 0, 0, MESG, 2 )
c
cmh end add

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx
 
94010   FORMAT( 10 ( A, :, I5, :, 2X ) )

      END
