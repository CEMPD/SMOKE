C Version "@(#)$Id$ $Source$ $Date$ 

      SUBROUTINE GETFYR(INERR,*)
C
C  GETFYR changes the value of PRTFLG to 4 and acknowledges that the user
C  has entered a PRTFLG=5 or 6, which means they are getting the effects
C  of non-FTP. Depending on how the user sets the PRTFLG, GETFYR does 
C  different things.  When PRTFLG=5 there are non-FTP effects all years
C  until the NFTPYR  (NONFTP=2).  But when the PRTFLG=6 the user
C  enters the first model year that they want begin the phase-in to turn off 
C  non-FTP effects (NONFTP=3) and they have the option to input the reductions
C  in non-FTP effects for each pollutant.
C
C  Called by ONESEC.
C
C  Calls QUITER.
C
C  Input on call:
C
C    common blocks:
C    /IOUCOM/ IOUERR
C    /NFTPDT/ NONFTP
C    /YEARS4/ IYSTART,IYEND
C
C  Output on return:
C
C    common blocks:
C    /NFTPDT/ NFTPYR,REDUCT,NFTPV,PHASSC
C    /FLAGS4/ PRTFLG
C
C    parameter list: INERR
C
C  Notes:
C
C  Aug-18-1994 @ CSC-bsg Request 462 GETFYR was added for MOBILE5b to get the
C  value of the variable NFTPYR from the user when the NONFTP Flag is 3, 
C  NONFTP is set in GETVER/HOR and depends on PRTFLG as read from the user
C  added calls to QUITER 157 and 158 for warnings related to non-FTP
C  added call to QUITER 159 and 160 for errors related to non-FTP
C  sets PRTFLG to all pollutants (PRTFLG=4)
C  Aug-24-1994 @ CSC-bsg Request 462 added COMMON YEARS4 to bounds check using
C  easier to change variables, IY1941 and IY2020.
C  Sep-15-1994 @ CSC-bsg Request 467 added default values for the in-use
C  driving behavior record.
C  Sep-29-1994 @ CSC-bsg Request 467 modified default values for the in-use
C  driving behavior record and added Methane as a user supplied rate.  Also
C  checks THC is larger than Methane
C  Sep-30-1994 @ CSC-bsg Request 467 Changed the IYEAR definitions and the
C  common variable PHASSC to go from 0 to 6 instead of 1 to 7.
C  Jun-20-1996 @ DynTel-yc Request 2-624  Non-FTP effects have been 
C                 nullified by setting non-FTP variables to zero.
C
      INCLUDE 'FLAGS4.I'
      INCLUDE 'IOUCOM.I'
      INCLUDE 'NFTPDT.I'
      INCLUDE 'YEARS4.I'
C
C  Jun-20-1996 @ DynTel-yc Request 2-624  Non-FTP effects have been 
C                nullified by setting non-FTP variables to zero.
C                To reactivate non-FTP effects, the variables have to 
C                be refilled with non-zero values that may be available 
C                from EPA.
C
      NFTPYR=IYEND
      PHASSC(1)=0.000
      PHASSC(2)=0.000
      PHASSC(3)=0.000
      PHASSC(4)=0.000
      PHASSC(5)=0.000
      PHASSC(6)=0.000
C
C  Not user supplied, but 0th element is needed to make program smoother
C
      PHASSC(0)=0.000
C
      REDUCT(1)=0.000
      REDUCT(2)=0.000
      REDUCT(3)=0.000
C
C
      NFTPV(1)=00.000
      NFTPV(2)=00.000
      NFTPV(3)=00.000
      NFTPV(4)=00.000
C
      PRTFLG=4
C
C  Sep-15-1994 @ CSC-bsg Request 467 added user input for the in-use driving
C  behavior record.
C  Sep-29-1994 @ CSC-bsg Request 467 modified default values for the in-use
C  driving behavior record and added Methane as a user supplied rate.
C
      IF(NONFTP.EQ.3) THEN
        READ(IOUGEN,100,ERR=99,END=99) NFTPV,
     *                                 REDUCT,NFTPYR,PHASSC(1),
     *                                 PHASSC(2),PHASSC(3),PHASSC(4),
     *                                 PHASSC(5),PHASSC(6)
  100   FORMAT(4(F6.3,1X),3(F4.3,1X),I4,6(1X,F4.3))
C
C  Jun-20-1996 @ DynTel-yc Request 2-624  Non-FTP effects are nullified 
C                by setting non-FTP variables to zero.
C
        DO 10 I=1,6
   10   PHASSC(I)=0.0
        DO 20 I=1,3
   20   REDUCT(I)=0.0
        DO 30 I=1,4
   30   NFTPV(I)=0.0
C
C  Bounds checking the first year that non-FTP is turned off
C  Aug-24-1994 @ CSC-bsg Request 462 added COMMON YEARS4 to bounds check using 
C  easier to change variables, IY1941 and IY2020.
C
        IF(NFTPYR.LT.IYSTART.OR.NFTPYR.GT.(IYEND+1)) 
     *    CALL QUITER(0.,NFTPYR,159,INERR)
      ENDIF
C
C  Sep-15-1994 @ CSC-bsg Request 2-467 New Non-FTP Effects message and the
C  elimination of error message 157
C
      CALL QUITER(0.,NFTPYR,158,INERR)
C
C  Sep-29-1994 @ CSC-bsg Request 467 modified default values for the in-use
C  driving behavior record and added Methane as a user supplied rate.  Also
C  checks THC is larger than Methane
C  Error check REDUCT, NFTPV, and PHASSC
C
      IF((REDUCT(1).LT.0.0).OR.(REDUCT(1).GT.1.0).OR.
     * (REDUCT(2).LT.0.0).OR.(REDUCT(2).GT.1.0).OR.
     * (REDUCT(3).LT.0.0).OR.(REDUCT(3).GT.1.0).OR.
     * (NFTPV(1).LT.NFTPV(4)).OR.
     * (NFTPV(1).LT.0.0).OR.(NFTPV(2).LT.0.0).OR.
     * (NFTPV(3).LT.0.0).OR.(NFTPV(4).LT.0.0).OR.
     * (PHASSC(1).LT.0.0).OR.(PHASSC(2).LT.0.0).OR.
     * (PHASSC(3).LT.0.0).OR.(PHASSC(4).LT.0.0).OR.
     * (PHASSC(5).LT.0.0).OR.(PHASSC(6).LT.0.0).OR.
     * (PHASSC(1).GT.1.0).OR.(PHASSC(2).GT.1.0).OR.
     * (PHASSC(3).GT.1.0).OR.(PHASSC(4).GT.1.0).OR.
     * (PHASSC(5).GT.1.0).OR.(PHASSC(6).GT.1.0).OR.
     * (PHASSC(6).LT.PHASSC(5)).OR.(PHASSC(5).LT.PHASSC(4)).OR.
     * (PHASSC(4).LT.PHASSC(3)).OR.(PHASSC(3).LT.PHASSC(2)).OR.
     * (PHASSC(2).LT.PHASSC(1))) GOTO 99
      RETURN
C
C  Error or EOF on any attempted read => take alternate
C  return 1 => run aborts.
C
   99 CALL QUITER(0.,0,160,INERR)
C
      RETURN1
      END
