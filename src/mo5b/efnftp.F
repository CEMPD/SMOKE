C Version "@(#)$Id$ $Source$ $Date$ 

      FUNCTION EFNFTP(MY,IV,IP)
C
C  EFNFTP gets the non-FTP effect which depends on several factors:
C  model year, vehicle type, pollutant, speed, type of hydrocarbon, NONFTP,
C  REDUCT, NFTPV, PHASSC, and NFTPYR which holds the first model year where
C  there are no more non-FTP effects.  NFTPYR gets its value depending on how
C  the user sets the PRTFLG, when PRTFLG=5 there are non-FTP effects all years
C  (NONFTP=2), and when the PRTFLG=6 the user sets the NFTPYR and the program
C  set NONFTP=3
C
C  Called by HCCALX and EFCALX
C
C  Calls PCTLEV, TOGFID, and VOCFID
C
C  Input on call:
C    parameter list: MY,IV,IP
C
C    common blocks:
C    /FLAGS1/ OXYFLG
C    /FLAGS4/ NMHFLG
C    /NFTPDT/ NONFTP,NFTPYR,PHASSC,REDUCT,NFTPV
C    /SCENE1/ SPD
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  REDLOC   R    Local reductions for pollutants depend on NFTPYR
C  RTRNVL   R    Local return value for function
C  FRNFTP   R    The interpolation fraction of non-FTP effects
C  IYEAR    I    Index into the PHASSC array
C
C
C  Notes:
C
C  Aug-18-1994 @ CSC-bsg Request 462 EFNFTP was added for MOBILE5b to get the 
C  value of non-FTP offset.  Uses the variable NFTPYR from the user when the 
C  NONFTP Flag is 3, NONFTP is set in GETVER/HOR and depends on PRTFLG as read 
C  from the user.  Values for non-FTP are preset here as given by the EPA
C  Non-Methane non-FTP emission value is here as a place holder to 
C  demonstrate the relationship between Non-Methane offset and THC.
C  Aug-24-1994 @ CSC-bsg Request 462 revised non-FTP offset numbers:
C  total HC, methane HC, CO and NOx.
C  Sep-20-1994 @ CSC-bsg Request 467 non-FTP offset number are stored in
C  COMMON.  In-use driving now depends on REDUCT, and PHASSC as well as
C  its original input parameters: vehicle speed etc...
C  Sep-30-1994 @ CSC-bsg Request 467 Discovered that the VOC and TOG were
C  not formulated correctly, they are now offset for Methane.  Changed the
C  IYEAR definitions and the common variable PHASSC to go from 0 to 6 instead
C  of 1 to 7.  Fixed formulation of reductions for new certification period.
C  Oct-6-1994 @ CSC-bsg Request 467 Changed reduction formula to be correct
C  Jan-11-1995 @CSC-zp  Request 0-448 Changed place where the value of the 
C  IYEAR was set.
C  Jun-20-1996 @ DynTel-yc Request 2-624  Base non-FTP emission factors in 
C                comment line was removed.
C
      INCLUDE 'FLAGS1.I'
      INCLUDE 'FLAGS4.I'
      INCLUDE 'NFTPDT.I'
      INCLUDE 'SCENE1.I'
C
      RTRNVL = 0.0
      FRNFTP = 1.0
C
C  Sep-30-1994 @ CSC-bsg Request 467 Changed the IYEAR definitions and the
C  common variable PHASSC to go from 0 to 6 instead of 1 to 7.
C  Set the IYEAR pointer
C
      IYEAR = MY - NFTPYR + 1
      IF(IYEAR.LE.0) THEN
        IYEAR = 0
      ELSEIF(IYEAR.GT.6) THEN
        IYEAR = 6
      ENDIF
C
      IF(NONFTP.EQ.1) GOTO 99
C
C  Heavy duty vehicles HDGV & HDDV (IV=4/7) do not receive a non-FTP offset.
C  Motorcycles MC (IV=8) also do not receive a non-FTP offset.
C  Light duty diesel cars and trucks LDDV & LDDT (IV=5/6) only get non-FTP
C  offsets for NOx emissions (IP=3).
C
      IF(IV.EQ.4) GOTO 99
      IF(IV.EQ.7) GOTO 99
      IF(IV.EQ.8) GOTO 99
      IF(IV.GE.5.AND.IP.NE.3) GOTO 99
C
C  Interpolate the fraction of non-FTP effects that are applicable
C
      IF(SPD(IV).LT.19.6) THEN
        FRNFTP = (SPD(IV) - 2.5)/17.1
      ELSEIF(SPD(IV).GT.55) THEN
        FRNFTP = (65.0 - SPD(IV))/10.0
      ELSE
        FRNFTP = 1.0
      ENDIF
C
C  Apply the non-FTP offset to the complement of the zero emitting vehicle
C  amount (1-ZEV).  The ZEV fraction, is model year (MY) and vehicle (IV) 
C  dependent.  PCTLEV is used to obtain the fraction.
C
      FRNFTP = FRNFTP * (1.0 - PCTLEV(MY,IV,8))
C
C  Aug-24-1994 @ CSC-bsg Request 462 revised non-FTP offset numbers: 
C  total HC, methane HC, CO and NOx.
C  Sep-20-1994 @ CSC-bsg Request 467 non-FTP offset number are stored in
C  COMMON.  In-use driving now depends on REDUCT, and PHASSC as well as
C  its original input parameters: vehicle speed etc...
C  Sep-30-1994 @ CSC-bsg Request 467 Discovered that the VOC and TOG were
C  not formulated correctly, they are now offset for Methane
C  Fixed formulation of reductions for new certification period.
C
      IF(IP.EQ.1) THEN
        IF(NMHFLG.EQ.1) THEN
          RTRNVL = NFTPV(1)
        ELSEIF(NMHFLG.EQ.2) THEN
          RTRNVL = NFTPV(1) - NFTPV(4)
        ELSEIF(NMHFLG.EQ.3) THEN
          RTRNVL = (NFTPV(1) - NFTPV(4)) * VOCFID(MY,IV)
        ELSEIF(NMHFLG.EQ.4) THEN
          RTRNVL = ((NFTPV(1) - NFTPV(4)) * TOGFID(MY,IV,OXYFLG)) + 
     +               NFTPV(4)
        ELSE
          RTRNVL = (NFTPV(1) - NFTPV(4)) * TOGFID(MY,IV,OXYFLG)
        ENDIF
        REDLOC = REDUCT(1)
      ELSEIF(IP.EQ.2) THEN
        RTRNVL = NFTPV(2)
        REDLOC = REDUCT(2)
      ELSE
        RTRNVL = NFTPV(3)
        REDLOC = REDUCT(3)
      ENDIF
C
C  Sep-30-1994 @ CSC-bsg Request 467 Discovered that the VOC and TOG were
C  not formulated correctly, they are now offset for Methane.  Changed the
C  IYEAR definitions and the common variable PHASSC to go from 0 to 6 instead
C  of 1 to 7.  Fixed formulation of reductions for new certification period.
C
      IF(IYEAR.EQ.0) REDLOC = 0.0
C
C  Oct-6-1994 @ CSC-bsg Request 467 Changed reduction formula to be correct
C  Calculate effect of phase-in of in-use driving by using IYEAR pointer
C  and reduction effect REDLOC
C
   99 EFNFTP = RTRNVL * FRNFTP * (1.0 - (PHASSC(IYEAR) * REDLOC))
      RETURN
      END
