C Version "@(#)$Id$ $Source$ $Date$ 

      SUBROUTINE OUTLAP
C
C  OUTLAP echoes the user supplied values for the fields of the local area
C  parameters ("LAP") record.  This record is echoed as the last echo of 1-time
C  data, if LOCFLG = 2, or after the 1st scenario record echo of each scenario,
C  if LOCFLG = 1.  The echo record length depends on OUTFMT: OUTFMT= 3 or 5
C  allows records up to 112 characters long, while OUTFMT = 4 limits records to
C  a maximum of 80 characters (the LAP record is only echoed by the descriptive
C  formats, 3 and 4).
C
C  Called by OUTHD3, OUTHD4, OUTDT3 and OUTDT4.
C
C  Input on call:
C
C    common blocks:
C    /CITPAR/ SCNAME
C    /CITRV1/ RVPBAS,RVPIUS,IUSESY
C    /FLAGS3/ OUTFMT
C    /IOUCOM/ IOUREP
C    /RFORM1/ IAFCNT,ASTMCL
C    /TEMPS/  TEMMIN,TEMMAX
C
C  Output on return:
C
C  None.
C
C  Notes:
C
C  Echo fields are vertically aligned by the colons preceding the values.
C  Note that the T formats reference a 113 and 81 column format, respectively.
C  Each output column 1 is reserved for carriage control, with the remainder
C  being printed.
C
C  The OUTFMT = 5 case added for In-House MOBILE4 Version 02 (BYMY) is
C  included in the 112 column format case of OUTLAP.
C
C
      INCLUDE 'CITPAR.I'
      INCLUDE 'CITRV1.I'
      INCLUDE 'FLAGS3.I'
      INCLUDE 'IOUCOM.I'
      INCLUDE 'RFORM1.I'
      INCLUDE 'TEMPS.I'
C
C  112 column format colons aligned at T13, T36, T60 and T91.
C
      IF(OUTFMT.EQ.3.OR.OUTFMT.EQ.5) WRITE(IOUREP,200)
     *  SCNAME,TEMMIN,TEMMAX,RVPBAS,RVPIUS,IUSESY
  200 FORMAT('0',T2,4A4,
     *           T48,'Minimum Temp: ',F4.0,' (F)',
     *           T79,'Maximum Temp: ',F4.0,' (F)'/
     *       ' ',T24,'Period 1 RVP: ',F4.1,T48,'Period 2 RVP: ',F4.1,
     *           T74,'Period 2 Start Yr: ',I4)
CC   *           T75,'Ethanol Fraction:   ',F4.2)
CC   *       ' ',T17,'Base fuel indicator: ',I3)
C
C  80 column format colons aligned at T34, T58 and T76.
C
      IF(OUTFMT.EQ.4) WRITE(IOUREP,210)
     *  SCNAME,TEMMIN,TEMMAX,RVPBAS,RVPIUS,IUSESY
  210 FORMAT('0',T2,4A4/
CC   *           T46,'Ethanol Frac:  ',F4.2/
     *       ' ',T22,'Minimum Temp:',F4.0,' (F)',
     *           T46,'Maximum Temp:',F4.0,' (F)'/
     *       ' ',T22,'Period 1 RVP: ',F4.1,T46,'Period 2 RVP: ',F4.1,
     *           T65,'Period 2 Yr: ',I4)
CC   *       ' ',T15,'Base fuel indicator: ',I3)
C
      RETURN
      END
