C Version "@(#)$Id$ $Source$ $Date$ 

      SUBROUTINE OUTTAM
C
C  OUTTAM prints the user supplied replacement tampering rates and the derived
C  TNKS = MISF - NCKS case in ID (instead of IDUSER) order.
C
C  Called by OUTHD3 and OUTHD4.
C
C  Input on call:
C
C    common blocks:
C    /FLAGS2/ IMFLAG
C    /FLAGS3/ OUTFMT
C    /IOUCOM/ IOUREP
C    /TAMEQ1/ TAMZML,TAMDR
C
C  Local array subscripts:
C
C  DISABL(2,9)  -  DISABL ( ICH , ID )
C  IGTYPE(3,3)  -  IGTYPE ( ITY , IGD )
C  IMTYPE(3,2)  -  IMTYPE ( ITY , IM )
C  USC(4)       -  USC ( ISC )
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  ADET     C*4  column heading for deterioration rates
C  AZML     C*4  column heading for zero mile levels
C  DISABL   C*8  names of disablement categories
C  IGTYPE   C*4  tampering rates table headings, by myg & vmt, as applies
C  IMTYPE   C*4  non-I/M and I/M tampering rates table headings
C  USC      C*4  underscore strings for table and column headings
C  ILOOP      I  either a one or two tells IM how many to loop.
C  IM         I  loop counting variable
C
C  Notes:
C
C  MOBILE4 prints two evaporative tampering categories (see ID definition).
C  and the separate MYG2 >50K set.  OUTTAM handles the new IGD dimension and
C  the new misfueling algorithm.  OUTTAM was commented out in Release MOBILE4
C
C  The calls to OUTTAM were commented out in Release MOBILE4 (the echo table
C  is large), but the calls are restored in I00 for diagnostic purposes.
C
C  The OUTFMT = 5 case added for In-House MOBILE4 Version 02 (BYMY) is
C  included in the 112 column format case of OUTTAM.
C
C  OUTTAM was modified in draft MOBILE4.1v06 to add the third MYG and remove
C  the over 50K rates.
C
C
      INCLUDE 'FLAGS2.I'
      INCLUDE 'FLAGS3.I'
      INCLUDE 'IOUCOM.I'
      INCLUDE 'TAMEQ1.I'
C
      CHARACTER*8 DISABL(2,9)
      CHARACTER*4 IMTYPE(3,2),IGTYPE(3,3),AZML,ADET,USC(4)
C
      DATA IMTYPE/'Non-','I/M ','Case','  I/','M Ca','se  '/
      DATA IGTYPE/'MYG1',':pre','1981',
     *            'MYG2',':198','1-83',
     *            'MYG2',': 19','84+ '/
      DATA AZML/'ZML '/,ADET/'DET '/
      DATA USC/'_   ','__  ','___ ','____'/
      DATA DISABL/'Air Pump','        ',
     *            'Catalyst','        ',
     *            'Fuel Inl','et      ',
     *            'Other Mi','sfueling',
     *            'EGR Syst','em      ',
     *            'Evap Can','ister   ',
     *            'PCV Syst','em      ',
     *            'Cap     ','        ',
     *            'All Misf','ueling  '/
C
      WRITE(IOUREP,200)
  200 FORMAT('0',
     *'Replacement Tampering and Misfueling Rates Input by User:')
C
      ILOOP=2
      IF(IMFLAG.EQ.1) ILOOP=1
      DO 20 IM=1,ILOOP
      IF(OUTFMT.EQ.3.OR.OUTFMT.EQ.5) WRITE(IOUREP,210)
     *  (IMTYPE(ITY,IM),ITY=1,3),(USC(4),ISC=1,3)
  210 FORMAT('0',T54,3A4/'+',T54,3A4)
      IF(OUTFMT.EQ.4) WRITE(IOUREP,220)
     *  (IMTYPE(ITY,IM),ITY=1,3),(USC(4),ISC=1,3)
  220 FORMAT('0',T45,3A4/'+',T45,3A4)
C
      DO 20 IGD=1,3
      IF(OUTFMT.EQ.4) GOTO 10
C
C  Output Format 3 or 5: wide (112 column) descriptive layout:
C
      WRITE(IOUREP,230)
     *  (IGTYPE(ITY,IGD),ITY=1,3),
     *  (USC(4),ISC=1,5),USC(1),USC(4),USC(1),USC(4)
C
  230 FORMAT('0',T54,3A4/'+',T54,3A4/
     *  '0',T29,'LDGV',T48,'LDGT1',T67,'LDGT2',T86,'HDGV'/
     *  '+',T29,A4,T48,2A4,T67,2A4,T86,A4)
C
      WRITE(IOUREP,240)
     *  (AZML,ADET,ITER=1,4),
     *  (USC(4),ISC=1,4),USC(2),(USC(4),USC(2),USC(4),USC(3),ISC=1,4),
     *  ((DISABL(ICH,ID),ICH=1,2),(TAMZML(ID,IVTAM,IGD,IM),
     *           TAMDR(ID,IVTAM,IGD,IM),IVTAM=1,4),ID=1,9)
C
  240 FORMAT('0',T6,'Component',T19,'|',4(5X,A4,6X,A4)/
     *  '+',T2,5A4,4(2X,2A4,' ',2A4)/
     *  (' ',2A8,T19,'|',T21,4(2X,F7.4,2X,F8.4)))
      GOTO 20
C
C  Output Format 4: narrow (80 column) descriptive format:
C
   10 WRITE(IOUREP,250)
     *  (IGTYPE(ITY,IGD),ITY=1,3),
     *  (USC(4),ISC=1,5),USC(1),USC(4),USC(1),USC(4)
C
  250 FORMAT('0',T45,3A4/'+',T45,3A4/
     *  '0',T26,'LDGV',T41,'LDGT1',T56,'LDGT2',T71,'HDGV'/
     *  '+',T26,A4,T41,2A4,T56,2A4,T71,A4)
      WRITE(IOUREP,260)
     *  (AZML,ADET,ITER=1,4),
     *  (USC(4),ISC=1,4),USC(2),(USC(4),USC(4),USC(3),ISC=1,4),
     *  ((DISABL(ICH,ID),ICH=1,2),(TAMZML(ID,IVTAM,IGD,IM),
     *           TAMDR(ID,IVTAM,IGD,IM),IVTAM=1,4),ID=1,9)
C
  260 FORMAT('0',T6,'Component',T19,'|',4(3X,A4,4X,A4)/
     *  '+',T2,5A4,4(A4,'__ ',2A4)/
     *  (' ',2A8,T19,'|',T21,4(F7.4,F8.4)))
C
   20 CONTINUE
C
      RETURN
      END
