      SUBROUTINE OUTTAM 
C 
C  OUTTAM prints the user supplied replacement tampering rates and the derived 
C  TNKS = MISF - NCKS case in ID (instead of IDUSER) order. 
C 
C  Called by OUTHD3. 
C 
C  Changes: (Last change first) 
C 
C  08 Jun 00 AIR Task 03: Removed nonblock DO constructs. 
C  01 May 00 @EPA- BG    Removed use of DATA statement initializations       
C  14 Dec 98 @DynTel-ZK 2-000 Explicit Typing 
C  11 Aug 98 @DynTel-MLA 2-000  Changed HDG and HDD names to conform to 
C            new naming conventions. 
C 
C  Input on call: 
C 
C    common blocks: 
C    /FLAGS2/ IMFLAG 
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
C  11 Jun 97 DynTel@MLA 2-622  The code was modified to print output for 
C            all of the new Heavy Duty and BUS vehicle types that were 
C            introduced in Mobile6. 
C 
C  6 March 98 Dyntel@MMS 2-659 statements processing OUTFMT=4 were removed. 
C 
C  20 April 98 Dyntel@BAG 2-661 Deleted dependencies on OUTFMT. 
C 
C  23 Jul 98 @DynTel-MLA 2-000  Changed LDGT names to conform to 
C            certification naming conventions. 
C 
      IMPLICIT NONE 
      INCLUDE 'FLAGS2.I' 
      INCLUDE 'FLAGS3.I' 
      INCLUDE 'IOUCOM.I' 
      INCLUDE 'TAMEQ1.I' 
C 
      INTEGER ILOOP, IM, IGD, ISC, ITY, ID, ITER, ICH, IVTAM  
      CHARACTER(LEN=4) :: AZML='ZML ' 
C 
      CHARACTER(LEN=4) :: ADET='DET ' 
C 
      CHARACTER(LEN=8), DIMENSION(2,9) :: DISABL=RESHAPE((/ 
     *            'Air Pump','        ', 
     *            'Catalyst','        ', 
     *            'Fuel Inl','et      ', 
     *            'Other Mi','sfueling', 
     *            'EGR Syst','em      ', 
     *            'Evap Can','ister   ', 
     *            'PCV Syst','em      ', 
     *            'Cap     ','        ', 
     *            'All Misf','ueling  '/),(/2,9/)) 
C 
      CHARACTER(LEN=4), DIMENSION(3,2) :: IMTYPE=RESHAPE((/ 
     * 'Non-','I/M ','Case','  I/','M Ca','se  '/),(/3,2/)) 
C 
      CHARACTER(LEN=4), DIMENSION(3,3) :: IGTYPE=RESHAPE((/ 
     *                'MYG1',':pre','1981', 
     *                'MYG2',':198','1-83', 
     *                'MYG2',': 19','84+ '/),(/3,3/)) 
C 
      CHARACTER(LEN=4), DIMENSION(4) :: USC=(/ 
     *                               '_   ','__  ','___ ','____'/) 
C 
      WRITE(IOUREP,200) 
  200 FORMAT('0', 
     *'Replacement Tampering and Misfueling Rates Input by User:') 
C 
      ILOOP=2 
      IF(IMFLAG.EQ.1) ILOOP=1 
      DO 30 IM=1,ILOOP 
      WRITE(IOUREP,210)  
     * (IMTYPE(ITY,IM),ITY=1,3),(USC(4),ISC=1,3) 
  210 FORMAT('0',T54,3A4/'+',T54,3A4) 
C 
      DO 20 IGD=1,3 
C 
C  Output Format 3 or 5: wide (112 column) descriptive layout: 
C 
C  Group 1: LDGV - - - - - - - - - - - - - - - - - - - - - - - - - - - 
C 
      WRITE(IOUREP,230) 
     *  (IGTYPE(ITY,IGD),ITY=1,3),(USC(4),ISC=1,3),USC(4),USC(1) 
C 
  230 FORMAT('0',T54,3A4/'+',T54,3A4/ 
     *  '0',T29,'LDGV',/T29,2A4) 
C 
      WRITE(IOUREP,232) AZML,ADET, 
     *  (USC(4),ISC=1,4),USC(2),USC(4),USC(2),USC(4),USC(2), 
     *  ((DISABL(ICH,ID),ICH=1,2),TAMZML(ID,1,IGD,IM), 
     *    TAMDR(ID,1,IGD,IM),ID=1,9) 
C 
  232 FORMAT('0',T6,'Component',T19,'|',5X,A4,6X,A4/ 
     *  '+',T2,5A4,2X,2A4,' ',2A4/ 
     *  (' ',2A8,T19,'|',T21,2X,F7.4,2X,F8.4)) 
C 
C Group 2: LDGT1, 2, 3, 4  - - - - - - - - - - - - - - - - - - - -  
C 
      WRITE(IOUREP,234) 
     *  (USC(4),USC(2),ISC=1,4) 
C 
  234 FORMAT( 
     *  '0',T29,'LDGT1',T48,'LDGT2',T67,'LDGT3',T86,'LDGT4'/ 
     *  '+',T29,2A4,T48,2A4,T67,2A4,T86,2A4) 
C 
      WRITE(IOUREP,236) 
     *  (AZML,ADET,ITER=1,4), 
     *  (USC(4),ISC=1,4),USC(2),(USC(4),USC(2),USC(4),USC(3),ISC=1,4), 
     *  ((DISABL(ICH,ID),ICH=1,2),(TAMZML(ID,IVTAM,IGD,IM), 
     *           TAMDR(ID,IVTAM,IGD,IM),IVTAM=2,5),ID=1,9) 
C 
  236 FORMAT('0',T6,'Component',T19,'|',4(5X,A4,6X,A4)/ 
     *  '+',T2,5A4,2X,4(2X,2A4,' ',2A4)/ 
     *  (' ',2A8,T19,'|',T21,4(2X,F7.4,2X,F8.4))) 
C 
C Group 3: HDGV2B-HDGV5  - - - - - - - - - - - - - - - - - - - - - - -  
C 
      WRITE(IOUREP,238) 
     *  USC(4),USC(2),(USC(4),USC(1),ISC=1,3) 
C 
  238 FORMAT( 
     *  '0',T29,'HDGV2B',T48,'HDGV3',T67,'HDGV4',T86,'HDGV5'/ 
     *  '+',T29,2A4,T48,2A4,T67,2A4,T86,2A4) 
C 
      WRITE(IOUREP,236) 
     *  (AZML,ADET,ITER=1,4), 
     *  (USC(4),ISC=1,4),USC(2),(USC(4),USC(2),USC(4),USC(3),ISC=1,4), 
     *  ((DISABL(ICH,ID),ICH=1,2),(TAMZML(ID,IVTAM,IGD,IM), 
     *           TAMDR(ID,IVTAM,IGD,IM),IVTAM=6,9),ID=1,9) 
C 
C Group 4: HDGV6-HDGV8B  - - - - - - - - - - - - - - - - - - - - - - -  
C 
      WRITE(IOUREP,240) 
     *  (USC(4),USC(1),ISC=1,4) 
C 
  240 FORMAT( 
     *  '0',T29,'HDGV6',T48,'HDGV7',T67,'HDGV8A',T86,'HDGV8B'/ 
     *  '+',T29,2A4,T48,2A4,T67,2A4,T86,2A4) 
C 
      WRITE(IOUREP,236) 
     *  (AZML,ADET,ITER=1,4), 
     *  (USC(4),ISC=1,4),USC(2),(USC(4),USC(2),USC(4),USC(3),ISC=1,4), 
     *  ((DISABL(ICH,ID),ICH=1,2),(TAMZML(ID,IVTAM,IGD,IM), 
     *           TAMDR(ID,IVTAM,IGD,IM),IVTAM=10,13),ID=1,9) 
C 
C  Group 5: GAS BUS  - - - - - - - - - - - - - - - - - - - - - - - - - 
C 
      WRITE(IOUREP,244) 
     *  USC(4),USC(3) 
C 
  244 FORMAT( 
     *  '0',T25,'GAS BUS',/T25,2A4) 
C 
      WRITE(IOUREP,246) AZML,ADET, 
     *  (USC(4),ISC=1,4),USC(2),USC(4),USC(2),USC(4),USC(2), 
     *  ((DISABL(ICH,ID),ICH=1,2),TAMZML(ID,14,IGD,IM), 
     *    TAMDR(ID,14,IGD,IM),ID=1,9) 
C 
  246 FORMAT('0',T6,'Component',T19,'|',5X,A4,6X,A4/ 
     *  '+',T2,5A4,2X,2A4,' ',2A4/ 
     *  (' ',2A8,T19,'|',T21,2X,F7.4,2X,F8.4)) 
  20  CONTINUE 
  30  CONTINUE 
C 
      RETURN 
      END 
