C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE SETEGS
C
C  SETEGS sets the effects group sizes by combining subsets of the
C  overlap category sizes and weighting by the technology group sizes.
C
C  Called by EFFGRP.
C
C  Input on call:
C
C    common blocks:
C    /LOOKUP/ IQG,IPG
C    /SIZCAL/ TGSUSE,CSAE
C
C  Output on return:
C
C    common block: /EGSCAL/ EGS
C
C  Local array subscripts:
C
C  INDXEG(6,7,2)  -  INDXEG ( IC, IEG, IPG )
C  INDXPG(2,2)    -  INDXPG ( ILU, IPG )
C  NEGC(7,2)      -  NEGC ( IEG, IPG )
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  INDXEG   I    secondary index array: pointers to category sizes
C                in IEG's size sum, given IPG
C  INDXPG   I    loop indices array: first and last nonzero effects group
C                sizes, given IPG
C  NC       I    loop test value: number of category sizes in IEG's size
C                equation, given IPG
C  NEGC     I    loop test values array: looks up NC
C
C  Notes:
C
C  None.
C
C
      DIMENSION INDXPG(2,2),NEGC(7,2),INDXEG(6,7,2)
C
      COMMON /EGSCAL/ AER(11,11,2,2),TGS(11,6,4),EGS(7,2)
      COMMON /LOOKUP/ IVTAM,IQG,IPG,JPGD,IHG,IGCSF
      COMMON /SIZCAL/ ATR(5,2),TGSUSE,CSIZE(12,2),CSAE(12)
C
      DATA INDXPG/1,4, 3,7/,NEGC/1,1,6,4,3*0, 2*0,3,2,1,3,2/
      DATA INDXEG/
     1  8,5*0, 8,5*0, 1,4,5,6,7,9, 2,3,10,11,2*0, 6*0, 6*0, 6*0,
     2  6*0, 6*0, 6,7,9,3*0, 10,11,4*0, 8,5*0, 1,4,5,3*0, 2,3,4*0/
C
C  Do not compute group sizes for cases known to be always 0 (5-7 for HC/CO
C  and 1-2 for NOX).
C
      IE1ST=INDXPG(1,IPG)
      IELST=INDXPG(2,IPG)
C
C  Equipment = air pump /catalyst => air pump disablement (air pump only
C                                    equipped) effects group increment = 0.0
C
      IF(IQG.EQ.2.AND.IPG.EQ.1) IE1ST=2
C
C  Converse air pump case:
C
C  Equipment = air pump only => only effects group 1 (ap only disabled) > 0.0
C
      IF(IQG.EQ.1.AND.IPG.EQ.1) IELST=1
C
C  For each effects group IEG, the size is the sum of the product of the
C  related category sizes times the technology group size for the technology
C  type IQG being evaluated on this pass.
C
      DO 20 IEG=IE1ST,IELST
C
C  Get the number of categories involved in group IEG.
C
      NC=NEGC(IEG,IPG)
C
C  Now do the sum.
C
      DO 10 IC=1,NC
      JCS=INDXEG(IC,IEG,IPG)
      EGS(IEG,IPG)=EGS(IEG,IPG)+TGSUSE*CSAE(JCS)
   10 CONTINUE
   20 CONTINUE
C
      RETURN
      END
