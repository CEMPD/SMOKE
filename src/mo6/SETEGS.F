      SUBROUTINE SETEGS 
C 
C  SETEGS sets the effects group sizes by combining subsets of the 
C  overlap category sizes and weighting by the technology group sizes. 
C 
C  Called by EFFGRP. 
C 
C  Changes: (Last change first) 
C 
C  01 May 00 @EPA- BG    Removed use of DATA statement initializations       
C  07 Jan 1999 @DynTel-ZK 2-000 Explicit Typing 
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
      IMPLICIT NONE 
C 
      INCLUDE 'EGSCAL.I' 
      INCLUDE 'LOOKUP.I' 
      INCLUDE 'SIZCAL.I' 
C 
      INTEGER IE1ST,IELST,IEG,NC,IC,JCS 
C              
      INTEGER, DIMENSION(2,2) :: INDXPG=RESHAPE((/1,4, 3,7/),(/2,2/)) 
      INTEGER, DIMENSION(7,2) :: NEGC=RESHAPE((/ 
     *             1,1,6,4,0,0,0,0,0,3,2,1,3,2/),(/7,2/)) 
      INTEGER, DIMENSION(6,7,2) :: INDXEG= RESHAPE((/ 
     * 8,0,0,0,0,0,  8,0,0,0,0,0,  1, 4,5,6,7,9,  2,3,10,11,0,0, 
     * 0,0,0,0,0,0,  0,0,0,0,0,0,  0, 0,0,0,0,0,  0,0, 0, 0,0,0, 
     * 0,0,0,0,0,0,  6,7,9,0,0,0, 10,11,0,0,0,0,  8,0, 0, 0,0,0, 
     * 1,4,5,0,0,0,  2,3,0,0,0,0/), 
     * (/6,7,2/)) 
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
