C Data is initialized in BD_IM1 and REINIT 
C 
C  03 Dec 01 @EPA-BAG        Remove AGEIM, VMTIM, and associated variables/ see 
C                            also IMPROG, FINDIM, AND INITIM 
C 25 Aug 00 AIR Task 08 Increased I/M related arrays from 5 to 7 to 
C           allow for EVAP I/M 
C 17 Jul 00 @EPA-elg  Added Array ASMRAT 
C 21 Sep 98 @DynTel-RJD 2-663 LDGIMC contains data used the 1981-1993 light 
C           duty gasoline vehicle and truck IM calculations 
C 
      REAL ASEBZML,IDR240,IDRID1,IDRID2, 
     *     IMCRED,WAVRDC,ASMCP1,ASMCP2,IMPEFF,ASMRAT 
      INTEGER IMSTRICT 
      COMMON /LDGIMC/ ASEBZML(13,3,5,3),IDR240(25,2,3,4,7), 
     *             IDRID1(3,2,3),IDRID2(3,2,3),IMCRED(3,25,5,2,3), 
     *             WAVRDC(7),ASMCP1(4,3,6), 
     *             ASMCP2(4,3,6),IMPEFF(7,3),IMSTRICT(7), 
     *             ASMRAT(3,18,25) 
 
