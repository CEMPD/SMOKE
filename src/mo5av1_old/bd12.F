C Version @(#)$Id$ $Source$ $Date$ 

      BLOCK DATA BD12
C
C  BLOCK DATA Subprogram 12: /SAVE01/,/SAVE02/,/SAVE03/,/SAVE04/,/SAVE05/,
C                            /SAVE06/
C
C  BLOCK DATA 12 defines the save areas and initializes flag JCALL.
C
C  Common block dictionary:
C
C   Name   Type               Description
C  ------  ----  --------------------------------------------
C  /SAVE01/:
C
C  JCALL    I    MOBILE flag used to control direction of storing loops:
C                1 = copy default values of input data areas into save areas
C                2 = copy save area values back into corresponding data areas
C  VA001    R    /IDLEQ1/ ZMLIDL
C  VA002    R             DRIDL
C
C  /SAVE02/:
C
C  VA003    R    /ICR4VA/ CRD4VA
C  VA004    R    /ICR4VB/ CRD4VA
C
C  /SAVE03/:
C
C  VA005    R    /IM12HC/ CR12HC
C  VA006    R    /IM12CO/ CR12CO
C
C  /SAVE04/:
C
C  VA007    R    /MYRSAV/ AMAR
C  VA008    R             JULMYR
C
C  /SAVE05/:
C
C  VA009    R    /MYRSAV/ TAMZML
C  VA010    R             TAMDR
C
C  /SAVE06/:
C
C  VA011    R    /VMXCOM/ VMTMIX
C  VA012    R    /ETPPAR/ ETP
C  IVA013   I    /ETPPAR/ IETPYR
C  VA014    R    /ETPPAR/ ETPPIV
C
C  Notes:
C
C  BD 12 was added for Version 01's conversion of MOBILE4 to a subroutine.
C
C  JCALL is used to initialize the save areas on the first call to MOBILE and
C  restore the default values to the input data areas at the outset of every
C  call thereafter.
C
C
      COMMON /SAVE01/ JCALL,VA001(15,3,8,2),VA002(15,3,8,2)
      COMMON /SAVE02/ VA003(24,3,16),VA004(24,3,16)
      COMMON /SAVE03/ VA005(19,20,5,2),VA006(19,20,5,2)
      COMMON /SAVE04/ VA007(25,8),VA008(25,8)
      COMMON /SAVE05/ VA009(9,4,3,2),VA010(9,4,3,2)
      COMMON /SAVE06/ VA011(8),VA012(7,3),IVA013(4),VA014(4)
C
      DATA JCALL/1/
C
      END
