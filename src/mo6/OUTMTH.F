      SUBROUTINE OUTMTH 
C 
C  OUTMTH writes evaporative HC message on unit IOUERR. 
C  ljn 3/5/93 MARK@HC=3 
C 
C  Called by OUTHD3 . 
C 
C  Changes: (Last change first) 
C 
C  01 May 00 @EPA- BG    Removed use of DATA statement initializations       
C  17 Apr 00 AIR Task 02: Removed spurious integer/real compare 
C  10 Nov 99 @DynTel-BG 1-006  Revised HCFLAG numbering options. 
C            HCFLAG =4 now does not have refueling included in evap. 
C            and evap. + exhaust 'composite' sums. Previously HCLFAG=3 
C            enabled the "no refueling" in these sums. 
C   3 Mar 99 @DynTel-MLA 2-659  Prints messages on IOUERR, not IOUREP. 
C  11 Dec 98 @DynTel-ZK 2-000 Explicit Typing  
C 
C  Input on call: 
C 
C    common blocks: 
C    /FLAGS4/ NMHFLG 
C    /IOUCOM/ IOUERR 
C 
C  Local array subscripts: 
C 
C  NONMTH(3)  -  NONMTH ( ICH ) 
C  REMARK(15) -  REMARK ( ICH ) 
C  TOTAL(2)   -  TOTAL ( ICH ) 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  NONMTH   C*4  character string 'Non-methane ' 
C  NMOG     C*4  character string 'NMOG' 
C  REMARK   C*4  character string comment that HC ef's include evap ef's 
C  TOG      C*4  character string 'TOG ' 
C  TOTAL    C*4  character string 'Total   ' 
C  VOC      C*4  character string 'VOC ' 
C 
C  Notes: 
C 
C  None. 
C 
      IMPLICIT NONE 
      INCLUDE 'FLAGS4.I' 
      INCLUDE 'IOUCOM.I' 
      INCLUDE 'FLAGS3.I' 
C 
C Added HC3 11 Jan 92, HC3NEW 18 Feb 93 - MLA 
C 
C 20 Apr 1998 Deleted code dependencies on OUTFMT. 
C 
      CHARACTER(LEN=6) :: TOTAL='Total ' 
      CHARACTER(LEN=12) :: NONMTH='Non-methane ' 
      CHARACTER(LEN=4) :: NMOG='NMOG' 
      CHARACTER(LEN=4) :: TOG='TOG ' 
      CHARACTER(LEN=4) :: VOC='VOC ' 
      CHARACTER(LEN=4), DIMENSION(15) :: REMARK=(/ 
     *'HC e','miss','ion ','fact','ors ','incl','ude ', 
     *'evap','orat','ive ','HC e','miss','ion ','fact','ors.'/) 
      CHARACTER(LEN=4), DIMENSION(24) :: HC3NEW=(/ 
     *'HC e','miss','ion ','fact','ors ','incl','ude ','all ', 
     *'evap','orat','ive ','HC e','miss','ion ','fact','ors,', 
     *' exc','ept ','for ','refu','elin','g em','issi','ons.'/) 
C 
C Changed all cases (NMHFLG 1, 2, 3,...) to print an extra message 
C when HCFLAG is equal to 3. 11 Jan 93 - MLA. 
C 
      IF(NMHFLG.EQ.1) THEN 
        IF(HCFLAG.EQ.4) THEN 
          WRITE (IOUERR,201) TOTAL,HC3NEW 
        ELSE 
          WRITE (IOUERR,200) TOTAL,REMARK 
        ENDIF 
      ENDIF 
C 
  200 FORMAT('0',A6,15A4) 
  201 FORMAT('0',A6,24A4) 
C 
      IF(NMHFLG.EQ.2) THEN 
        IF(HCFLAG.EQ.4) THEN 
          WRITE (IOUERR,211) NONMTH,HC3NEW 
        ELSE 
          WRITE (IOUERR,210) NONMTH,REMARK 
        ENDIF 
      ENDIF 
C 
  210 FORMAT('0',A12,15A4) 
  211 FORMAT('0',A12,24A4) 
C 
      IF(NMHFLG.EQ.3) THEN 
        IF(HCFLAG.EQ.4) THEN 
          WRITE (IOUERR,221) VOC,HC3NEW 
        ELSE 
          WRITE (IOUERR,220) VOC,REMARK 
        ENDIF 
      ENDIF 
C 
  220 FORMAT('0',A4,15A4) 
  221 FORMAT('0',A4,24A4) 
C 
      IF(NMHFLG.EQ.4) THEN 
        IF(HCFLAG.EQ.4) THEN 
          WRITE (IOUERR,221) TOG, HC3NEW 
        ELSE 
          WRITE (IOUERR,220) TOG,REMARK 
        ENDIF 
      ENDIF 
C 
      IF(NMHFLG.EQ.5) THEN 
        IF(HCFLAG.EQ.4) THEN 
          WRITE (IOUERR,231) NMOG,HC3NEW 
        ELSE 
          WRITE (IOUERR,230) NMOG,REMARK 
        ENDIF 
      ENDIF 
C 
  230 FORMAT('0',A4,' ',15A4) 
  231 FORMAT('0',A4,' ',24A4) 
C 
C 
      RETURN 
      END 
