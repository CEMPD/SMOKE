C Version "@(#)$Id$ $Source$ $Date$ 

      SUBROUTINE OUTMTH
C
C  OUTMTH writes evaporative HC message on unit IOUREP.
C  ljn 3/5/93 MARK@HC=3
C
C  Called by OUTHD3 and OUTHD4.
C
C  Input on call:
C
C    common blocks:
C    /FLAGS4/ NMHFLG
C    /IOUCOM/ IOUREP
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
C
      INCLUDE 'FLAGS4.I'
      INCLUDE 'IOUCOM.I'
      INCLUDE 'FLAGS3.I'
C
C Added HC3 11 Jan 92, HC3NEW 18 Feb 93 - MLA
      CHARACTER*6  TOTAL
      CHARACTER*12 NONMTH
      CHARACTER*4  NMOG,TOG,VOC,REMARK(15)
      CHARACTER*4  HC3NEW(24)
C
      DATA TOTAL/'Total '/
C
      DATA NONMTH/'Non-methane '/
C
      DATA NMOG,TOG,VOC/'NMOG','TOG ','VOC '/
C
      DATA REMARK/
     *'HC e','miss','ion ','fact','ors ','incl','ude ',
     *'evap','orat','ive ','HC e','miss','ion ','fact','ors.'/
C
C Added variable HC3NEW - 18 Feb 93 - MLA
C
      DATA HC3NEW/
     *'HC e','miss','ion ','fact','ors ','incl','ude ','all ',
     *'evap','orat','ive ','HC e','miss','ion ','fact','ors,',
     *' exc','ept ','for ','refu','elin','g em','issi','ons.'/
C
C Changed all cases (NMHFLG 1, 2, 3,...) to print an extra message
C when HCFLAG is equal to 3. 11 Jan 93 - MLA.
C
      IF(NMHFLG.EQ.1) THEN
        IF((HCFLAG.EQ.3.AND.OUTFMT.EQ.3).OR.
     *     (HCFLAG.EQ.3.AND.OUTFMT.EQ.5)) THEN
          WRITE (IOUREP,201) TOTAL,HC3NEW
        ELSE
          WRITE (IOUREP,200) TOTAL,REMARK
        ENDIF
      ENDIF
C
  200 FORMAT('0',A6,15A4)
  201 FORMAT('0',A6,24A4)
C
      IF(NMHFLG.EQ.2) THEN
        IF((HCFLAG.EQ.3.AND.OUTFMT.EQ.3).OR.
     *     (HCFLAG.EQ.3.AND.OUTFMT.EQ.5)) THEN
          WRITE (IOUREP,211) NONMTH,HC3NEW
        ELSE
          WRITE (IOUREP,210) NONMTH,REMARK
        ENDIF
      ENDIF
C
  210 FORMAT('0',A12,15A4)
  211 FORMAT('0',A12,24A4)
C
      IF(NMHFLG.EQ.3) THEN
        IF((HCFLAG.EQ.3.AND.OUTFMT.EQ.3).OR.
     *     (HCFLAG.EQ.3.AND.OUTFMT.EQ.5)) THEN
          WRITE (IOUREP,221) VOC,HC3NEW
        ELSE
          WRITE (IOUREP,220) VOC,REMARK
        ENDIF
      ENDIF
C
  220 FORMAT('0',A4,15A4)
  221 FORMAT('0',A4,24A4)
C
      IF(NMHFLG.EQ.4) THEN
        IF((HCFLAG.EQ.3.AND.OUTFMT.EQ.3).OR.
     *     (HCFLAG.EQ.3.AND.OUTFMT.EQ.5)) THEN
          WRITE (IOUREP,221) TOG, HC3NEW
        ELSE
          WRITE (IOUREP,220) TOG,REMARK
        ENDIF
      ENDIF
C
      IF(NMHFLG.EQ.5) THEN
        IF((HCFLAG.EQ.3.AND.OUTFMT.EQ.3).OR.
     *     (HCFLAG.EQ.3.AND.OUTFMT.EQ.5)) THEN
          WRITE (IOUREP,231) NMOG,HC3NEW
        ELSE
          WRITE (IOUREP,230) NMOG,REMARK
        ENDIF
      ENDIF
C
  230 FORMAT('0',A4,' ',15A4)
  231 FORMAT('0',A4,' ',24A4)
C
C
      RETURN
      END
