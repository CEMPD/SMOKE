      SUBROUTINE VEHCOUNT(ICY,FRYA,FRYB) 
C 
C  VEHCOUNT computes the vehicle counts by vehicle class. 
C 
C  Called by REGMOD. 
C 
C  Calls. 
C 
C  Changes: (Last change first) 
C 
C  01 May 00 @EPA- BG    Removed use of DATA statement initialization 
C  12 Apr 00 EPA@djb New for MOBILE6.  Calculates vehicle class VMT 
C            from default mileage, age dist, diesel fractions, and 
C            vehicle counts or user inputs of those values. 
C 
C  Input on call: 
C 
C    Constants:  
C    MAXVEH from include file IVTYPE.I 
C 
C    argument list: ICY,FRYA,FRYB 
C 
C    common blocks: 
C    /EVAL/   MEVAL 
C    /MAXIMA/ MAXYRS 
C    /MYRCAL/ AGED, DAF 
C    /VMXCOM/ VCOUNT 
C     
C 
C  Output on return: 
C 
C    argument list: ICY,FRYA,FRYB 
C 
C    common blocks: 
C    /VMXCOM/ VBCLASS 
C 
C  Local array subscripts: 
C 
C    AGEMILE(30) - AGEMILE(IVD) 
C    RATIO(16)   - RATIO(IVC) 
C    VCLASS(16)  - VCLASS(IVC) 
C    VMTBC(30)   - VMTBC(IVD) 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  AGEMILE  R    average mileage accumulation by vehicle class ( AGED * TFMYM ) 
C  FRYA     R    fraction of a year that is remaining until the next July 1 
C  FRYB     R    fraction of a year that has elapsed since the last is July 1 
C  I81      I    number of whole years that have elapsed since July 1 1981 
C  ICY      I    calendar year 
C  IV       I    standard veh class groupings (28) 
C  IVC      I    combined gas/diesel veh class groupings (16) 
C  IVD      I    complete veh class groupings (30) 
C  IVX      I    maps IV to IVD 
C  JDX      I    vehicle age index; JDX=25-IDX+1 
C  RATIO    R    ratio of the user supplied VMT compared to default values 
C  T81      R    total time that has elapsed in years since July 1 1981 
C  TOTCNT   R    total number of vehicles in this calendar year 
C  TOTMILE  R    sum of the product of the annual mileage and vehicle count 
C  VCLASS   R    interpolated vehicle count for the current year ICY and 
C                and month MEVAL (16) 
C  VMTBC    R    default VMT distribution 
C 
C 
C  Notes: 
C 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'EVAL.I' 
      INCLUDE 'FLAGS1.I' 
      INCLUDE 'IVTYPE.I' 
      INCLUDE 'MAXIMA.I' 
      INCLUDE 'MYRCAL.I' 
      INCLUDE 'VMXCOM.I' 
C 
      REAL  FRYA 
      REAL  FRYB 
      REAL  T81 
      REAL  TOTCNT 
      REAL  TOTMILE 
      REAL, DIMENSION(16) :: RATIO 
      REAL, DIMENSION(30) :: AGEMILE 
      REAL, DIMENSION(16) :: VCLASS 
      REAL, DIMENSION(30) :: VMTBC 
C 
      INTEGER   I81 
      INTEGER   ICY 
      INTEGER   IV 
      INTEGER   IVC 
      INTEGER   IVD 
      INTEGER   IVX 
      INTEGER   JDX 
C 
      VCLASS  = 0.0 
      AGEMILE = 0.0 
      VMTBC = 0.0 
      RATIO = 0.0 
C 
C  Determine the array pointer value for the VCOUNT array. 
C 
      IF(MEVAL.GE.7) THEN  
       I81=ICY-1981 
      ELSE IF (MEVAL.LT.7) THEN 
       I81=ICY-1982 
      END IF 
      T81=ICY-1981+((MEVAL-7)/12.0) 
C 
C IVC indicates the composite (gas/disel) vehicle class. 
C For all IV that are pre-1982 ICY, use the 1982 VCOUNT data. 
C For all IV that are post-2020 ICY, use the 2020 VCOUNT data. 
C ICY indexes into VCOUNT, ICY requires to be offset by "1981" 
C the start of the calander years + 1  = "1982" 
C 
      DO 25 IVC=1,16 
       IF(T81.LE.1.0) THEN 
          VCLASS(IVC)=VCOUNT(1,IVC) 
       ELSEIF(T81.LE.38) THEN 
          VCLASS(IVC)=FRYA*VCOUNT(I81,IVC)+FRYB*VCOUNT(I81+1,IVC) 
       ELSEIF(T81.GT.38) THEN 
          VCLASS(IVC)=VCOUNT(39,IVC) 
       ENDIF 
   25 CONTINUE 
C 
C  Split the composite vehicle counts into their gas/diesel portions 
C  based on the sum of the dieselized age distributions from REGMOD. 
C 
      TOTCNT=0.0 
      DO IVC=1,16 
        TOTCNT=TOTCNT+VCLASS(IVC) 
        VALL(IVC)=VCLASS(IVC)*DAF(IVC) 
        IF(IVC.LE.14) THEN 
          VALL(IVC+16)=VCLASS(IVC)*DAF(IVC+16) 
        ENDIF 
      END DO !IVC 
C 
C  Calculate distribution of vehicles by class VBLCASS. 
C  Calculate the vehicle counts by class VTEMP. 
C  LDDT1 and LDDT2 are combined and LDDT3 and LDDT4 are combined. 
C 
      DO IV=1,MAXVEH 
        IF(IV.LE.13) THEN 
          REGMIX(IV)=VALL(IV)/TOTCNT 
          VTEMP(IV)=VALL(IV) 
        ELSEIF(IV.EQ.14) THEN 
          REGMIX(IV)=VALL(17)/TOTCNT 
          VTEMP(IV)=VALL(17) 
        ELSEIF(IV.EQ.15) THEN 
          REGMIX(IV)=(VALL(18)+VALL(19))/TOTCNT 
          VTEMP(IV)=VALL(18)+VALL(19) 
        ELSEIF(IV.GE.16.AND.IV.LE.23) THEN 
          REGMIX(IV)=VALL(IV+6)/TOTCNT 
          VTEMP(IV)=VALL(IV+6) 
        ELSEIF(IV.EQ.24) THEN 
          REGMIX(IV)=VALL(16)/TOTCNT 
          VTEMP(IV)=VALL(16) 
        ELSEIF(IV.EQ.25) THEN 
          REGMIX(IV)=VALL(14)/TOTCNT 
          VTEMP(IV)=VALL(14) 
        ELSEIF(IV.EQ.26) THEN 
          REGMIX(IV)=VALL(15)/TOTCNT 
          VTEMP(IV)=VALL(15) 
        ELSEIF(IV.EQ.27) THEN 
          REGMIX(IV)=VALL(30)/TOTCNT 
          VTEMP(IV)=VALL(30) 
        ELSEIF(IV.EQ.28) THEN 
          REGMIX(IV)=(VALL(20)+VALL(21))/TOTCNT 
          VTEMP(IV)=VALL(20)+VALL(21) 
        ELSE 
          REGMIX(IV)=0.0 
          VTEMP(IV)=0.0 
        ENDIF 
      END DO !IV 
C 
C  Calculate normalized VMTMIX 
C 
      DO IVD=1,30 
        IF(IVD.LE.13) THEN 
          IVX=IVD 
        ELSEIF(IVD.GE.14.AND.IVD.LE.15) THEN 
          IVX=IVD+11 
        ELSEIF(IVD.EQ.16) THEN 
          IVX=24 
        ELSEIF(IVD.EQ.17) THEN 
          IVX=14 
        ELSEIF(IVD.GE.18.AND.IVD.LE.19) THEN 
          IVX=15 
        ELSEIF(IVD.GE.20.AND.IVD.LE.21) THEN 
          IVX=28 
        ELSEIF(IVD.GE.22.AND.IVD.LE.29) THEN 
          IVX=IVD-6 
        ELSEIF(IVD.EQ.30) THEN 
          IVX=27 
        ELSE 
          IVX=1 
        ENDIF 
        DO JDX=1,MAXYRS 
          AGEMILE(IVD)=AGEMILE(IVD) 
     *                +(AGED(JDX,IVD)*TFMYM(JDX,IVX)) 
        END DO !JDX 
      END DO   !IVD 
C 
      TOTMILE=0.0 
      DO IVD=1,30 
        VMTBC(IVD)=AGEMILE(IVD)*VALL(IVD) 
        TOTMILE=TOTMILE+VMTBC(IVD) 
      END DO !IVD 
C 
      IF(TOTMILE.GT.0.0) THEN 
        DO IVD=1,30 
          VMTBC(IVD)=VMTBC(IVD)/TOTMILE 
        END DO !IVD 
      ENDIF 
C 
C  If the user inputs the VMT distribution, then recalculate the 
C  vehicle counts and distribute the VMT by gas/diesel. 
C 
      IF(VMFLAG.EQ.2) THEN 
        DO IVC=1,16 
          IF(IVC.LE.14) THEN 
            RATIO(IVC)=VMTUSR(IVC)/(VMTBC(IVC)+VMTBC(IVC+16)) 
          ELSE 
            RATIO(IVC)=VMTUSR(IVC)/VMTBC(IVC) 
          ENDIF 
        END DO !IVC 
        DO IV=1,MAXVEH 
          IF(IV.LE.13) THEN 
            VMTMIX(IV)=RATIO(IV)*VMTBC(IV) 
            VTEMP(IV)=RATIO(IV)*VALL(IV) 
          ELSEIF(IV.EQ.14) THEN 
            VMTMIX(IV)=RATIO(1)*VMTBC(17) 
            VTEMP(IV)=RATIO(1)*VALL(17) 
          ELSEIF(IV.EQ.15) THEN 
            VMTMIX(IV)=RATIO(2)*VMTBC(18)+RATIO(3)*VMTBC(19) 
            VTEMP(IV)=RATIO(2)*VALL(18)+RATIO(3)*VALL(19) 
          ELSEIF(IV.GE.16.AND.IV.LE.23) THEN 
            VMTMIX(IV)=RATIO(IV-10)*VMTBC(IV+6) 
            VTEMP(IV)=RATIO(IV-10)*VALL(IV+6) 
          ELSEIF(IV.EQ.24) THEN 
            VMTMIX(IV)=RATIO(16)*VMTBC(16) 
            VTEMP(IV)=RATIO(16)*VALL(16) 
          ELSEIF(IV.EQ.25) THEN 
            VMTMIX(IV)=RATIO(14)*VMTBC(14) 
            VTEMP(IV)=RATIO(14)*VALL(14) 
          ELSEIF(IV.EQ.26) THEN 
            VMTMIX(IV)=RATIO(15)*VMTBC(15) 
            VTEMP(IV)=RATIO(15)*VALL(15) 
          ELSEIF(IV.EQ.27) THEN 
            VMTMIX(IV)=RATIO(14)*VMTBC(30) 
            VTEMP(IV)=RATIO(14)*VALL(30) 
          ELSEIF(IV.EQ.28) THEN 
            VMTMIX(IV)=RATIO(4)*VMTBC(20)+RATIO(5)*VMTBC(21) 
            VTEMP(IV)=RATIO(4)*VALL(20)+RATIO(5)*VALL(21) 
          ELSE 
            VMTMIX(IV)=0.0 
            VTEMP(IV)=0.0 
          ENDIF 
        END DO !IV 
      ENDIF 
C 
      RETURN 
      END 
