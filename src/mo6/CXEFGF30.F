      SUBROUTINE CXEFGF30(IV,MY,ICY,CXNOLD,CXYSLD,INERR) 
C 
C  Calculates carbon emissions for gasoline vehicles at PM30. 
C 
C Changes: 
C 
C 03 Jan 02 AIR Task 15-29: Corrected a long-standing PART5 bug in the 
C           formula used to calculate CXNOLD. 
C 23 OCT 01 AIR Task 25: Values returned are at PM30. 
C 16 Oct 01 AIR Task 25: New MOBILE6.2T routine (based on 
C           PART5 routine), extended to include deterioration. 
C 13 Feb 01 AIR Task 11: Change code to use MOBILE6 vehicle 
C           classes and technology fractions. 
C 14 Nov 00 AIR Task 06: Change DO loop maximum to NCAT. 
C 10 Nov 00 AIR Task 06: Updated code to MOBILE6 standards. 
C 
C  Called by: CEFPM 
C 
C  Calls: PMFRACTS 
C 
C  Output on Return: leaded and unleaded carbon emissions for gasoline 
C  vehicles 
C                     
C  Common Blocks:  
C  /EFPMC/ PMGZML,PMGDR1,PMGDR2,NCAT 
C  /PART4/ PMDRAGE 
C  /SULFR/ EXPECV 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ----------------------------------------------- 
C  CXICAT   R    Emission value of PMGZML for a given ICAT 
C  FCATFR   R    Fraction catalyst type for ICAT 1/2/3/4 depends 
C                on IV 
C  ICAT     I    1=LEADED,, 2=UNLEADED CAT NO AIR, 
C                3=UNLEADED CAT AIR, 4=UNLEADED NO CAT 
C  MYPM     I    model year index into the PMGZML, PMGDR1 and 
C                PMGDR2 arrays, 1=1950-, 71=2020+ 
C  CAT      R    Total catalyst fraction 
C   
C  Notes: 
C  May-12-1994 @ CSC-bsg Request 430 created for revisions of emission 
C  factor for gasoline vehicles. 
C  May-24-1994 @ CSC-bsg Request 434 integration and regression testing 
C  Jun-17-1994 @ CSC-bsg Request 430 MC use lead instead of gas carbon 
C  Jun-24-1994 @ CSC-bsg Request 430 MC use lead and gas carbonn 
C 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'CUMCOM.I'   ! CUMMIL 
      INCLUDE 'EFGAS.I'    ! PMGZML,PMGDR1,PMGDR2,NCAT 
      INCLUDE 'FRACT.I'    ! EXPECV 
      INCLUDE 'IVTYPE.I'   ! VTMC 
      INCLUDE 'PART4.I'    ! PMDRAGE 
C 
C  Declare external functions. 
C 
      REAL, EXTERNAL :: PMFRACTS 
C 
C  Declare parameter list. 
C 
      INTEGER ICY 
      INTEGER INERR 
      INTEGER IV 
      INTEGER MY 
      REAL    CXNOLD 
      REAL    CXYSLD       
C 
C  Declare local variables/arrays. 
C 
      INTEGER IAGE 
      INTEGER ICAT 
      INTEGER MYPM 
      REAL    BODOM 
      REAL    CAT 
      REAL    CXICAT(4) 
      REAL    FCATFR(4) 
      REAL    KODOM 
C 
C  Determine the age 
C 
      IAGE = ICY - MY + 1 
C 
C  Determine the base and kink mileage, in 10K. 
C 
      IF(IAGE.LE.PMDRAGE(IV)) THEN 
        BODOM=CUMMIL(IAGE,IV)/10000. 
        KODOM=0. 
      ELSE 
        IF(PMDRAGE(IV).EQ.0) THEN 
          BODOM=0. 
        ELSE 
          BODOM=CUMMIL(PMDRAGE(IV),IV)/10000. 
        END IF 
        KODOM=CUMMIL(IAGE,IV)/10000. - BODOM 
      END IF       
C 
C  Determine the model year pointer into the PM_BER arrays. 
C 
      MYPM=MY-1950+1 
      IF(MYPM.LT. 1) MYPM= 1 
      IF(MYPM.GT.71) MYPM=71 
C 
C  Motorcyles only have leaded fuel carbon emissions 
C 
      IF(IV.EQ.VTMC) THEN 
        CXNOLD=0.0 
        CXYSLD=(PMGZML(MYPM,1,IV) 
     *        + PMGDR1(MYPM,1,IV) * BODOM 
     *        + PMGDR2(MYPM,1,IV) * KODOM) 
        RETURN 
      ENDIF 
C 
C  Assign fraction using unleaded cat no air result of FCTNS 
C 
      FCATFR(1)=1.0 
      FCATFR(2)=PMFRACTS(MY,ICY,IV,4,INERR) 
      FCATFR(3)=PMFRACTS(MY,ICY,IV,5,INERR) 
      FCATFR(4)=PMFRACTS(MY,ICY,IV,3,INERR) 
      CAT=FCATFR(2)+FCATFR(3) 
C 
C  Assign Carbon array to element of MYPTRC array model year pointer 
C 
      DO 30 ICAT=1,NCAT 
        CXICAT(ICAT)= PMGZML(MYPM,ICAT,IV) 
     *              + PMGDR1(MYPM,ICAT,IV) * BODOM 
     *              + PMGDR2(MYPM,ICAT,IV) * KODOM 
  30  CONTINUE 
C 
      CXYSLD=CXICAT(1)*FCATFR(1) 
C 
      IF(CAT.EQ.0.) THEN 
        CXNOLD=CXICAT(4)*FCATFR(4) 
      ELSE 
        CXNOLD=(CXICAT(2)*FCATFR(2)+CXICAT(3)*FCATFR(3))* 
     *          PMFRACTS(MY,ICY,IV,2,INERR)*EXPECV/CAT+ 
     *          CXICAT(4)*FCATFR(4) 
      END IF 
C 
      RETURN 
      END 
