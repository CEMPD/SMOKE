      SUBROUTINE CALCUM() 
C 
C  CALCUM adjusts model year mileage arrays. 
C 
C  Called by REGMOD and CHKINP. 
C 
C  Calls (None) 
C 
C  Changes :  
C 
C  07 Jul 00 AIR Task 05: Removed unused ICY from parameter list. 
C  07 Apr 00 @EPA-djb Change comments and 
C            indentation of code. Call to 
C            HDDMYM removed. 
C  25 Aug 98 @DYNTEL-ZK 2-000 Explicit Typing 
C 
C  Input on call: 
C 
C    Constants:  
C    MAXVEH from include file IVTYPE.I 
C 
C    argument list: NONE 
C    common blocks: 
C    /EVAL/   MEVAL 
C    /FLAGS2/ MYMRFG 
C    /MAXIMA/ MAXYRS 
C    /MYRSAV/ AMAR 
C    /VVNAMS/ VVHEAVY  (IVTYPE.I) 
C    /VVTYPS/ VTMC  (IVTYPE.I) 
C 
C  Output on return: 
C 
C    argument list: INERR 
C    common blocks: 
C    /EVAL/   FRC,FRN  
C    /MYRCAL/ XMYM 
C 
C  Local array subscripts: 
C 
C  CUM(25)  -  CUM ( JDX ) 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  CUM      R    cumulative mileage distribution as of January 1 or July 1 
C  FRC      R    fraction of January 1 or July 1 registrations in Ith  
C                year of operation with model year index I.  
C  FRN      R    fraction of January 1 or July 1 registrations in Ith  
C                year of operation with model year index I+1 (ie. the  
C                vehicles were sold before the model year change over 
C                date). 
C 
C  Notes: 
C 
C  CALCUM was created for MOBILE4.1v6. 
C  11 Jun 97 DynTel@MLA 2-622  Added IVTYPE.I for Mobile6 and modified the 
C            code to use the symbolic vehicle type names and vehicle type 
C            vectors that were introduced in Mobile6. 
C 
C  02 Mar 98 DynTel@RJD 2-204 Elevated FRC and FRN to the status of /EVAL/          
C            common block variables so that REGMOD.FOR can evaluate FRC,  
C            FRN and XMYM via a call to CALCUM.FOR 
C 
      IMPLICIT NONE  
      INCLUDE 'EVAL.I' 
      INCLUDE 'FLAGS2.I' 
      INCLUDE 'IVTYPE.I' 
      INCLUDE 'MAXIMA.I' 
      INCLUDE 'MYRCAL.I' 
      INCLUDE 'MYRSAV.I' 
C 
C  Generate the mileage accrual rates. 
C 
C  If MOBILE4 defaults used, determine the HDDV calendar year mileage 
C  by registration weighting. Otherwise, use user input mileage. 
 
      REAL AMC, AMN, CUM(25) 
      INTEGER IV, JDX 
C 
C  Compute the annual mileage / 100,000  distributions, XMYM. 
C 
C  Evaluate FRC and FRN in terms of MEVAL. Note: this module has only been 
C  tested for MEVAL = 1 and 7. 
C 
      DO IV=1,MAXVEH 
        IF (MEVAL.LT.10)THEN 
          FRC(IV)=.25+((MEVAL-1)/12.0) 
        ELSE 
          FRC(IV)=(MEVAL-10)/12.0 
        ENDIF 
        IF(VVHEAVY(IV).EQ.1.OR.IV.EQ.VTMC)THEN 
          FRC(IV)=(MEVAL-1)/12.0 
        ENDIF 
        FRN(IV) = 1.0-FRC(IV) 
      END DO  !IV 
C 
C     Evaluate cummulative mileage CUM(JDX) for vehicles with an average 
C     age of JDX-1+0.5*FRC(IV) years.  
C 
      DO 30 IV=1,MAXVEH 
C 
       AMN=AMAR(1,IV) 
       AMC=AMAR(2,IV) 
       CUM(1)=.5*FRC(IV)*AMN 
       CUM(2)=FRC(IV)*(.5*FRC(IV)*AMC+AMN)+FRN(IV) 
     *   *((1.+FRC(IV))*.5*AMN) 
       AMC=AMC+AMN 
C 
        DO 10 JDX=3,MAXYRS 
          CUM(JDX)=FRC(IV)*(.5*FRC(IV)*AMAR(JDX,IV)+AMC)+FRN(IV) 
     *     *(.5*(1.+FRC(IV))*AMAR(JDX-1,IV)+AMN) 
          AMN=AMC 
          AMC=AMC+AMAR(JDX,IV) 
   10   CONTINUE 
C 
C  Subroutine GETCUM produces the cumulative mileage partial sums array, 
C  CUMMIL, from XMYM. 
C 
        XMYM(1,IV)=CUM(1) 
        DO 20 JDX=2,MAXYRS 
          XMYM(JDX,IV)=CUM(JDX)-CUM(JDX-1) 
   20   CONTINUE 
C 
   30 CONTINUE !IV 
C 
      RETURN 
      END 
