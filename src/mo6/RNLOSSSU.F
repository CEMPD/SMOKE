      SUBROUTINE RNLOSSSUM(IH,IVTL,MY,ACTUAL_MY,JDX,RNGLOS) 
C 
C 
C  RNLOSSSUM computes the an hourly running loss in grams/mile. There is 
C  a running loss value for each model year and gas vehicle class at the 
C  given RVP and running Loss temperature using 4-point linear interpolation. 
C  Note: This routine is not called for any Diesel vehicles or for motorcycles. 
C 
C    Called by EVP_HRLY 
C 
C    Calls: RNLOSS 
C 
C  Changes: 
C 
C  13 Feb 01 @EPA-elg increased the size of FVMT arrays to include vehicle class 
C  13 Dec 00 @EPA-bag t-x6 MY and ACTUAL_MY added to parameter list for CAA switch 
C 01 Aug 00 @EPA-djb Bug 210 - Added effect from NO SPEED CORRECTION label 
C           and added and altered comments. 
C 18 Jan 00 @DynTel-MLA 1-011  Added ICY to the parameter list. (Required 
C           for the call to RNLOSS.) 
C 13 Jan 00 @DynTel-ddj  1-013  Placed AVSPR and AVSPL into common 
C           block (SPEED9.I). 
C 13 Dec 99 @DynTel-bg 1-006 Added parameters for road types; running loss is now passed 
C           through the subroutine argument list; running loss is also passed through the 
C           RNLOSS subroutine's parameter list 
C  5 May 99 @DynTel-HQ 2-664 Add this routine to sum up the run loss 
C                              from rnloss subroutine 
C 
C    Constants: 
C    MAXVEH from include file IVTYPE.I 
C 
C     common blocks: 
C     /FLAGS5/ SPD_FLAG 
C     /SPEED9/ FVMT, SVMT  
C 
      IMPLICIT NONE 
      INCLUDE 'FLAGS5.I' 
      INCLUDE 'SPEED9.I' 
C 
C 
C     Declare parameters. 
C 
      INTEGER, INTENT(IN) :: IH 
      INTEGER, INTENT(IN) :: IVTL 
      INTEGER, INTENT(IN) :: MY 
      INTEGER, INTENT(IN) :: ACTUAL_MY 
      INTEGER, INTENT(IN) :: JDX 
      REAL,    INTENT(OUT), DIMENSION(5) :: RNGLOS 
C 
C     Declare local variables. 
C 
      INTEGER I 
      REAL    SPD 
      REAL    RNLOSSTEM 
      INTEGER, PARAMETER :: IFREEWAY=1, IARTERIAL=2        ! roadway type indices 
      INTEGER, PARAMETER :: ILOCAL=3, IRAMP=4, IALLROADS=5 ! roadway type indices 
C 
C  initialize 
      RNGLOS(IFREEWAY)=0.0 
      RNGLOS(IARTERIAL)=0.0 
      RNGLOS(ILOCAL)=0.0 
      RNGLOS(IRAMP)=0.0 
      RNGLOS(IALLROADS)=0.0 
C 
      IF(SPD_FLAG.EQ.1) THEN 
      CALL RNLOSS(IH,IVTL,MY,ACTUAL_MY,JDX,2.5,RNLOSSTEM) 
      RNGLOS(IFREEWAY)= RNGLOS(IFREEWAY)+RNLOSSTEM*svmt(IFREEWAY,1,IH)   !freeway 
      RNGLOS(IARTERIAL)=RNGLOS(IARTERIAL)+RNLOSSTEM*svmt(IARTERIAL,1,IH) !arterial 
      SPD=0.0 
      DO I=2,14 
      SPD=SPD+5.0 
      CALL RNLOSS(IH,IVTL,MY,ACTUAL_MY,JDX,SPD,RNLOSSTEM) 
      RNGLOS(IFREEWAY)= RNGLOS(IFREEWAY)+RNLOSSTEM*svmt(IFREEWAY,I,IH)   !freeway 
      RNGLOS(IARTERIAL)=RNGLOS(IARTERIAL)+RNLOSSTEM*svmt(IARTERIAL,I,IH) !arterial 
      END DO 
      CALL RNLOSS(IH,IVTL,MY,ACTUAL_MY,JDX,AVSPL,RNLOSSTEM) 
      RNGLOS(ILOCAL)=RNLOSSTEM                       !local 
      CALL RNLOSS(IH,IVTL,MY,ACTUAL_MY,JDX,AVSPR,RNLOSSTEM) 
      RNGLOS(IRAMP)=RNLOSSTEM                        !ramp 
C 
      RNGLOS(IALLROADS)=RNGLOS(IFREEWAY)*FVMT(IFREEWAY,IH,IVTL)+ 
     *         RNGLOS(IARTERIAL)*FVMT(IARTERIAL,IH,IVTL)+ 
     *         RNGLOS(ILOCAL)*FVMT(ILOCAL,IH,IVTL)+ 
     *         RNGLOS(IRAMP)*FVMT(IRAMP,IH,IVTL) 
      ELSE     !Speed corrections turned off 
        CALL RNLOSS(IH,IVTL,MY,ACTUAL_MY,JDX,19.6,RNLOSSTEM) 
        DO I=1,5 
          RNGLOS(I)=RNLOSSTEM 
        END DO 
      ENDIF 
C 
      RETURN 
      END 
