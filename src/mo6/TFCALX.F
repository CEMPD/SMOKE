      SUBROUTINE TFCALX(INERR) 
C 
C  TFCALX computes the travel fractions and composite class VMT mix. 
C 
C  Called by REGMOD. 
C 
C  Calls QUITER. 
C 
C  Changes: (Last change first) 
C 
C  07 Jul 00 AIR Task 05: Removed unused ICY from parameter list. 
C  07 Apr 00 @EPA-djb Major mods in logic and arrays. 
C  28 Feb 00 @DynTel-MLA 1-011  Added DBVCOUNT to collect database output. 
C  24 Jan 00 @DynTel-MLA 1-011  Added VMTLDDT. 
C  12 Jan 99 @DynTel-ZK 2-000 Explicit Typing 
C  9 Nov. 98 @DynTel-HXQ 2-680 to delete VMTMIX(28). 
C  06 Nov 98 @Dyntel-BAG 2-000 Renamed VMTGT1 and VMTGT2 to VMTGT12 and VMTGT34, respectively 
C  11 Jun 97 DynTel@MLA 2-622  Added IVTYPE.I for Mobile6 and modified the 
C            code to use the symbolic vehicle type names and vehicle type 
C            vectors that were introduced in Mobile6. Adjusted the database 
C            indexes to m atch the new vehicle type numbering in Mobile6. 
C 
C  02 Mar 98 DynTel@RJD 2-204 Changed module to interpolate between vehicle  
C            count data stored in BD20.FOR so that identical vehicle counts  
C            are not used for January 1st and July 1st calculations in any  
C            given year.                   
C 
C  03 Mar 98 @DynTel-RJD 2-204  Changed the name of the array storing the Model 
C            Year Registration fractions from JANMYR to MEVMYR to reflect the 
C            fact that January is no longer the only month of evaluation that 
C            the Mobile model can handle. 
C 
C  10 Mar 98 @DynTel-RJD 2-622  Moved calculation of VMTMIX for composite  
C            vehicles from the OUTDT modules to this one. 
C  
C  12 Mar 98 @DynTel-RJD 2-622  Moved calculation of VMTMIX for composite  
C            LDGV vehicles from the EFCALX module to this one. 
C 
C  23 Jul 98 @DynTel-NH  2-000  Removed declaration of I81 as a REAL. 
C   
C  23 Jul 98 @DynTel-MLA 2-000  Changed LDGT names to conform to 
C            certification naming conventions. 
C 
C  Input on call: 
C 
C    Constants:  
C    MAXVEH from include file IVTYPE.I 
C 
C    argument list: INERR 
C 
C    common blocks: 
C    /EVAL/   MEVAL 
C    /FLAGS1/ VMFLAG 
C    /MAXIMA/ MAXYRS 
C    /MYRCAL/ MEVMYR,TFMYM 
C    /REGISF/ GSFVCT 
C    /VVTYPS/ VTLDGT1, VTLDGT2, VTLDGT3, VTLDGT4  (IVTYPE.I) 
C 
C  Output on return: 
C 
C    argument list: INERR 
C    common blocks: 
C    /MYRCAL/ TF 
C    /VMXCOM/ DBVCOUNT,REGMIX,TFNORM,VMTMIX,VMTGT12,VMTGT34,VMTHDG, 
C             VMTLDDT,VMTHDD,VMTLDGT 
C 
C  Local array subscripts: 
C 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  IDX      I    vehicle age loop index(1-25 years) 
C  INERR    I    count of errors 
C  IV       I    standard veh class groupings (28) 
C  IVC      I    combined gas/diesel veh class groupings (16) 
C  IVD      I    separate gas/diesel veh class groupings (30) 
C  JDX      I    vehicle age index; JDX=25-IDX+1 
C  TOTVMT   R    total vmt: used to normalize the default vmt mix 
C 
C  Notes: 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'BASEQ9.I' 
      INCLUDE 'EVAL.I' 
      INCLUDE 'FLAGS1.I' 
      INCLUDE 'IVTYPE.I' 
      INCLUDE 'MAXIMA.I' 
      INCLUDE 'MYRCAL.I' 
      INCLUDE 'REGISF.I' 
      INCLUDE 'VMXCOM.I' 
C 
      REAL    TOTVMT 
      INTEGER IV, JDX, INERR, IDX 
C 
C  Calculate normalized VMTMIX 
C 
      DO IV=1,MAXVEH 
        TFNORM(IV)=0.0 
        DO JDX=1,MAXYRS 
          TFNORM(IV)=TFNORM(IV)+TFMYM(JDX,IV)*MEVMYR(JDX,IV) 
C 
C  Initialize the travel fraction array 
C 
          TF(JDX,IV)=0.0 
        END DO !JDX 
      END DO   !IV 
C 
      TOTVMT=0.0 
      DO IV=1,MAXVEH 
        DBVCOUNT(IV)=VTEMP(IV) 
        VMTMIX(IV)=VTEMP(IV)*TFNORM(IV) 
        TOTVMT=TOTVMT+VMTMIX(IV) 
      END DO   !IV 
C 
      IF(TOTVMT.GT.0.0) THEN 
        DO IV=1,MAXVEH 
          VMTMIX(IV)=VMTMIX(IV)/TOTVMT 
        END DO !IV 
      ENDIF 
C 
C  Check for veh classes with no mileage accumulation before 
C  calculating the travel fraction by age. (TF) 
C 
      DO IV=1,MAXVEH 
        IF(TFNORM(IV).LE.0.0.AND.VMTMIX(IV).GT.0.0) THEN 
          CALL QUITER(0.,IV,48,INERR) 
        ELSEIF(TFNORM(IV).LE.0.0) THEN 
          CALL QUITER(0.,IV,48,INERR) 
        ELSE 
          DO JDX=1,MAXYRS 
            IDX=MAXYRS-JDX+1 
            TF(IDX,IV)=(TFMYM(JDX,IV)*MEVMYR(JDX,IV))/TFNORM(IV) 
          END DO !JDX 
        ENDIF 
      END DO     !IV 
C 
C Calculate VMTMIX for composite vehicles used in descriptive output. 
C 
      VMTGT12 = VMTMIX(2)  + VMTMIX(3) 
      VMTGT34 = VMTMIX(4)  + VMTMIX(5) 
      VMTHDG  = VMTMIX(6)  + VMTMIX(7)  + VMTMIX(8)  + VMTMIX(9) + 
     *          VMTMIX(10) + VMTMIX(11) + VMTMIX(12) + VMTMIX(13) + 
     *          VMTMIX(25) 
      VMTHDD  = VMTMIX(16) + VMTMIX(17) + VMTMIX(18) + VMTMIX(19) + 
     *          VMTMIX(20) + VMTMIX(21) + VMTMIX(22) + VMTMIX(23) + 
     *          VMTMIX(26) + VMTMIX(27) 
      VMLDGT  = VMTMIX(2)  + VMTMIX(3)  + VMTMIX(4)  + VMTMIX(5) 
      VMTLDDT = VMTMIX(15) + VMTMIX(28) 
C 
      RETURN 
      END 
