C
C  Changes: (Last change first)
C
C  11 Jul 02 AIR Task 36: Hourly relative humidity correct made EMFAC
C            array hour-dependent
C  18 Jan 01 @EPA-bag 2007 rule for heavy duty vehicles; increased
C         size of ERUZML, ERUDR, ESUZML and ESUDR arrays
C  08 Jun 00 @EPA-djb Bug113 Removed HIGHER array.
C  16 May 00 @EPA-ddj BUG27 collapsed parameters for array EMFAC to
C            (ISR,ISTD,IFDS,INH).
C   3 March 00 @EPA ELG   Added EMFAC(IP,IV,JDX,ISR,ISTD,IFDS,INH)
C                         Added HIGHER(INH,IV,IP,JDX,IFDS)
C   9 July 99 @ Dyntel JWR 2-705 Increased all 1st dimensions to 33 to
C             allow user inputs from 1950 - present.
C
C   6 JAN 99 @DynTel-HXQ 2-668 Pre-1981 LDGV & LDGT basic exhaust
C             engine start and running; removed ERU50K array
C 13 Apr 98  @DynTel-RJD  Added arrays ESUZML and ESUDR to accommodate
C             starts emissions.
C
C  2 Apr 97 - DynTel@MLA  Expanded arrays ERUZML and ERUDR to accomodate 
C             the new BUS vehicle categories.
C
C  The initialization of the data could be found in BD04.
C
      REAL ERUZML,ERUDR,ESUZML,ESUDR,EMFAC
      COMMON /BASEQ5/ ERUZML(35,3,28,2),ERUDR(35,3,28,2),
     &                ESUZML(35,3,28,2),ESUDR(35,3,28,2),
     &                EMFAC(2,10,3,2,24)
