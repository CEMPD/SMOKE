C 
C  Changes: 
C 
C  18 Jan 01 @EPA-bag 2007 rule for heavy duty vehicles; increased 
C         size of ERBZML,ERBDR,SEBDR,SEBZML,REBDR, and REBZML arrays 
C   6 JAN 99 @DynTel-HXQ 2-668 Pre-1981 LDGV & LDGT basic exhaust 
C             engine start and running; removed ERB50K array 
C  22 Oct. 98 - DynTel@HXQ  2-680 Expanded arrays ERBZML, ERBDR, SEBDR, 
C 		SEBZML, REBDR, REBZML to accomodate the new model years. 
C  2 Apr 97 - DynTel@MLA  Expanded arrays ERBZML and ERBDR to accomodate  
C             the new BUS vehicle categories. 
C  9 Apr 98  @DynTel-RJD  2-654  Added additional arrays to accommodate 
C            seperate starting and running emissions. 
C  29 June 98 @dyntel -MMS expanded start and running arrays  
C 
C  The  initialization of the data could be found in BD04. 
C 
       
      REAL ERBZML,ERBDR,SEBDR,SEBZML, 
     &     REBDR,REBZML 
 
      COMMON /BASEQ1/ ERBZML(35,3,28,2),ERBDR(35,3,28,2), 
     &                SEBDR(35,3,28,2),SEBZML(35,3,28,2), 
     &                REBDR(35,3,28,2),REBZML(35,3,28,2) 
