C 
C The data is initialized in BD30. 
C  
C  Changes: (Last change first) 
C 
C 20 Jun 00 @EPA-elg Changed RFGRVP array to 18 and removed GLRFGRVP. 
C 16 May 00 @EPA-ddj RFG BUG, Removed RFG_S array. 
C  4 May 00 @EPS-ddj D-28 Changed dimension on RFG_S from 2 to 5. 
C 28 Mar 00 @DynTel-rwk 1-028, Added global RFG RVP variable for 
C           reinitialization 
C 25 Feb 00 @DynTel-ddj 1-028, Changed dimensions on RFGSYR for 3 to 2 
C           and on RFGRVP from (5,3) to (2,2).  Added array RFG_CF. 
C 19 May 98 @DynTel-MLA 2-659  Added GBLRFGSYR for Mobile6. This variable 
C           holds the global value for RFGSYR(2), which might be reset, 
C           depending on the value of RFGFLG. 
C 
      INTEGER RFGSYR, GLRFGSYR 
C 
      REAL RFGRVP, RFG_CF 
C 
      COMMON /RFORM2/ RFGSYR(2), RFGRVP(2,18), RFG_CF(3,2,18,5), 
     &              GLRFGSYR(2) 
