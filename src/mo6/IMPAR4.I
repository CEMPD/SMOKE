C 
C  IMPAR4 contains the variables that holds the values from  
C         the 'Smart Line' in the I/M credit file.  
C 
C  Changes: (Last change first) 
C 
C 25 Aug 00 @AIR    Expanded arrays CUTHC, CUTCO, CUTNO to 
C           accomodate 7 IM programs. 
C 22 Jun 99 @DynTel-MLA 2-000  Removed SMARTL (unused). 
C 17 Feb 99 @DynTel-MLA 2-663  Expanded the cutpoint arrays to hold 
C           information by model year for 4 different vehicle classes 
C           for each I/M program. 
C 
C 18 Dec 98 @DynTel-MLA 2-663  Added IVCUT to map IV into an index for 
C           CUTHC, CUTCO, and CUTNOX. 
C 
C  Common Block Dictionary: 
C 
C   Name    Type                  Description 
C  ------   ----  ----------------------------------------------------- 
C  CUTHC     R     Cutpoint for HC 
C  CUTCO     R     Cutpoint for CO 
C  CUTNO     R     Cutpoint for NOx 
C  IVCUTIDX  I     Maps IV into an index for CUTHC, CUTCO, CUTNOX. 
C  SMRTL1    C     Keeps information from the 
C                  'smart line Tech1&2 credit description'          
C 
C 
C  Notes: 
C 
C  The data is initialized in BD_IM2. 
C 
      REAL CUTHC,CUTCO,CUTNO 
      CHARACTER*80 SMRTL1 
      INTEGER IVCUTIDX 
      COMMON /IMPAR4/ CUTHC(25,4,7),CUTCO(25,4,7),CUTNO(25,4,7), 
     *                SMRTL1,IVCUTIDX(27) 
