C Version "@(#)$Id$ $Source$ $Date$ 

C
C  IMPAR4 contains the variables that holds the values from 
C         the 'Smart Line' in the I/M credit file. 
C
C  Common Block Dictionary:
C
C   Name    Type                  Description
C  ------   ----  -----------------------------------------------------
C  CUTHC     R     Cutpoint for HC
C  CUTCO     R     Cutpoint for CO
C  CUTNO     R     Cutpoint for NOx
C  SMARTL    C     Keeps information from the
C                  'smart line I/M credit description'         
C  SMRTL1    C     Keeps information from the
C                  'smart line Tech1&2 credit description'         
C
C
C  Notes:
C
C
      REAL CUTHC,CUTCO,CUTNO
      CHARACTER*80 SMARTL,SMRTL1
      COMMON /IMPAR4/ CUTHC(2),CUTCO(2),CUTNO(2),SMARTL(2),
     *                SMRTL1                                                    
