      BLOCK DATA BD48 
C 
C  BLOCK DATA Subprogram 48: /ATOX1/,/ATOX3/,/ATOX4/,ATOX5/ 
C 
C  BLOCK DATA 43 initializes the Air Toxics commons. 
C 
C  Changes: (Last change first) 
C 
C  06 Sep 01 AIR Task 15: New Block Data. 
C 
C  Common block array subscripts: 
C 
C  ATOXY(3,5)   -  ATOXY(I,IEVAP) 
C  TOXCHK(6)    -  TOXCHK (ICHK) 
C 
C  Common block dictionary: 
C 
C   Name    Type                  Description 
C  ------   ----  ----------------------------------------------------- 
C  /ATOX1/ 
C  AROM      R    percentage of aromatic compounds in gasoline on a Volume Basis (required input) 
C  BENZ      R    percentage of benzene in gasoline on a Volume Basis (required input) 
C  E200      R    percentage of vapor of gasoline fuel at 200F (required input) 
C  E300      R    percentage of vapor of gasoline fuel at 300F (required input) 
C  OLEF      R    percentage of olefin in gasoline on a Volume Basis (required input) 
C  TOXCHK    I    checklist of required air toxics run inputs (1=not input, 2=input) for: 
C                 1: E200, 2: E300, 3: Aromatics, 4: Olefins, 5: Benzene, 6: Oxygenate 
C  TOXFLG    L    switch to run the air toxics subroutine 
C  ATOXY     R    Oxygenate content, wt%, for 1:MTBE, 2:ETBE, 3:ETHANOL, 4:TAME, 5:NoOxy 
C  MAXIAT    I    Number of AT pollutant types 
C  ATSSFLAG  L    switch to print out the air toxics spreadsheet 
C   
C 
      IMPLICIT NONE 
C 
      INCLUDE 'ATOX1.I' 
      INCLUDE 'ATOX3.I' 
      INCLUDE 'ATOX4.I' 
      INCLUDE 'ATOX5.I' 
C 
C  /ATOX1/ 
C 
      DATA AROM    /0./ 
      DATA ATOXY   /15*0./ 
      DATA ATSSFLAG/.FALSE./ 
      DATA BENZ    /0./ 
      DATA E200    /0./ 
      DATA E300    /0./ 
      DATA MAXIAT  /6/ 
      DATA OLEF    /0./ 
      DATA TOXCHK  /6*1/ 
      DATA TOXFLG  /.FALSE./ 
C 
C  /ATOX3/ 
C 
      DATA ACETRATIO  /10*0./ 
      DATA BENZRATIO  /10*0./ 
      DATA BUTARATIO  /10*0./ 
      DATA FORMRATIO  /10*0./ 
      DATA GROUP      /10*0./ 
      DATA MTBERATIO  /10*0./ 
      DATA VOCCOMPLEX /10*0./ 
C 
C  /ATOX4/ 
C 
      DATA ATCONV  /24*0./ 
      DATA ATNHWT  /2*0./ 
      DATA ATHEF   /240*0./ 
      DATA ATEMFAC /240*0./ 
      DATA ATHEV   /600*0./ 
      DATA DBAT    /3024000*0./ 
      DATA EFAT    /1008*0./ 
      DATA VAT     /36*0./ 
C 
C  /ATOX5/ 
C 
      DATA ATHAPS    /.FALSE./ 
      DATA ATHBEGMY  /50*0/ 
      DATA ATHBER    /560000*0./ 
      DATA ATHDET    /50*0./ 
      DATA ATHEKIND  /50*' '/ 
      DATA ATHENDMY  /50*0/ 
      DATA ATHPMAP   /50*0/ 
      DATA ATHPNAME  /50*' '/ 
      DATA ATHPOL    /50*0/ 
      DATA ATHRATIO  /50*' '/ 
      DATA ATHRATZML /50*0./ 
      DATA ATHVTYPE  /50*0/ 
      DATA NATH      /0/         
C 
      END 
       
