      BLOCK DATA BDEV01 
C 
C  BLOCK DATA for the evap I/M and OBD data 
C 
C 
C  Changes: (Last change first) 
C 
C  15 Sep 00 AIR Task 08: New block data. 
C 
C 
C  Common block array subscripts: 
C 
C  TEST_GC(16)     - TEST_GC   ( MY - MYB_GC ) 
C  TEST_FP(18)     - TEST_GC   ( MY - MYB_FP ) 
C  FAIL_IL(3,26)   - FAIL_IL   ( I, JDX) 
C  FAIL_AZFP(3,26) - FAIL_AZFP ( I, JDX) 
C  EVNOIM(4,26,2)  - EVNOIM    (IPF, JDX, IOBD) 
C  EVIM(4,26,2)    - EVIM      (IPF, JDX, IOBD) 
C  EVSAW(4,25,2)   - EVSAW     (IPF, JDX, IOBD) 
C 
C 
C  Common block dictionary: 
C 
C   Name    Type                  Description 
C  ------   ----  ----------------------------------------------------- 
C  TEST_GC   R    Gas cap testability rates. 
C  TEST_FP   R    Fill-pipe testability rates. 
C  FAIL_IL   R    Gas cap failure rates - Illinois data, 
C                  I=1=MY95-, I=2=MY99+, I=3=ETP. 
C  FAIL_AZFP R    Fill-pip failure rates - Arizona data, 
C                  I=1=MY95-, I=2=MY99+, I=3=ETP. 
C  FAIL_AZGC R    Gas cap failure with fill-pipe pass rates - Arizona data, 
C                  I=1=MY95-, I=2=MY99+. 
C  EVNOIM    R    Evap pass/fail weighting for no I/M case, basic values,  
C                  IOBD=1=without OBD, IOBD=2=with OBD  
C                  IPF=1=pass pressure/pass purge, 
C                  IPF=2=pass pressure/fail purge, 
C                  IPF=3=fail pressure/pass purge,  
C                  IPF=4=fail pressure/fail purge. 
C  EVIM      R    Evap pass/fail weighting for with I/M case, basic values  
C                  IOBD=1=without OBD, IOBD=2=with OBD  
C                  IPF=1=pass pressure/pass purge, 
C                  IPF=2=pass pressure/fail purge, 
C                  IPF=3=fail pressure/pass purge,  
C                  IPF=4=fail pressure/fail purge. 
C  MYB_GC    I    First model year for gas cap testability. 
C  MYB_FP    I    First model year for fill-pipe testability. 
C  EVSAW     R    Evap pass/fail weighting, no I/M and I/M sawtooth weighted,  
C                  IOBD=1=without OBD, IOBD=2=with OBD  
C                  IPF=1=pass pressure/pass purge, 
C                  IPF=2=pass pressure/fail purge, 
C                  IPF=3=fail pressure/pass purge,  
C                  IPF=4=fail pressure/fail purge. 
C  GLLRDC    R    Gross liquid leaker fraction reduction 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'BDEVIM.I' 
      INCLUDE 'EVIM.I' 
C 
C  Gas cap testability model year start and fractions 
C 
      DATA MYB_GC/1971/ 
C 
      DATA TEST_GC/ 
C        1971-  1972   1973   1974  
     * 0.8841,0.8698,0.8479,0.8146, 
C        1975   1976   1977   1978  
     * 0.7834,0.7522,0.7878,0.7522, 
C        1979   1980   1981   1982  
     * 0.9329,0.9426,0.9610,0.9689, 
C        1983   1984   1985   1986+  
     * 0.9745,0.9804,0.9999,1.0000/ 
C 
C  Fill pipe testability model year start and fractions 
C 
      DATA MYB_FP/1981/ 
C 
      DATA TEST_FP/ 
C        1981-  1982   1983   1984  
     * 0.6093,0.6292,0.6811,0.7115, 
C        1985   1986   1987   1988  
     * 0.7315,0.7718,0.7862,0.7771, 
C        1989   1990   1991   1992  
     * 0.7666,0.7563,0.7936,0.7965, 
C        1993   1994   1995   1996  
     * 0.7879,0.7616,0.6401,0.4391, 
C        1997   1998+  
     * 0.3519,0.1191/ 
C 
C  Initialize failure rates 
C 
      DATA FAIL_IL,FAIL_AZFP,FAIL_AZGC/78*0.,78*0.,52*0./ 
C 
C  Initialize the basic evap pass/fail rates, no I/M, with I/M 
C 
      DATA EVNOIM,EVIM/208*0.,208*0./ 
C 
C  Initialize sawtooth weighted evap no I/M and I/M pass/fail arrays 
C 
      DATA EVSAW/200*0./ 
C 
C  Initialize the gross liquid leaker reduction variable 
C 
      DATA GLLRDC/0./ 
C 
      END 
  
