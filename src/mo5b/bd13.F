C Version "@(#)$Id$ $Source$ $Date$ 

      BLOCK DATA BD13
C
C  BLOCK DATA Subprogram 13: /EVADU1/,/EVADU2/
C
C  BLOCK DATA 13 initializes the parameters and rates of the diurnal
C  component of evaporative HC.
C
C  Common block array subscripts:
C
C  DUEQ(3,2,3)   -  DUEQ  ( IFDS ,  IVV , IEVEQ )
C  EFDU (3,2,5,5)-  EFDU  ( KEYEQ , IRVP , ISTD , IVGAS )
C  HIDU(3,2)     -  HIDU  ( KEYEQ, IRVP )
C  R2DU(2,2)     -  R2DU  ( KSTD , IREJN )
C
C  Common block dictionary:
C
C   Name   Type              Description
C  ------  ----  ---------------------------------------------
C  /EVADU1/:
C
C  DUEQ     R    1981+ LDGV and LDGT1/2 and Cert DU equations parameters
C  EFDU     R    EF DU rates for IVGAS vehicles for Pass/Fail and 2 RVP fuels
C
C  /EVADU2/:
C
C  HIDU     R    High alt, 1972-76, LDGV & LDGT1, DU rates on 2 RVP fuels
C  R2DU     R    values for high RVP diurnal based on altitude and standard
C                for In-Use condition (IDUAF==> 1=FTP, 2=IN-USE)
C
C  Local array subscripts:
C
C  LDDU(3,2,5)  -  LDDU    ( KEYEQ, IRVP , ISTD  )
C  T1DU(3,2,5)  -  T1DU    ( KEYEQ, IRVP , ISTD  )
C  T2DU(3,2,5)  -  T2DU    ( KEYEQ, IRVP , ISTD  )
C  HDDU(3,2,5)  -  HDDU    ( KEYEQ, IRVP , ISTD  )
C  CCLU(3,2,5)  -  LDDU    ( KEYEQ, IRVP , ISTD  )
C
C  Local array dictionary:
C
C   Name   Type                Description
C  ------  ----  -----------------------------------------------
C  LDDU     R    EF DU rates for pre-81 LDGV
C  T1DU     R    EF DU rates for pre-81 LDGT1
C  T2DU     R    EF DU rates for pre-81 LDGT2
C  HDDU     R    EF DU rates for pre-85 HDGV
C  CLDU     R    EF DU rates for MC
C
C
C  Notes:
C
C  These data structures are in support of the Release MOBILE4
C  evaporative emissions calculation methodology.
C  BLOCK DATA 13 was modified to reflect the changes in diurnal,
C  the changes account for the difference between a 1 hour and
C  an 8 hour measurement of diurnal emissions.
C
C
      INCLUDE 'EVADU1.I'
      INCLUDE 'EVADU2.I'
C
      REAL LDDU
C
      DIMENSION LDDU(3,2,5),T1DU(3,2,5),
     *          T2DU(3,2,5),HDDU(3,2,5),CLDU(3,2,5)
C
      EQUIVALENCE
     *     (EFDU(1,1,1,1),LDDU(1,1,1)),
     *     (EFDU(1,1,1,2),T1DU(1,1,1)),
     *     (EFDU(1,1,1,3),T2DU(1,1,1)),
     *     (EFDU(1,1,1,4),HDDU(1,1,1)),
     *     (EFDU(1,1,1,5),CLDU(1,1,1))
C
C  /EVADU1/:
C
C  DUEQ provides the following coefficients and constants:
C
C  For MOBILE4.1, lines a to c are changed to reflect
C  LDGV and LDGT diurnal rates for pass vehicles.
C
C  DU81 = a + b*UDI + c*UDI**2
C
C                         LDGV                    LDGT1/2
C                CARB     TBI      PFI     CARB     TBI      PFI
      DATA DUEQ/
     a          -4.2904,-0.7233,-0.0679,  -7.6803,-3.8537,-3.8537,
     b           5.9904, 0.0   , 0.0   ,   9.7703, 4.9437, 4.9437,
     c           0.0   , 2.0433, 1.4179,   0.0   , 0.0   , 0.0   /
C
C  Continuation line code for LDDU, T1DU, T2DU, HDDU and CLDU
C  are MY groups.
C
C                      At 9.0 RVP                 At 11.5 RVP
C                       Failed   Failed            Failed   Failed
C                Pass    Purge Pressure     Pass    Purge Pressure
      DATA LDDU/
     1          34.27,   34.27,   34.27,   56.18,   56.18,   56.18,
     2          23.74,   34.27,   34.27,   46.07,   56.18,   56.18,
     3          16.81,   17.71,   24.12,   31.28,   34.96,   38.88,
     4          12.99,   13.77,   20.18,   22.29,   24.98,   28.89,
     *           0.00,    0.00,    0.00,    0.00,    0.00,    0.00/
C
      DATA T1DU/
     1          40.10,   40.10,   40.10,   65.73,   65.73,   65.73,
     2          27.77,   40.10,   40.10,   53.90,   65.73,   65.73,
     3          19.66,   20.72,   28.22,   36.60,   40.91,   45.49,
     4          15.20,   16.11,   23.61,   26.08,   29.22,   33.80,
     *           0.00,    0.00,    0.00,    0.00,    0.00,    0.00/
C
      DATA T2DU/
     1          50.52,   50.52,   50.52,   86.08,   86.08,   86.08,
     2          15.20,   16.11,   23.61,   26.08,   29.22,   33.80,
     *           0.00,    0.00,    0.00,    0.00,    0.00,    0.00,
     *           0.00,    0.00,    0.00,    0.00,    0.00,    0.00,
     *           0.00,    0.00,    0.00,    0.00,    0.00,    0.00/
C
      DATA HDDU/
     1          50.52,   50.52,   50.52,   86.08,   86.08,   86.08,
     *           0.00,    0.00,    0.00,    0.00,    0.00,    0.00,
     *           0.00,    0.00,    0.00,    0.00,    0.00,    0.00,
     *           0.00,    0.00,    0.00,    0.00,    0.00,    0.00,
     *           0.00,    0.00,    0.00,    0.00,    0.00,    0.00/
C
      DATA CLDU/
     1          14.72,    0.00,    0.00,   20.21,    0.00,    0.00,
     2          16.98,    0.00,    0.00,   24.36,    0.00,    0.00,
     3          17.32,    0.00,    0.00,   24.99,    0.00,    0.00,
     4          17.43,    0.00,    0.00,   25.19,    0.00,    0.00,
     5          17.54,    0.00,    0.00,   25.39,    0.00,    0.00/
C
C  /EVADU2/:
C
C                      At 9.0 RVP                 At 11.5 RVP
C                       Failed   Failed            Failed   Failed
C                Pass    Purge Pressure     Pass    Purge Pressure
C
C
      DATA HIDU/
     *          20.81,   24.47,   38.11,   43.78,   55.71,   66.80/
C
C                    Low      |     High
      DATA R2DU/1.9581,2.0677, 2.2086,2.3667/
C
C
      END
