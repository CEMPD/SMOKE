C Version @(#)$Id$ $Source$ $Date$ 

      BLOCK DATA BD35
C
C  BLOCK DATA Subprogram 35: /BYMYC1/,/BYMYC2/,/BYMYC3/,/BYMYC4/,/BYMYC5/
C
C  BLOCK DATA 35 defines the By Model Year Table generation common blocks.
C
C  Common block array subscripts:
C
C  BYBEF4(3,25,8)  -  BYBEF4 ( IP, JDX, IV )
C  BYCUM(25,8)     -  BYCUM  (     JDX, IV )
C  BYEVAP(25,8)    -  BYEVAP (     JDX, IV )
C  BYFER(3,20,8)   -  BYFER  ( IP, JDX, IV )
C  BYIMCR(3,25,4)  -  BYIMCR ( IPIM, JDX, IVIM )
C  BYIMRE(3,25,4)  -  BYIMRE ( IPIM, JDX, IVIM )
C  BYTAM(3,20,4)   -  BYTAM  ( IP, JDX, IVTAM )
C  BYREFL(25,4)    -  BYREFL (     JDX, IVTAM )
C  BYRSTL(25,8)    -  BYRSTL (     JDX, IV )
C  BYRUNL(25,4)    -  BYRUNL (     JDX, IVTAM )
C  BYVTN(8)        -  BYVTN  (          IV )
C  IVPICK(8)       -  IVPICK (          IV )
C  PN(10)          -  PN     ( ICH )
C  PSEC(3)         -  PSEC   ( IP )
C  REDSUM(3,4)     -  REDSUM ( IPIM, IVIM )
C
C  Common block dictionary:
C
C   Name   Type               Description
C  ------  ----  --------------------------------------------
C  /BYMYC1/:
C
C  IVPICK   I    By Model Year (BYMY) Table vehicle class case selection vector.
C                For each IV,   1 = skip the Table(s)   2 = print the Table(s)
C  IMPICK   I    flag allows user to skip or pick printing of I/M BYMY Table :
C                    1 = skip    2 = print
C                components BYMY Table:
C                    1 = skip    2 = print
C  BYCUM    R    BYMY version of CUMMIL: same except zeroed out where the
C                corresponding TF is zero.  BYCUM is used as a "clean"
C                proxy for CUMMIL in the BYMY tables.  CUMMIL could be used
C                directly if it was recomputed every scenario TF was redone,
C                which is not true when the user enters AMAR.
C
C  /BYMYC2/:
C
C  BYBEF4   R    by model year basic emission rate for given IDX x IV x IV,
C                corrected for RVP, temperature, operating mode and, if
C                applicable, methane and / or I/M: BEF4 = BEF - OMTTAM
C  BYTAM    R    by model year final tampering rates = OMTTAM * RVOPCF
C                (tampering offset corrected for op mode, temp & all RVP cases)
C  BYFER    R    by model year final emission rate for IDX x IV x IP = COMPEF
C
C  /BYMYC3/:
C
C  BYEVAP   R    by model year evaporative HC component = CCEVRT
C  BYRUNL   R    by model year running loss HC component = RULOSS
C  BYREFL   R    by model year refueling loss HC component = RLRATE
C  BYRSTL   R    by model year resting loss HC component = RSTLOS
C
C  /BYMYC4/:
C
C  BYVTN   C*26  full BYMY Table Vehicle Type names
C  YFM     C*17  column titles for model year, travel fraction & cumulative
C                mileage
C  PN      C*7   pollutant names for column headings of BYMY tables
C  U7      C*7   underline of column heading title in BYMY table
C  PSEC    C*3   pollutant names for section headings of BYMY table
C  VB      C*1   vertical bar
C
C  /BYMYC5/:
C
C  BYIMCR   R    by model year I/M credit percentage for given IPIM x JDX x IVIM
C                BYIMCR is applied after the corrections for operating mode,
C                temperature, RVP, and, if applies, methane.
C
C  BYIMRE   R    by model year reduction in gpm in the BEF gained by applying
C                BYIMCR.  BYIMRE is subsequently weighted by SALHCF, TF, and
C                RVOPCF.  I/M FER + BYIMRE = no I/M FER, for each MY.
C
C  REDSUM   R    I/M reduction BYIMRE summed over ICY's 25 MYs, in gpm, with
C                all correction factors applied & weighted by TF.
C
C  Notes:
C
C  BD35 was added for MOBILE4 In-House Version 02 (By Model Year Tables).
C  BD 35 provide the parameters and by model year results save areas.
C
C
      CHARACTER*26 BYVTN
      CHARACTER*17 YFM
      CHARACTER*7 PN,U7
      CHARACTER*3 PSEC
      CHARACTER*1 VB
C
      COMMON /BYMYC1/ IVPICK(8),IMPICK,BYCUM(25,8)
      COMMON /BYMYC2/ BYBEF4(3,25,8),BYTAM(3,25,4),BYFER(3,25,8)
      COMMON /BYMYC3/BYEVAP(25,8),BYRUNL(25,4),BYREFL(25,4),BYRSTL(25,8)
      COMMON /BYMYC4/ BYVTN(8),YFM,PN(10),U7,PSEC(3),VB
      COMMON /BYMYC5/ BYIMCR(3,25,4),BYIMRE(3,25,4),REDSUM(3,4)
C
C  /BYMYC1/:
C
      DATA IVPICK/8*1/,IMPICK/1/,BYCUM/200*0.0/
C
C  /BYMYC2/:
C
      DATA BYBEF4/600*0.0/,BYTAM/300*0.0/,BYFER/600*0.0/
C
C  /BYMYC3/:
C
      DATA BYEVAP,BYRUNL,BYREFL,BYRSTL/600*0.0/
C
C
C  /BYMYC4/:
C
C
      DATA BYVTN/
     1            ' Light Duty Gas Vehicles  ',
     2            ' Light Duty Gas Trucks 1  ',
     3            ' Light Duty Gas Trucks 2  ',
     4            ' Heavy Duty Gas Vehicles  ',
     5            'Light Duty Diesel Vehicles',
     6            ' Light Duty Diesel Trucks ',
     7            'Heavy Duty Diesel Vehicles',
     8            '       Motorcycles        '/
C
      DATA YFM/'Year   TF   Miles'/
      DATA PN/' BEF4  ','Tamper ','SALHCF ','Evapor ','Refuel ',
     *        'Runnin ','Restin ','  FER  ','Credit ','Reduc  '/
      DATA U7/' ______'/
      DATA PSEC/'HC ','CO ','NOX'/
      DATA VB/'|'/
C
C  /BYMYC5/:
C
      DATA BYIMCR/300*0.0/,BYIMRE/300*0.0/,REDSUM/12*0.0/
C
      END
