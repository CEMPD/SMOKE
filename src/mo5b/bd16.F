C Version "@(#)$Id$ $Source$ $Date$ 

      BLOCK DATA BD16
C
C  BLOCK DATA Subprogram 16: /IOUCOM/,/IVPCOM/,/MAXIMA/,/QUITXQ/,/YEARS4/
C
C  BLOCK DATA 16 initializes miscellaneous integer-valued parameters.
C
C  Common block array subscripts:
C
C  IVPTRL(8)  -  IVPTRL ( IVLM )
C  IVPTRT(8)  -  IVPTRT ( IVTEMP )
C  IVPTRA(8)  -  IVPTRA ( IVA )
C  IVPTRB(8)  -  IVPTRB ( IVB )
C
C  Common block dictionary:
C
C  /IOUCOM/: a special format is utilized for the i/o unit number (#)
C          COMMON block:
C
C  Unit#   Usage                     Description
C  ------  ------  ----------------------------------------------------------
C  IOUIMD  input   user supplied I/M credits (when IMFLAG=2)
C  IOUGEN  input   user supplied general input data
C  IOUREP  output  formatted report file
C  IOUERR  output  error/warning messages file
C  IOUASK  output  prompting messages file
C
C   Name   Type               Description
C  ------  ----  --------------------------------------------
C  /IVPCOM/:
C  IVPTRL   I    IV subgroup pointer into MYGSP1 & BFRCOF (see IVLM subscript
C                definition)
C  IVPTRT   I    IV subgroup pointer into MYGTCF (see IVTEMP subscript
C                definition)
C  IVPTRA   I    IV subgroup pointer into ALHRET (see IVA
C                subscript definition)
C  IVPTRB   I    IV subgroup pointer into GPBSCO (see IVB subscript
C                definition)
C
C  /QUITXQ/:
C
C
C  N1QUIT   I    switch: 1 = first QUITER call has been made, 0 = it has not
C                On first call, extra vertical spacing is inserted between the
C                run title and the QUITER message.
C
C  /MAXIMA/:
C
C  MAXVEH   I    maximum number of vehicle types - additional categories
C                may be created by combining 2 or more of the base 8 types
C                (LDGT = LDGT1 + LDGT2)
C  MAXLTW   I    maximum number of vehicle types that can be affected by
C                BIGALH's optional correction factors for a/c, extra load
C                and trailer towing ( = 3 = LDGV, LDGT1 & LDGT2 )
C  MAXPOL   I    maximum number of pollutants covered by MOBILE4
C                ( = 3 = HC, CO & NOx )
C  MAXREG   I    maximum number of regions allowed in MOBILE4
C                ( = 3 = low, high & mid )
C  MAXYRS   I    maximum number of model years in ICY's fleet (1 thru 20+)
C
C  /YEARS4/:
C
C  IYSTART  I    4-digit representation of start year
C  IY1960   I    4-digit representation of year 1960
C  IY1994   I    4-digit representation of year 1994
C  IYEND    I    4-digit representation of end year 
C
C  Notes:
C
C  None.
C
C
      INCLUDE 'IOUCOM.I'
      INCLUDE 'IVPCOM.I'
      INCLUDE 'MAXIMA.I'
      INCLUDE 'PARAME.I'
      INCLUDE 'PARAMS.I'
      INCLUDE 'QUITXQ.I'
      INCLUDE 'YEARS4.I'
C
C Change for different computers
C on the MAC /4,9,2*6,9/ on the IBM /4,5,2*6,0/
C and on MTS /4,1,2*7,6/
C
      DATA IOUIMD,IOUGEN,IOUREP,IOUERR,IOUASK/4,5,2*6,0/
C
      DATA IVPTRL/1,2,3,0,4,5,0,6/
C
      DATA IVPTRT/1,2,3,4,0,0,0,5/
C
      DATA IVPTRA/1,2,3,0,0,0,0,4/
C
      DATA IVPTRB/0,0,0,1,2,3,4,0/
C
      DATA MAXVEH,MAXLTW,MAXPOL,MAXREG,MAXYRS/8,3,3,2,25/
C
      DATA N1QUIT/1/
C
      DATA IYSTART,IY1960,IY1994,IYEND/NSYR,1960,1994,NEYR/
C
      END

