C Version @(#)$Id$ $Source$ $Date$ 

      BLOCK DATA BD32
C
C  BLOCK DATA Subprogram 32: /IMPAR4/,/IMPAR5/,/IMPAR6/,
C  /IMPAR7/, and /IMPAR8/
C
C  BLOCK DATA 32 initializes I/M program parameters.
C
C  Common block array subscripts:
C
C  CRHDGV(2,2)    -  CRHDGV ( IP, IFREQ )
C  DISCNT(3)      -  DISCNT ( INTYP )
C  MYGIM2(3,3)    -  MYGIM2 ( IG, IVLDG )
C  CUTHC(2)       -  CUTHC  ( IMPNUM )
C  CUTCO(2)       -  CUTCO  ( IMPNUM )
C  CUTNO(2)       -  CUTNO  ( IMPNUM )
C  IFREQ(2)       -  IFREQ  ( IMPNUM )
C  INTYP(2)       -  INTYP  ( IMPNUM )
C  IFILET(1)      -  IFILET ( <TECH1&2> )
C  IFILET(2)      -  IFILET ( <TECH4+> )
C
C  Common block dictionary:
C
C   Name   Type              Description
C  ------  ----  --------------------------------------------------
C
C  /IMPAR4/
C
C  CUTHC    R    Cut point to be found in the I/M credit data file or use
C                default HC
C  CUTCO    R    Cut point to be found in the I/M credit data file or use
C                default CO
C  CUTNO    R    Cut point to be found in the I/M credit data file or use
C                default NOx
C
C  /IMPAR5/:
C
C  MYGIM2   I    mygs for IM technologies 1 and 2 (3 is treated as 2, for I/M)
C  LBIM4P   I    lb of last myg of the tech 4+ class (its ub is 2020)
C  CRHDGV   R    BEF I/M credit adjustment for HDGV for HC & CO, annu/biennial
C  DISCNT   R    I/M program efficiency discount as a f(INTYP)
C
C  /IMPAR6/:
C
C  IFREQ    I    I/M Program vehicle inspection frequency: 1=annual,2=biennial
C  INTYP    I    I/M Program type: 1 = Test only,
C                2 = Test and Repair (Computerized),
C                3 = Test and Repair (manual)
C  IFILET   I    Tells MOBILE5 whether or not to use default I/M credit data
C                filenames where 1=default 2=user inputs
C
C  /IMPAR7/
C
C  IREAD    I    Tells MOBILE5 that 1=low altitude has been read or 2=high
C                and 0=no data have been read
C
C /IMPAR8/
C
C  ICALC    I    Tells ATPAER which rates to use with or with out IM240.
C                1=with out an IM240 program and 2=with IM240.
C
C
C  Notes:
C
C  IFILET was added for version 2 of MOBILE5 to allow
C  default filenames in I/M calculations
C  and they are stored in arrays CUTHC, CUTCO, and CUTNO
C  Also /IMPAR7/ was added to indicate whether MOBILE5 has read
C  any I/M credits yet.
C
C
      CHARACTER*40 TNAME
      COMMON /IMPAR4/ CUTHC(2),CUTCO(2),CUTNO(2)
      COMMON /IMPAR5/ MYGIM2(3,3),LBIM4P,CRHDGV(2,2),DISCNT(3)
      COMMON /IMPAR6/ IFREQ(2),INTYP(2),IFILET(2)
      COMMON /IMPAR7/ TNAME,IREAD
      COMMON /IMPAR8/ ICALC
C
C  /IMPAR4/
C
      DATA CUTHC/2*220./, CUTCO/2*1.2/, CUTNO/2*999./
C
C  /IMPAR5/
C
C  MYGIM2 continuation code is IVLDG.
C                       TECH
C                   1     2     3
      DATA MYGIM2/
     1            1974, 1979, 1980,
     2            1974, 1980, 1983,
     3            1978, 1980, 1983/
C
      DATA LBIM4P/1998/,CRHDGV/.1250,.1575, .0625,.1050/
C
      DATA DISCNT/1.0, 0.50, 0.50/
C
C  /IMPAR8/
C
      DATA ICALC/1/
C
C  /IMPAR7/
C
      DATA IREAD/0/
C
C  /IMPAR6/
C
      DATA IFREQ,INTYP,IFILET/6*1/
C
      END
