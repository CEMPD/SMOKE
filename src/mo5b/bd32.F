C Version "@(#)$Id$ $Source$ $Date$ 

      BLOCK DATA BD32
C
C  BLOCK DATA Subprogram 32: /IMPAR4/,/IMPAR5/,/IMPAR6/,
C  /IMPAR7/,/PPEFF1/
C
C  BLOCK DATA 32 initializes I/M program parameters.
C
C  Common block array subscripts:
C
C  CRHDGV(3,2)    -  CRHDGV ( IP, IFREQ )
C  DCNTSV(3, 3)   -  DCNTSV ( IP, INTYP )
C  DISCNT(3, 3)   -  DISCNT ( IP, INTYP )
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
C  LBIM4P   I    lb of last myg of the tech 4+ class (its ub is 2051)
C  CRHDGV   R    HDGV BEF I/M credit adjustment (HC, CO, NOx) annual/biennial
C  DCNTSV   R    Saves the default values for DISCNT.
C  DISCNT   R    I/M program efficiency discount as a f(IP, INTYP)
C
C  /IMPAR6/:
C
C  IFREQ    I    I/M Program vehicle inspection frequency: 1=annual,2=biennial
C  INTYP    I    I/M Program type: 1 = Test only,
C                2 = Test and Repair (Computerized),
C                3 = Test and Repair (manual)
C                4 = Retest Based Hybrid Program 
C  IFILET   I    Tells MOBILE5 whether or not to use default I/M credit data
C                filenames where 1=default 2=user inputs
C
C  /IMPAR7/
C
C  IREAD    I    Tells MOBILE5 that 1=low altitude has been read or 2=high
C                and 0=no data have been read
C  SVALTC   I    0=No I/M credits have been read. 1=Standard I/M credits
C                have been read. 2=Alternate I/M credits have been read.
C
C  /PPEFFD/
C  PPEFF1    R    holds the default values of the credits for pressure and 
C                 purge test effectiveness
C  PPEFF     R    holds the values that have to be for the current run.
C
C
C  Notes:
C
C  IFILET was added for version 2 of MOBILE5 to allow
C  default filenames in I/M calculations
C  and they are stored in arrays CUTHC, CUTCO, and CUTNO
C  Also /IMPAR7/ was added to indicate whether MOBILE5 has read
C  any I/M credits yet.
C  May 18, 1993 @ ARC-bsg Subtask 238, the IMPAR8 common was
C  removed, which contained the variable ICALC.
C  June-15-1993 @ ARC-bk Subtask 244 the common block IMPAR5
C  was modified, the HDGV I/M credit was expanded to include
C  an initial value for NOx.
C
CC  May-1995 @ CSC-zp Req. 2-495 "Add MOB5a-H to MOBILE5b"
C                Initializes values of the credits for pressure/purge
C                test effectiveness PPEFF1(25,2)
C
C  1 Feb 96 @DynTel-MLA 2-611 Added IP dimension to DISCNT. Added a new
C     array, DCNTSV, with the same dismensions as DISCNT. (See IMPAR5.)
C
C  5 Jun 96 @DynTel-MLA 2-620 Added variable SVALTC to IMPAR7.
C
      INCLUDE 'IMPAR4.I'
      INCLUDE 'IMPAR5.I'
      INCLUDE 'IMPAR6.I'
      INCLUDE 'IMPAR7.I'
      INCLUDE 'PPEFFD.I'
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
      DATA LBIM4P/1998/,CRHDGV/.1250,.1575,0.00, .0625,.1050,0.00/
C
C     We now need separate I/M effectiveness rates for each I/M test
C     type (1=Test only, 2=Test and repair-computerized, and 3=Test
C     and repair-manual) and for each pollutant. Users are now allowed
C     to enter effectiveness rates for each pollutant for test-and-repair
C     programs.
C
C     DCNTSV holds the default effectiveness rates. This array is used to
C     reinitialize DISCNT. In each row, there are values for each pollutant.
C     The continuation codes for each table is the I/M test type.
CC
      DATA DCNTSV/1.00, 1.00, 1.00,
     2            0.50, 0.50, 0.50,
     3            0.50, 0.50, 0.50/
C
      DATA DISCNT/1.00, 1.00, 1.00,
     2            0.50, 0.50, 0.50,
     3            0.50, 0.50, 0.50/
C
C  /IMPAR7/
C
      DATA IREAD,SVALTC/2*0/
C
C  /IMPAR6/
C
      DATA IFREQ,INTYP,IFILET/6*1/
C
C  /PPEFFD/
C 
      DATA PPEFF1 /
     * .187,.412,.947,.916,.740,.724,.787,.831,.862,.882,
     * .898,.911,.958,.937,.920,.925,.933,.939,.943,.948,
     * .952,.980,.964,.961,.960,
     * .187,.412,.947,.916,.727,.723,.575,.822,.711,.875,
     * .788,.906,.875,.933,.908,.888,.814,.899,.834,.914,
     * .858,.951,.926,.975,.953 /
C  
      END
