      BLOCK DATA BD_IM2 
C 
C  BLOCK DATA Subprogram IM2: /IMPAR4/,/IMPAR5/,/IMPAR6/, 
C  /IMPAR7/,/PPEFF1/ 
C 
C  BLOCK DATA IM2 initializes I/M program parameters. 
C 
C  Changes: (Last change first) 
C   
C 
C  15 Sep 00 AIR Task 08: Increased the size of CUTHC, CUTCO, CUTNO, IFREQ, 
C            INTYP to 7 to accomodate EVAP I/M. 
C  14 Sep 00 @EPA-djb  Bug 251, Initialize ATPUEF and PPGUEF. 
C  29 Aug 00 @EPA-elg  bug241 Initialize FUELCF from new IMBLOCK.I 
C  07 Jul 00 @EPA-elg  Initialize CUTHC, CUTCO and CUTNO to final IM240 cutpoints. 
C  26 Jun 00 @EPA-elg  Changed default I/M discount rates to 100% - No discount. 
C  01 Oct 99 @Dyntel-JWR 2-705 Removed data statement that incorrectly 
C            initialized cutpoints. The initialization is carried out 
C            in REINIT.FOR  
C  03 Aug 99 @DynTel-BG  2-663  Renamed BD32.FOR to BD_IM2.FOR 
C  16 Dec 98 @DynTel-MLA 2-663  Added IVCUT to map IV into an index for 
C            the cutpoint arrays CUTHC, CUTCO, and CUTNOX. 
C  19 Nov 98 @DynTel-MLA 2-663  Removed references to Retest Based Hybrid 
C            programs, since they are no longer supported in Mobile6. 
C 
C  Common block array subscripts: 
C 
C  CRHDGV(3,2)    -  CRHDGV ( IP, IFREQ ) 
C  DCNTSV(3, 3)   -  DCNTSV ( IP, INTYP ) 
C  DISCNT(3, 3)   -  DISCNT ( IP, INTYP ) 
C  MYGIM2(3,5)    -  MYGIM2 ( IG, IVLDG ) 
C  CUTHC(25,4,7)  -  CUTHC  ( IMPNUM, IVCUT, IMPGM ) 
C  CUTCO(25,4,7)  -  CUTCO  ( IMPNUM, IVCUT, IMPGM ) 
C  CUTNO(25,4,7)  -  CUTNO  ( IMPNUM, IVCUT, IMPGM ) 
C  IFREQ(7)       -  IFREQ  ( IMPNUM ) 
C  INTYP(7)       -  INTYP  ( IMPNUM ) 
C  IFILET(5)      -  IFILET ( <TECH1&2> ) 
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
C  IVCUT    I    Pointer array maps IV to index for CUTHC, CUTCO, CUTNOX. 
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
C  IFREQ    I    I/M Program vehicle inspection frequency: 1=annual,2=biennial,3=triennial 
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
C  SVALTC   I    0=No I/M credits have been read. 1=Standard I/M credits 
C                have been read. 2=Alternate I/M credits have been read. 
C 
C  /PPEFFD/ 
C  PPEFF1    R    holds the default values of the credits for pressure and  
C                 purge test effectiveness 
C  PPEFF     R    holds the values that have to be for the current run. 
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
C  1 Feb 96 DynTel@MLA 2-611 Added IP dimension to DISCNT. Added a new 
C     array, DCNTSV, with the same dismensions as DISCNT. (See IMPAR5.) 
C 
C  5 Jun 96 DynTel@MLA 2-620 Added variable SVALTC to IMPAR7. 
C 
C  3 Dec 96 DynTel@MLA 2-622  Expanded the arrays to accomodate the new 
C           LDGT vehicle categories. 
C 
C 23 Jul 98 @DynTel-MLA 2-000  Changed LDGT names to conform to 
C           certification naming conventions. 
C 
      IMPLICIT NONE 
      INCLUDE 'IMPAR4.I' 
      INCLUDE 'IMPAR5.I' 
      INCLUDE 'IMPAR6.I' 
      INCLUDE 'IMPAR7.I' 
      INCLUDE 'PPEFFD.I' 
      INCLUDE 'IMBLOCK.I' 
C 
C  /IMPAR4/ 
C 
 
      DATA IVCUTIDX /1,1,2,2,3,4,4,4,4,4,4,4,4, 
     &               0,0,0,0,0,0,0,0,0,0,0,4,0,0/ 
C 
C  /IMPAR5/ 
C 
C  MYGIM2 continuation code is the old Mobile5 IVLDG code (1=LDGV, 
C  2=LDGT12, 3=LDGT34). 
C                       TECH 
C                   1     2     3 
      DATA MYGIM2/ 
     1            1974, 1979, 1980, 
     2            1974, 1980, 1983, 
     2            1974, 1980, 1983, 
     3            1978, 1980, 1983, 
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
C     The continuation codes for each table are the I/M test types. 
CC 
      DATA DCNTSV/1.00, 1.00, 1.00, 
     2            1.00, 1.00, 1.00, 
     3            1.00, 1.00, 1.00/ 
C 
      DATA DISCNT/1.00, 1.00, 1.00, 
     2            1.00, 1.00, 1.00, 
     3            1.00, 1.00, 1.00/ 
C 
      DATA ATPUEF / 0.0 / 
      DATA PPGUEF / 0.5 / 
C 
C  /IMPAR7/ 
C 
      DATA IREAD,SVALTC/2*0/ 
C 
C  /IMPAR6/ 
C 
      DATA IFREQ,INTYP,IFILET/15*1/ 
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
C  /LDGIMC/ 
C 
      DATA CUTHC / 700 * 0.8/ 
      DATA CUTCO / 700 * 15.0/ 
      DATA CUTNO / 700 * 2.0/ 
C 
C  /IMBLOCK1/ 
C 
      DATA FUELCF /2880 * 1.0/ 
C 
      END 
