C 
C  Changes: (Last change first) 
C 
C 15 Sep 00 AIR Task 08  Increase MAXIMPAR to 28. 
C 25 Aug 00 AIR Task 08: Increased MAXIMPGM to 7 to handle EVAP I/M and 
C                        MAXIMPAR to handle Evap Fractional Benefits. 
C 13 Apr 98 AIR Task 02: Removed ALTITUDE from CHKLIST and set MAXREQ 
C                        to 3. ALTITUDE is by default set to Low (1) 
C 04 Apr 00 AIR Task 02: Changed REGION to ALTITUDE, removed ASTM CLASS 
C                        and set MAXREQ from 5 to 4. 
C 26 May 99 @DynTel-HXQ 2-698  Change MAXREQ from 7 to 5 
C 08 Apr 98 @DynTel-bag 2-000  Removed ambient temperature from CHKLIST. 
C 10 Dec 98 @DynTel-MLA 2-663  Expanded IMCHECK (added SIM ONLY case). 
C 
C 27 Oct 98 @DynTel-MLA 2-663  Added IMCHECK, MAXIMPAR, and MAXIMPGM 
C           to support I/M programs in Mobile6. 
C 
C The variables in this file keep track of required user inputs. 
C 
C This include file is new for Mobile6. 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------ 
C  /CHKLST/: 
C 
C  CHKLIST   I   Checklist of required user inputs (general inputs) 
C  GLCHECK   I   'Global data' copy of CHKLIST 
C  IMCHECK   I   Checklist of inputs for I/M programs 
C  MAXIMPAR  I   Number of items on the I/M checklist 
C  MAXIMPGM  I   Maximum number of I/M programs allowed 
C  MAXREQ    I   Number of required user inputs (general inputs) 
C 
C  CHKLST(1) = Min/max temperatures 
C  CHKLST(2) = Fuel RVP 
C  CHKLST(3) = Calendar year 
C 
C  IMCHECK(1,n)  = Was any user data entered for I/M program n? 
C  IMCHECK(2,n)  = Program start year 
C  IMCHECK(3,n)  = Program end year 
C  IMCHECK(4,n)  = Program frequency  
C  IMCHECK(5,n)  = Program type (test and repair or test only,n) 
C  IMCHECK(6,n)  = Test type (Idle, IM240, ...,n) 
C  IMCHECK(7,n)  = Program cutpoints 
C  IMCHECK(8,n)  = First and last model years covered 
C  IMCHECK(9,n)  = Vehicles covered by program 
C  IMCHECK(10,n) = Program stringency 
C  IMCHECK(11,n) = Program compliance rate 
C  IMCHECK(12,n) = Program waiver rates 
C  IMCHECK(13,n) = Number of failed RSD tests required before I/M 
C  IMCHECK(14,n) = Number of RSD tests administered 
C  IMCHECK(15,n) = Number of vehicles in RSD catchment area. 
C  IMCHECK(16,n) = RSD coverage  
C  IMCHECK(17,n) = Number of vehicles failing RSD 
C  IMCHECK(18,n) = Number of vehicles in RSD catchment area 
C  IMCHECK(19,n) = CO cutpoint for RSD test (used by dirty screening) 
C  IMCHECK(20,n) = NOX cutpoint for RSD test (used by clean screening) 
C  IMCHECK(21,n) = RSD clean screening losses 
C  IMCHECK(22,n) = COIM (either COIM PCT or INCLUDE COIM) 
C  IMCHECK(23,n) = Selective I/M -- RSD ONLY 
C  IMCHECK(24,n) = Selective I/M -- COIM ONLY 
C  IMCHECK(25,n) = I/M grace period (for warning message only) 
C  IMCHECK(26,n) = Evap I/M Fractional Benefits 
C  IMCHECK(27,n) = Waiver Reduction 
C  IMCHECK(28,n) = Exempt Age 
C   
C 
      INTEGER CHKLIST,GLCHECK,IMCHECK,MAXIMPGM,MAXIMPAR,MAXREQ 
C 
      PARAMETER (MAXREQ = 3, MAXIMPAR = 28, MAXIMPGM = 7) 
C 
      COMMON /CHKLST/ CHKLIST(MAXREQ),GLCHECK(MAXREQ), 
     &                IMCHECK(MAXIMPAR,MAXIMPGM) 
 
