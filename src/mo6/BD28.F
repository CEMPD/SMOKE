      BLOCK DATA BD28 
C 
C  BLOCK DATA Subprogram 28: /RLCOM1/,/RLCOM2/,/RLCOM3/ 
C 
C  BLOCK DATA 28 defines the refueling loss HC emission factor common blocks. 
C 
C  Changes: (Last change first) 
C 
C 26 Feb 01 @EPA-djb Added HDGV2b to the ORVR arrays (OBIPY & OBPIP). 
C 19 Oct 00 - @EPA-bag    Removed obsolete IOBMY and IVOB 
C   2 Jun 99 @DynTel-MLA 2-694  Removed obsolete variable OBDF. 
C  19 Feb 99 @DynTel-ddj 2-696 Added hourly indices to RLRATE. 
C  11 Aug 98 @DynTel-MLA 2-000  Changed HDG and HDD names to conform to 
C            new naming conventions. 
C 
C  Common block array subscripts: 
C 
C  RLRATE(25,14,24) -  RLRATE ( IDX, IVTAM, IH ) 
C  ROADFE(32,14)    -  ROADFE ( MDX, IVTAM ) 
C  S2EFF(14)        -  S2EFF  ( IVTAM ) 
C 
C  Common block dictionary: 
C 
C   Name   Tlied parameters for refueling loss calculation 
C 
C  ROADFE   R    road fuel economy rates for converting refueling loss HC rates 
C                from grams per gallon (GPG) to grams per mile (GPM). Obtained 
C                from the MOBILE4 Fuel Consumption Model, Default Scenario. ROADFE 
C                has units of miles/gallon. 
C  DISPL    R    HC displacement component of refueling loss, in GPG. 
C  SPILL    R    HC spillage component of refueling loss, in GPG. 
C  OBED     R    Onboard VRS displacement loss reduction efficiency percentage 
C  OBES     R    Onboard VRS spillage loss reduction efficiency percentage 
C  OBPIY    I    Onboard VRS phase-in years 
C  OBPIP    R    Onboard VRS phase-in percentages 
C 
C  /RLCOM2/ user supplied parameters for refueling loss calculation: 
C 
C  IS2SY    I    Stage II vrs controls start year. 
C  NPHASE   I    Number of years in the phase-in period of Stage II vrs 
C                controls. 
C  S2EFF    R    Stage II effciency = percentage of refueling losses removed 
C                by S2 installation, assuming 100% installation. 
C 
C  /RLCOM3/ calculated refueling loss rates for current scenario: 
C 
C  RLRATE   R    refueling loss HC rates for current scenario, as determined by 
C                subprogram REFUEL.  The loss rates are weighted by ICY's 
C                corresponding TF's and summed to get a scenario's REFL_GM_MI. 
C 
C  Notes: 
C 
C    5-May-1994 @ CSC-tmm request 2-411 
C       Added data items to handle phase-in of onboard fuel recovery systems 
C 
C  11 Jun 97 DynTel@MLA 2-622  Modified the arrays to accomodate the new 
C            vehicle types that were introduced in Mobile6. 
C 
C 23 Jul 98 @DynTel-MLA 2-000  Changed LDGT names to conform to 
C           certification naming conventions. 
C 
C 
C  The MDX mapping to MY for ROADFE is: 
C 
C   MDX     MY 
C    1     1952-1969 
C    2-31  1970-1999 (i.e., MDX=MY-1969+1) 
C    32    2000-2051 
C 
C 
      IMPLICIT NONE 
      INCLUDE 'PARAME.I' 
      INCLUDE 'RLCOM1.I' 
      INCLUDE 'RLCOM2.I' 
      INCLUDE 'RLCOM3.I' 
C 
      REAL RDFE1(32,4),RDFE2(32,4),RDFE3(32,4),RDFE4(32,2) 
      EQUIVALENCE (ROADFE(1,1),RDFE1(1,1)),(ROADFE(1,5),RDFE2(1,1)), 
     &            (ROADFE(1,9),RDFE3(1,1)),(ROADFE(1,13),RDFE4(1,1)) 
C 
C  /RLCOM1/: 
C 
C The continuation codes for all the arrays below are the old Mobile5 
C IVTAM codes (1=LDGV, 2=LDGT12, 3=LDGT34, 4=HDGV). The code for the 
C new GAS BUS vehicle type is B. 
C 
      DATA RDFE1/ 
C  LDGV 
     1 12.46, 12.46, 12.06, 12.15, 11.95, 11.95, 13.20, 14.58, 
     1 15.32, 16.55, 16.83, 19.42, 20.74, 21.55, 21.52, 21.85, 
     1 22.47, 23.21, 23.42, 23.74, 23.31, 23.03, 22.73, 22.71, 
     1 22.68, 22.66, 22.64, 22.64, 22.64, 22.64, 22.64, 22.64, 
C  LDGT1 
     2 10.64, 10.64, 10.64, 10.64, 10.45, 10.45, 11.45, 12.05, 
     2 13.09, 12.78, 12.31, 15.30, 16.43, 16.77, 17.31, 17.02, 
     2 17.16, 17.88, 18.05, 17.58, 17.28, 17.29, 17.10, 17.05, 
     2 16.99, 16.93, 16.87, 16.87, 16.87, 16.87, 16.86, 16.86, 
C  LDGT2 
     2 10.64, 10.64, 10.64, 10.64, 10.45, 10.45, 11.45, 12.05, 
     2 13.09, 12.78, 12.31, 15.30, 16.43, 16.77, 17.31, 17.02, 
     2 17.16, 17.88, 18.05, 17.58, 17.28, 17.29, 17.10, 17.05, 
     2 16.99, 16.93, 16.87, 16.87, 16.87, 16.87, 16.86, 16.86, 
C  LDGT3 
     3 10.65, 10.65, 10.65, 10.65, 10.46, 10.46, 11.45, 12.06, 
     3 13.10, 12.78, 12.34, 15.30, 16.37, 16.63, 17.23, 16.98, 
     3 17.15, 17.84, 18.03, 17.54, 17.21, 17.20, 16.95, 16.86, 
     3 16.77, 16.68, 16.58, 16.58, 16.58, 16.58, 16.58, 16.58/ 
      DATA RDFE2/ 
C  LDGT4 
     3 10.65, 10.65, 10.65, 10.65, 10.46, 10.46, 11.45, 12.06, 
     3 13.10, 12.78, 12.34, 15.30, 16.37, 16.63, 17.23, 16.98, 
     3 17.15, 17.84, 18.03, 17.54, 17.21, 17.20, 16.95, 16.86, 
     3 16.77, 16.68, 16.58, 16.58, 16.58, 16.58, 16.58, 16.58, 
C  HDGV2B 
     4  7.76,  7.88,  7.88,  7.88,  7.90,  8.06,  8.22,  8.46, 
     4  8.68,  8.87,  9.01,  9.56,  9.73, 10.07, 10.13, 10.15, 
     4 10.05, 10.22, 10.25, 10.27, 10.27, 10.20, 10.48, 10.49, 
     4 10.50, 10.51, 10.52, 10.52, 10.53, 10.53, 10.53, 10.53, 
C  HDGV3 
     4  7.76,  7.88,  7.88,  7.88,  7.90,  8.06,  8.22,  8.46, 
     4  8.68,  8.87,  9.01,  9.56,  9.73, 10.07, 10.13, 10.15, 
     4 10.05, 10.22, 10.25, 10.27, 10.27, 10.20, 10.48, 10.49, 
     4 10.50, 10.51, 10.52, 10.52, 10.53, 10.53, 10.53, 10.53, 
C  HDGV4 
     4  7.76,  7.88,  7.88,  7.88,  7.90,  8.06,  8.22,  8.46, 
     4  8.68,  8.87,  9.01,  9.56,  9.73, 10.07, 10.13, 10.15, 
     4 10.05, 10.22, 10.25, 10.27, 10.27, 10.20, 10.48, 10.49, 
     4 10.50, 10.51, 10.52, 10.52, 10.53, 10.53, 10.53, 10.53/ 
      DATA RDFE3/ 
C  HDGV5 
     4  7.76,  7.88,  7.88,  7.88,  7.90,  8.06,  8.22,  8.46, 
     4  8.68,  8.87,  9.01,  9.56,  9.73, 10.07, 10.13, 10.15, 
     4 10.05, 10.22, 10.25, 10.27, 10.27, 10.20, 10.48, 10.49, 
     4 10.50, 10.51, 10.52, 10.52, 10.53, 10.53, 10.53, 10.53, 
C  HDGV6 
     4  7.76,  7.88,  7.88,  7.88,  7.90,  8.06,  8.22,  8.46, 
     4  8.68,  8.87,  9.01,  9.56,  9.73, 10.07, 10.13, 10.15, 
     4 10.05, 10.22, 10.25, 10.27, 10.27, 10.20, 10.48, 10.49, 
     4 10.50, 10.51, 10.52, 10.52, 10.53, 10.53, 10.53, 10.53, 
C  HDGV7 
     4  7.76,  7.88,  7.88,  7.88,  7.90,  8.06,  8.22,  8.46, 
     4  8.68,  8.87,  9.01,  9.56,  9.73, 10.07, 10.13, 10.15, 
     4 10.05, 10.22, 10.25, 10.27, 10.27, 10.20, 10.48, 10.49, 
     4 10.50, 10.51, 10.52, 10.52, 10.53, 10.53, 10.53, 10.53, 
C  HDGV8A 
     4  7.76,  7.88,  7.88,  7.88,  7.90,  8.06,  8.22,  8.46, 
     4  8.68,  8.87,  9.01,  9.56,  9.73, 10.07, 10.13, 10.15, 
     4 10.05, 10.22, 10.25, 10.27, 10.27, 10.20, 10.48, 10.49, 
     4 10.50, 10.51, 10.52, 10.52, 10.53, 10.53, 10.53, 10.53/ 
      DATA RDFE4/ 
C  HDGV8B 
     4  7.76,  7.88,  7.88,  7.88,  7.90,  8.06,  8.22,  8.46, 
     4  8.68,  8.87,  9.01,  9.56,  9.73, 10.07, 10.13, 10.15, 
     4 10.05, 10.22, 10.25, 10.27, 10.27, 10.20, 10.48, 10.49, 
     4 10.50, 10.51, 10.52, 10.52, 10.53, 10.53, 10.53, 10.53, 
C  GAS BUS 
     B  7.76,  7.88,  7.88,  7.88,  7.90,  8.06,  8.22,  8.46, 
     B  8.68,  8.87,  9.01,  9.56,  9.73, 10.07, 10.13, 10.15, 
     B 10.05, 10.22, 10.25, 10.27, 10.27, 10.20, 10.48, 10.49, 
     B 10.50, 10.51, 10.52, 10.52, 10.53, 10.53, 10.53, 10.53/ 
C 
      DATA DISPL/5.46/,SPILL/.31/,OBED/.98/,OBES/.50/ 
C 
C    5-May-1994 @ CSC-tmm request 2-411 
C       Added data items to handle phase-in of onboard fuel recovery systems 
C 
      DATA OBPIY/ 
C  LDGV 
     1 1998, 1999, 2000, 
C  LDGT1 
     2 2001, 2002, 2003, 
C  LDGT2 
     3 2001, 2002, 2003, 
C  LDGT3 
     4 2004, 2005, 2006, 
C  LDGT4 
     5 2004, 2005, 2006, 
C  HDGV2b 
     6 2004, 2005, 2006/ 
C 
      DATA OBPIP/ 
C  LDGV 
     1 0.40, 0.80, 1.00, 
C  LDGT1 
     2 0.40, 0.80, 1.00, 
C  LDGT2 
     3 0.40, 0.80, 1.00, 
C  LDGT3 
     4 0.40, 0.80, 1.00, 
C  LDGT4 
     5 0.40, 0.80, 1.00, 
C  HDGV2b 
     6 0.00, 0.80, 1.00/ 
C 
C 
C  /RLCOM2/: 
C 
      DATA IS2SY/NEYR/,NPHASE/4/,S2EFF/14*0.0/ 
C 
C  /RLCOM3/: 
C 
      DATA RLRATE/8400*0.0/ 
C 
      END 
