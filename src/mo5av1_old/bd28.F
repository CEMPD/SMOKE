C Version @(#)$Id$ $Source$ $Date$ 

      BLOCK DATA BD28
C
C  BLOCK DATA Subprogram 28: /RLCOM1/,/RLCOM2/,/RLCOM3/
C
C  BLOCK DATA 28 defines the refueling loss HC emission factor common blocks.
C
C  Common block array subscripts:
C
C  IVOB(4)       -  IVOB   ( IVTAM )
C  OBDF(4)       -  OBDF   ( IVTAM )
C  RLRATE(25,4)  -  RLRATE ( IDX, IVTAM )
C  ROADFE(32,4)  -  ROADFE ( MDX, IVTAM )
C
C  Common block dictionary:
C
C   Name   Tlied parameters for refueling loss calculation
C
C  ROADFE   R    road fuel economy rates for converting refueling loss HC rates
C                from grams per gallon (GPG) to grams per mile (GPM). Obtained
C                from the MOBILE4 Fuel Consumption Model, Default Scenario.
C  DISPL    R    HC displacement component of refueling loss, in GPG.
C  SPILL    R    HC spillage component of refueling loss, in GPG.
C  OBED     R    Onboard VRS displacement loss reduction efficiency percentage
C  OBES     R    Onboard VRS spillage loss reduction efficiency percentage
C  OBDF     R    Onboard VRS deterioration factor
C
C  /RLCOM2/ user supplied parameters for refueling loss calculation:
C
C  IOBMY    I    first model year of onboard vrs installation.  All gas vehicles
C                produced in this and all subsequent years are assumed to have
C                onboard vapor recovery systems installed.
C  IVOB     I    switch vector controlling which tampering vehicle classes
C                are to be affected by ob vrs's this run:
C                  IVOB(IVTAM) = 1 (no) or 2 (yes), IVTAM=1,2,3,4
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
C                corresponding TF's and summed to get a scenario's EFREFL.
C
C  Notes:
C
C  The MDX mapping to MY for ROADFE is:
C
C   MDX     MY
C    1     1941-1969
C    2-31  1970-1999 (i.e., MDX=MY-1969+1)
C    32    2000-2020
C
C
      COMMON /RLCOM1/ ROADFE(32,4),DISPL,SPILL,OBED,OBES,OBDF(4)
      COMMON /RLCOM2/ IOBMY,IVOB(4),IS2SY,NPHASE,S2EFF(4)
      COMMON /RLCOM3/ RLRATE(25,4)
C
C  /RLCOM1/:
C
      DATA ROADFE/
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
     3 10.65, 10.65, 10.65, 10.65, 10.46, 10.46, 11.45, 12.06,
     3 13.10, 12.78, 12.34, 15.30, 16.37, 16.63, 17.23, 16.98,
     3 17.15, 17.84, 18.03, 17.54, 17.21, 17.20, 16.95, 16.86,
     3 16.77, 16.68, 16.58, 16.58, 16.58, 16.58, 16.58, 16.58,
C  HDGV
     4  7.76,  7.88,  7.88,  7.88,  7.90,  8.06,  8.22,  8.46,
     4  8.68,  8.87,  9.01,  9.56,  9.73, 10.07, 10.13, 10.15,
     4 10.05, 10.22, 10.25, 10.27, 10.27, 10.20, 10.48, 10.49,
     4 10.50, 10.51, 10.52, 10.52, 10.53, 10.53, 10.53, 10.53/
C
C  Suggested nonzero OBDF values are .0002 for LDGV and .0005 LDGT1/2 & HDGV.
C
      DATA DISPL/5.46/,SPILL/.31/,OBED/.98/,OBES/.50/,OBDF/4*0.0/
C
C  /RLCOM2/:
C
      DATA IOBMY/2020/,IVOB/4*1/,IS2SY/2020/,NPHASE/4/,S2EFF/4*0.0/
C
C  /RLCOM3/:
C
      DATA RLRATE/100*0.0/
C
      END
