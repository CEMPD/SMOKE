C Version "@(#)$Id$ $Source$ $Date$ 

C
C  DTGGAS contains the variables that hold the values of dtergent gasoline
C  effects.
C
C  Common Block Dictionary:
C
C   Name    Type                  Description
C  ------   ----  -----------------------------------------------------
C
C  DGPISY    I    Phase-in start year.
C  DGPIEY    I    Phase-in end year.
C  DGPHSF    R    Phase-in percent.
C  FIDFRC    R    Percent of fuel injected vehicles with deposits.
C  FIDEFF    R    Fuel injector deposit effects on HC, CO, and NOx.
C  IVDFRC    R    Percent of vehicles with intake valve deposits.
C  IVDEFF    R    Intake valve deposit effects on HC, CO, and NOx.
C
C  Notes:
C
C  13-Sep-1996 @DynTel-yc Request #2-630 The detergent gas phase-in effect
C               was updated by refining phase-in period
C
      INTEGER DGPISY,DGPIEY
      REAL DGPHSF,FIDFRC,FIDEFF(3),IVDFRC,IVDEFF(3)
C
      COMMON /DTGGAS/ DGPISY,DGPIEY,DGPHSF,FIDFRC,FIDEFF,IVDFRC,IVDEFF
