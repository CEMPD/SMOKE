C 
C  DTGGAS contains the variables that hold the values of detergent gasoline 
C  effects. 
C 
C  Changes: (Last change first) 
C 
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
C  The data is initialized in GETDTG.FOR 
C 
C  Notes: 
C 
C  25 Mar 98 @DynTel-MLA 2-659  In Mobile6, the initialization of this 
C            data was moved from RDCITY to INISCN. Alternate user input 
C            is now read in GETDTG instead of RDCITY. 
C 
C  18 May 98 @DynTel-MLA 2-659  Added global parameters for DTGGAS. 
C 
C 
      INTEGER DGPISY,DGPIEY,GLDGPEY,GLDGPSY 
      REAL DGPHSF,FIDFRC,FIDEFF(3),IVDFRC,IVDEFF(3) 
      REAL GLFIDFC,GLFIDEF(3),GLIVDFC,GLIVDEF(3),GLDGPHSF 
C 
      COMMON /DTGGAS/ DGPISY,DGPIEY,DGPHSF,FIDFRC,FIDEFF,IVDFRC, 
     &       IVDEFF,GLDGPSY,GLDGPEY,GLDGPHSF,GLFIDEF,GLFIDFC, 
     &       GLIVDEF,GLIVDFC 
