      SUBROUTINE PMGASSL(SPD,SO4,SO2) 
C 
C  Called by: EFCALXP5 
C 
C Changes: 
C 
C 13 Feb 01 AIR Task 11: New MOBILE6 routine. 
C 13 Feb 01 AIR Task 11: Updated to use MOBILE6 vehicle classes, 
C           fuel economies, and tech fractiosn. Also added fuel 
C           sulfur correction to calculations. 
C 17 Nov 00 AIR Task 06: Removed Dead Code. 
C 10 Nov 00 AIR Task 06: Updated code to MOBILE6 standards. 
C 
C  Calls PMTAMPER, PMLEDPTR 
C 
C  Input on call: 
C    parameter list: SPD 
C 
C    common blocks: See Includes below 
C 
C  Output on return: SO4,SO2 
C 
C  Local variable dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------------- 
C  DSLFCH   R    34.8 MPH or faster D-sulf for catalyst unleaded vehicles 
C  DSLFCL   R    19.6 MPH or slower D-sulf for catalyst unleaded vehicles 
C  DSULFC   R    Direct sulfate emission factor for catalyst vehicles 
C  DSULFN   R    Direct sulfate emission factor for noncatalyst vehicles 
C  FCNVRC   R    Estimated percent of sulfur that is converted to sulfate in a catalyst vehicle 
C  FCNVRN   R    Estimated percent of sulfur that is converted to sulfate in a noncatalyst vehicle 
C  SO2      R    Gaseous sulfur emission factor 
C  SO2C     R    Gaseous sulfur emissions from catalyst vehicles 
C  SO2N     R    Gaseous sulfur emissions from noncatalyst vehicles 
C  SO4      R    Direct sulfate emission factor for fleet average gasoline vehicle 
C  SPD      R    Speed, mph 
C 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'FRACT.I'  ! WATER,FDNSTY 
      INCLUDE 'PART2.I'  ! BSWGTR,CTLFRC,FECONM,FRACA,GSWGTR,SULFCF 
C 
C  Declare parameter list. 
C 
      REAL    SO4 
      REAL    SO2 
      REAL    SPD 
C 
C  Declare local variables/arrays. 
C 
      REAL    DSLFCH 
      REAL    DSLFCL 
      REAL    DSULFC 
      REAL    DSULFN 
      REAL    FCNVRC 
      REAL    FCNVRN 
      REAL    SO2C 
      REAL    SO2N 
C 
C     Compute the grams/mile emission factor of the direct sulfate material: 
C  
C     13.6078 = (453.592 * 3.)/100. 
C     where: 453.592 = number of grams in a pound 
C                  3 = weight ratio of SO4 to sulfur 
C     and the division by 100 is to correct for GSWGTR the gasoline weight  
C     ratio of sulfur content in fuel. 
C 
C     WATER  = weight ratio of seven water molecules to sulfate  
C              (see p. 6-25 of "Regulatory Impact Analysis, 
C              Control of Sulfur and Aromatics Contents of On- 
C              Highway Diesel Fuel") 
C     FDNSTY = fuel density in lbs/gallon, 1=gasoline, 2=diesel 
C     FECONM = fuel economy 
C     GSWGTR = gasoline weight ratio of sulfur content in fuel 
C     DCNVRT  = percent of sulfur exhausted that is directly converted 
C               into sulfuric acid. 
C 
C   
C  Equation 1:  This is the equation that is used to calculate direct  
C      sulfate as a function of the sulfur content and fuel economy:  
C      (used for calculating the fraction of sulfur in the fuel that 
C      gets converted to sulfate, DCNVRT. 
C      DCSULF = 13.6078*(1.+ WATER)*FDNSTY(1)*GSWGTU*DCNVRT/FECONMC 
C 
C  Programming Note: PART5 had a special method to compute SO4 for HDGV. 
C  This method was replaced with the PART5 LDGV/LDGT method. 
C 
        IF (SPD.LE.19.6) THEN 
C 
C 
          DSULFN = .002 
          DSULFC =  
     2           ((FRACA(1)+FRACA(3))*.005 
     3           +(FRACA(2)+FRACA(4))*.016)     
C 
        ELSE IF (SPD.GE.34.8) THEN 
          DSULFN = .001 
          DSULFC =   
     2           (FRACA(1)* .005 
     3           +FRACA(3)* .001 
     4           +FRACA(2)* .020 
     5           +FRACA(4)* .025) 
        ELSE 
C 
C  Interpolate for speeds between 19.6 and 34.8 
C 
          DSULFN = ((.002-.001)*(SPD-19.6)/(19.6-34.8)+.002) 
          DSLFCL =  
     2           ((FRACA(1)+FRACA(3))*.005 
     3           +(FRACA(2)+FRACA(4))*.016) 
          DSLFCH =   
     2           (FRACA(1)* .005 
     3           +FRACA(3)* .001  
     4           +FRACA(2)* .020  
     5           +FRACA(4)* .025) 
          DSULFC = (DSLFCL-DSLFCH)*(SPD-19.6)/(19.6-34.8)+DSLFCL 
C 
C  end of interpolation section 
C 
        ENDIF 
C 
      SO4 = CTLFRC*DSULFC + (1. - CTLFRC)*DSULFN 
C 
C  Apply fuel sulfur correction to SO4 
C 
      SO4 = SO4 * SULFCF 
C 
C     Calculate the fraction of sulfur in the fuel that was converted 
C     to sulfate (based on Equation 1 above) for catalyst and non- 
C     catalyst vehicles. Since the DSULFC and DSULFN values are based 
C     on 340 ppm (0.034 wt%) fuel, that sulfur content must be used here. 
C 
      FCNVRC = DSULFC*FECONM/( 13.6078*(1.+WATER)*FDNSTY(1)*BSWGTR ) 
      FCNVRN = DSULFN*FECONM/( 13.6078*(1.+WATER)*FDNSTY(1)*BSWGTR ) 
C 
C     Compute the exhausted SO2: 
C 
C     9.072 = (453.592 * 2.)/100. 
C     where: 453.592 = number of grams in a pound 
C                 2. = weight ratio of SO2 to sulfur 
C     and the division by 100 is to correct for BSWGTR the gasoline weight 
C     ratio of sulfur content in the fuel. 
C 
      SO2C = 9.072*FDNSTY(1)*BSWGTR*(1.-FCNVRC)/FECONM 
      SO2N = 9.072*FDNSTY(1)*BSWGTR*(1.-FCNVRN)/FECONM 
C 
      SO2 = CTLFRC*SO2C + (1.-CTLFRC)*SO2N 
C 
C  Apply fuel sulfur correction to SO2 
C 
      SO2 = SO2 * SULFCF 
C 
      RETURN 
      END 
       
