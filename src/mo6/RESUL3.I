C 
C  Changes: (Last change first) 
C 
C  12 Mar 02 AIR Task Bug384: Print diesel evap and All Vehicles 
C  03 Dec 01 AIR Task 27: Increased VFTP size from 3 to 4 to include CO2. 
C 
C The data is initialized in BD27. 
C 
      REAL VFTP, VEXH, VEVAP, VLOSS, VRUNLS, VRSTLS, VIDLE, 
     *     VDRNL, VHTSK, VCRNK           
      COMMON /RESUL3/ VFTP(5),VEXH,VEVAP,VLOSS,VRUNLS,VRSTLS,VIDLE(3), 
     *                VDRNL,VHTSK,VCRNK 
