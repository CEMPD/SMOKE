C 
C  Changes: (Last change first) 
C 
C  12 Sep 00 @EPA-djb Bug 237, Changed GLVMTMIX array size to 16 same as VMTUSR. 
C  07 Apr 00 @EPA-djb Added VALL, VTEMP and VMTUSR. 
C  28 Feb 00 @DynTel-MLA 1-011  Added DBVCOUNT to collect database output. 
C  24 Jan 00 @DynTel-MLA 1-011 Added VMTLDDT. 
C  06 Nov 98 @Dyntel-BAG 2-000 Renamed VMTGT1 and VMTGT2 to VMTGT12 and VMTGT34, respectively 
C 10 Mar 98 @DynTel-RJD 2-622 Added VMTMIX fractions for composite vehicles.   
C  7 Apr 97 - DynTel@MLA  Expanded arrays REGMIX, TFNORM, VCOUNT, and VMTMIX  
C             to accomodate the new BUS vehicle types. 
C 20 Nov 96 DynTel@MLA Expanded the arrays to accomodate the new Heavy Duty 
C           Gasoline, Heavy Duty Diesel, and LDGT vehicle categories. 
C 
C 
C The data is initialized in BD20. 
C 
      REAL REGMIX, TFNORM, VCOUNT, VMTMIX, VMTGT12, VMTGT34, VMTHDD,  
     *     VMTHDG, GLVMTMIX, VMLDGT, VMTLDDT, DBVCOUNT,VMTUSR, 
     *     VALL,VTEMP 
      COMMON /VMXCOM/ REGMIX(28),TFNORM(28),VMTMIX(28),VCOUNT(39,28), 
     *                VMTGT12, VMTGT34, VMTHDD, VMTHDG, GLVMTMIX(16), 
     *                VMLDGT,VMTLDDT,DBVCOUNT(28),VMTUSR(16), 
     *                VALL(30),VTEMP(30) 
