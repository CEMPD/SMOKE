C 
C  IOUCOM contains the variables that hold the values of the unit number 
C   for INPUT and OUTPUT files. 
C 
C  Changes: (Last change first) 
C 
C  31 Aug 01 AIR Task 22-24: Added IOUSS 
C  11 Jul 01 AIR Task 15: added I/O for Air Toxic output. 
C  20 Feb 01 AIR Task 11: added I/O for PM output. 
C  02 Oct 00 @EPA-djb, Bug230 Removed unused I/O units. 
C   15 Jun 00 @EPA-bag,  Bug120 Removed unused SSACTION and SSOPEN variables 
C  11 Apr 00 AIR Task 01: Replaced all occurences IOUASK with IOUIN and IOUOUT 
C            so that Lahey Standard for keyboard/monitor I/O is followed. 
C  20 Sep 99 @DynTel-mla 2-000 Added variables for new file input procedure 
C  18 Jun 99 @DynTel-ddj 2-000 Added FTP1, FTP2, and FTP3 to allow writing 
C            out of FTP emission #'s; these are currently commented out 
C   17 Feb 99 @Dyntel-MLA 2-663 Added IOUXIM to read I/M program descriptions 
C             from an external file. 
C 
C   10 Sep 98 @DynTel-BAG 2-661 WMOUT now 80 characters long and added 
C             logical variable READOUTFILE 
C 
C  Common Block Dictionary: 
C 
C   Name    Type                  Description 
C  ------   ----  ----------------------------------------------------- 
C 
C 
C   WMOUT   C*80  name of the user supplied output file 
C   IOUALT   I    unit number that connected with alternate input file 
C   IOUBAT   I    unit number that connected with batch input file 
C   IOUERR   I    unit number that connected with error   file 
C   IOUGEN   I    unit number that connected with MOBILE5 inputs 
C   IOUIMD   I    unit number that connected with I/M credit file 
C   IOUIN    I    unit number that connected with input console (keyboard) 
C   IOUNULL  I    Dummy IO unit - output written here disappears! 
C   IOUOUT   I    unit number that connected with output console (monitor) 
C   IOUPM    I    unit number that connected with PM output 
C   IOUREP   I    unit number that connected with report  file 
C   IOUXIM   I    unit number that connected with external I/M description 
C   IOUAT    I    unit number that connected with Air Toxic output 
C   IOUSS    I    unit number that connected with spreadsheet output 
C 
C  Notes: 
C 
C   25 Mar 98 @DynTel-MLA 2-659 Added IOUALT for Mobile6. 
C 
C   20-Apr-98 @DynTel-BAG 2-661 Added the spreadsheet output file  
C             unit number variable IOUSSH 
C 
C  The data is initialized in BD16 
C 
      CHARACTER*80 WMOUT 
      INTEGER IOUIMD,IOUGEN,IOUREP,IOUALT,IOUERR, IOUIN, IOUOUT, 
     &        IOUBAT, IOUNULL, IOUXIM, IOUPM, IOUAT, IOUSS 
      COMMON /IOUCOM/ IOUIMD,IOUGEN,IOUREP,IOUALT,IOUERR,IOUIN,IOUOUT, 
     &                IOUBAT,IOUNULL,IOUXIM,IOUPM,IOUAT,IOUSS, 
     &                WMOUT 
