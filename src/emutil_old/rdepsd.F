 
        SUBROUTINE RDEPSD( IERR, IOUNIT, ESTRING )

C***************************************************************************
C  program body starts at line 70
C
C  DESCRIPTION:
C    READS THE USER INPUT FILE AND STORES THE DATA IN THE /EPISODE/
C    PACKET IN COMMON VARIABLES
C
C    This program was written because the standard EPS program to 
C    read the EPISODE packet (RDCNTL) did not save all of the 
C    information required for this program.
C
C***************************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1998, MCNC--North Carolina Supercomputing Center
C All Rights Reserved
C
C See file COPYRIGHT for conditions of use.
C
C Environmental Programs Group
C MCNC--North Carolina Supercomputing Center
C P.O. Box 12889
C Research Triangle Park, NC  27709-2889
C
C env_progs@mcnc.org
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C***************************************************************************

        IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE "param.inc"
        INCLUDE "uamcntl.inc"

C...........   LOCAL VARIABLES and their descriptions:
C...........   Point Sources input and output arrays

        INTEGER        IERR        ! Main routine error flag
        INTEGER        JERR        ! Error flag
        INTEGER        IOUNIT      ! File unit number of USERIN file
        INTEGER        IDXSPC      ! Species counter
        INTEGER        NHOURS      ! Number of hours from USERIN file

        CHARACTER*(*)  ESTRING(*)  ! String to store CONTROL info
        CHARACTER*100  BUFFER      ! Character buffer 
        CHARACTER*256  MESG        ! Buffer for M3EXIT() messages
 
C...........   Data statements

        CHARACTER*(KEYLEN) KEYWRD
        INTEGER        IREC, IOS

        DATA IREC   /0/
        DATA IOS    /0/
        DATA KEYWRD /'/EPISODE/'/

C***********************************************************************
C   begin body of program RDEPSD

        IERR = IFAIL

        CALL FNDKEY( JERR, IOUNIT, KEYWRD)   ! EPS subroutine
        IF( JERR .EQ. IRDERR ) THEN
            WRITE( MESG, 94000 )
     &            'Error: failure to correctly read the '//KEYWRD//
     &            ' packet.'
            CALL M3EXIT( 'FNDKEY', 0, 0, MESG, 1 )
        ENDIF 

        IF( JERR .NE. ISUCES ) THEN
            WRITE( MESG, 94000 )
     &            'Error: this program requires the '//KEYWRD//
     &            ' packet of USERIN.'
            CALL M3EXIT( 'RUAMPKS', 0, 0, MESG, 1 )
        ENDIF 

        IREC = IREC + 1
        READ( IOUNIT,93000,ERR=9,IOSTAT=IOS ) ESTRING(1)  ! Do not save
        IREC = IREC + 1
        READ( IOUNIT,93000,ERR=9,IOSTAT=IOS ) ESTRING(1)  ! Do not save
        IREC = IREC + 1
        READ( IOUNIT,93010,ERR=9,IOSTAT=IOS ) FNOTE       ! Save file descriptor
        IREC = IREC + 1
        READ( IOUNIT,93000,ERR=9,IOSTAT=IOS ) ESTRING(1)  ! Save
        IREC = IREC + 1
        READ( IOUNIT,93000,ERR=9,IOSTAT=IOS ) ESTRING(2)  ! Save
        IREC = IREC + 1
        READ( IOUNIT,93000,ERR=9,IOSTAT=IOS ) ESTRING(3)  ! Save
        IREC = IREC + 1
        READ( IOUNIT,93000,ERR=9,IOSTAT=IOS ) ESTRING(4)  ! Save
        IREC = IREC + 1
        READ( IOUNIT,93000,ERR=9,IOSTAT=IOS ) ESTRING(5)  ! Save

    9   IF( IOS .NE. 0 ) THEN
            WRITE( MESG, 94010 )
     &            'Error: reading '//KEYWRD//' packet of userin file '//
     &            'at line ',IREC
            CALL M3EXIT( 'RUAMPKS', 0, 0, MESG, 1 )
        ENDIF 

C.........  Set number of species in userin file

        BUFFER = ESTRING(1)
        READ( BUFFER(1:10), 94030 ) NSPECS

        IF( NSPECS .LE. 0 .OR. NSPECS .GT. MXSPEC ) THEN
             WRITE( MESG, 94020 )
     &             'Error: number of species in USERIN file: ',NSPECS,
     &             ' exceeds max: ',MXSPEC
            CALL M3EXIT( 'RUAMPKS', 0, 0, MESG, 1 )
        ENDIF

C.........  Read species names

        DO 10 IDXSPC = 1, NSPECS
            IREC = IREC + 1
            READ( IOUNIT,93020,IOSTAT=IOS ) SPCNAM( IDXSPC )
   10   CONTINUE

        IF( IOS .NE. 0 ) THEN
            WRITE( MESG, 94010 )
     &            'Error: reading '//KEYWRD//' packet of userin file '//
     &            'at line ',IREC
            CALL M3EXIT( 'RUAMPKS', 0, 0, MESG, 1 )
        ENDIF 

C.........  Read time data (checking for valid data) ---
C
        IREC = IREC + 1
        READ( IOUNIT,93030,IOSTAT=IOS ) IBGDAT,BEGTIM,IENDAT,ENDTIM

        IF( IOS .NE. 0 ) THEN
            WRITE( MESG, 94010 )
     &            'Error: reading '//KEYWRD//' packet of userin file '//
     &            'at line ',IREC
            CALL M3EXIT( 'RUAMPKS', 0, 0, MESG, 1 )
        ENDIF

        IF( IBGDAT .GT. IENDAT ) THEN
            WRITE( MESG, 94020 )
     &            'Error: invalid julian dates: beginning: ',IBGDAT,
     &            ' ending: ',IENDAT
            CALL M3EXIT( 'RUAMPKS', 0, 0, MESG, 1 )
        ENDIF

        NHOURS = (IENDAT-IBGDAT)*24 + NINT( (ENDTIM-BEGTIM)/100. )
        IF( NHOURS .GT. MXHOUR ) THEN
            WRITE( MESG, 94020 )
     &            'Error: number of hours: ',NHOURS,' exceeds max: ',
     &            MXHOUR
            CALL M3EXIT( 'RUAMPKS', 0, 0, MESG, 1 )
        ENDIF

        IERR = ISUCES
  100   RETURN

C******************  FORMAT  STATEMENTS   *****************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )
93010   FORMAT( A60 )
93020   FORMAT( A10 )
93030   FORMAT( I10,F10.0,I10,F10.0 )

C...........   Internal buffering formats.............94xxx

94000   FORMAT( 10A )
94010   FORMAT( 12( A, :, I8, 1X ) )
94020   FORMAT( A, I8, A, I8 )
94030   FORMAT( I10 )

        END
