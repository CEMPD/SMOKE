
        PROGRAM MVCONDNS

C***********************************************************************
C  program body starts at line
C
C  DESCRIPTION:
C    Program MVSETUP uses the SMOKE mobile source file to create an condensed
C    list of all sources in the inventory by FIPS code, roadtype, vehicle type,
C    and including the speed from the inventory (if any).  A default speed
C    is determined for each county and road type, which is not associated 
C    with a vehicle type.
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C    Original by M. Houyoux 9/99
C
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE)
C                
C File: @(#)$Id$
C
C COPYRIGHT (C) 2001, MCNC--North Carolina Supercomputing Center
C All Rights Reserved
C
C See file COPYRIGHT for conditions of use.
C
C Environmental Programs Group
C MCNC--North Carolina Supercomputing Center
C P.O. Box 12889
C Research Triangle Park, NC  27709-2889
C
C smoke@ncsc.org
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C***********************************************************************
 
C...........   MODULES for public variables
C...........   This module is the source inventory arrays
        USE MODSOURC

C.........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE
 
C...........   INCLUDES:
        INCLUDE 'EMCNST3.EXT'
        INCLUDE 'PARMS3.EXT'
        INCLUDE 'IODECL3.EXT'
        INCLUDE 'FDESC3.EXT'
      
C...........   EXTERNAL FUNCTIONS 
        CHARACTER*2   CRLF
        INTEGER       FINDR2
        INTEGER       FINDR3
        INTEGER       INDEX1
        INTEGER       PROMPTFFILE
        CHARACTER*16  PROMPTMFILE

        EXTERNAL      CRLF, FINDR2, FINDR3, INDEX1, PROMPTFFILE, 
     &                PROMPTMFILE

C...........   PARAMETERS and their descriptions:
        CHARACTER*50, PARAMETER :: CVSW = '$Name$' ! CVS release tag

C.........  Array that contains the names of the inventory variables needed 
C           for this program
        CHARACTER(LEN=IOVLEN3) IVARNAMS( MXINVARR )

C.........  Local allocatable arrays
        INTEGER, ALLOCATABLE :: INDX  ( : ) ! Sorting index
        INTEGER, ALLOCATABLE :: DEFFIP( : ) ! Cy/St/Co codes for default out
        INTEGER, ALLOCATABLE :: DEFRCL( : ) ! Road class codes for default out
        REAL   , ALLOCATABLE :: DEFSPD( : ) ! Speeds for default outputs
        REAL   , ALLOCATABLE :: RFIP  ( : ) ! FIPS codes converted to reals
        REAL   , ALLOCATABLE :: RRCLAS( : ) ! road classes converted to reals
        REAL   , ALLOCATABLE :: RVTYPE( : ) ! road classes converted to reals

C...........   File units and logical/physical names
        INTEGER         LDEV    !  log-device
        INTEGER      :: SDEV = 0!  ASCII part of inventory unit no.
        INTEGER      :: ODEV = 0!  Output file unit no.

        CHARACTER*16    ANAME   !  logical name for ASCII inventory input file
        CHARACTER*16    ENAME   !  logical name for i/o api inventory input file

C...........   LOCAL VARIABLES and their descriptions:
        INTEGER         J, K, N, S

        INTEGER         ENLEN   ! length of the emissions inven name
        INTEGER         IOS     ! i/o status
        INTEGER         FIP     ! tmp cy/st/co code
        INTEGER         LFIP    ! cy/st/co code from previous iteration
        INTEGER         LRCL    ! road class code from previous iteration
        INTEGER         NINVARR ! no. of inventory characteristics
        INTEGER         NDEFLT  ! no. default entries
        INTEGER         MAXSCNT ! max value of SCNT for cy/st/co-roadclass
        INTEGER         RCL     ! tmp roadclass code
        INTEGER         SCNT    ! count of speeds of same value
        INTEGER         VTP     ! tmp vehicle type number

        REAL         :: SPD = 0.    ! tmp speed
        REAL            LSPD    ! speed from previous iteration
        REAL            SPDOFMAX! speed associated with MAXSCNT

        LOGICAL    :: EFLAG   = .FALSE.          ! error flag
        LOGICAL    :: NFLAG   = .TRUE.           ! true: new cy/st/co-roadclass
        LOGICAL    :: SFLAG   = .FALSE.          ! true: speeds in file

        CHARACTER*300             MESG           ! message field

        CHARACTER*16, PARAMETER :: PROGNAME = 'MVCONDNS'

C***********************************************************************
C   begin body of program MVCONDNS
 
        LDEV = INIT3()
 
C.........  Write out copywrite, version, web address, header info, and prompt
C           to continue running the program.
        CALL INITEM( LDEV, CVSW, PROGNAME )
 
C.........  Limit source category to mobile sources
        CATEGORY = 'MOBILE'
        CATDESC  = 'Mobile'
        CATLEN   = LEN_TRIM( CATEGORY )

C.........  Get inventory file names given source category
        CALL GETINAME( CATEGORY, ENAME, ANAME )

C.........  Get file names and open files
        ENAME = PROMPTMFILE( 
     &          'Enter logical name for the I/O API INVENTORY file',
     &          FSREAD3, ENAME, PROGNAME )
        ENLEN = LEN_TRIM( ENAME )

C.........  Get header description of inventory file, error if problem
        IF( .NOT. DESC3( ENAME ) ) THEN
            MESG = 'Could not get description of file "' //
     &             ENAME( 1:ENLEN ) // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C.........  Otherwise, store source-category-specific header information, 
        ELSE
       
            CALL GETSINFO

        END IF

C.........  Read source characteristics from mobile source inventory files
        NINVARR = 3
        IVARNAMS( 1 ) = 'IFIP'
        IVARNAMS( 2 ) = 'IRCLAS'
        IVARNAMS( 3 ) = 'IVTYPE'

        J = INDEX1( 'SPEED', NIACT, ACTVTY )
        IF ( J .GT. 0 ) THEN
            NINVARR = 4
            IVARNAMS( 4 ) = 'SPEED'
            SFLAG = .TRUE.
        END IF

C.........  Allocate memory for and read in required inventory characteristics
        CALL RDINVCHR( CATEGORY, ENAME, SDEV, NSRC, NINVARR, IVARNAMS )

C.........  Allocate memory for the sorting index and for store the default
C           cy/st/co-roadclass-speed combinations
        ALLOCATE( INDX( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INDX', PROGNAME )
        ALLOCATE( DEFFIP( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'DEFFIP', PROGNAME )
        ALLOCATE( DEFRCL( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'DEFRCL', PROGNAME )
        ALLOCATE( DEFSPD( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'DEFSPD', PROGNAME )
 
C.........  Initialize sorting index
        DO S = 1, NSRC
            INDX( S ) = S
        END DO

C.........  Copy integer codes to real arrays for sorting
        ALLOCATE( RFIP( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'RFIP', PROGNAME )
        ALLOCATE( RRCLAS( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'RRCLAS', PROGNAME )
        ALLOCATE( RVTYPE( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'RVTYPE', PROGNAME )
        
        RFIP   = REAL( IFIP   ) ! array
        RRCLAS = REAL( IRCLAS ) ! array
        RVTYPE = REAL( IVTYPE ) ! array

C.........  Sort sources by cy/st/co code, roadclass, speed, and vehicle type
        IF ( SFLAG ) THEN
            CALL SORTR4( NSRC, INDX, RFIP, RRCLAS, SPEED, RVTYPE )
        ELSE
            CALL SORTI3( NSRC, INDX, IFIP, IRCLAS, IVTYPE )
        END IF

C.........  Loop through sorted sources and count the number of different
C           cy/st/co code-roadclass-speed combinations, and the total for each.
C           Store the total count of each combination, and their associated
C           source characteristics.
        LFIP    = IFIP( 1 )
        LRCL    = IRCLAS( 1 )
        LSPD    = IMISS3
        N       = 0
        SCNT    = 0
        MAXSCNT = 0 
        DO S = 1, NSRC

           J   = INDX  ( S )
           FIP = IFIP  ( J )
           RCL = IRCLAS( J )
           IF( SFLAG ) SPD = SPEED ( J )

C............  If new cy/st/co and/or road class group, store info from
C              previous group and set flag that indicates current source
C              is new
           IF( FIP .NE. LFIP .OR. RCL .NE. LRCL .OR. S .EQ. NSRC ) THEN
               N = N + 1

               DEFFIP( N ) = LFIP
               DEFRCL( N ) = LRCL
               DEFSPD( N ) = SPDOFMAX
               LFIP   = FIP
               LRCL   = RCL

               SCNT = 0      ! Initialize count for new source
               MAXSCNT = 0   ! Initialize maximum for new source

           END IF

           IF( SPD .NE. LSPD ) SCNT = 0  ! Initialize count for new speed

           SCNT = SCNT + 1               ! Increment speed counter

           IF( SCNT .GT. MAXSCNT ) THEN  ! Store max values for current source

               MAXSCNT = SCNT
               SPDOFMAX = SPD

           END IF

           LSPD = SPD                    ! Update previous speed

        END DO

        NDEFLT = N

        ODEV = PROMPTFFILE(  
     &        'Enter name for CONDENSED MOBILE SOURCES output file',
     &        .FALSE., .TRUE., 'OUTFILE', PROGNAME )

C.........  Loop through all sources.  In each group of cy/st/co codes and 
C           road classes, output the default and non-default entries for
C           speeds by vehicle type.
        LFIP = IMISS3
        LRCL = IMISS3
        DO S = 1, NSRC

            FIP = IFIP( S )
            RCL = IRCLAS( S )
            VTP = IVTYPE( S )
            IF( SFLAG ) SPD = SPEED( S )

            IF ( SFLAG ) THEN
                K = FINDR3( FIP, RCL, SPD, NDEFLT, 
     &                      DEFFIP, DEFRCL, DEFSPD )
            ELSE
                K = FINDR2( FIP, RCL, NDEFLT, DEFFIP, DEFRCL )
            END IF

            IF( FIP .NE. LFIP .OR. RCL .NE. LRCL ) NFLAG = .TRUE.

            IF( K .GT. 0 ) THEN

                IF( NFLAG ) WRITE( ODEV,93300 ) FIP, RCL, 0, SPD
                NFLAG = .FALSE.

            ELSE

                WRITE( ODEV,93300 ) FIP, RCL, VTP, SPD

            END IF

            LFIP = FIP
            LRCL = RCL

        END DO

C.........  Normal completion
        CALL M3EXIT( PROGNAME, 0, 0, ' ', 0 )
	
C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93300   FORMAT( I6.6, ',', I8.2, ',', I5.4, ',', F5.1 )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I7, :, 1X ) )
 
	END PROGRAM MVCONDNS
