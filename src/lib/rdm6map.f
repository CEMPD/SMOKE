
      SUBROUTINE RDM6MAP( FDEV )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C     Reads the MOBILE6-to-SMOKE vehicle type mapping file and stores
C     data
C
C  PRECONDITIONS REQUIRED:
C     File opened on unit FDEV
C     Unique inventory lists created by GENUSLST
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created 9/2005 by C. Seppanen
C
C**************************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2004, Environmental Modeling for Policy Development
C All Rights Reserved
C 
C Carolina Environmental Program
C University of North Carolina at Chapel Hill
C 137 E. Franklin St., CB# 6116
C Chapel Hill, NC 27599-6116
C 
C smoke@unc.edu
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C***************************************************************************

C.......  MODULES for public variables
C.......  This module contains the lists of unique source characteristics
      USE MODLISTS, ONLY: NINVVTYP, INVVTYP

C.......  This module contains emission factor related variables      
      USE MODEMFAC, ONLY: M6VEHMAP, SMKVEH2EF, NVTYPE
      
      IMPLICIT NONE
      
C.......  INCLUDE files
      INCLUDE 'IODECL3.EXT'   ! I/O API function declarations
      INCLUDE 'EMCNST3.EXT'   ! emissions constant parameters
      INCLUDE 'M6CNST3.EXT'   ! MOBILE6 parameters

C.......  EXTERNAL functions
      LOGICAL  BLKORCMT
      INTEGER  FINDC
      INTEGER  INDEX1
      
      EXTERNAL BLKORCMT, FINDC, INDEX1

C.......  Subroutine arguments
      INTEGER, INTENT(IN) :: FDEV               ! M6MAP file unit number

C.......  Local parameters
      INTEGER, PARAMETER :: NSEG = 2            ! number of segments in line

C.......  Local fixed arrays
      CHARACTER(VTPLEN3) SEGMENT( NSEG )        ! segments of line

C.......  Local allocatable arrays
      LOGICAL, ALLOCATABLE :: TMPINVCHK( : )    ! true: inventory vehicle type found
      
C.......  Local variables
      INTEGER            EPR                    ! emission process counter
      INTEGER            IOS                    ! i/o status
      INTEGER            IREC                   ! record counter
      INTEGER            M6IDX                  ! index into master list of M6 types
      INTEGER            POS                    ! tmp position counter
      INTEGER            SMKIDX                 ! index into inv. list of SMOKE types
      INTEGER            VTYP                   ! vehicle type counter
      
      LOGICAL         :: EFLAG = .FALSE.        ! true: error detected
      
      CHARACTER(256)     LINE                   ! line of file
      CHARACTER(256)     MESG                   ! message buffer
      CHARACTER(VTPLEN3) M6TYPE                 ! MOBILE6 vehicle type
      CHARACTER(VTPLEN3) SMKTYPE                ! SMOKE vehicle type
      
      CHARACTER(16)   :: PROGNAME = 'RDM6MAP'   ! program name
      
C***********************************************************************
C   begin body of subroutine RDVMIX

C.......  Allocate and initialize arrays
      ALLOCATE( TMPINVCHK( NINVVTYP ), STAT=IOS )
      CALL CHECKMEM( IOS, 'TMPINVCHK', PROGNAME )
      ALLOCATE( M6VEHMAP( MXM6VTYP ), STAT=IOS )
      CALL CHECKMEM( IOS, 'M6VEHMAP', PROGNAME )
      ALLOCATE( SMKVEH2EF( MXM6EPR, NINVVTYP ), STAT=IOS )
      CALL CHECKMEM( IOS, 'SMKVEH2EF', PROGNAME )

      TMPINVCHK = .FALSE.   ! array
      M6VEHMAP  = 0         ! array
      SMKVEH2EF = -1        ! array

C.......  Read and process lines of file
      IREC = 0
      DO
          READ( FDEV, '(A)', IOSTAT=IOS ) LINE

C...........  Exit if we've reached the end of the file          
          IF( IOS < 0 ) EXIT
          
          IREC = IREC + 1
          
C...........  Check for I/O errors
          IF( IOS /= 0 ) THEN
              EFLAG = .TRUE.
              WRITE( MESG, 94010 ) 'I/O error', IOS, 'reading ' //
     &            'MOBILE6 mapping file at line', IREC
              CALL M3MESG( MESG )
              CYCLE
          END IF
      
C...........  Skip blank or comment lines
          IF( BLKORCMT( LINE ) ) CYCLE
          
C...........  Parse line into segments
          CALL PARSLINE( LINE, NSEG, SEGMENT )
          
          M6TYPE  = SEGMENT( 1 )
          SMKTYPE = SEGMENT( 2 )

C...........  Check if MOBILE6 type is in master list
          M6IDX = INDEX1( TRIM( M6TYPE ), MXM6VTYP, M6VTYPES )
          IF( M6IDX < 1 ) THEN
              EFLAG = .TRUE.
              WRITE( MESG, 94010 ) 'ERROR: Invalid MOBILE6 ' //
     &            'vehicle type ' // TRIM( M6TYPE ) // 
     &            ' at line', IREC
              CALL M3MESG( MESG )
              CYCLE
          END IF

C...........  Check if SMOKE type is in inventory
          SMKIDX = FINDC( TRIM( SMKTYPE ), NINVVTYP, INVVTYP )
          IF( SMKIDX < 1 ) THEN
              WRITE( MESG, 94010 ) 'WARNING: Vehicle type ' //
     &            TRIM( SMKTYPE ) // ' at line', IREC,
     &            'is not in the mobile source inventory'
              CALL M3MESG( MESG )
              CYCLE
          END IF

C...........  Store SMOKE vehicle type for MOBILE6 type

C...........  Check if this MOBILE6 type is a duplicate
          IF( M6VEHMAP( M6IDX ) /= 0 ) THEN
              EFLAG = .TRUE.
              WRITE( MESG, 94010 ) 'ERROR: Duplicate MOBILE6 ' //
     &            'vehicle type ' // TRIM( M6TYPE ) // ' at line', IREC
              CALL M3MESG( MESG )
              CYCLE
          END IF

          M6VEHMAP ( M6IDX )  = SMKIDX
          TMPINVCHK( SMKIDX ) = .TRUE.
      
      END DO

C.......  Check for any errors
      IF( EFLAG ) THEN
          MESG = 'Problem reading M6MAP file'
          CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      END IF
      
C.......  Check that all MOBILE6 types have been assigned
      DO VTYP = 1, MXM6VTYP
          IF( M6VEHMAP( VTYP ) == 0 ) THEN
              MESG = 'WARNING: MOBILE6 vehicle type ' // 
     &               TRIM( M6VTYPES( VTYP ) ) // ' is not mapped ' //
     &               'to an inventory vehicle type'
              CALL M3MESG( MESG )
          END IF
      END DO

C.......  Check that all inventory types have been assigned
      DO VTYP = 1, NINVVTYP
          IF( .NOT. TMPINVCHK( VTYP ) ) THEN
              MESG = 'WARNING: Inventory vehicle type ' //
     &               TRIM( INVVTYP( VTYP ) ) // ' is not mapped ' //
     &               'to any MOBILE6 vehicle types'
              CALL M3MESG( MESG )
          END IF
      END DO
      
      DEALLOCATE( TMPINVCHK )
      
C.......  Build SMOKE vehicle type / emission process mapping array
      DO EPR = 1, MXM6EPR
          DO VTYP = 1, MXM6VTYP
              IF( M6VEH2EF( EPR, VTYP ) /= -1 ) THEN
                  SMKVEH2EF( EPR, M6VEHMAP( VTYP ) ) = 1
              END IF
          END DO
          
          POS = 0
          DO VTYP = 1, NINVVTYP
              IF( SMKVEH2EF( EPR, VTYP ) == 1 ) THEN
                  POS = POS + 1
                  SMKVEH2EF( EPR, VTYP ) = POS
              END IF
          END DO
      END DO

C.......  Store number of vehicle types in MODEMFAC so MOBILE6 code doesn't
C         have to use MODINFO
      NVTYPE = NINVVTYP

      RETURN

C******************  FORMAT  STATEMENTS   ******************************

94010 FORMAT( 10( A, :, I8, :, 1X ) )

      END SUBROUTINE RDM6MAP
