
       LOGICAL FUNCTION DESCSET( ROOTNAME, FILENUM )

!***********************************************************************
!  Function body starts at line 40
!
!  DESCRIPTION:
!     Get description of file set
!
!  PRECONDITIONS REQUIRED:
!     File set has been opened with OPENSET
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!     DESC3 - get file description for individual file
!
!  REVISION HISTORY:
!     Created 6/02 by C. Seppanen
!     09/2025 by HT UNC-IE: Use M3UTILIO
!
!***************************************************************************
!
! Project Title: FileSetAPI
! File: @(#)$Id$
!
! COPYRIGHT (C) 2004, Environmental Modeling for Policy Development
! All Rights Reserved
!
! Carolina Environmental Program
! University of North Carolina at Chapel Hill
! 137 E. Franklin St., CB# 6116
! Chapel Hill, NC 27599-6116
!
! smoke@unc.edu
!
! Pathname: $Source$
! Last updated: $Date$ 
!
!*************************************************************************
       USE M3UTILIO

!........  Modules for public variables
       USE MODFILESET

       IMPLICIT NONE
       
!........  Include files
c      INCLUDE 'IODECL3.EXT'  ! I/O API function declarations
       
!........  External functions
c      INTEGER,       EXTERNAL :: INDEX1
       CHARACTER(10), EXTERNAL :: GETCFDSC
c      INTEGER,       EXTERNAL :: STR2INT
       
!........  Function arguments
       CHARACTER(*), INTENT(IN) :: ROOTNAME  ! logical file name for file set
       INTEGER,      INTENT(IN) :: FILENUM   ! position of file in set or ALLFILES
  
       INTEGER, PARAMETER :: ALLFILES = -1

!........  Local variables
       INTEGER            I           ! counter
       INTEGER            IOS         ! I/O status
       INTEGER            FILEIDX     ! file index
       INTEGER            VARPOS      ! position in variable information arrays
       INTEGER            NFILEINT    ! number of files as an integer
       INTEGER            NVARINT     ! number of variables as an integer

       CHARACTER(2)   NFILESTR        ! number of files as a string
       CHARACTER(4)   NVARSTR         ! number of variables as a string       
       CHARACTER(16)  ROOTNAME16      ! fixed length root file name
       CHARACTER(16)  LNAME           ! temporary logical name
       CHARACTER(256) MESG            ! message buffer

       CHARACTER(16) :: FUNCNAME = 'DESCSET'  ! function name
       
!---------------------------------
!  Begin body of function DESCSET
!---------------------------------

!........  Check length of file name
       IF( LEN( ROOTNAME ) > 16 ) THEN
           MESG = 'Max file name length (16) exceeded for "' // 
     &            ROOTNAME // '"'
           CALL M3MSG2( MESG )
           DESCSET = .FALSE.
           RETURN
       END IF

!........  Get file index
       ROOTNAME16 = ROOTNAME
       FILEIDX = INDEX1( ROOTNAME16, MXFILE3, RNAMES )

!........  If file is not open, exit with error
       IF( FILEIDX == 0 ) THEN
           MESG = 'File set "' // TRIM( ROOTNAME ) // '" is not ' //
     &            'currently open'
           CALL M3MSG2( MESG )
           DESCSET = .FALSE.
           RETURN
       END IF

!........  Deallocate variable information arrays if necessary
       IF( ALLOCATED( VARS_PER_FILE ) ) THEN
           DEALLOCATE( VARS_PER_FILE, STAT=IOS )
           IF( IOS > 0 ) THEN
               MESG = 'Failure deallocating memory for ' //
     &                '"VARS_PER_FILE" variable'
               CALL M3EXIT( FUNCNAME, 0, 0, MESG, 2 )
           END IF
       END IF
       
       IF( ALLOCATED( VTYPESET ) ) THEN
           DEALLOCATE( VTYPESET, STAT=IOS )
           IF( IOS > 0 ) THEN
               MESG = 'Failure deallocating memory for ' //
     &                '"VTYPESET" variable'
               CALL M3EXIT( FUNCNAME, 0, 0, MESG, 2 )
           END IF
       END IF

       IF( ALLOCATED( VNAMESET ) ) THEN
           DEALLOCATE( VNAMESET, STAT=IOS )
           IF( IOS > 0 ) THEN
               MESG = 'Failure deallocating memory for ' //
     &                '"VNAMESET" variable'
               CALL M3EXIT( FUNCNAME, 0, 0, MESG, 2 )
           END IF
       END IF
       
       IF( ALLOCATED( VUNITSET ) ) THEN
           DEALLOCATE( VUNITSET, STAT=IOS )
           IF( IOS > 0 ) THEN
               MESG = 'Failure deallocating memory for ' //
     &                '"VUNITSET" variable'
               CALL M3EXIT( FUNCNAME, 0, 0, MESG, 2 )
           END IF
       END IF
       
       IF( ALLOCATED( VDESCSET ) ) THEN
           DEALLOCATE( VDESCSET, STAT=IOS )
           IF( IOS > 0 ) THEN
               MESG = 'Failure deallocating memory for ' //
     &                '"VDESCSET" variable'
               CALL M3EXIT( FUNCNAME, 0, 0, MESG, 2 )
           END IF
       END IF

!........  Check if a single file is requested
       IF( FILENUM /= ALLFILES ) THEN
 
!............  Make sure file number is positive
           IF( FILENUM < 1 ) THEN
               MESG = 'File number must be positive or "ALLFILES"'
               CALL M3MSG2( MESG )
               DESCSET = .FALSE.
               RETURN
               
!............  Check that number is not greater than total number of files
           ELSE IF( FILENUM > SIZE( FILE_INFO( FILEIDX )%LNAMES ) ) THEN
               WRITE( MESG,92010 ) 'Invalid file number requested; ' //
     &                'file set "' // TRIM( ROOTNAME ) // '" contains ',
     &                SIZE( FILE_INFO( FILEIDX )%LNAMES ), ' files '
               CALL M3MSG2( MESG )
               DESCSET = .FALSE.
               RETURN
           ELSE

!................  Set logical file name
               LNAME = FILE_INFO( FILEIDX )%LNAMES( FILENUM )

!................  Try to get file description
               IF( .NOT. DESC3( LNAME ) ) THEN
                   DESCSET = .FALSE.
                   RETURN
               END IF

!................  Store file set information               
               NFILESET = 1
               NVARSET = NVARS3D
               
!................  Allocate variable information arrays           
               ALLOCATE( VARS_PER_FILE( NFILESET ), STAT=IOS )
               CALL CHECKMEM( IOS, 'VARS_PER_FILE', FUNCNAME )
               ALLOCATE( VTYPESET( NVARSET ), STAT=IOS )
               CALL CHECKMEM( IOS, 'VTYPESET', FUNCNAME )
               ALLOCATE( VNAMESET( NVARSET ), STAT=IOS )
               CALL CHECKMEM( IOS, 'VNAMESET', FUNCNAME )
               ALLOCATE( VUNITSET( NVARSET ), STAT=IOS )
               CALL CHECKMEM( IOS, 'VUNITSET', FUNCNAME )
               ALLOCATE( VDESCSET( NVARSET ), STAT=IOS )
               CALL CHECKMEM( IOS, 'VDESCSET', FUNCNAME )

!................  Store info for this file
               VARS_PER_FILE( 1 ) = NVARS3D
               VTYPESET( 1:NVARS3D ) = VTYPE3D( 1:NVARS3D )
               VNAMESET( 1:NVARS3D ) = VNAME3D( 1:NVARS3D )
               VUNITSET( 1:NVARS3D ) = UNITS3D( 1:NVARS3D )
               VDESCSET( 1:NVARS3D ) = VDESC3D( 1:NVARS3D )           
           END IF
       ELSE

!............  Open first file of set
           LNAME = FILE_INFO( FILEIDX )%LNAMES( 1 )
           IF( .NOT. DESC3( LNAME ) ) THEN
               DESCSET = .FALSE.
               RETURN
           END IF

!............  Get total number of files from file header                   
           NFILESTR = GETCFDSC( FDESC3D, '/NUMBER OF FILES/', 
     &                          .TRUE. )
           NFILEINT = STR2INT( NFILESTR )
           IF( NFILEINT == IMISS3 ) THEN
               MESG = 'Invalid number of files in header of file ' //
     &                'set "' // TRIM( ROOTNAME ) // '"'
               CALL M3MSG2( MESG ) 
               DESCSET = .FALSE.
               RETURN
           END IF

!............  Get total number of variables from file header                
           NVARSTR  = GETCFDSC( FDESC3D, '/NUMBER OF VARIABLES/', 
     &                          .TRUE. )
           NVARINT  = STR2INT( NVARSTR )
           IF( NVARINT == IMISS3 ) THEN
               MESG = 'Invalid number of variables in header of ' //
     &                'file set "' // TRIM( ROOTNAME ) // '"'
               CALL M3MSG2( MESG )
               DESCSET = .FALSE.
               RETURN
           END IF
           
!............  Store file set information from header
           NFILESET = NFILEINT
           NVARSET = NVARINT

!............  Allocate variable information arrays           
           ALLOCATE( VARS_PER_FILE( NFILESET ), STAT=IOS )
           CALL CHECKMEM( IOS, 'VARS_PER_FILE', FUNCNAME )
           ALLOCATE( VTYPESET( NVARSET ), STAT=IOS )
           CALL CHECKMEM( IOS, 'VTYPESET', FUNCNAME )
           ALLOCATE( VNAMESET( NVARSET ), STAT=IOS )
           CALL CHECKMEM( IOS, 'VNAMESET', FUNCNAME )
           ALLOCATE( VUNITSET( NVARSET ), STAT=IOS )
           CALL CHECKMEM( IOS, 'VUNITSET', FUNCNAME )
           ALLOCATE( VDESCSET( NVARSET ), STAT=IOS )
           CALL CHECKMEM( IOS, 'VDESCSET', FUNCNAME )

!............  Store info for this file
           VARS_PER_FILE( 1 ) = NVARS3D
           VTYPESET( 1:NVARS3D ) = VTYPE3D( 1:NVARS3D )
           VNAMESET( 1:NVARS3D ) = VNAME3D( 1:NVARS3D )
           VUNITSET( 1:NVARS3D ) = UNITS3D( 1:NVARS3D )
           VDESCSET( 1:NVARS3D ) = VDESC3D( 1:NVARS3D )
           VARPOS = NVARS3D + 1
           
!............  Loop through remaining files
           DO I = 2, NFILESET

!................  Set logical file name
               LNAME = FILE_INFO( FILEIDX )%LNAMES( I )

!................  Try to get file description
               IF( .NOT. DESC3( LNAME ) ) THEN
                   DESCSET = .FALSE.
                   RETURN
               END IF

!................  Store current file description                
               VARS_PER_FILE( I ) = NVARS3D
               VTYPESET( VARPOS:VARPOS + NVARS3D - 1 ) = 
     &                                   VTYPE3D( 1:NVARS3D )
               VNAMESET( VARPOS:VARPOS + NVARS3D - 1 ) = 
     &                                   VNAME3D( 1:NVARS3D )
               VUNITSET( VARPOS:VARPOS + NVARS3D - 1 ) = 
     &                                   UNITS3D( 1:NVARS3D )
               VDESCSET( VARPOS:VARPOS + NVARS3D - 1 ) = 
     &                                   VDESC3D( 1:NVARS3D )
               VARPOS = VARPOS + NVARS3D
           END DO

       END IF
       
       DESCSET = .TRUE.

!---------- Format statements --------------

92010   FORMAT ( A, :, I3, :, A )

       END FUNCTION DESCSET
