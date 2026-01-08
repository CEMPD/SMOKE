
       LOGICAL FUNCTION CREATESET( FIDX, RNAME, PHYSNAME, FSTATUS )

!***********************************************************************
!  Function body starts at line 41
!
!  DESCRIPTION:
!     Creates a new file set based on the file description information,
!     opens individual files and sets up internal memory structures
!
!  PRECONDITIONS REQUIRED:
!     File description in FDESC3.EXT and MODFILESET 
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!     OPEN3 - opens an individual file
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
c      INCLUDE 'IODECL3.EXT' ! I/O API function declarations

!........  External functions
c      LOGICAL, EXTERNAL :: SETENVVAR

!........  Subroutine arguments
       INTEGER,        INTENT(IN) :: FIDX      ! file index
       CHARACTER(16),  INTENT(IN) :: RNAME     ! root logical file name
       CHARACTER(256), INTENT(IN) :: PHYSNAME  ! physical file name
       INTEGER,        INTENT(IN) :: FSTATUS   ! file mode (read, read/write)
       
!........  Local variables
       INTEGER I, J, K              ! counters
       INTEGER IOS                  ! I/O status
       INTEGER VARPOS               ! variable position
       INTEGER BLANKS               ! number of blank spaces in FDESC3D
       INTEGER FDESCPOS             ! position in FDESC3D array
       
       LOGICAL EFLAG                ! true: error has occurred
       
       CHARACTER(2)   INTBUF        ! integer string buffer
       CHARACTER(16)  LNAME         ! temporary logical file name
       CHARACTER(256) TEMPNAME      ! temporary physical file name
       CHARACTER(300) MESG          ! message buffer

       CHARACTER(16) :: PROGNAME = 'CREATESET'  ! program name

!-------------------------------------
!  Begin body of function CREATESET
!-------------------------------------

!........  Allocate memory for logical file names and variables                      
       ALLOCATE( FILE_INFO( FIDX )%LNAMES( NFILESET ), STAT=IOS )
       CALL CHECKMEM( IOS, 'FILE_INFO%LNAMES', PROGNAME )
       ALLOCATE( FILE_INFO( FIDX )%VARS( NVARSET,2 ), STAT=IOS )
       CALL CHECKMEM( IOS, 'FILE_INFO%VARS', PROGNAME )

!........  Loop through individual files in the set
       VARPOS = 1
       DO I = 1, NFILESET
               
!............  Set file description values                   
           NVARS3D = VARS_PER_FILE( I )
           VTYPE3D( 1:NVARS3D ) = 
     &              VTYPESET( VARPOS:VARPOS + NVARS3D - 1 )
           VNAME3D( 1:NVARS3D ) = 
     &              VNAMESET( VARPOS:VARPOS + NVARS3D - 1 )
           UNITS3D( 1:NVARS3D ) = 
     &              VUNITSET( VARPOS:VARPOS + NVARS3D - 1 )
           VDESC3D( 1:NVARS3D ) = 
     &              VDESCSET( VARPOS:VARPOS + NVARS3D - 1 )

!............  On first time through, find space in FDESC3D array           
           IF( I == 1 ) THEN
               BLANKS = 0
               DO J = 1, MXDESC3
                   IF( FDESC3D( J ) /= '' ) THEN
                       BLANKS = 0
                   ELSE
                       BLANKS = BLANKS + 1
                   END IF

!....................  If we have 3 consecutive open spots, add file info                   
                   IF( BLANKS == 3 ) THEN
                       WRITE( FDESC3D( J - 2 ),93010 ) 
     &                        '/NUMBER OF FILES/', NFILESET
                       WRITE( FDESC3D( J - 1 ),93010 ) 
     &                        '/FILE POSITION/', I
                       WRITE( FDESC3D( J ),93010 ) 
     &                        '/NUMBER OF VARIABLES/', NVARSET
     
!........................  Save position for writing subsequent info
                       FDESCPOS = J - 1
                       EXIT
                   END IF

!....................  Couldn't find space
                   IF( J == MXDESC3 ) THEN
                       MESG = 'No spaces available in FDESC3D array'
                       CALL M3WARN( PROGNAME, 0, 0, MESG )
                       CREATESET = .FALSE.
                       RETURN
                   END IF
               END DO
           ELSE
!................  Add current file position number
               WRITE( FDESC3D( FDESCPOS ),93010 )
     &                '/FILE POSITION/', I
           END IF

!................  Create new logical and physical file names if needed
           IF( NFILESET > 1 ) THEN
               WRITE( INTBUF, '(I2)' ) I
               INTBUF = ADJUSTL( INTBUF )
               LNAME = TRIM( RNAME ) // TRIM( INTBUF )
               CALL APPENDNAME( PHYSNAME, I, TEMPNAME, EFLAG )
               IF( EFLAG ) THEN
                   CREATESET = .FALSE.
                   RETURN
               END IF

!................  Set new env. variable                       
               IF( .NOT. SETENVVAR( LNAME, TEMPNAME ) ) THEN
                   CREATESET = .FALSE.
                   RETURN
               END IF
           ELSE
               LNAME = RNAME
           END IF

!............  Try to open file           
           IF( .NOT. OPEN3( LNAME, FSTATUS, PROGNAME ) ) THEN
               CREATESET = .FALSE.
               RETURN
           END IF

!............  Store logical file names and variables                        
           FILE_INFO( FIDX )%LNAMES( I ) = LNAME
           FILE_INFO( FIDX )%VARS( VARPOS:VARPOS+NVARS3D-1,1 ) = 
     &                                           VNAME3D( 1:NVARS3D )
           FILE_INFO( FIDX )%VARS( VARPOS:VARPOS+NVARS3D-1,2 ) = LNAME
           VARPOS = VARPOS + NVARS3D

       END DO  ! loop over individual files
       
!........  Reset FDESC3D array back to original 
!..        (remove file number and position info)
       FDESC3D( FDESCPOS-1:FDESCPOS+1 ) = ''
       
       CREATESET = .TRUE.
       RETURN

!---------- Format statements --------------

93010   FORMAT ( A, I4 )
       
       END FUNCTION CREATESET
