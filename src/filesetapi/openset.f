
       LOGICAL FUNCTION OPENSET( ROOTNAME, FSTATUS, PGNAME )

!***********************************************************************
!  Function body starts at line 48
!
!  DESCRIPTION:
!     Opens a file set; can be already opened, existing on disk, or new
!
!  PRECONDITIONS REQUIRED:
!     If opening a new file, file description information must be set
!     in FDESC3.EXT and MODFILESET
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!     APPENDNAME - appends a number to a file name
!     DESC3 - get file description for individual file
!     CHKFILESET - checks open file against description
!     CHKSETDESC - checks file description
!     CLEANUP - cleans up internal memory structures
!     CLOSESET - closes a file set
!     CREATESET - creates a new file set
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
c      INCLUDE 'IODECL3.EXT'  ! I/O API function declarations

!........  External functions
c      LOGICAL,       EXTERNAL :: SETENVVAR
       INTEGER,       EXTERNAL :: RMFILE
c      INTEGER,       EXTERNAL :: INDEX1
       CHARACTER(10), EXTERNAL :: GETCFDSC
c      INTEGER,       EXTERNAL :: STR2INT
       LOGICAL,       EXTERNAL :: CLOSESET
       LOGICAL,       EXTERNAL :: CHKSETDESC
       LOGICAL,       EXTERNAL :: CHKFILESET
       LOGICAL,       EXTERNAL :: CREATESET
        
!........  Function arguments
       CHARACTER(*), INTENT(IN) :: ROOTNAME  ! logical file name for file set
       INTEGER,      INTENT(IN) :: FSTATUS   ! file mode/status (read/write, unknown)
       CHARACTER(*), INTENT(IN) :: PGNAME    ! name of calling program

!........  Local variables
       INTEGER            I           ! counter
       INTEGER            IOS         ! I/O status
       INTEGER            FILEIDX     ! file index
       INTEGER            VARPOS      ! position in variable information arrays
       INTEGER            NFILEINT    ! number of files as an integer
       INTEGER            NVARINT     ! number of variables as an integer

       LOGICAL, SAVE :: INITIAL = .TRUE. ! true: first time through function
       LOGICAL            FILEEXIST   ! true: individual file exists on disk
       LOGICAL         :: SETEXIST = .FALSE.   ! true: file set exists on disk
       LOGICAL            TEMPEXIST   ! true: temporary file name exists on disk
       LOGICAL            EFLAG       ! true: error occured

       CHARACTER(2)       NFILESTR    ! number of files as a string
       CHARACTER(4)       NVARSTR     ! number of variables as a string
       CHARACTER(2)       INTBUF      ! integer string buffer
       CHARACTER(16)      ROOTNAME16  ! fixed length root file name
       CHARACTER(16)      PGNAME16    ! fixed length calling program name
       CHARACTER(16)      LNAME       ! temporary logical name
       CHARACTER(256)     ENVNAME     ! environment variable value for ROOTNAME
       CHARACTER(256)     TEMPNAME    ! temporary physical file name
       CHARACTER(256)     MESG        ! message buffer
       
       CHARACTER(16) :: FUNCNAME = 'OPENSET'  ! function name

!---------------------------------
!  Begin body of function OPENSET
!---------------------------------

!........  Initialize arrays first time 
       IF( INITIAL ) THEN
           RNAMES = CMISS3
           
           DO I = 1, MXFILE3
               NULLIFY( FILE_INFO( I )%VARS, FILE_INFO( I )%LNAMES )
           END DO
 
           INITIAL = .FALSE.
       END IF

!........  Check length of file name
       IF( LEN( ROOTNAME ) > 16 ) THEN
           MESG = 'Max file name length (16) exceeded for "' // 
     &            ROOTNAME // '"'
           CALL M3MSG2( MESG )
           OPENSET = .FALSE.
           RETURN
       END IF

!........  Check if file is already open
       ROOTNAME16 = ROOTNAME
       FILEIDX = INDEX1( ROOTNAME16, MXFILE3, RNAMES )

!........  If file is already open, check status
       IF( FILEIDX /= 0 ) THEN
        
!............  Trying to open as a new file
           IF( FSTATUS == FSNEW3 ) THEN
               MESG = 'File ' // ROOTNAME16 // ' already opened; ' //
     &                'cannot subsequently create in "NEW".'
               CALL M3WARN( FUNCNAME, 0, 0, MESG )
               OPENSET = .FALSE.
               RETURN

!............  Trying to open as read/write or unknown when already readonly
           ELSE IF( FILE_INFO( FILEIDX )%RDONLY .AND.  
     &            ( FSTATUS == FSRDWR3 .OR. FSTATUS == FSUNKN3 ) ) THEN
               MESG = 'File ' // ROOTNAME16 // ' already opened ' //
     &                'READONLY; cannot subsequently open it ' //
     &                'for READ/WRITE.'
               CALL M3WARN( FUNCNAME, 0, 0, MESG )
               OPENSET = .FALSE.
               RETURN

!............  Trying to open as unknown (already checked that status is not readonly)
           ELSE IF( FSTATUS == FSUNKN3 ) THEN

!................  Check consistency of file set description
               IF( .NOT. CHKSETDESC( ROOTNAME16 ) ) THEN
                   MESG = 'Bad file description for file set ' //
     &                    ROOTNAME16
                   CALL M3WARN( FUNCNAME, 0, 0, MESG )
                   OPENSET = .FALSE.
                   RETURN
               END IF

!................  Check description against file info
               IF( .NOT. CHKFILESET( FILEIDX ) ) THEN
                   MESG = 'File description does not match for ' //
     &                    'file set ' // ROOTNAME16
                   CALL M3WARN( FUNCNAME, 0, 0, MESG )
                   OPENSET = .FALSE.
                   RETURN
               END IF

!................  Loop through individual files
               VARPOS = 1
               DO I = 1, NFILESET
               
!....................  Set necessary file description values                   
                   NVARS3D = VARS_PER_FILE( I )
                   VTYPE3D( 1:NVARS3D ) = 
     &                      VTYPESET( VARPOS:VARPOS + NVARS3D - 1 )
                   VNAME3D( 1:NVARS3D ) = 
     &                      VNAMESET( VARPOS:VARPOS + NVARS3D - 1 )

!....................  Try to open file
                   IF( .NOT. OPEN3( FILE_INFO( FILEIDX )%LNAMES( I ), 
     &                              FSTATUS, FUNCNAME ) ) THEN
                       OPENSET = .FALSE.
                       RETURN
                   END IF

                   VARPOS = VARPOS + NVARS3D
               END DO

!................  Successfully opened files                    
               OPENSET = .TRUE.
               RETURN

!............  Trying to open as create               
           ELSE IF( FSTATUS == FSCREA3 ) THEN
                
!................  Check consistency of file set description
               IF( .NOT. CHKSETDESC( ROOTNAME16 ) ) THEN
                   MESG = 'Bad file description for file set ' //
     &                    ROOTNAME16
                   CALL M3WARN( FUNCNAME, 0, 0, MESG )
                   OPENSET = .FALSE.
                   RETURN
               END IF

!................  Close already open file
               IF( CLOSESET( ROOTNAME16 ) ) THEN
                   MESG = 'File ' // ROOTNAME16 // ' already ' //
     &                    'opened. Closing, deleting, and ' //
     &                    're-opening it'
                   CALL M3WARN( FUNCNAME, 0, 0, MESG )
               ELSE
                   MESG = 'File ' // ROOTNAME16 // ' already ' //
     &                    'opened. Could not close to reopen ' //
     &                    'with status FSCREA3'
                   CALL M3WARN( FUNCNAME, 0, 0, MESG )
                   OPENSET = .FALSE.
                   RETURN
               END IF

!............  Checked all problematic combinations
           ELSE
           
               OPENSET = .TRUE.
               RETURN
           
           END IF

!........  Otherwise, file is not open
       ELSE
           IF( FSTATUS == FSNEW3 .OR. FSTATUS == FSUNKN3 .OR.
     &         FSTATUS == FSCREA3 ) THEN
        
!................  Check consistency of file set description
               IF( .NOT. CHKSETDESC( ROOTNAME16 ) ) THEN
                   MESG = 'Bad file description for file set ' //
     &                    ROOTNAME16
                   CALL M3WARN( FUNCNAME, 0, 0, MESG )
                   OPENSET = .FALSE.
                   RETURN
               END IF
           END IF
       END IF

!........  Find a place for the new file
       FILEIDX = INDEX1( CMISS3, MXFILE3, RNAMES )
       IF( FILEIDX == 0 ) THEN
           MESG = 'Could not open ' // ROOTNAME16 // 
     &            'Maximum number of files already have been opened.'
           CALL M3WARN( FUNCNAME, 0, 0, MESG )
           OPENSET = .FALSE.
           RETURN
       END IF

!........  Get env. variable setting for ROOTNAME
       CALL NAMEVAL( ROOTNAME16, ENVNAME )

!........  Find out if file already exists on disk
!..        First try original env. variable setting
       INQUIRE( FILE = ENVNAME, EXIST = FILEEXIST )

!........  If original file doesn't exist, try first file of set
       IF( .NOT. FILEEXIST ) THEN
           CALL APPENDNAME( ENVNAME, 1, TEMPNAME, EFLAG )
           IF( EFLAG ) THEN
               MESG = 'Could not generate individual file names ' //
     &                'for file set "' // TRIM( ROOTNAME ) // '"'
               CALL M3WARN( FUNCNAME, 0, 0, MESG )
               CALL CLEANUP( FILEIDX )
               OPENSET = .FALSE.
               RETURN
           END IF
           INQUIRE( FILE = TEMPNAME, EXIST = SETEXIST )
       END IF
       
!........  Store logical file name and readonly status
       RNAMES( FILEIDX ) = ROOTNAME16
       FILE_INFO( FILEIDX )%RDONLY = ( FSTATUS .EQ. FSREAD3 )

!........  Open file based on status and if file or set exists
       IF( FSTATUS == FSREAD3 .OR. FSTATUS == FSRDWR3 ) THEN

!............  If single file exists, we don't have a file set
           IF( FILEEXIST ) THEN
                
!................  Try to open single file
               IF( .NOT. OPEN3( ROOTNAME16, FSTATUS, FUNCNAME ) ) THEN
                   CALL CLEANUP( FILEIDX )
                   OPENSET = .FALSE.
                   RETURN
               END IF
               
!................  Get file description
               IF( .NOT. DESC3( ROOTNAME16 ) ) THEN
                   CALL CLEANUP( FILEIDX )
                   OPENSET = .FALSE.
                   RETURN
               END IF

!................  Allocate memory for logical file names and variables                 
               ALLOCATE( FILE_INFO( FILEIDX )%LNAMES( 1 ), STAT=IOS )
               CALL CHECKMEM( IOS, 'FILE_INFO%LNAMES', FUNCNAME )
               ALLOCATE( FILE_INFO( FILEIDX )%VARS( NVARS3D,2 ), 
     &                   STAT=IOS )
               CALL CHECKMEM( IOS, 'FILE_INFO%VARS', FUNCNAME )

!................  Store logical file names and variables                   
               FILE_INFO( FILEIDX )%LNAMES( 1 ) = ROOTNAME16
               FILE_INFO( FILEIDX )%VARS( 1:NVARS3D,1 ) = 
     &                                          VNAME3D( 1:NVARS3D )
               FILE_INFO( FILEIDX )%VARS( 1:NVARS3D,2 ) = ROOTNAME16

!............  Otherwise, if the file set exists, open the individual files
           ELSE IF( SETEXIST ) THEN

!................  Create new logical file name and set the env. variable
               LNAME = TRIM( ROOTNAME16 ) // '1'
               IF( .NOT. SETENVVAR( LNAME, TEMPNAME ) ) THEN
                   CALL CLEANUP( FILEIDX )
                   OPENSET = .FALSE.
                   RETURN
               END IF

!................  Try to open the first file in the set
               IF( .NOT. OPEN3( LNAME, FSTATUS, FUNCNAME ) ) THEN
                   CALL CLEANUP( FILEIDX )
                   OPENSET = .FALSE.
                   RETURN
               END IF

!................  Get file description               
               IF( .NOT. DESC3( LNAME ) ) THEN
                   CALL CLEANUP( FILEIDX )
                   OPENSET = .FALSE.
                   RETURN
               END IF

!................  Get total number of files from file header                   
               NFILESTR = GETCFDSC( FDESC3D, '/NUMBER OF FILES/', 
     &                              .TRUE. )
               NFILEINT = STR2INT( NFILESTR )
               IF( NFILEINT == IMISS3 ) THEN
                   MESG = 'Invalid number of files in header of ' //
     &                    'file set "' // TRIM( ROOTNAME ) // '"'
                   CALL M3MSG2( MESG ) 
                   CALL CLEANUP( FILEIDX )
                   OPENSET = .FALSE.
                   RETURN
               END IF

!................  Get total number of variables from file header                
               NVARSTR  = GETCFDSC( FDESC3D, '/NUMBER OF VARIABLES/', 
     &                              .TRUE. )
               NVARINT  = STR2INT( NVARSTR )
               IF( NVARINT == IMISS3 ) THEN
                   MESG = 'Invalid number of variables in header of ' //
     &                    'file set "' // TRIM( ROOTNAME ) // '"'
                   CALL M3MSG2( MESG )
                   CALL CLEANUP( FILEIDX )
                   OPENSET = .FALSE.
                   RETURN
               END IF

!................  Allocate memory for logical file names and variables                      
               ALLOCATE( FILE_INFO( FILEIDX )%LNAMES( NFILEINT ), 
     &                   STAT=IOS )
               CALL CHECKMEM( IOS, 'FILE_INFO%LNAMES', FUNCNAME )
               ALLOCATE( FILE_INFO( FILEIDX )%VARS( NVARINT,2 ), 
     &                   STAT=IOS )
               CALL CHECKMEM( IOS, 'FILE_INFO%VARS', FUNCNAME )

!................  Store logical file names and variables                        
               FILE_INFO( FILEIDX )%LNAMES( 1 ) = LNAME
               FILE_INFO( FILEIDX )%VARS( 1:NVARS3D,1 ) = 
     &                                     VNAME3D( 1:NVARS3D )
               FILE_INFO( FILEIDX )%VARS( 1:NVARS3D,2 ) = LNAME
               VARPOS = NVARS3D + 1

!................  Loop through remaining files in set                   
               DO I = 2, NFILEINT

!....................  Create new logical and physical file names
                   WRITE( INTBUF, '(I2)' ) I
                   INTBUF = ADJUSTL( INTBUF )
                   LNAME = TRIM( ROOTNAME16 ) // TRIM( INTBUF )
                   CALL APPENDNAME( ENVNAME, I, TEMPNAME, EFLAG )
                   IF( EFLAG ) THEN
                       CALL CLEANUP( FILEIDX )
                       OPENSET = .FALSE.
                       RETURN
                   END IF

!....................  Set new env. variable                       
                   IF( .NOT. SETENVVAR( LNAME, TEMPNAME ) ) THEN
                       CALL CLEANUP( FILEIDX )
                       OPENSET = .FALSE.
                       RETURN
                   END IF

!....................  Try to open file
                   IF( .NOT. OPEN3( LNAME, FSTATUS, FUNCNAME ) ) THEN
                       CALL CLEANUP( FILEIDX )
                       OPENSET = .FALSE.
                       RETURN
                   END IF

!....................  Get file description                       
                   IF( .NOT. DESC3( LNAME ) ) THEN
                       CALL CLEANUP( FILEIDX )
                       OPENSET = .FALSE.
                       RETURN
                   END IF
                     
!....................  Store logical file names and variables                           
                   FILE_INFO(FILEIDX)%LNAMES( I ) = LNAME
                   FILE_INFO(FILEIDX)%VARS( VARPOS:VARPOS+NVARS3D-1,1 )
     &                                           = VNAME3D( 1:NVARS3D )
                   FILE_INFO(FILEIDX)%VARS( VARPOS:VARPOS+NVARS3D-1,2 )
     &                                           = LNAME
                   VARPOS = VARPOS + NVARS3D
                      
               END DO  ! loop over remaining files in file set

           ELSE
!................  File doesn't exist on disk
               MESG = ROOTNAME16 // ':' // TRIM( ENVNAME )
               CALL M3MSG2( MESG )
               CALL M3WARN( FUNCNAME, 0, 0, 'File not available.' )
               
               CALL CLEANUP( FILEIDX )
               OPENSET = .FALSE.
               RETURN
           END IF
       
       ELSE IF( FSTATUS == FSNEW3 ) THEN

!............  For status new, neither file should exist
           IF( FILEEXIST .OR. SETEXIST ) THEN
               MESG = ROOTNAME16 // ':' // TRIM( ENVNAME )
               CALL M3MSG2( MESG )
               MESG = 'File already exists on disk, cannot open ' //
     &                'as "NEW".'
               CALL M3WARN( FUNCNAME, 0, 0, MESG )
               CALL CLEANUP( FILEIDX )
               OPENSET = .FALSE.
               RETURN

!............  Otherwise, try to create new files
           ELSE
               IF( .NOT. CREATESET( FILEIDX, ROOTNAME16, ENVNAME, 
     &                              FSTATUS ) ) THEN
                   CALL CLEANUP( FILEIDX )
                   OPENSET = .FALSE.
                   RETURN
               END IF
           END IF

       ELSE IF( FSTATUS == FSUNKN3 ) THEN

!............  If file or set exists, try to open it 
           IF( FILEEXIST .OR. SETEXIST ) THEN

!................  Allocate memory for logical file names and variables                      
               ALLOCATE( FILE_INFO( FILEIDX )%LNAMES( NFILESET ), 
     &                   STAT=IOS )
               CALL CHECKMEM( IOS, 'FILE_INFO%LNAMES', FUNCNAME )
               ALLOCATE( FILE_INFO( FILEIDX )%VARS( NVARSET,2 ), 
     &                   STAT=IOS )
               CALL CHECKMEM( IOS, 'FILE_INFO%VARS', FUNCNAME )

!................  Loop through individual files in the set
               VARPOS = 1
               DO I = 1, NFILESET
               
!....................  Set necessary file description values                   
                   NVARS3D = VARS_PER_FILE( I )
                   VTYPE3D( 1:NVARS3D ) = 
     &                      VTYPESET( VARPOS:VARPOS + NVARS3D - 1 )
                   VNAME3D( 1:NVARS3D ) = 
     &                      VNAMESET( VARPOS:VARPOS + NVARS3D - 1 )

!........................  Create new logical and physical file names if needed
                   IF( NFILESET > 1 ) THEN
                       WRITE( INTBUF, '(I2)' ) I
                       INTBUF = ADJUSTL( INTBUF )
                       LNAME = TRIM( ROOTNAME16 ) // TRIM( INTBUF )
                       CALL APPENDNAME( ENVNAME, I, TEMPNAME, EFLAG )
                       IF( EFLAG ) THEN
                           CALL CLEANUP( FILEIDX )
                           OPENSET = .FALSE.
                           RETURN
                       END IF

!........................  Set new env. variable                       
                       IF( .NOT. SETENVVAR( LNAME, TEMPNAME ) ) THEN
                           CALL CLEANUP( FILEIDX )
                           OPENSET = .FALSE.
                           RETURN
                       END IF
                   ELSE
                       LNAME = ROOTNAME16
                   END IF

!....................  Try to open file                   
                   IF( .NOT. OPEN3( LNAME, FSTATUS, FUNCNAME ) ) THEN
                       CALL CLEANUP( FILEIDX )
                       OPENSET = .FALSE.
                       RETURN
                   END IF

!....................  Store logical file names and variables                        
                   FILE_INFO(FILEIDX)%LNAMES( I ) = LNAME
                   FILE_INFO(FILEIDX)%VARS( VARPOS:VARPOS+NVARS3D-1,1 )
     &                                           = VNAME3D( 1:NVARS3D )
                   FILE_INFO(FILEIDX)%VARS( VARPOS:VARPOS+NVARS3D-1,2 )
     &                                           = LNAME
                   VARPOS = VARPOS + NVARS3D

               END DO  ! loop over individual files

!............  Otherwise, try to create new files
           ELSE
               IF( .NOT. CREATESET( FILEIDX, ROOTNAME16, ENVNAME, 
     &                              FSTATUS ) ) THEN
                   CALL CLEANUP( FILEIDX )
                   OPENSET = .FALSE.
                   RETURN
               END IF
           END IF
       
       ELSE IF( FSTATUS == FSCREA3 ) THEN

!............  If file or set exists, delete it
           IF( FILEEXIST .OR. SETEXIST ) THEN

!................  Loop through files in set
               DO I = 1, NFILESET

!....................  If more than one file, generate subsequent file names
                   IF( NFILESET > 1 ) THEN
                       CALL APPENDNAME( ENVNAME, I, TEMPNAME, EFLAG )
                       IF( EFLAG ) THEN
                           CALL CLEANUP( FILEIDX )
                           OPENSET = .FALSE.
                           RETURN
                       END IF
                   ELSE
                       TEMPNAME = ENVNAME
                   END IF

!....................  Check that file exists
                   INQUIRE( FILE = TEMPNAME, EXIST = TEMPEXIST )

!....................  Try to delete files
                   IF( TEMPEXIST ) THEN
                       IOS = RMFILE( TEMPNAME )
                       IF( IOS /= 0 ) THEN
                           WRITE( MESG, 93010 ) 'Error number', IOS,
     &                          'removing file ' // TEMPNAME
                           CALL M3WARN( FUNCNAME, 0, 0, MESG )
                           CALL CLEANUP( FILEIDX )
                           OPENSET = .FALSE.
                           RETURN
                       END IF
                   END IF

               END DO
           END IF

!............  Try to create new files
           IF( .NOT. CREATESET( FILEIDX, ROOTNAME16, ENVNAME, 
     &                          FSTATUS ) ) THEN
               CALL CLEANUP( FILEIDX )
               OPENSET = .FALSE.
               RETURN
           END IF

!........  Illegal FSTATUS value
       ELSE

           MESG = 'File opening error: illegal FSTATUS argument.'
           CALL M3WARN( FUNCNAME, 0, 0, MESG )
           MESG = 'Legal values: 1-READONLY, 2-READ/WRITE, 3-NEW, ' //
     &            '4 -UNKNOWN'
           CALL M3MSG2( MESG )
           WRITE( MESG, 93010 ) 'Value supplied by caller:', FSTATUS
           CALL M3MSG2( MESG )

           CALL CLEANUP( FILEIDX )
           OPENSET = .FALSE.
           RETURN

       END IF

!........  Update number of open file sets
       NOPENSETS = NOPENSETS + 1
       
       OPENSET = .TRUE.
       
       RETURN

!---------- Format statements --------------

93010   FORMAT ( 5 ( A, :, I9, :, 2X ) )
                           
       END FUNCTION OPENSET
