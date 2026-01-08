
       LOGICAL FUNCTION READSET( ROOTNAME, VNAME, LAYER, FILENUM, 
     &                           JDATE, JTIME, BUFFER )

!***********************************************************************
!  Function body starts at line 41
!
!  DESCRIPTION:
!     Reads information from a file set
!
!  PRECONDITIONS REQUIRED:
!     File set opened by OPENSET
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!     READ3 - reads from an individual file
!
!  REVISION HISTORY:
!     Created 6/02 by C. Seppanen
!     09/2025 by HT UNC-IE: USe M3UTILIO     
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
c      INTEGER, EXTERNAL :: INDEX1
       
!........  Function arguments
       CHARACTER(*), INTENT(IN)  :: ROOTNAME  ! logical file name for file set
       CHARACTER(*), INTENT(IN)  :: VNAME     ! variable name or ALLVAR3
       INTEGER,      INTENT(IN)  :: LAYER     ! layer number or ALLAYS3
       INTEGER,      INTENT(IN)  :: FILENUM   ! position of file in set or ALLFILES
       INTEGER,      INTENT(IN)  :: JDATE     ! date (YYYYDDD)
       INTEGER,      INTENT(IN)  :: JTIME     ! time (HHMMSS)
       REAL,         INTENT(OUT) :: BUFFER(*) ! array to hold data

       INTEGER, PARAMETER :: ALLFILES = -1

!........  Local arrays
       INTEGER            TYPSIZE( 6 )           !  sizeof( variable ) / sizeof( real )
#if _CRAY || REAL8
       DATA               TYPSIZE / 1, 1, 1, 1, 1, 1 /
#endif
#if ! ( _CRAY || REAL8 )
       DATA               TYPSIZE / 1, 1, 1, 1, 1, 2 /
#endif

!........  Local variables
       INTEGER            I, J        ! counters
       INTEGER            FILEIDX     ! file index
       INTEGER            MEMIDX      ! index into BUFFER array
       INTEGER            NLAYERS     ! number of layers to read (either 1 or all)
       INTEGER            DELTA       ! increment amount

       CHARACTER(16)  ROOTNAME16      ! fixed length root file name
       CHARACTER(16)  VNAME16         ! fixed length variable name
       CHARACTER(16)  LNAME           ! temporary logical name
       CHARACTER(16)  SIGNAME         ! logical name for specific file
       CHARACTER(256) MESG            ! message buffer
       
       CHARACTER(16) :: FUNCNAME = 'READSET'  ! function name

!---------------------------------
!  Begin body of function READSET
!---------------------------------

!........  Check length of file name
       IF( LEN( ROOTNAME ) > 16 ) THEN
           MESG = 'Max file name length (16) exceeded for "' // 
     &            ROOTNAME // '"'
           CALL M3MSG2( MESG )
           READSET = .FALSE.
           RETURN
       END IF

!........  Check length of variable name
       IF( LEN( VNAME ) > 16 ) THEN
           MESG = 'Max variable name length (16) exceeded for "' // 
     &            VNAME // '"'
           CALL M3MSG2( MESG )
           READSET = .FALSE.
           RETURN
       END IF
       VNAME16 = VNAME

!........  Get file index
       ROOTNAME16 = ROOTNAME
       FILEIDX = INDEX1( ROOTNAME16, MXFILE3, RNAMES )
       
!........  If file is not open, exit with error
       IF( FILEIDX == 0 ) THEN
           MESG = 'File set "' // TRIM( ROOTNAME ) // '" is not ' //
     &            'currently open'
           CALL M3MSG2( MESG )
           READSET = .FALSE.
           RETURN
       END IF

!........  Check that layer number is valid
       IF( LAYER /= ALLAYS3 ) THEN
           IF( LAYER < 1 ) THEN
               MESG = 'Layer number must be positive or "ALLAYS3"'
               CALL M3MSG2( MESG )
               READSET = .FALSE.
               RETURN
           END IF
       END IF

!........  Check if a single file is requested
       IF( FILENUM /= ALLFILES ) THEN

!............  Make sure file number is positive
           IF( FILENUM < 1 ) THEN
               MESG = 'File number must be positive or "ALLFILES"'
               CALL M3MSG2( MESG )
               READSET = .FALSE.
               RETURN
               
!............  Check that number is not greater than total number of files
           ELSE IF( FILENUM > SIZE( FILE_INFO( FILEIDX )%LNAMES ) ) THEN
               WRITE( MESG,92010 ) 'Invalid file number requested; ' //
     &                'file set "' // TRIM( ROOTNAME ) // '" contains ',
     &                SIZE( FILE_INFO( FILEIDX )%LNAMES ), ' files '
               CALL M3MSG2( MESG )
               READSET = .FALSE.
               RETURN
               
!............  Otherwise, set logical file name
           ELSE
               SIGNAME = FILE_INFO( FILEIDX )%LNAMES( FILENUM )
           END IF
           
!............  If all variables are requested, read file
           IF( VNAME16 == ALLVAR3 ) THEN
               IF( .NOT. READ3( SIGNAME, VNAME, LAYER, 
     &                          JDATE, JTIME, BUFFER ) ) THEN
                   READSET = .FALSE.
                   RETURN
               ELSE
                   READSET = .TRUE.
                   RETURN
               END IF
           END IF           
       END IF

!........  Check that requested variable is in file set
       IF( VNAME16 /= ALLVAR3 ) THEN
           LNAME = ' '

!............  Loop over variables in file set
           DO I = 1, SIZE( FILE_INFO( FILEIDX )%VARS, 1 )

!................  If variable matches, set logical file name
               IF( FILE_INFO( FILEIDX )%VARS( I,1 ) == VNAME16 ) THEN
                   LNAME = FILE_INFO( FILEIDX )%VARS( I,2 )
                   EXIT
               END IF
           END DO

!............  Couldn't find variable name, exit with error
           IF( LNAME == ' ' ) THEN
               MESG = 'Requested variable "' // VNAME16 // '" not ' //
     &                'available in file set "' // ROOTNAME16 // '"'
               CALL M3WARN( FUNCNAME, JDATE, JTIME, MESG )
               READSET = .FALSE.
               RETURN
           END IF
           
!............  If a single file is requested, make sure it matches
           IF( FILENUM /= ALLFILES ) THEN
               IF( LNAME /= SIGNAME ) THEN
                   MESG = 'Requested variable "' // VNAME16 // '" ' //
     &                    'not in specified file of file set "' //
     &                    ROOTNAME16 // '"'
                   CALL M3MSG2( MESG )
                   READSET = .FALSE.
                   RETURN
               END IF
           END IF
           
!............  Try to read file
           IF( .NOT. READ3( LNAME, VNAME16, LAYER, 
     &                      JDATE, JTIME, BUFFER ) ) THEN
               READSET = .FALSE.
               RETURN
           END IF
           
       ELSE

           MEMIDX = 1

!............  Loop over individual files
           DO I = 1, SIZE( FILE_INFO( FILEIDX )%LNAMES )

!................  Set logical file name
               LNAME = FILE_INFO( FILEIDX )%LNAMES( I )

!................  Get description of file
               IF( .NOT. DESC3( LNAME ) ) THEN
                   READSET = .FALSE.
                   RETURN
               END IF
               
!................  Try to read from file
               IF( .NOT. READ3( LNAME, VNAME16, LAYER, 
     &                          JDATE, JTIME, BUFFER( MEMIDX ) ) ) THEN
                   READSET = .FALSE.
                   RETURN
               END IF

!................  Increment BUFFER index based on total amount of memory this file takes
               IF( LAYER == ALLAYS3 ) THEN
                   NLAYERS = NLAYS3D
               ELSE
                   NLAYERS = 1
               END IF

!................  For each file type, calculate amount of space one variable takes
               SELECT CASE( FTYPE3D )
               CASE( CUSTOM3 )
                   DELTA = NCOLS3D * NLAYERS
               CASE( GRDDED3, TSRIES3, PTRFLY3 )
                   DELTA = NCOLS3D * NROWS3D * NLAYERS
               CASE( BNDARY3 )
                   DELTA = 2*NTHIK3D * ( NCOLS3D + NROWS3D + 2*NTHIK3D )
                   DELTA = DELTA * NLAYERS
               CASE( IDDATA3 )               
                   MEMIDX = MEMIDX + 1 + NROWS3D    ! account for space used by number of sites
                                                    !  and site IDs 
                   DELTA = NROWS3D * NLAYERS
               CASE( PROFIL3 )
                   MEMIDX = MEMIDX + 1 + 5*NROWS3D  ! number of sites, site IDs,
                                                    !  profile-count list, and location lists
                   DELTA = NCOLS3D * NROWS3D * NLAYERS
               CASE( GRNEST3 )
                   MEMIDX = MEMIDX + 1 + 6*NROWS3D  ! no. sites, site IDs, profile-count list,
                                                    ! X,Y location lists, and cellsize lists
                   DELTA = NCOLS3D * NROWS3D * NLAYERS
               CASE( SMATRX3 )
                   MEMIDX = MEMIDX + NROWS3D + NCOLS3D  ! no. columns/row, column index
                   DELTA = NCOLS3D
               CASE DEFAULT
                   MESG = 'Cannot read all variables from file set ' //
     &                    TRIM( ROOTNAME ) // ' due to the file type'
                   CALL M3MSG2( MESG )
                   READSET = .FALSE.
                   RETURN
               END SELECT

!................  Loop through variables and increment index               
               DO J = 1, NVARS3D
                   MEMIDX = MEMIDX + ( DELTA * TYPSIZE( VTYPE3D( J ) ) )
               END DO
           END DO

       END IF

       READSET = .TRUE.

!---------- Format statements --------------

92010   FORMAT ( A, :, I3, :, A )

       END FUNCTION READSET
