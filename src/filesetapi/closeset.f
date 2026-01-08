
       LOGICAL FUNCTION CLOSESET( ROOTNAME )

!***********************************************************************
!  Function body starts at line 41
!
!  DESCRIPTION:
!     Closes a file set
!
!  PRECONDITIONS REQUIRED:
!     File set has been opened with OPENSET
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!     CLOSE3 - closes an individual file
!     CLEANUP - cleans up internal memory structures
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
c      INTEGER, EXTERNAL :: INDEX1
       
!........  Function arguments
       CHARACTER(*), INTENT(IN) :: ROOTNAME  ! logical file name for file set

!........  Local variables
       INTEGER            I           ! counter
       INTEGER            FILEIDX     ! file index
       
       CHARACTER(16)  ROOTNAME16  ! fixed length root file name
       CHARACTER(256) MESG        ! message buffer

!---------------------------------
!  Begin body of function CLOSESET
!---------------------------------

!........  Check length of file name
       IF( LEN( ROOTNAME ) > 16 ) THEN
           MESG = 'Max file name length (16) exceeded for "' // 
     &            TRIM( ROOTNAME ) // '"'
           CALL M3MSG2( MESG )
           CLOSESET = .FALSE.
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
           CLOSESET = .FALSE.
           RETURN
       END IF

!........  Loop through individual files
       DO I = 1, SIZE( FILE_INFO( FILEIDX )%LNAMES )
           IF( .NOT. CLOSE3( FILE_INFO( FILEIDX )%LNAMES( I ) ) ) THEN
               CLOSESET = .FALSE.
               RETURN
           END IF
       END DO

!........  Write message only if file set contains more than one file       
       IF( SIZE( FILE_INFO( FILEIDX )%LNAMES ) > 1 ) THEN       
           MESG = 'Closing file set "' // TRIM( ROOTNAME ) // '"'
           CALL M3MSG2( MESG )
       END IF      
 
       CALL CLEANUP( FILEIDX )
       NOPENSETS = NOPENSETS - 1

       CLOSESET = .TRUE.

       END FUNCTION CLOSESET
