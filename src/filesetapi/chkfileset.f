
       LOGICAL FUNCTION CHKFILESET( FIDX )

!***********************************************************************
!  Function body starts at line 41
!
!  DESCRIPTION:
!     Checks a file description against an already opened file.
!
!  PRECONDITIONS REQUIRED:
!     File description information has been filled out in FDESC3.EXT 
!     and MODFILESET
!     File set has already been opened by OPENSET
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION HISTORY:
!     Created 6/02 by C. Seppanen
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

!........  Modules for public variables
       USE MODFILESET
       
       IMPLICIT NONE

!........  Function arguments
       INTEGER, INTENT(IN) :: FIDX  ! file index
              
!........  Local variables
       INTEGER TOTALFILES  ! total number of files
       INTEGER TOTALVARS   ! total number of variables
       
       CHARACTER(16)  FNAME        ! logical file name
       CHARACTER(300) MESG         ! message buffer
       
       LOGICAL :: EFLAG = .FALSE.  ! true: error found
       
!------------------------------------
!  Begin body of function CHKFILESET
!------------------------------------
       
       EFLAG = .FALSE.
       FNAME = RNAMES( FIDX )
       
!........  Check that total number of files is consistent
       TOTALFILES = SIZE( FILE_INFO( FIDX )%LNAMES )
       IF( TOTALFILES /= NFILESET ) THEN
           WRITE( MESG,91010 )
           CALL M3MSG2( MESG )
           WRITE( MESG,91020 ) 'Inconsistent number of files ' //
     &                         'for file set ' // TRIM( FNAME )
           CALL M3MSG2( MESG )
           WRITE( MESG,91030 ) 'Value from file:  ', TOTALFILES
           CALL M3MSG2( MESG )
           WRITE( MESG,91030 ) 'Value from caller:', NFILESET
           CALL M3MSG2( MESG )
           EFLAG = .TRUE.
       END IF

!........  Check that total number of vars is consistent
       TOTALVARS = SIZE( FILE_INFO( FIDX )%VARS,1 )
       IF( TOTALVARS /= NVARSET ) THEN
           WRITE( MESG,91010 )
           CALL M3MSG2( MESG )
           WRITE( MESG,91020 ) 'Inconsistent number of variables ' //
     &                         'for file set ' // TRIM( FNAME )
           CALL M3MSG2( MESG )
           WRITE( MESG,91030 ) 'Value from file:  ', TOTALVARS
           CALL M3MSG2( MESG )
           WRITE( MESG,91030 ) 'Value from caller:', NVARSET
           CALL M3MSG2( MESG )
           EFLAG = .TRUE. 
       END IF

!........  Set function value and return
       IF( EFLAG ) THEN
           CHKFILESET = .FALSE.
       ELSE
           CHKFILESET = .TRUE.
       END IF

       RETURN

!---------- Format statements --------------

91010   FORMAT ( 5X, '>>> WARNING in function CHKFILESET <<<' )

91020   FORMAT ( 5X , A )

91030   FORMAT ( 5X , A , I4 )
       
       END FUNCTION CHKFILESET
