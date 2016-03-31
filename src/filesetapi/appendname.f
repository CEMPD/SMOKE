
       SUBROUTINE APPENDNAME( ORIGNAME, FILENUM, NEWNAME, EFLAG )

!***********************************************************************
!  Subroutine body starts at line 38
!
!  DESCRIPTION:
!     Appends an integer to the given file name. Inserts number before
!     the file extension; for example, asmat_l.test.ncf -> asmat_l.test.1.ncf
!
!  PRECONDITIONS REQUIRED:
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
       
       IMPLICIT NONE

!........  Subroutine arguments
       CHARACTER(*), INTENT(IN)  :: ORIGNAME  ! original file name
       INTEGER,      INTENT(IN)  :: FILENUM   ! number to append to name
       CHARACTER(*), INTENT(OUT) :: NEWNAME   ! new file name
       LOGICAL,      INTENT(OUT) :: EFLAG     ! true if error occurs

!........  Local arguments
       INTEGER          IDX     ! string index
       INTEGER          INTLEN  ! length of integer string
       INTEGER          ORGLEN  ! length of original name
       INTEGER          NEWLEN  ! length of new name
       INTEGER          TMPLEN  ! length of temporary name
       CHARACTER(2)     INTBUF  ! buffer for integer value

!--------------------------------------
!  Begin body of subroutine APPENDNAME
!--------------------------------------

       EFLAG = .FALSE.

!........  Find .ncf extension by looking from the end of the string
       IDX = INDEX( ORIGNAME, '.', .TRUE. )

!........  Check that a location was found       
       IF( IDX == 0 ) THEN
           EFLAG = .TRUE.
           RETURN
       END IF

!........  Make sure NEWNAME is as long as ORIGNAME
       IF( LEN( NEWNAME ) < LEN( ORIGNAME ) ) THEN
           EFLAG = .TRUE.
           RETURN
       END IF

!........  Check that FILENUM is two characters or less
       IF( FILENUM < 1 .OR. FILENUM >= 100 ) THEN
           EFLAG = .TRUE.
           RETURN
       END IF

!........  Write integer to string
       WRITE( INTBUF, '(I2)' ) FILENUM
       INTBUF = ADJUSTL( INTBUF )
       INTLEN = LEN_TRIM( INTBUF )

!........  Get length of names
       ORGLEN = LEN( ORIGNAME )
       NEWLEN = LEN( NEWNAME )
 
!........  Create new name from original name      
       NEWNAME( 1:IDX-1 ) = ORIGNAME( 1:IDX-1 )
       NEWNAME( IDX:IDX+INTLEN ) = '.' // TRIM( INTBUF )
       
       TMPLEN = IDX+INTLEN
       NEWNAME( TMPLEN+1:NEWLEN ) = ORIGNAME( IDX:ORGLEN-INTLEN-1 )

       END SUBROUTINE APPENDNAME
        
