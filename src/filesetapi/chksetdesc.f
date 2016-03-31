       
       LOGICAL FUNCTION CHKSETDESC( FNAME )

!***********************************************************************
!  Function body starts at line 39
!
!  DESCRIPTION:
!     Checks a file description for internal consistency.
!
!  PRECONDITIONS REQUIRED:
!     File description information in MODFILESET has been entered
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
       CHARACTER(16), INTENT(IN) :: FNAME  ! logical file name
              
!........  Local variables
       INTEGER              I                        ! counter
       INTEGER              IOS                      ! I/O status
       INTEGER              TOTALVARS                ! total number of variables

       REAL                 MXVARS                   ! real MXVARS3 from ioapi param
 
       LOGICAL ::           EFLAG = .FALSE.          ! true: error found
       
       CHARACTER(300)   MESG                     ! message buffer
       
       CHARACTER(16) :: FUNCNAME = 'CHKSETDESC'  ! function name
       
!------------------------------------
!  Begin body of function CHKSETDESC
!------------------------------------
       
       EFLAG = .FALSE.
       MXVARS = MXVARS3

!........  Check that total number of variables is valid
       IF( NVARSET <= 0 ) THEN
           MESG = 'Illegal number of variable for file set ' //
     &            TRIM( FNAME )
           CALL M3MSG2( MESG )
           CHKSETDESC = .FALSE.
           RETURN
       END IF

!........  Check that number of variables per file is set
       IF( .NOT. ALLOCATED( VARS_PER_FILE ) ) THEN
           MESG = 'Number of variables per file array is not ' //
     &            'allocated for file set ' // TRIM( FNAME ) // ';'
           CALL M3MSG2( MESG )
           WRITE( MESG,'( 10(A,I5) )' ) 
     &         'using default of ',MXVARS3,' variables per file'
           CALL M3MSG2( MESG )
           
           NFILESET = CEILING( NVARSET / MXVARS )  ! real division, convert back to int
           ALLOCATE( VARS_PER_FILE( NFILESET ), STAT=IOS )
           CALL CHECKMEM( IOS, 'VARS_PER_FILE', FUNCNAME )

           TOTALVARS = NVARSET
           I = 1
           DO
               IF( TOTALVARS - MXVARS3 <= 0 ) THEN
                   VARS_PER_FILE( I ) = TOTALVARS
                   EXIT
               ELSE
                   VARS_PER_FILE( I ) = MXVARS3
               END IF
               
               TOTALVARS = TOTALVARS - MXVARS3 
               I = I + 1
           END DO

       END IF
       
!........  Check that total number of files is consistent
       IF( SIZE( VARS_PER_FILE ) /= NFILESET ) THEN
           MESG = 'Inconsistent number of files in file set ' // 
     &            TRIM( FNAME )
           CALL M3MSG2( MESG )
           EFLAG = .TRUE.
       END IF

!........  Count total number of variables in VARS_PER_FILE;
!..        also check that number of variables per file is not more than MXVARS3
       TOTALVARS = 0
       DO I = 1, NFILESET
           TOTALVARS = TOTALVARS + VARS_PER_FILE( I )
           IF( VARS_PER_FILE( I ) > MXVARS3 ) THEN
               WRITE( MESG,'( 10(A,I5) )' ) 
     &             'More than',MXVARS3,' variables per file '//
     &             'in file set ' // TRIM( FNAME )
               CALL M3MSG2( MESG )
               EFLAG = .TRUE.
           END IF
       END DO

!........  Check that total number of vars is consistent
       IF( TOTALVARS /= NVARSET ) THEN
           MESG = 'Inconsistent number of variables in file set ' //
     &            TRIM( FNAME )
           CALL M3MSG2( MESG )
           EFLAG = .TRUE. 
       END IF

!........  Make sure variable information is set
       IF( .NOT. ALLOCATED( VTYPESET ) .OR. 
     &     .NOT. ALLOCATED( VNAMESET ) .OR.
     &     .NOT. ALLOCATED( VUNITSET ) .OR.
     &     .NOT. ALLOCATED( VDESCSET ) ) THEN
           MESG = 'One or more variable information arrays are ' //
     &            'not allocated for file set ' // TRIM( FNAME )
           CALL M3MSG2( MESG )
           CHKSETDESC = .FALSE.
           RETURN
       END IF
       
!........  Check size of variable information arrays
       IF( SIZE( VTYPESET ) < NVARSET .OR.
     &     SIZE( VNAMESET ) < NVARSET .OR.
     &     SIZE( VUNITSET ) < NVARSET .OR.
     &     SIZE( VDESCSET ) < NVARSET ) THEN
           MESG = 'One or more variable information arrays are ' //
     &            'not as large as the number of variables in ' //
     &            'file set ' // TRIM( FNAME )
           CALL M3MSG2( MESG )
           EFLAG = .TRUE.
       END IF

!........  Set function value and return
       IF( EFLAG ) THEN
           CHKSETDESC = .FALSE.
       ELSE
           CHKSETDESC = .TRUE.
       END IF
       
       END FUNCTION CHKSETDESC
