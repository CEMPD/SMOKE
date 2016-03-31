
       CHARACTER(16) FUNCTION PROMPTSET( PROMPT, FMODE, DEFAULT,
     &                                   CALLER )

!***********************************************************************
!  Function body starts at line 39
!
!  DESCRIPTION:
!     Prompts the user for a logical file name, then tries to open
!     specified file set
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!     OPENSET - opens a file set
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

!........  External functions
       LOGICAL, EXTERNAL :: ENVYN
       LOGICAL, EXTERNAL :: GETYN
       LOGICAL, EXTERNAL :: OPENSET
       
!........  Function arguments
       CHARACTER(*), INTENT(IN) :: PROMPT   ! prompt for user
       INTEGER,      INTENT(IN) :: FMODE    ! file opening mode
       CHARACTER(*), INTENT(IN) :: DEFAULT  ! default logical file name
       CHARACTER(*), INTENT(IN) :: CALLER   ! name of calling program

!........  Local parameters
       CHARACTER(16), PARAMETER :: NONE16 = 'NONE'
       
!........  Local variables
       INTEGER            IOS                ! I/O status
       INTEGER            IDX                ! index in logical file name

       LOGICAL, SAVE ::   PROMPTON           ! true: prompt for input
       LOGICAL, SAVE ::   INITIAL = .TRUE.   ! true: first time
       LOGICAL            NFLAG              ! true: "NONE" is in the prompt 
       
       CHARACTER(16)  LNAME                  ! logical file name
       CHARACTER(300) MESG                   ! message buffer
       CHARACTER(512) BUFFER                 ! prompt buffer
       
       CHARACTER(16) :: FUNCNAME = 'PROMPTSET'  ! function name
       
!------------------------------------
!  Begin body of function PROMPTSET
!------------------------------------

!........  On first time, check if prompt should be shown
       IF( INITIAL ) THEN
           PROMPTON = ENVYN( 'PROMPTFLAG', 'Prompt for input flag',
     &                       .TRUE., IOS )
           INITIAL = .FALSE.
       END IF

!........  Decide if 'NONE' is a valid response
       NFLAG = ( INDEX( PROMPT, '"NONE"' ) > 0 )
       
       IF( PROMPTON ) THEN

!............  Construct actual prompt
           WRITE( BUFFER,94000 ) TRIM( PROMPT ), ' [',
     &                           TRIM( DEFAULT ), '] >>'

!............  Loop until valid name is given and file set is opened
           DO
               LNAME = ' '
               WRITE( *,95000 ) TRIM( BUFFER ) // ' '
               READ ( *,'(A16)',IOSTAT=IOS ) LNAME

!................  Check if response was read               
               IF( IOS > 0 ) THEN
                   MESG = 'Could not read your response'
                   CALL M3MSG2( MESG )
                   IF( GETYN( 'Try again?', .TRUE. ) ) THEN
                       CYCLE
                   ELSE
                       MESG = 'Could not read logical name for file set'
                       CALL M3EXIT( FUNCNAME, 0, 0, MESG, 2 )
                   END IF
               ELSE IF( IOS < 0 ) THEN
                   MESG = 'Ending program "' // TRIM( CALLER ) // '".'
                   CALL M3EXIT( CALLER, 0, 0, MESG, 2 )
               END IF

!................  Check logical name               
               IDX = INDEX( LNAME, '!' )
               IF( IDX > 0 ) THEN
                   LNAME( IDX:LEN( LNAME ) ) = ' '
               END IF
               
               IF( LNAME == ' ' ) THEN
                   LNAME = DEFAULT
               END IF
               
               IF( NFLAG .AND. ( LNAME == NONE16 ) ) THEN
                   PROMPTSET = NONE16
                   RETURN
               END IF
               
!................  Try to open file set
               IF( .NOT. OPENSET( LNAME, FMODE, CALLER ) ) THEN 
                   MESG = 'Could not open file set "' //
     &                    TRIM( LNAME ) // '".'
                   CALL M3MSG2( MESG )
                   IF( GETYN( 'Try again?', .TRUE. ) ) THEN
                       CYCLE
                   ELSE
                       MESG = 'Ending program "' // 
     &                 TRIM( CALLER ) // '".'
                       CALL M3EXIT( CALLER, 0, 0, MESG, 2 )
                   END IF
               ELSE
                   EXIT
               END IF
           END DO

!........  Otherwise, don't prompt for input
       ELSE
           LNAME = DEFAULT

!............  Check if NONE is valid
           IF( NFLAG ) THEN
               IF( LNAME == NONE16 ) THEN
                   PROMPTSET = NONE16
                   RETURN
               END IF

!................  Check if logical file name is set               
               CALL ENVSTR( LNAME, 'Input file name', ' ', BUFFER, IOS )

!................  If not set (-2) or empty (-1)               
               IF( IOS < 0 ) THEN
                   PROMPTSET = NONE16
                   RETURN
               END IF
           END IF

!............  Try to open file set           
           IF( .NOT. OPENSET( LNAME, FMODE, CALLER ) ) THEN
               MESG = 'Could not open file set "' // 
     &                TRIM( LNAME ) // '".'
               CALL M3MSG2( MESG )
               CALL M3EXIT( CALLER, 0, 0, MESG, 2 )
           END IF
       END IF
       
       PROMPTSET = LNAME
       RETURN

!---------- Format statements --------------

94000  FORMAT( 12( A, : ) )

95000  FORMAT ( /5X, A )          !  generic prompt format.

       END FUNCTION PROMPTSET
