
       INTEGER FUNCTION CVTRDTYPE( SMKROAD, RLASAFLAG, ULASAFLAG )

C***********************************************************************
C  function body starts at line 68
C
C  DESCRIPTION:
C       Converts inventory road types into MOBILE6 type
C
C  PRECONDITIONS REQUIRED: none
C
C  SUBROUTINES AND FUNCTIONS CALLED: none
C
C  REVISION  HISTORY:
C     10/01: Created by C. Seppanen
C
C***********************************************************************
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
C***********************************************************************
       
       IMPLICIT NONE

C........  Includes
       INCLUDE 'M6CNST3.EXT'  ! MOBILE6 constants

C........  Function arguments
       INTEGER, INTENT (IN) :: SMKROAD   ! road type in SMOKE code ( 1 - 19 )
       LOGICAL, INTENT (IN) :: RLASAFLAG ! true: treat rural local roads as arterial
       LOGICAL, INTENT (IN) :: ULASAFLAG ! true: treat urbal local roads as arterial

C........  Road type parameters
       INTEGER, PARAMETER :: RURALINTERSTATE = 1   ! rural interstate
       INTEGER, PARAMETER :: RURALPRINCART   = 2   ! rural principle arterial
       INTEGER, PARAMETER :: RURALMINORART   = 6   ! rural minor arterial
       INTEGER, PARAMETER :: RURALMAJORCOLL  = 7   ! rural major collector
       INTEGER, PARAMETER :: RURALMINORCOLL  = 8   ! rural minor collector
       INTEGER, PARAMETER :: RURALLOCAL      = 9   ! rural local
       INTEGER, PARAMETER :: URBANINTERSTATE = 11  ! urban interstate
       INTEGER, PARAMETER :: URBANFREEWAY    = 12  ! urban freeway
       INTEGER, PARAMETER :: URBANPRINCART   = 14  ! urban principle arterial
       INTEGER, PARAMETER :: URBANMINORART   = 16  ! urban minor arterial
       INTEGER, PARAMETER :: URBANCOLL       = 17  ! urban collector
       INTEGER, PARAMETER :: URBANLOCAL      = 19  ! urban local

C........  LOCAL VARIABLES and their descriptions:
       LOGICAL :: EFLAG      = .FALSE.   ! true: error found

       CHARACTER(300)          MESG      !  message buffer 
        
       CHARACTER(16) :: PROGNAME = 'CVTRDTYPE' ! program name
        
C***********************************************************************
C   begin body of function CVTRDTYPE

       SELECT CASE( SMKROAD )
       CASE( RURALINTERSTATE ) 
           CVTRDTYPE = M6FREEWAY
       CASE( RURALPRINCART:RURALMINORCOLL )    ! rural arterials through collectors
           CVTRDTYPE = M6ARTERIAL
       CASE( RURALLOCAL )
           IF( RLASAFLAG ) THEN
               CVTRDTYPE = M6ARTERIAL
           ELSE
               CVTRDTYPE = M6LOCAL
           END IF
       CASE( URBANINTERSTATE:URBANFREEWAY )    ! urban interstate and freeway
           CVTRDTYPE = M6FREEWAY
       CASE( URBANPRINCART:URBANCOLL )         ! urban arterials through collectors
           CVTRDTYPE = M6ARTERIAL
       CASE( URBANLOCAL )
           IF( ULASAFLAG ) THEN
               CVTRDTYPE = M6ARTERIAL
           ELSE
               CVTRDTYPE = M6LOCAL
           END IF
       CASE DEFAULT                   
           EFLAG = .TRUE.
 
           WRITE( MESG, 94010 ) 'ERROR: Road type ', SMKROAD,
     &                          'is not recognized.'
           CALL M3MESG( MESG )
       END SELECT

       IF( EFLAG ) THEN       
           MESG = 'Problem converting road type'
           CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )       
       END IF

       RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )  
      
C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )
       
       END FUNCTION CVTRDTYPE
