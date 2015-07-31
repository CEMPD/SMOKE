
       INTEGER FUNCTION CVTVEHTYPE( SMKVEH )

C***********************************************************************
C  subroutine body starts at line 53
C
C  DESCRIPTION:
C       Converts inventory vehicle type to MOBILE6 type
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

C.........  Function arguments
       INTEGER, INTENT (IN) :: SMKVEH   ! vehicle type in SMOKE code

C...........   LOCAL VARIABLES and their descriptions:
       LOGICAL :: EFLAG      = .FALSE.   ! true: error found

       CHARACTER(300)          MESG      !  message buffer 
        
       CHARACTER(16) :: PROGNAME = 'CVTVEHTYPE' ! program name
        
C***********************************************************************
C   begin body of function CVTVEHTYPE

       SELECT CASE( SMKVEH )
       CASE( 100 )              ! light duty gasoline vehicles
           CVTVEHTYPE = 1
       CASE( 102 )              ! light duty gasoline trucks 1
           CVTVEHTYPE = 2
       CASE( 104 )              ! light duty gasoline trucks 2
           CVTVEHTYPE = 3
       CASE( 107 )              ! heavy duty gasoline vehicles
           CVTVEHTYPE = 4
       CASE( 3000 )             ! light duty diesel vehicles
           CVTVEHTYPE = 5
       CASE( 3006 )             ! light duty diesel trucks
           CVTVEHTYPE = 6
       CASE( 3007 )             ! heavy duty diesel vehicles
           CVTVEHTYPE = 7
       CASE( 108 )              ! motorcycles
           CVTVEHTYPE = 8
       CASE DEFAULT                   
           EFLAG = .TRUE.
 
           WRITE( MESG, 94010 ) 'ERROR: Vehicle type ', SMKVEH,
     &                          'is not recognized.'
           CALL M3MESG( MESG )
       END SELECT

       IF( EFLAG ) THEN       
           MESG = 'Problem converting vehicle type'
           CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )       
       END IF

       RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )  
      
C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )
       
       END FUNCTION CVTVEHTYPE
