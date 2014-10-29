
        LOGICAL FUNCTION USEEXPGEO()

C***********************************************************************
C  function USEEXPGEO body starts at line 50
C
C  DESCRIPTION:
C      The first time this function is called it will get the value
C      of the environment variable USE_EXP_GEO_CODES to indicate if 
C      expanded geographic codes should be used. Subsequent calls
C      will return the saved setting.
C
C  REVISION  HISTORY:
C     Created 12/13 by C. Seppanen
C
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2005, Environmental Modeling for Policy Development
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
C****************************************************************************

        IMPLICIT NONE

C...........   EXTERNAL FUNCTIONS and their descriptions
        LOGICAL        ENVYN
        
        EXTERNAL       ENVYN

C...........   Other local variables
        INTEGER        IOS    ! i/o status
        
        LOGICAL, SAVE :: FIRSTTIME = .TRUE.
        LOGICAL, SAVE :: EXPGEOFLAG = .FALSE.

        CHARACTER(256) MESG        ! message buffer

        CHARACTER(16) :: PROGNAME = 'USEEXPGEO' ! program name

C****************************************************************************
C   begin body of function USEEXPGEO

        IF ( FIRSTTIME ) THEN
            MESG = 'Use expanded geographic codes'
            EXPGEOFLAG = ENVYN( 'USE_EXP_GEO_CODES', MESG, .FALSE., IOS )
            
            FIRSTTIME = .FALSE.
        END IF
        
        USEEXPGEO = EXPGEOFLAG
        
        RETURN
 
        END FUNCTION USEEXPGEO
