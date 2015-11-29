
        SUBROUTINE GETBASYR( NSRC, BASEYEAR )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION: 
C     This subroutine examines the years for all of the sources, and determines
C     the base year from the most frequent of these base years.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:  M3EXIT
C
C  REVISION  HISTORY:
C       Created 4/99 by M Houyoux
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
C****************************************************************************

C...........   MODULES for public variables
C...........   This module is the inventory arrays
        USE MODSOURC, ONLY: INVYR

        IMPLICIT NONE

C...........   INCLUDES:
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   ARGUMENTS and their descriptions:
        INTEGER, INTENT (IN) :: NSRC        ! number of sources
        INTEGER, INTENT(OUT) :: BASEYEAR    ! base inventory year

C...........   EXTERNAL FUNCTIONS:
        CHARACTER(2)    CRLF
        EXTERNAL        CRLF

C...........   Argument dimensioned source arrrays
        INTEGER, ALLOCATABLE :: INDX ( : )  ! sorting index for processing years
        INTEGER, ALLOCATABLE :: YRGRP( : )  ! group number of year

C...........   LOCAL VARIABLES their descriptions:
        INTEGER         J, S       ! indices
        INTEGER         ICNT       ! counter
        INTEGER         IOS        ! i/o status
        INTEGER         GRPMAX     ! count of most prevalent year
        INTEGER         NYEARS     ! number of different years in inventory
        INTEGER         PGRP       ! year group from previous loop iteration
        INTEGER         PYEAR      ! year from previous loop iteration

        CHARACTER(300)  MESG       ! Message buffer
  
        CHARACTER(16) :: PROGNAME = 'GETBASYR'    ! Program name

C***********************************************************************
C   begin body of subroutine GETBASYR

C.........  Allocate memory for local arrays
        ALLOCATE( INDX( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INDX', PROGNAME )
        ALLOCATE( YRGRP( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'YRGRP', PROGNAME )

C.........  Create sorting index for years
        DO S = 1, NSRC

            INDX( S ) = S

        END DO

C.........  Sort years (from MODSOURC)
        CALL SORTI1( NSRC, INDX, INVYR )        

C.........  Count the number of different years in inventory
        PYEAR = 0
        ICNT  = 0
        DO S = 1, NSRC
         
            J = INDX( S )

            IF( INVYR( J ) .NE. PYEAR ) THEN

                ICNT = ICNT + 1

                PYEAR      = INVYR( J )

            END IF

            YRGRP( S ) = ICNT

        END DO

C.........  Post-process list of year groups to find most prevalent year
        PGRP   = 1
        ICNT   = 0
        GRPMAX = 0
        NYEARS = 0
	IF ( NSRC .EQ. 1) THEN

            BASEYEAR = PYEAR
	    NYEARS = 1
	    PGRP   = YRGRP( 1 )
            PYEAR = INVYR( INDX( 1 ) )	   	

	ELSE
            DO S = 1, NSRC

                IF( YRGRP( S ) .NE. PGRP .OR. S .EQ. NSRC ) THEN

                    IF( ICNT .GT. GRPMAX ) THEN
                        GRPMAX = ICNT
                        BASEYEAR = PYEAR
                    END IF

                    NYEARS = NYEARS + 1
                    ICNT   = 0
                    PGRP   = YRGRP( S ) 

                END IF

                PYEAR = INVYR( INDX( S ) )
                ICNT  = ICNT + 1

            END DO
        ENDIF
	
        IF( NYEARS .EQ. 1 ) THEN

            WRITE( MESG,94010 ) 
     &             'NOTE: Inventory base year set to ', BASEYEAR

        ELSE

            WRITE( MESG,94010 ) 
     &             'NOTE: The number of inventory years was ', NYEARS, 
     &             'and the inventory base' // CRLF() // BLANK10 //
     &             'year was set to the most frequent year, ', BASEYEAR

        END IF

        CALL M3MSG2( MESG )

C.........  Deallocate local memory
        DEALLOCATE( INDX, YRGRP )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats.............94xxx
 
94010   FORMAT( 10( A, :, I6, :, 2X ) )

        END SUBROUTINE GETBASYR
