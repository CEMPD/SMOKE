
        SUBROUTINE GETINAME( CATEGORY, ENAME, ANAME )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION: 
C     This subroutine sets the inventory file names based on the source category
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C       Created 3/99 by M Houyoux
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

        IMPLICIT NONE

C...........   ARGUMENTS and their descriptions:
        CHARACTER(*), INTENT (IN) :: CATEGORY  ! category name
        CHARACTER(*), INTENT(OUT) :: ENAME     ! i/o api inventory file name
        CHARACTER(*), INTENT(OUT) :: ANAME     ! ASCII inventory file name

C...........   LOCAL VARIABLES their descriptions:

        CHARACTER(300)  MESG                   ! Message buffer

        CHARACTER(16) :: PROGNAME = 'GETINAME'  ! Program name

C***********************************************************************
C   begin body of subroutine GETINAME

        SELECT CASE( CATEGORY )

        CASE( 'AREA' )

            ENAME = 'AREA'
            ANAME = 'ASRC'

        CASE( 'BIOGEN' )
 
            ENAME = 'BGRD'
            ANAME = '    '

        CASE( 'MOBILE' )

            ENAME = 'MOBL'
            ANAME = 'MSRC'

        CASE( 'POINT' )

            ENAME = 'PNTS'
            ANAME = 'PSRC'

        CASE DEFAULT
            MESG = 'Category ' // CATEGORY // ' not known in program.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END SELECT

        RETURN

        END SUBROUTINE GETINAME
