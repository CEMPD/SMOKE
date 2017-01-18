        LOGICAL FUNCTION FLTERR( P, Q )

C***********************************************************************
C  function body starts at line 52
C
C  DESCRIPTION:
C      This function determines if two floating point numbers are
C      "definitely unequal". It compares the square of the
C      normalized difference against the square of the tolerance and
C      returns TRUE if and only if the numbers are significantly
C      different.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created by C. Seppanen 6/04 from FLTERR.EXT
C
C**************************************************************************
C
C Project Title: EDSS Tools Library
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
C Pathname: $Source: /afs/isis/depts/cep/emc/apps/archive/edss_tools/edss_tools/src/lib/ioapi_grd_size.f
C Last updated: $Date$ 
C
C***************************************************************************

        IMPLICIT NONE

C.........  Function arguments        
        REAL, INTENT(IN) :: P
        REAL, INTENT(IN) :: Q

C***********************************************************************
C   begin body of function FLTERR
        
        FLTERR = ((P - Q)**2 .GT. 1.0E-11*( P*P + Q*Q + 1.0E-5 ))

        RETURN
        
        END FUNCTION FLTERR

C......................................................................        
        
        LOGICAL FUNCTION DBLERR( PD, QD )
        
C***********************************************************************
C  function body starts at line 85
C
C  DESCRIPTION:
C      This function is the same as FLTERR but works on double
C      floating point numbers.
C
C  PRECONDITIONS REQUIRED:

C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created by C. Seppanen 6/04 from FLTERR.EXT
C
C**************************************************************************
        
        IMPLICIT NONE

C.........  Function arguments        
        REAL(8), INTENT(IN) :: PD
        REAL(8), INTENT(IN) :: QD

C***********************************************************************
C   begin body of function FLTERR

        DBLERR = ((PD - QD)**2 .GT. 1.0D-12*( PD*PD + QD*QD + 1.0D-5 ))

        RETURN
        
        END FUNCTION DBLERR
