 
      LOGICAL FUNCTION CHKVAR( OUTVAR, VARLIST, NVARDIM, IPTR )

C***************************************************************************
C  function body starts at line 53
C
C  DESCRIPTION:
C     LOGICAL FUNCTION CHKVAR checks to see if OUTVAR is in VARLIST, and if so
C     it returns a .TRUE. and a pointer value.
C  SUBROUTINES AND FUNCTIONS CALLED:
C       Models-3 I/O
C
C.... Input variables
C       OUTVAR        ! desired UAM output variable from main routine
C       VARLIST       ! list of NetCDF input variables
C       NVARDIM       ! actual number of input variables
C.... Output variables
C       IPTR          ! pointer from output variable to input variable list
C
C***************************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1996, MCNC--North Carolina Supercomputing Center
C All Rights Reserved
C
C See file COPYRIGHT for conditions of use.
C
C Environmental Programs Group
C MCNC--North Carolina Supercomputing Center
C P.O. Box 12889
C Research Triangle Park, NC  27709-2889
C
C env_progs@mcnc.org
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C***************************************************************************

      IMPLICIT NONE

      INTEGER*4 I, IPTR, NVARDIM

      CHARACTER*10 OUTVAR
      CHARACTER*16 VARLIST(NVARDIM)

C***********************************************************************
C   begin body of function CHKVAR

      CHKVAR = .FALSE.
      IPTR = 0
      I = 1
      DOWHILE ( .NOT. CHKVAR .AND. I .LE. NVARDIM )
       IF(OUTVAR//'      ' .EQ. VARLIST(I)) THEN
        CHKVAR = .TRUE.
        IPTR = I
       ENDIF

       I = I + 1 
      ENDDO

      RETURN
      END
