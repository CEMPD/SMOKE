 
C Version @(#)$Id$ $Source$ $Date$ 
 
C.... LOGICAL FUNCTION CHKVAR checks to see if OUTVAR is in VARLIST, and if so
C     it returns a .TRUE. and a pointer value.

      LOGICAL FUNCTION CHKVAR( OUTVAR, VARLIST, NVARDIM, IPTR )

C.... Input variables
C       OUTVAR        ! desired UAM output variable from main routine
C       VARLIST       ! list of NetCDF input variables
C       NVARDIM       ! actual number of input variables
C.... Output variables
C       IPTR          ! pointer from output variable to input variable list

      IMPLICIT NONE

      INTEGER*4 I, IPTR, NVARDIM

      CHARACTER*10 OUTVAR
      CHARACTER*16 VARLIST(NVARDIM)

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
