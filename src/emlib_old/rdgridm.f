
        SUBROUTINE RDGRIDM( CATEGORY, FNAME, 
     &                      RDIM, CDIM, TDIM, NA, IA, MA )

C***********************************************************************
C  subroutine body starts at line 75
C
C  DESCRIPTION:
C       This subroutine reads the gridding matrix for area, mobile, or point
C       sources, and compares the dimensions of the matrix with the 
C       compiled dimensions.
C
C  PRECONDITIONS REQUIRED:
C       NA must be dimensioned by RDIM and IA/MA must be dimensioned by CDIM
C       and put in COMMON statement in calling program.  File name FNAME
C       must be defined and opened.
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       TRIMLEN, DESC3
C
C  REVISION  HISTORY:
C       Started 1/98 by M Houyoux
C
C***********************************************************************
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
C************************************************************************

        IMPLICIT NONE

C.........  Include files
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  0-based I/O API file desc. data structures

C.........  External functions
        INTEGER         TRIMLEN

        EXTERNAL        TRIMLEN

C.........  Subroutine arguments and their descriptions

        CHARACTER*6     CATEGORY     ! 'AREA', 'MOBILE', or 'POINT'
        CHARACTER*16    FNAME        ! gridding matrix file name
        INTEGER         RDIM         ! matrix row dimension
        INTEGER         CDIM         ! matrix column dimension
        INTEGER         TDIM         ! source number check (stored in NTHIK3D)
        INTEGER         NA( RDIM )   ! sparse gridding matrix
        INTEGER         IA( CDIM ) 
        REAL            MA( CDIM )

C.........  Other local variables

        CHARACTER*256   MESG
 
C***********************************************************************
C   begin body of program RDGRIDM

C.........  Read gridding matrix

        IF( .NOT. READ3( FNAME, 'ALL', 1, 0, 0, NA ) ) THEN

            MESG = 'Could not read ' // CATEGORY // 'GRIDDING MATIX' //
     &             'from file "' // FNAME( 1:TRIMLEN( FNAME ) ) // '"'
            CALL M3EXIT( 'RDGRIDM', 0, 0, MESG, 2 )

        ENDIF

C.........  Check consistency of gridding matrix with dimensions

        IF( .NOT. DESC3( FNAME ) ) THEN
            MESG = 'Could not get description of file "' //
     &              FNAME( 1:TRIMLEN( FNAME ) ) // '"'

        ELSEIF( NROWS3D .NE. RDIM ) THEN
            WRITE( MESG, 94010 )
     &           'Row dimension mismatch. ' // CATEGORY //
     &           'GRIDDING MATRIX:', NROWS3D, 'program:', RDIM
            CALL M3EXIT( 'RDGRIDM', 0, 0, MESG, 2 )

        ELSEIF( NCOLS3D .NE. CDIM ) THEN
            WRITE( MESG, 94010 )
     &           'Column dimension mismatch. ' // CATEGORY //
     &           'GRIDDING MATRIX:', NCOLS3D, 'program:', CDIM
            CALL M3EXIT( 'RDGRIDM', 0, 0, MESG, 2 )

        ELSEIF( NTHIK3D .NE. TDIM ) THEN
            WRITE( MESG, 94010 )
     &           'Source number dimension mismatch. ' // CATEGORY //
     &           'GRIDDING MATRIX:', NTHIK3D, 'program:', TDIM
            CALL M3EXIT( 'RDGRIDM', 0, 0, MESG, 2 )

        ENDIF 

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

94010   FORMAT ( 10 ( A, :, I10, :, 2X ) )

        END

