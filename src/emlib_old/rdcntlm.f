
        SUBROUTINE RDCNTLM( CATEGORY, FNAME, MXSRC, NSRCS, 
     &                      MXIPOL, NINVP, CTLFACS, CINV )

C***********************************************************************
C  subroutine body starts at line 75
C
C  DESCRIPTION:
C       This subroutine reads the control matrix for area, mobile, or point
C       sources, and compares the dimensions of the matrix with the 
C       compiled dimensions.
C
C  PRECONDITIONS REQUIRED:
C       
C       File name FNAME must be defined and opened.
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       TRIMLEN, DESC3
C
C  REVISION  HISTORY:
C       Started 7/98 by M Houyoux
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
        INTEGER         MXSRC        ! maximum number of sources (in)
        INTEGER         NSRCS        ! actual number of sources (out)
        INTEGER         MXIPOL       ! maximum number of pollutants (in)
        INTEGER         NINVP        ! actual number of pollutants (out)
        REAL            CTLFACS( MXSRC, MXIPOL ) ! spec matrix
        CHARACTER*16    CINV( MXIPOL )! names of spc matrix variables

C.........  Other local variables

        INTEGER         V
        LOGICAL         EFLAG
        CHARACTER*256   MESG
 
C***********************************************************************
C   begin body of program RDCNTLM

C.........  Get description of speciation matrix
 
        EFLAG = .FALSE.
        IF( .NOT. DESC3( FNAME ) ) THEN
            MESG = 'Could not get description of file "' //
     &              FNAME( 1:TRIMLEN( FNAME ) ) // '"'
            CALL M3EXIT( 'RDCNTLM', 0, 0, MESG, 2 )

        ELSEIF( NROWS3D .GT. MXSRC ) THEN
            EFLAG = .TRUE.
            WRITE( MESG, 94010 )
     &           'Source dimension mismatch. ' // CATEGORY //
     &           'CONTROL MATRIX:', NROWS3D, 'program:', MXSRC
            CALL M3MSG2( MESG )
        ENDIF

        IF( NVARS3D .GT. MXIPOL ) THEN
            EFLAG = .TRUE.
            WRITE( MESG, 94010 )
     &           'Variables dimension mismatch. ' // CATEGORY //
     &           'CONTROL MATRIX:', NVARS3D, 'program:', MXIPOL            
            CALL M3MSG2( MESG )
        ENDIF

        IF( EFLAG ) THEN
            MESG = 'Program - file dimension mismatch'
            CALL M3EXIT( 'RDCNTLM', 0, 0, MESG, 2 )

        ELSE
            NSRCS = NROWS3D
            NINVP = NVARS3D

        ENDIF

C.........  Read speciation matrix for each species
        DO 101 V = 1, NINVP
             CINV( V ) =  VNAME3D( V )(9:16)
             IF( .NOT. READ3( FNAME, VNAME3D( V ), 1, 0, 0, 
     &                        CTLFACS( 1,V )               ) ) THEN

                MESG= 'Could not read ' //CATEGORY// ' CONTROL MATRIX'//
     &                'from file "' // FNAME( 1:TRIMLEN(FNAME) ) // '"'
                CALL M3EXIT( 'RDCNTLM', 0, 0, MESG, 2 )

            ENDIF
101     CONTINUE

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

94010   FORMAT ( 10 ( A, :, I10, :, 2X ) )

        END

