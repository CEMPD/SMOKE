
        SUBROUTINE OPENMRGOUT

C***********************************************************************
C  subroutine OPENMRGOUT body starts at line
C
C  DESCRIPTION:
C      The purpose of this subroutine is to open all of the necessary
C      output files for the merge routine (both I/O API and ASCII files)
C
C  PRECONDITIONS REQUIRED:  
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C       Created 2/99 by M. Houyoux
C
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1998, MCNC--North Carolina Supercomputing Center
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
C****************************************************************************

C.........  MODULES for public variables
C.........  This module contains the major data structure and control flags
        USE MODMERGE

        IMPLICIT NONE

C.........  INCLUDES:
        
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file desc. data structures

C.........  EXTERNAL FUNCTIONS and their descriptions:
        
        CHARACTER*2     CRLF
        INTEGER         ENVINT  
c        INTEGER         ENVYN  
        INTEGER         GETDATE  
        INTEGER         GETNUM  
        CHARACTER*14    MMDDYY
        INTEGER         PROMPTFFILE  
        CHARACTER*16    PROMPTMFILE  

        EXTERNAL  CRLF, ENVINT, ENVYN, GETDATE, GETNUM, PROMPTFFILE, 
     &            PROMPTMFILE

C.........  Other local variables

        INTEGER         L, L1, L2, V

        CHARACTER*50    UNIT         ! tmp units buffer
        CHARACTER*300   BUFFER       ! tmp buffer
        CHARACTER*300   MESG         ! message buffer

        CHARACTER*16 :: PROGNAME = 'OPENMRGOUT' ! program name

C***********************************************************************
C   begin body of subroutine OPENMRGOUT

C.........  Set default output file names
        CALL MRGONAMS

C.........  Set up and open I/O API output file
        IF( LGRDOUT ) THEN

            FTYPE3D = GRDDED3
            P_ALP3D = DBLE( P_ALP )
            P_BET3D = DBLE( P_BET )
            P_GAM3D = DBLE( P_GAM )
            XCENT3D = DBLE( XCENT )
            YCENT3D = DBLE( YCENT )
            XORIG3D = DBLE( XORIG )
            YORIG3D = DBLE( YORIG )
            XCELL3D = DBLE( XCELL )
            YCELL3D = DBLE( YCELL )
            SDATE3D = SDATE
            STIME3D = STIME
            TSTEP3D = TSTEP
            NCOLS3D = NCOLS
            NROWS3D = NROWS
            NLAYS3D = EMLAYS
            NTHIK3D = 1
            GDTYP3D = GDTYP
            VGTYP3D = VGTYP
            VGTOP3D = VGTOP
            GDNAM3D = GRDNM

! NOTE: Must insert FDDESC, including met header information, version of
! programs for all input files

C.............  Set constants number and values for variables
            UNIT = INVUNIT
            IF( SFLAG ) THEN
                NVARS3D = NMSPC

                L  = INDEX( INVUNIT, '/' )
                L1 = LEN_TRIM( INVUNIT )
                L2 = INDEX( SPCUNIT, '/' )
                UNIT = SPCUNIT( 1:L2 ) // INVUNIT( L+1:L1 )

                DO V = 1, NVARS3D
                    VNAME3D( V ) = EMNAM( V )
                    UNITS3D( V ) = UNIT
                    VDESC3D( V ) = 'Species ' // EMNAM( V )
                    VTYPE3D( V ) = M3REAL
                END DO

            ELSEIF( VFLAG ) THEN
                NVARS3D = 1

                VNAME3D( 1 ) = 'VMT'
                UNITS3D( 1 ) = INVUNIT
                VDESC3D( 1 ) = 'Vehicle Miles Traveled, VMT'
                VTYPE3D( 1 ) = M3REAL

            ELSE
                NVARS3D = NIPOL

                DO V = 1, NVARS3D
                    VNAME3D( V ) = EINAM( V )
                    UNITS3D( V ) = INVUNIT
                    VDESC3D( V ) = 'Pollutant ' // EINAM( V )
                    VTYPE3D( V ) = M3REAL
                END DO

            ENDIF
          
C.............  Prompt for and open file
            IF( AFLAG ) THEN
                AONAME = PROMPTMFILE(  
     &            'Enter name for AREA-SOURCE GRIDDED OUTPUT file',
     &            FSUNKN3, AONAME, PROGNAME )
            END IF 

            IF( MFLAG ) THEN
                MONAME = PROMPTMFILE(  
     &            'Enter name for MOBILE-SOURCE GRIDDED OUTPUT file',
     &            FSUNKN3, MONAME, PROGNAME )
            END IF 

            IF( PFLAG ) THEN
                PONAME = PROMPTMFILE(  
     &            'Enter name for POINT-SOURCE GRIDDED OUTPUT file',
     &            FSUNKN3, PONAME, PROGNAME )
            END IF 

            IF( XFLAG ) THEN
                TONAME = PROMPTMFILE(  
     &            'Enter name for MULTI-SOURCE GRIDDED OUTPUT file',
     &            FSUNKN3, TONAME, PROGNAME )
            END IF 

        END IF  ! End of gridded output

C.........  Open report file(s) and output header(s)

        IF( LREPSTA .OR. LREPCNY ) THEN

! NOTE: Fill this in later...

        END IF  ! End of state and/or county output

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats.............94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94050   FORMAT( A, 1X, I2.2, A, 1X, A, 1X, I6.6, 1X,
     &          A, 1X, I3.3, 1X, A, 1X, I3.3, 1X, A   )

        END SUBROUTINE OPENMRGOUT
