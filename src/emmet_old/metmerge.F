C ****************************************************************************
C *
C * Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C *                System
C * File Version: @(#)$Id$
C *
C * COPYRIGHT (C) 1996, MCNC--North Carolina Supercomputing Center
C * All Rights Reserved
C *
C * See file COPYRIGHT for conditions of use.
C *
C * Environmental Programs Group
C * MCNC--North Carolina Supercomputing Center
C * P.O. Box 12889
C * Research Triangle Park, NC  27709-2889
C *
C * env_progs@mcnc.org
C *
C * Pathname: $Source$
C * Last updated: $Date$ 
C *
C
****************************************************************************/

        PROGRAM METMERGE
 
C***********************************************************************
C  program body starts at line 86
C
C  DESCRIPTION:
C       This program converts RGRND (Solar Radiation Reaching Ground) from
C       MET_CRO_2D to PAR (Photosynthetic Active Radiation), and merges it
C       with the Surface Temperature from a netCDF file (e.g. jul_kuo.cdf)
C       to generate a MET INPUT FILE for TMPBIO. Modifly NTMAX's value in
C       case the time step of the source netCDF file exceeds NTMAX.
C
C  REVISION  HISTORY:
C       Prototype  8/96 by C. SUNNY LO
C       Modified  12/96 by C. SUNNY LO
C       
C***********************************************************************
 
       IMPLICIT NONE
 
       INCLUDE 'GRDIMS3.EXT'   !  grid parameters
       INCLUDE 'PARMS3.EXT'    !  I/O API parameters
       INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
       INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.
 
C........  External functions
 
       LOGICAL         DSCGRID
       INTEGER         PROMPTFFILE
       CHARACTER*16    PROMPTMFILE
       REAL            GETREAL
       LOGICAL         GETYN
 
       EXTERNAL        DSCGRID, PROMPTFFILE,
     &                 PROMPTMFILE, GETREAL, GETYN
 
C........  Local variables
       INTEGER         LDEV      !  UNIT NUMBER FOR LOG FILE 
       CHARACTER*16    MC2NAME   !  logical name for MET_CRO_3D
       CHARACTER*16    SFCNAME   !  logical name for Surface Met File
       CHARACTER*16    ONAME     !  logical name for METMERGE_OUT
       CHARACTER*256   MESG      !  message buffer 
       CHARACTER*16    GNAME     !  scratch copy of grid name
       CHARACTER*16    CNAME     !  coordinate system name from DSCGRID()

       INTEGER         I, J, T, JDATE, JTIME, TSTEP

C........  Modify NTMAX if the source netCDF file has timestep > NTMAX
       INTEGER         NTMAX
       PARAMETER       ( NTMAX = 125 )

       REAL            WM2LY, LY2UE, SOL2PAR
       PARAMETER       ( WM2LY  = 0.001433,
     &                   LY2UE  = 2916.0
     &                 )

C........   Arrays for I/O API

       REAL            RGRND    ( NCOLS, NROWS, NTMAX)
       REAL            TEMP     ( NCOLS, NROWS, NTMAX)

C***********************************************************************
C   begin body of program METMERGE
     
       LDEV =INIT3()

C........  Process Input File (MCIP output) 

       MC2NAME = PROMPTMFILE(
     &        'Enter logical name for MET_CRO_2D file',
     &         FSREAD3, 'MC2NAME', 'METMERGE' )

       SFCNAME = PROMPTMFILE(
     &        'Enter logical name for Surface Met file',
     &         FSREAD3, 'SFCNAME', 'METMERGE' )
 
       IF ( .NOT. DESC3( SFCNAME ) ) THEN
           MESG = 'Could not get description of file  ' //
     &            SFCNAME
           CALL M3EXIT( 'METMERGE', 0, 0, MESG, 2 )
       END IF
 
       IF ( .NOT. DSCGRID( GRDNM, CNAME, GDTYP3D,
     &            P_ALP3D, P_BET3D, P_GAM3D, XCENT3D, YCENT3D,
     &            XORIG3D, YORIG3D, XCELL3D, YCELL3D,
     &            NCOLS3D, NROWS3D, NTHIK3D ) ) THEN
 
          GNAME = GRDNM
          MESG  = 'COULD NOT GET DESCRIPTION OF GRID'
          CALL M3EXIT( 'METCDF', 0, 0, MESG, 2 )
 
       END IF          !  IF DSCGRID() FAILED

       JDATE = SDATE3D
       JTIME = STIME3D
       TSTEP = 10000

       GDNAM3D = GRDNM
       FTYPE3D = GRDDED3
       NVARS3D = 2
       NLAYS3D = 1
       NTHIK3D = 1
       VGTYP3D = IMISS3
       VGTOP3D = AMISS3
       FDESC3D( 1 ) = 'MET input file for TMPBIO'

C........ FDESC3D( 2 ) and FDESC3D( 3 ) for scenario and scheme 
C........ will be imported from SFCNAME

       DO 30 I = 4, MXDESC3
          FDESC3D( I ) = ' '
30     END DO
 
       VNAME3D( 1 ) = 'PAR'
       VDESC3D( 1 ) = 'PAR'
       UNITS3D( 1 ) = 'uE/m2/sec'
       VTYPE3D( 1 ) = M3REAL
 
       VNAME3D( 2 ) = 'TA'
       VDESC3D( 2 ) = 'Temperature at 1.5 m'
       UNITS3D( 2 ) = 'K'
       VTYPE3D( 2 ) = M3REAL

       ONAME = PROMPTMFILE(
     &        'Enter logical name for MET output file',
     &         FSUNKN3, 'METMERGE_OUT', 'METMERGE' )

C........    Read Variables one timestep at a time

       DO 50 T = 1, MXREC3D 

          IF ( .NOT. READ3( MC2NAME, 'RGRND', ALLAYS3, JDATE, 
     &         JTIME, RGRND(1,1,T) ) ) THEN
               MESG = 'Could not read from file ' // MC2NAME
              CALL M3EXIT( 'METMERGE', 0, 0, MESG, 2 )
          END IF

          IF ( .NOT. READ3( SFCNAME, 'TEMP_00001', 1, JDATE, 
     &         JTIME, TEMP(1,1,T) ) ) THEN
               MESG = 'Could not read from file ' // SFCNAME
              CALL M3EXIT( 'METMERGE', 0, 0, MESG, 2 )
          END IF

          CALL NEXTIME( JDATE, JTIME, TSTEP )

50     END DO    ! T = 1, MXREC3D

C........    Calculation

       SOL2PAR = 0.5 * WM2LY * LY2UE

       DO 80 T = 1, MXREC3D
          DO 70 I = 1, NCOLS
             DO 60 J = 1, NROWS
                RGRND( I, J, T ) = SOL2PAR * RGRND( I, J, T )
60           END DO 
70        END DO 
80     END DO 

C........    Write the interpolated field into netCDF file

       JDATE = SDATE3D
       JTIME = STIME3D

       DO 100 T = 1, MXREC3D
          IF ( .NOT. WRITE3( ONAME, VNAME3D(1), 
     &         JDATE, JTIME, RGRND(1,1,T) ) ) THEN
             MESG = 'Could not write ' // VNAME3D(1)
             CALL M3EXIT( 'METCDF', 0, 0, MESG, 2 )
          END IF

          IF ( .NOT. WRITE3( ONAME, VNAME3D(2), 
     &         JDATE, JTIME, TEMP(1,1,T) ) ) THEN
             MESG = 'Could not write ' // VNAME3D(1)
             CALL M3EXIT( 'METCDF', 0, 0, MESG, 2 )
          END IF

          CALL NEXTIME( JDATE, JTIME, TSTEP )

100    END DO

C........    Finished processing one time step, call nextime

       CALL M3EXIT( 'METMERGE', 0, 0,
     &              'Normal completion of PROGRAM METMERGE', 0 )

       STOP
       END
