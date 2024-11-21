
      PROGRAM GEN_AFDUST_TFRAC

C***********************************************************************
C
C  DESCRIPTION: Takes gridded SMOKE emissions file and applies
C               a user-supplied factor for each individual species
C               and applies it for a certain geographical region
C               obtained from input mask file..
C               The resulting emissions are output to a netCDF file.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C             10/00 : Prototype by JMV
C
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id: geofac.f,v 1.10 2007/07/11 19:18:32 bbaek Exp $
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
C Pathname: $Source: /nas01/depts/ie/cempd/apps/SMOKE_archive/smoke_archive/smoke/smoke/src/emutil/geofac.f,v $
C Last updated: $Date: 2007/07/11 19:18:32 $ 
C
C*************************************************************************

      USE MODFILESET
c     USE M3UTILIO

      IMPLICIT NONE

C...........   INCLUDES:

         INCLUDE 'IODECL3.EXT'     ! I/O API function declarations
         INCLUDE 'EMCNST3.EXT'     ! Emissions constants
         INCLUDE 'SETDECL.EXT'   !  FileSetAPI function declarations

C...........   EXTERNAL FUNCTIONS and their descriptions:

        INTEGER         GETFLINE
        INTEGER         PROMPTFFILE
        CHARACTER(16)   PROMPTMFILE
        INTEGER         TRIMLEN
        LOGICAL         ENVYN
        EXTERNAL  ENVYN, GETFLINE, PROMPTFFILE, PROMPTMFILE,
     &            TRIMLEN
           
C...........   PARAMETERS and their descriptions:

        CHARACTER(50), PARAMETER ::
     &  CVSW = '$Name SMOKEv5.1_July2024$'  ! CVS release tag

        INTEGER, PARAMETER :: NCLASS = 6

C...........  LOCAL VARIABLES

        REAL  TMPFRAC  
        REAL, ALLOCATABLE :: LUFRAC( :, : )       ! Landuse fraction 
        REAL, ALLOCATABLE :: TSUM( :, :, : )      ! species factors
        REAL, ALLOCATABLE :: CFRAC( : )           ! species factors
        REAL, ALLOCATABLE :: TFRAC( :, : )        ! fraction

        CHARACTER*16, ALLOCATABLE :: B4VAR( : )     ! Landuse fraction var
        CHARACTER*16, ALLOCATABLE :: CAPVAR( : )    ! species factors
        CHARACTER*16, ALLOCATABLE :: CCVAR( : )     ! species factors
        CHARACTER*16, ALLOCATABLE :: OUTNMVAR( : )  ! Output variable name

c       INTEGER  TSTEP                           ! time step
        INTEGER  I, J, K , L, M, N, Z            ! counters
        INTEGER  LDEV                            ! log file unit number
        INTEGER  RDEV, CDEV                      ! species factors unit number
        INTEGER  IFOUND, ICLASS
c       INTEGER  NSPECS                          ! number of species
        INTEGER  HR                              ! hour loop counter
        INTEGER  NLINES, MLINES                  ! number of species factors 
c       INTEGER  NSTEPS                          ! number of time steps
        INTEGER  IOS                             ! iostat
c       INTEGER  SDATE                           ! start date
c       INTEGER  STIME                           ! start time
        LOGICAL  LU_PCT_YN
        LOGICAL :: GET_LUFRAC = .FALSE.


c GRIDCRO2D
        INTEGER  FTYPE1, SDATE1, STIME1, TSTEP1, NTHIK1, NCOLS1, NROWS1
        INTEGER  NSTEPS1, NLAYS1, NVARS1, NSPECS1, GDTYP1
        REAL*8   P_ALP1, P_BET1, P_GAM1, XCENT1, YCENT1, XORIG1, YORIG1
        REAL*8   XCELL1, YCELL1
        INTEGER  VGTYP1
        REAL     VGTOP1
        INTEGER  mxrec1

c LUFRAC_CRO
        INTEGER  FTYPE2, SDATE2, STIME2, TSTEP2, NTHIK2, NCOLS2, NROWS2
        INTEGER  NSTEPS2, NLAYS2, NVARS2, NSPECS2, GDTYP2
        REAL*8   P_ALP2, P_BET2, P_GAM2, XCENT2, YCENT2, XORIG2, YORIG2
        REAL*8   XCELL2, YCELL2
        INTEGER  VGTYP2
        REAL     VGTOP2
        INTEGER  mxrec2
        INTEGER  NLU
        CHARACTER*2 chnfrac
        CHARACTER*16, ALLOCATABLE :: LUVAR(:)

c ----- 
        CHARACTER(16)  ENAME         ! logical name for GRIDCRO2D input file
        CHARACTER(16)  FNAME         ! logical name for LUFRAC_CRO input file
        CHARACTER(16)  ONAME         ! logical name for output file
        CHARACTER(300) MESG          ! message buffer for M3EXIT()

        CHARACTER(16) :: PROGNAME = 'GEN_AFDUST_TFRAC'   !  program name

C***********************************************************************
C   begin body of program

        LDEV = INIT3()

C.........  Write out copyright, version, web address, header info, and prompt
C           to continue running the program.

        CALL INITEM( LDEV, CVSW, PROGNAME )

C.........  Read Mapping of LUFRAC to landuse capture class mapping file
        RDEV = PROMPTFFILE(
     &           'Enter logical name for BELD4 to capture class file',
     &           .TRUE., .TRUE., 'BELD4TOCAPTURE', PROGNAME )

        NLINES = GETFLINE( RDEV, 'BELD4TOCAPTURE file' )
        WRITE( 6, * ) 'Number of lines BELD4TOCAPTURE = ',NLINES

        ALLOCATE( B4VAR( NLINES-1 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'B4VAR', PROGNAME )

        ALLOCATE( CAPVAR( NLINES-1 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CAPVAR', PROGNAME )

        READ( RDEV, * )

        DO N = 1, NLINES-1
          READ( RDEV, * )B4VAR( N ), CAPVAR( N ) 
          write( 6, * ) n, B4VAR( N ), CAPVAR( N )
        ENDDO

C.........  Read Fractions values of Capture Class file
        CDEV = PROMPTFFILE(
     &           'Enter logical name for capture class fractions file',
     &           .TRUE., .TRUE., 'CAPFRACS', PROGNAME )

        MLINES = GETFLINE( CDEV, 'CAPFRACS file' )
        WRITE( 6, * ) 'Number of lines CAPFRACS = ',MLINES

        ALLOCATE( CCVAR( MLINES-1 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CCVAR', PROGNAME )

        ALLOCATE( CFRAC( MLINES-1 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CFRAC', PROGNAME )

        READ( CDEV, * )

        DO N = 1, MLINES-1
          READ( CDEV, * )CCVAR( N ), CFRAC( N )
          write( 6, * ) n, ccvar( N ), cfrac( N )
        ENDDO

        ALLOCATE( OUTNMVAR( MLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'OUTNMVAR', PROGNAME )

C.........  Prompt for name of NetCDF input file

        ENAME = PROMPTSET(
     &       'Enter logical name for GRIDCRO2D file',
     &        FSREAD3, 'GRIDCRO2D', PROGNAME )

        IF ( .NOT. DESC3( ENAME ) ) THEN
            MESG = 'Could not get description of file "' //
     &             TRIM( ENAME ) // '"'
            CALL M3EXIT( PROGNAME , 0, 0, MESG, 2 )
        END IF

C.........  Assign local variables
        TSTEP1  = TSTEP3D
        SDATE1  = SDATE3D
        STIME1  = STIME3D
        NSPECS1 = NVARS3D
        NSTEPS1 = MXREC3D
        NCOLS1 = NCOLS3D
        NROWS1 = NROWS3D
        NLAYS1 = NLAYS3D
        NVARS1 = NVARS3D
        GDTYP1 = GDTYP3D
        P_ALP1 = P_ALP3D
        P_BET1 = P_BET3D
        P_GAM1 = P_GAM3D
        XCENT1 = XCENT3D
        YCENT1 = YCENT3D
        XORIG1 = XORIG3D
        YORIG1 = YORIG3D
        XCELL1 = XCELL3D
        YCELL1 = YCELL3D
        VGTYP1 = VGTYP3D
        VGTOP1 = VGTOP3D

C.........  Allocate memory for Landuse fraction

        ALLOCATE( LUFRAC( NCOLS3D, NROWS3D ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LUFRC', PROGNAME )

        ALLOCATE( TFRAC( NCOLS3D, NROWS3D ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TFRAC', PROGNAME )

        ALLOCATE( TSUM( NCOLS3D, NROWS3D, MLINES-1 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TSUM', PROGNAME )

        TSUM  = 0.0000 ! array
        TFRAC = 0.0000 ! array

C........ Whether Landuse fractions are percentage
        LU_PCT_YN  = ENVYN ( 'LU_PCT_YN', MESG, .TRUE., IOS )

C........ Prompt for seperate LUFRAC files (output from MCIP post v4.5)

        CALL   ENVSTR(  'LUFRACRO', 'Get LUFRAC_CRO path', '', FNAME, IOS )

        IF ( IOS .NE. 0 ) THEN
            MESG = 'LUFRAC_CRO input file not defined. Assuming LUFRAC variables in ' // TRIM( ENAME )
            CALL M3MESG(MESG)
        ELSE
            FNAME = PROMPTSET(
     &       'Enter logical name for LUFRAC_CRO file',
     &        FSREAD3, 'LUFRACRO', PROGNAME )

            IF ( .NOT. DESC3( FNAME ) ) THEN
               MESG = 'Could not get description of file "' //
     &             TRIM( FNAME ) // '"'
               CALL M3EXIT( PROGNAME , 0, 0, MESG, 2 )
            END IF

C.........  Assign local variables

            TSTEP2  = TSTEP3D
            SDATE2  = SDATE3D
            STIME2  = STIME3D
c           NSPECS2 = NVARS3D
            NSTEPS2 = MXREC3D
            NCOLS2 = NCOLS3D
            NROWS2 = NROWS3D
            NLAYS2 = NLAYS3D
c           NVARS2 = NVARS3D
            GDTYP2 = GDTYP3D
            P_ALP2 = P_ALP3D
            P_BET2 = P_BET3D
            P_GAM2 = P_GAM3D
            XCENT2 = XCENT3D
            YCENT2 = YCENT3D
            XORIG2 = XORIG3D
            YORIG2 = YORIG3D
            XCELL2 = XCELL3D
            YCELL2 = YCELL3D
c           VGTYP2 = VGTYP3D
c           VGTOP2 = VGTOP3D

            if ( NCOLS1 .ne. NCOLS2 ) stop 'NCOLS3D inconsistent'
            if ( NROWS1 .ne. NROWS2 ) stop 'NROWS3D inconsistent'
            if ( GDTYP1 .ne. GDTYP2 ) stop 'GDTYP3D inconsistent'
            if ( P_ALP1 .ne. P_ALP2 ) stop 'P_ALP3D inconsistent'
            if ( P_BET1 .ne. P_BET2 ) stop 'P_BET3D inconsistent'
            if ( P_GAM1 .ne. P_GAM2 ) stop 'P_GAM3D inconsistent'
            if ( XCENT1 .ne. XCENT2 ) stop 'XCENT3D inconsistent'
            if ( YCENT1 .ne. YCENT2 ) stop 'YCENT3D inconsistent'
            if ( XORIG1 .ne. XORIG2 ) stop 'XORIG3D inconsistent'
            if ( YORIG1 .ne. YORIG2 ) stop 'YORIG3D inconsistent'
            if ( XCELL1 .ne. XCELL2 ) stop 'XCELL3D inconsistent'
            if ( YCELL1 .ne. YCELL2 ) stop 'YCELL3D inconsistent'

            NLU = NLAYS2

            ALLOCATE( LUVAR( NLU ), STAT=IOS )
            CALL CHECKMEM( IOS, 'LUVAR', PROGNAME )

            GET_LUFRAC = .TRUE.

C........ Read LUFRAC_CRO file
c           CALL READ_LUFRAC (NCOLS3D, NROWS3D, NLINES, B4VAR, CAPVAR, 
c    &           MLINES, CCVAR, CFRAC, LU_PCT_YN, TSUM, TFRAC)

        END IF
              
C.............  Write to screen because WRITE3 only writes to LDEV
        IF ( GET_LUFRAC ) THEN

           DO L = 1, NLU

              WRITE(chnfrac,'(i2)')L                       ! should range from 0-99
              IF ( L .lt. 10 ) chnfrac = '0'//chnfrac(2:2)
              LUVAR(L) = 'LUFRAC_'//chnfrac
c             print*, LUVAR(L)

              LUFRAC = 0.0   !  array

              IF ( .not. read3(FNAME, 'LUFRAC',L,sdate2,stime2,LUFRAC)) THEN
                 MESG = 'Could not read variable LUFRAC' //
     &                  ' from file ' // TRIM( FNAME )
                 CALL M3EXIT( PROGNAME, SDATE2, STIME2, MESG, 2 )
              ENDIF

              IFOUND = -9
              DO N = 1, NLINES-1
                IF ( LUVAR( L ) .EQ. B4VAR( N ) ) THEN
                  write( 6, * ) LUVAR( L )," CLASS FOUND ",CAPVAR( N )
                  IFOUND = L
                  DO M = 1, MLINES-1
                     IF ( CAPVAR( N ) .EQ. CCVAR( M ) ) THEN
                        ICLASS = M
                     ENDIF
                  ENDDO
                ENDIF
              ENDDO

              IF ( IFOUND .LT. 0 ) THEN
                 write( 67, * ) LUVAR( L )," CLASS NOT FOUND " 
              ELSE
                 DO I = 1, NCOLS2
                 DO J = 1, NROWS2
                    TSUM( I, J, ICLASS ) = TSUM( I,J, ICLASS ) + LUFRAC( I, J )
                    IF ( LU_PCT_YN ) THEN
                       TMPFRAC = LUFRAC( I, J ) * 0.01 * CFRAC( ICLASS )
                    ELSE
                       TMPFRAC = LUFRAC( I, J ) * CFRAC( ICLASS )
                    ENDIF
                    TFRAC( I, J ) = TFRAC( I, J ) + TMPFRAC
                 ENDDO
                 ENDDO
              ENDIF

           ENDDO

        ELSE
               
           DO  L = 1, NSPECS1

              LUFRAC = 0.0   !  array

C.....................  Read input file for time and species of interest

              IF( .NOT. READ3( ENAME, VNAME3D( L ), 1,
     &                       SDATE1, STIME1, LUFRAC ) ) THEN
                 MESG = 'Could not read variable ' //
     &                 TRIM( VNAME3D( L ) ) // ' from file ' //
     &                 TRIM( ENAME )
                 CALL M3EXIT( PROGNAME, SDATE1, STIME1, MESG, 2 )
              END IF

              IFOUND = -9
              DO N = 1, NLINES-1

                 IF ( VNAME3D( L ) .EQ. B4VAR( N ) ) THEN
                    write( 6, * ) VNAME3D( L )," CLASS FOUND ",CAPVAR( N )
                    IFOUND = L
                    DO M = 1, MLINES-1
                       IF ( CAPVAR( N ) .EQ. CCVAR( M ) ) THEN
                           ICLASS = M
                       ENDIF
                    ENDDO
                 ENDIF

              ENDDO

              IF ( IFOUND .LT. 0 ) THEN
                 write( 67, * ) VNAME3D( L )," ACK CLASS NOT FOUND " 
              ELSE

                 DO I = 1, NCOLS3D
                 DO J = 1, NROWS3D
                    TSUM( I, J, ICLASS ) = TSUM( I,J, ICLASS ) + LUFRAC( I, J )
                    IF ( LU_PCT_YN ) THEN
                       TMPFRAC = LUFRAC( I, J ) * 0.01 * CFRAC( ICLASS )
                    ELSE
                       TMPFRAC = LUFRAC( I, J ) * CFRAC( ICLASS )
                    ENDIF
                    TFRAC( I, J ) = TFRAC( I, J ) + TMPFRAC
                 ENDDO
                 ENDDO

              ENDIF

           ENDDO

        ENDIF


c       write(*,*)TFRAC(10,10)

C.........  Open output file

c set up the output file by resetting some parameters of the ingrid file to outgrid specs.
        IF ( .NOT. DESC3( ENAME ) ) THEN
            MESG = 'Could not get description of file "' //
     &             TRIM( ENAME ) // '"'
            CALL M3EXIT( PROGNAME , 0, 0, MESG, 2 )
        END IF
       
        NVARS3D = MLINES
        DO M = 1, MLINES - 1
          VNAME3D ( M ) = CCVAR( M )
          VDESC3D ( M ) = "fractional land use"
          VTYPE3D ( M ) = 5
          UNITS3D ( MLINES ) = ''
          IF (LU_PCT_YN) UNITS3D ( MLINES ) = 'percent'
        ENDDO

        VNAME3D ( MLINES ) = 'xportfrac'
        UNITS3D ( MLINES ) = ''
        VDESC3D ( MLINES ) = 'capture fraction'
        VTYPE3D ( MLINES ) = 5

        MESG  = 'Output file'
        ONAME = 'OUTFILE'
        ONAME = PROMPTMFILE( MESG, FSUNKN3, ONAME, PROGNAME )

        DO M = 1, MLINES - 1

           IF( .NOT. WRITE3( ONAME, VNAME3D( M ),
     &                    SDATE1, STIME1, TSUM(:,:,M)) ) THEN

              MESG = 'Could not write variable ' //
     &               TRIM( VNAME3D( M ) ) // ' to file ' //
     &               TRIM( ONAME )
             CALL M3EXIT( PROGNAME, SDATE1, STIME1, MESG, 2 )

           END IF

        ENDDO


        IF( .NOT. WRITE3( ONAME, VNAME3D( MLINES ),
     &                    SDATE1, STIME1, TFRAC ) ) THEN

             MESG = 'Could not write variable ' //
     &               TRIM( VNAME3D( MLINES ) ) // ' to file ' //
     &               TRIM( ONAME )
             CALL M3EXIT( PROGNAME, SDATE1, STIME1, MESG, 2 )

        END IF


C.........   End of program:

        CALL M3EXIT( PROGNAME, 0, 0, ' ', 0 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Error and warning message formats..... 91xxx

C...........   Informational (LOG) message formats... 92xxx

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A, F10.5 ) 

C...........   Internal buffering formats............ 94xxx

      END PROGRAM  GEN_AFDUST_TFRAC
