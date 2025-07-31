
        PROGRAM UAM2NCF

C***********************************************************************
C
C  DESCRIPTION: Takes UAM EMISSIONS gridded files and creates
C               netCDF gridded file which can be merge with other
C               gridded SMOKE netCDF files.
C
C  PRECONDITIONS REQUIRED:  
C               Tested only for lat-lon projection
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C             4/00 : Prototype by JMV
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
C*************************************************************************

      IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'PARMS3.EXT'      ! I/O API constants
        INCLUDE 'FDESC3.EXT'      ! I/O API file description data structure
        INCLUDE 'IODECL3.EXT'     ! I/O API function declarations
        INCLUDE 'EMCNST3.EXT'     !

C...........   EXTERNAL FUNCTIONS and their descriptions:

        LOGICAL         DSCM3GRD
        CHARACTER(10)   HHMMSS
        CHARACTER(14)   MMDDYY
        INTEGER         PROMPTFFILE 
        CHARACTER(16)   PROMPTMFILE
        INTEGER         TRIMLEN
        INTEGER         WKDAY
        
        EXTERNAL        DSCM3GRD, HHMMSS, MMDDYY, 
     &                  PROMPTFFILE, PROMPTMFILE, TRIMLEN, 
     &                  WKDAY

C...........   PARAMETER

        CHARACTER(50), PARAMETER :: 
     &  CVSW = '$Name SMOKEv5.2_Jul2025$' ! CVS release tag

C.......  Input and output emissions

        REAL, ALLOCATABLE ::  EMOB( :, : )         ! UAM input emissions 
        REAL, ALLOCATABLE ::  EMISG( : )        ! netCDF output emissions 

C...........   UAM header variables

        INTEGER  IFILE(10), NOTE(60)
        INTEGER  IX, IY, NXCLL, NYCLL, NX, NY, NZ
        INTEGER  NZLOWR, NZUPPR
        INTEGER  IBGDAT, IENDAT, IZONE
        INTEGER  NSEG, NSPECS, IDATE
        INTEGER         JDATE, JTIME        
        INTEGER, ALLOCATABLE ::  MSPEC ( :, : )  ! species names

        REAL     HTSUR, HTLOW, HTUPP
        REAL     ORGX, ORGY, UTMX, UTMY, DELTAX, DELTAY
        REAL     BEGTIM, ENDTIM
 
C...........   Other local variables

        INTEGER         ITIME
        INTEGER         BTIME                !  beginning time of day in run
        INTEGER         I, J, IJ, M, T, V    ! counters
        INTEGER         KDATE, LDATE, LTIME  ! previous date and time in loop
        INTEGER         LDEV                 ! logfile unit number
        INTEGER         NSTEPS               ! number of time steps
 
        INTEGER         ISEG                 ! segments
        INTEGER         IOS                  ! iostat error
        INTEGER         NGRID                ! ncols by nrows
        INTEGER         NCOLS, NROWS         ! number of cols and rows

        INTEGER         IDEV           !  unit number for emission input file
        LOGICAL         WFLAG, EFLAG   !  error flag

        CHARACTER(16)    ENAME   !  logical name for emission output file
        CHARACTER(16), ALLOCATABLE ::  KSPEC( : )  ! netCDF output species
        CHARACTER(16)            COORD    !  coordinate system name 
        CHARACTER(16)            COORUNIT !  coordinate system projection units
        CHARACTER(16)            GRDNM    !  grid name
        CHARACTER(80)            GDESC    !  grid description

        CHARACTER(256)  MESG    !  message buffer

        CHARACTER(16) :: PROGNAME = 'UAM2NCF'   !  program name

C***********************************************************************
C   begin body of program  
        
        LDEV = INIT3()

C.........  Write out copyright, version, web address, header info, and prompt
C           to continue running the program.

        CALL INITEM( LDEV, CVSW, PROGNAME )

C.......   Get UAM emissions file name 

        IDEV = PROMPTFFILE( 
     &          'Enter logical name for UAM EMISSIONS file',
     &         .TRUE., .FALSE., 'UAMEMIS', PROGNAME )

C.........  Get grid name from the environment and read grid parameters
        IF( .NOT. DSCM3GRD( GRDNM, GDESC, COORD, GDTYP3D, COORUNIT,
     &                          P_ALP3D, P_BET3D, P_GAM3D, XCENT3D,
     &                          YCENT3D, XORIG3D, YORIG3D, XCELL3D,
     &                          YCELL3D, NCOLS, NROWS, NTHIK3D ) ) THEN

                MESG = 'Could not get Models-3 grid description.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        WRITE( MESG,93000 )
     &              'Reading UAM emissions file header'
        CALL M3MESG( MESG )

C.......   Read UAM emissions file header

        READ(IDEV) IFILE,NOTE,NSEG,NSPECS,IDATE,BEGTIM,JDATE,ENDTIM
        READ(IDEV) ORGX,ORGY,IZONE,UTMX,UTMY,DELTAX,DELTAY,NX,NY,
     &             NZ,NZLOWR,NZUPPR,HTSUR,HTLOW,HTUPP
        READ(IDEV) IX,IY,NXCLL,NYCLL

        WRITE( MESG,93000 )
     &              'Checking grid parameters from header'
        CALL M3MESG( MESG )

C.......   Checking grid name parameters versus UAM header info

        IF ( NX .NE. NCOLS .OR. NY .NE. NROWS ) THEN
             WRITE( MESG,94010 ) 'Grid mismatch : ',
     &       'UAM File dimensions ', NX, NY, 
     &       'Grid name dimension ', NCOLS, NROWS 
             CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ENDIF

        IF ( XCELL3D .NE. DELTAX .OR. YCELL3D .NE. DELTAY ) THEN
             WRITE( MESG,94020 ) 'WARNING possible grid mismatch : ', 
     &       ' UAM File resolution(X,Y)  ', DELTAX, DELTAY,  
     &       ' Grid name resolution(X,Y) ', XCELL3D, YCELL3D 
             CALL M3MSG2( MESG )
        ENDIF

C.......  Allocate memory for emissions species variables

        ALLOCATE( MSPEC ( 10, NSPECS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MSPEC', PROGNAME )

        ALLOCATE( KSPEC ( NSPECS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'KSPEC', PROGNAME )

        WRITE( MESG,93000 )
     &              'Reading UAM-V species emitted from header'
        CALL M3MESG( MESG )

C.......  Read species from UAM-V header

        READ(IDEV) ((MSPEC(I,J),I=1,10),J=1,NSPECS)

        NSTEPS = ( JDATE - IDATE ) * 24 + NINT( ENDTIM - BEGTIM ) + 1


C.......  Write UAM file description to FDESC3

        WRITE( FDESC3D( 2 ),93015 )( NOTE( M ),M=1,60  )

        DO V = 1, NSPECS

          WRITE( KSPEC(V),93005 )( MSPEC( M,V ),M=1,10 )

        ENDDO


        IF ( IDATE .LT. 60000 ) THEN
          JDATE = IDATE + 2000000
        ELSE
          JDATE = IDATE + 1900000
        ENDIF
  
        JTIME = BEGTIM
        SDATE3D = JDATE
        STIME3D = NINT( BEGTIM ) 
        TSTEP3D = 10000
        NVARS3D = NSPECS
        FTYPE3D = GRDDED3   !  shares most of file-description with input file.
        NCOLS3D = NCOLS
        NROWS3D = NROWS
        NLAYS3D = 1
        GDNAM3D = GRDNM
        VGTYP3D = IMISS3
        VGTOP3D = BADVAL3

        DO  33 V = 1, NSPECS
            VNAME3D( V ) = KSPEC( V )
            UNITS3D( V ) = 'mole/hr'
            VDESC3D( V ) = 'UAM gridded speciated emissions'
            VTYPE3D( V ) = M3REAL
33      CONTINUE
        FDESC3D( 1 ) = 'UAM speciated gridded emissions.'

   
        WRITE( MESG,93000 )
     &              'Opening output file'
        CALL M3MESG( MESG )
 
C.........  Open output file
    
        ENAME = PROMPTMFILE( 
     &          'Enter logical name for OUTPUT EMIS file ' 
     &          , FSUNKN3, 'E2DNCF',  PROGNAME      )

C........  Allocate memory for emission arrays 

        ALLOCATE( EMOB( NX, NY  ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMOB', PROGNAME )

        NGRID = NCOLS3D * NROWS3D

        ALLOCATE( EMISG( NGRID  ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMISG', PROGNAME )

C.......   Write out emissions values

        BTIME = JTIME
        LDATE = 0
        DO  199  T = 1, NSTEPS

            READ( IDEV , END=999) IBGDAT, BEGTIM, IENDAT, ENDTIM 

            IF ( IBGDAT .LT. 60000 ) THEN
             KDATE = IBGDAT + 2000000
            ELSE
             KDATE = IBGDAT + 1900000
            ENDIF

C...............   If this is a new month, or new day, write message
            IF ( LDATE .NE. KDATE ) THEN
 
                MESG = 'Processing ' //
     &                 DAYS( WKDAY( KDATE ) ) // MMDDYY( KDATE )
                CALL M3MSG2( MESG( 1:TRIMLEN( MESG ) ) )

            END IF

            ITIME = NINT ( BEGTIM ) * TSTEP3D
  
C.............  Write to screen because WRITE3 only writes to LDEV
            WRITE( *, 93020 ) HHMMSS( ITIME )


            DO  122  V = 1,  NSPECS

C.............  Read UAM emissions

              READ ( IDEV ) ISEG, (MSPEC( M, V),M = 1,10), 
     &                      ((EMOB(I,J), I=1,NX),J=1,NY)


C.............  2d to 1d array

              DO J = 1, NY
                   DO I = 1, NX
                    IJ = ( J - 1 ) * NX + I
                    EMISG( IJ ) = EMOB( I, J )
                   ENDDO
              ENDDO


              IF ( .NOT. WRITE3( ENAME, KSPEC( V ), JDATE, JTIME,
     &                     EMISG ) ) THEN

                   CALL M3EXIT( PROGNAME , JDATE, JTIME,
     &                          'Could not write "' //
     &                           KSPEC( V ) //
     &                          '" to ' // ENAME, 2 )

              END IF         !  if write3() failed

 
122         CONTINUE    !  end loop on output vbles for this input vble


            LDATE = JDATE
            LTIME = JTIME

            CALL NEXTIME( JDATE, JTIME, TSTEP3D )


199     CONTINUE          !  end loop on time steps

999     CONTINUE          !  exit program
      

C.........   End of program:

        CALL M3EXIT( PROGNAME, 0, 0, ' ', 0 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Error and warning message formats..... 91xxx

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

93005   FORMAT( 10A1 )

93010   FORMAT( A16 )

93015   FORMAT( 60A1 )

93020   FORMAT( 8X, 'at time ', A8 )
 
C...........   Internal buffering formats............ 94xxx

94010   FORMAT( A, 10( A, 2I10)  )

94020   FORMAT( A, 10( A, 2F15.7)  ) 

        END PROGRAM  UAM2NCF
