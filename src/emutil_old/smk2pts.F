C
      PROGRAM SMK2PTS

C***********************************************************************
C  program body starts at line 149
C
C  DESCRIPTION:
C       Uses smoke generated NetCDF ETSC file of elevated point sources
C       and UAM's USERIN file to generate the input file to the UAM PTSRCE 
C       elevated point sources pre-processor.
C
C  PRECONDITIONS REQUIRED:
C       Code requires the following input information
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       Models-3 I/O
C       GETSTR, SORTI3, NEXTIME
C       Specific      : RDEPSD, CHKTMP 
C       EPS functions : RDREGN
C
C  REVISION  HISTORY:
C       Prototype  9/95 by MRH
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

        IMPLICIT NONE

C...........   INCLUDES:
        INCLUDE 'PTDIMS3.EXT'
        INCLUDE 'PARMS3.EXT'
        INCLUDE 'IODECL3.EXT'
        INCLUDE 'FDESC3.EXT'

C...........   UAM includes (in /200p_scratch/uamguides/eps2.0b/src/include)
        INCLUDE 'param.inc'
        INCLUDE 'uamcntl.inc'
        INCLUDE 'uamregn.inc'
        INCLUDE 'mrguam.inc' 

C...........   EXTERNAL FUNCTIONS and their descriptions:

        INTEGER       GETNUM
        LOGICAL       GETYN
        INTEGER       PROMPTFFILE
        CHARACTER*16  PROMPTMFILE
        INTEGER       TRIMLEN
   
        EXTERNAL      GETYN, GETNUM, PROMPTFFILE, PROMPTMFILE, TRIMLEN

C...........   LOCAL FUNCTIONS and their descriptions:

        LOGICAL       CHKVAR  ! Checks for the output vars in the input var list

C...........   LOCAL PARAMETERS and their descriptions
        INTEGER       ISL              !  Length of string for reading CONTROL
        INTEGER       MXS              !  Number of lines in CONTROL
        PARAMETER   ( ISL=70, MXS=5 )

C...........   LOCAL VARIABLES and their descriptions:
        INTEGER   IFIP  ( NPSRC )  ! Point source FIP values
        INTEGER   IPLANT( NPSRC )  ! Point source plant numbers
        INTEGER   ISTACK( NPSRC )  ! Point source stack numbers
        REAL      STKHT ( NPSRC )  ! Point source stack heights [m]
        REAL      STKDM ( NPSRC )  ! Point source stack diameters [m]
        REAL      STKTK ( NPSRC )  ! Point source stack exit temperatures [K]
        REAL      STKVE ( NPSRC )  ! Point source stack exit velocities [m/hr]
        REAL      XLOCA ( NPSRC )  ! Point source x-location [km]
        REAL      YLOCA ( NPSRC )  ! Point source y-location [km]

        INTEGER   NELVIN                  ! Current source count
        INTEGER   INDXIN( NPELV )         ! Subscript into NetCDF inventory file
        REAL      EMIS  ( NPELV, MXSPEC ) ! Emissions values ( gm-mol/hr )

        COMMON / ELEVIN / NELVIN, INDXIN, EMIS

        REAL      ETEMP ( NPELV )  ! Temporary emissions array
        INTEGER   EFIP  ( NPELV )  ! Elevated source FIP 
        INTEGER   EPLANT( NPELV )  ! Elevated source plant numbers
        INTEGER   ESTACK( NPELV )  ! Elevated source stack numbers
        INTEGER   ENDX  ( NPELV )  ! Elevated source index for sorting
        INTEGER   IEND  ( NPELV )  ! Pointer to end of sum list for same sources
    
        INTEGER   IBD           ! Start date for each output time period (YYDDD)
        INTEGER   IBT           ! Beginning time for a time step period
        INTEGER   IED           ! End date for output time period (YYDDD)
        INTEGER   IET           ! Ending time for a time step period
        INTEGER   IERR          ! Error flag, 0 = no error
        INTEGER   ISPTR(MXSPEC) ! location of output vars in input vars list
        INTEGER   JDATE         ! Date variable (YYYYDDD)
        INTEGER   JDATEP1       ! Date variable plus one time step (YYYYDDD)
        INTEGER   JTIME         ! Time variable (HHMMSS)
        INTEGER   JTIMEP1       ! Time variable plus one time step (HHMMSS)
        INTEGER   NELPT         ! Actual count of elevated sources
        INTEGER   NPOUT         ! Count of number of output sources
        INTEGER   NSTEPS        ! User specified number of output time steps
        INTEGER   ODEV          ! UAM file desigator for output file
        INTEGER   TSTEP         ! Time step (HHMMSS)
        INTEGER   UDEV          ! Userin file desigator for applicable UAM run

        CHARACTER*10  VTYPE     ! User specified vertical method type

C...........   Other local variables :
        INTEGER         BDATE, BTIME, EDDATE, EDTIME
        INTEGER         I, ICNT, IDX, IOS, ISPC
        INTEGER         J, JDX, LDEV, M, NELVPREV, T
        LOGICAL         IGO

        REAL            ESUM, VH, XMAX, YMAX
        CHARACTER*10    BSPEC
        CHARACTER*10    DEFAULT
        CHARACTER*(ISL) ESTRING( MXS ), BUFFER
        CHARACTER*16    ENAME   !  logical name for point-source input file
        CHARACTER*16    UNAME   !  " (upper-level) output file, or "NONE"
        CHARACTER*256   MESG    !  buffer for M3EXIT() messages
        
C.... Packets
        INTEGER      MXPKTS
        PARAMETER  ( MXPKTS = 7 )
        CHARACTER*60 PNAME ( MXPKTS )
        DATA         PNAME 
     &             / 'CONTROL', 'REGION', 'POINT SOURCES', 
     &               'TIME INTERVAL', 'METHOD', 'VERTICAL METHOD', 
     &               'EMISSIONS VALUES' /

C***********************************************************************
C   begin body of program SMK2PTS

        LDEV = INIT3()

        CALL INITEM( LDEV )
    
C.........  Get file names; open input files

        ENAME = PROMPTMFILE(
     &          'Enter name for POINT SOURCE input file',
     &          FSREAD3, 'PNTS', 'SMK2PTS' )

        UNAME = PROMPTMFILE(
     &          'Enter name for ELEVATED EMISSIONS input file',
     &          FSREAD3, 'ETSC', 'SMK2PTS' )
      
        UDEV = PROMPTFFILE(
     &        'Enter logical name for USERIN file applying to UAM run',
     &        .TRUE., .TRUE., 'USERIN', 'SMK2PTS' )

        NSTEPS = GETNUM( 24, 999999, NSTEPS,
     &                   'Enter output duration (hours)' )

C.........  Get type of plume rise calculation

        WRITE( MESG,94010 ) 
     &         'Enter vertical method: "PLUMERISE" or "STACKHGT"'
        DEFAULT = 'PLUMERISE'
        CALL GETSTR( MESG, DEFAULT, VTYPE )
        IF( VTYPE .EQ. 'PLUMERISE' ) VTYPE = 'PLUME RISE'

C.........  Read UNAME header information

        IF( .NOT. DESC3( UNAME ) ) THEN
        CALL M3ERR( 'SMK2PTS', 0, 0,
     &              'Could not get description of file ' // UNAME,
     &              .TRUE. )
        END IF

C.........  Read the EPISODE and UAMREGN packets from the USERIN file

        CALL RDEPSD( IERR, UDEV, ESTRING )
        IF ( IERR .NE. 0 )
     &     CALL M3ERR( 'SMK2PTS', 0, 0,
     &                 'Error reading EPISODE packet' ,
     &                 .TRUE. )
 
C.........  Read the UAMREGN packet from the USERIN file

        CALL RDREGN( IERR, UDEV )                    ! Note: EPS function
        IF (IERR .NE. 0)
     &       CALL M3ERR( 'SMK2PTS', 0, 0,
     &                   'Error reading UAMREGN packet' ,
     &                   .TRUE. )


C.........  Compare the NetCDF header information with the USERIN information

C.... a) Check variable names and set up pointer for correct order

        DO 11 ISPC = 1, NSPECS
           
           BUFFER = SPCNAM( ISPC )
           M = TRIMLEN( BUFFER )
           BSPEC = BUFFER( 1:M )
           IF( .NOT. 
     &         CHKVAR( BSPEC, VNAME3D, NVARS3D, ISPTR( ISPC ) ) ) THEN
        
              WRITE( MESG,94010 ) BSPEC // ' not found ' //
     &                            'in input.  Zero and continue?'
              IF( .NOT. GETYN( MESG, .TRUE. ) ) THEN
                  CALL M3EXIT( 'SMK2PTS', 0, 0, 'Ending program', 2 )
              ENDIF

           ENDIF

   11   CONTINUE

C.... b) Check temporal parameters and initialize JTIME, JDATE, and TSTEP

        CALL CHKTMP( IERR, JDATE, EDDATE, JTIME, EDTIME, TSTEP )
        IF ( IERR .NE. 0) 
     &  CALL M3ERR( 'SMK2PTS', 0, 0,
     &              'Fatal difference between input and output ' //
     &              'time file specifications', .TRUE. )

C.........  Read time independent point source data from PNTS file

        IF ( .NOT. READ3( ENAME, 'XLOCA', ALLAYS3,0,0, XLOCA ) ) THEN
            CALL M3EXIT( 'SMK2PTS', 0, 0,
     &                   'Error reading XLOCA from file' // ENAME, 2)
        END IF
 
        IF ( .NOT. READ3( ENAME, 'YLOCA', ALLAYS3,0,0, YLOCA ) ) THEN
            CALL M3EXIT( 'SMK2PTS', 0, 0,
     &                   'Error reading YLOCA from file' // ENAME, 2)
        END IF
 
        IF ( .NOT. READ3( ENAME, 'IPLANT', ALLAYS3,0,0, IPLANT ) ) THEN
            CALL M3EXIT( 'SMK2PTS', 0, 0,
     &                   'Error reading IPLANT from file' // ENAME, 2 )
        END IF
 
        IF ( .NOT. READ3( ENAME, 'ISTACK', ALLAYS3,0,0, ISTACK ) ) THEN
            CALL M3EXIT( 'SMK2PTS', 0, 0,
     &                   'Error reading ISTACK from file' // ENAME, 2 )
        END IF

        IF ( .NOT. READ3( ENAME, 'IFIP', ALLAYS3,0,0, IFIP ) ) THEN
            CALL M3EXIT( 'SMK2PTS', 0, 0,
     &                   'Error reading FIP from file' // ENAME, 2)
        END IF
 
        IF ( .NOT. READ3( ENAME, 'STKHT', ALLAYS3,0,0, STKHT ) ) THEN
            CALL M3EXIT( 'SMK2PTS', 0, 0,
     &                   'Error reading STKHT from file' // ENAME, 2)
        END IF
 
        IF ( .NOT. READ3( ENAME, 'STKDM', ALLAYS3,0,0, STKDM ) ) THEN
            CALL M3EXIT( 'SMK2PTS', 0, 0,
     &                   'Error reading STKDM from file' // ENAME, 2)
        END IF
 
        IF ( .NOT. READ3( ENAME, 'STKTK', ALLAYS3,0,0, STKTK ) ) THEN
            CALL M3EXIT( 'SMK2PTS', 0, 0,
     &                   'Error reading STKTK from file' // ENAME, 2)
        END IF
 
        IF ( .NOT. READ3( ENAME, 'STKVE', ALLAYS3,0,0, STKVE ) ) THEN
            CALL M3EXIT( 'SMK2PTS', 0, 0,
     &                   'Error reading STKVE from file' // ENAME, 2)
        END IF

C.........  Read the first time step of the emissions file to get a copy 
C   ......  of the index.  Must assume all time steps have same number of
C   ......  sources

        IF( .NOT. READ3( UNAME, ALLVAR3, ALLAYS3, 
     &                   JDATE, JTIME, NELVIN ) ) THEN
            CALL M3EXIT( 'SMK2PTS', JDATE, JTIME,
     &                   'Could not read time step from ' // UNAME, 2 )
        ENDIF
        NELPT = NELVIN

C.........  Copy IPLANT and ISTACK to elevated-only arrays
        
        DO 16 I = 1, NELPT
        
            EFIP  ( I ) = IFIP  ( INDXIN( I ) )
            EPLANT( I ) = IPLANT( INDXIN( I ) )
            ESTACK( I ) = ISTACK( INDXIN( I ) )
            ENDX  ( I ) = I

   16   CONTINUE

C.........  Sort elevated pointer by EFIP, EPLANT, ESTACK

        CALL SORTI3( NELPT, ENDX, EFIP, EPLANT, ESTACK )

C.........  Save beginning date and time values

        BDATE = JDATE
        BTIME = JTIME

C.........  Initiallize ETEMP and IEND
        DO 17 I = 1, NELPT
            ETEMP( I ) = 0.
            IEND ( I ) = 1
   17   CONTINUE

C.........  Decide which sources have at least one non-zero emission

        DO 20 T = 1, NSTEPS

            IF( .NOT. READ3( UNAME, ALLVAR3, ALLAYS3, 
     &                       JDATE, JTIME, NELVIN ) ) THEN
                CALL M3EXIT( 'SMK2PTS', JDATE, JTIME,
     &                       'Could not read time step from ' 
     &                        // UNAME, 2 )
            ENDIF

            DO 19 I = 1, NELPT

                DO 18 ISPC = 1, MXSPEC
                    ETEMP( I ) = ETEMP( I ) + EMIS( ENDX( I ), ISPC )
   18           CONTINUE

   19       CONTINUE

            CALL NEXTIME( JDATE, JTIME, TSTEP )

   20   CONTINUE

        DO 21 I = 1, NELPT
            IF( ETEMP( I ) .EQ. 0. ) THEN
                ENDX( I ) = 0
            ENDIF               
   21   CONTINUE

C.........  Eliminate sources that are not in domain

        XMAX = XORIG + NXCELL * DXCELL
        YMAX = YORIG + NYCELL * DYCELL
        DO 22 I = 1, NELPT

C.............  Based on PTSRCE's expections, use <= and >= below
            IF( ENDX( I ) .NE. 0 ) THEN
                IDX = INDXIN( ENDX( I ) )
                IF( XLOCA( IDX ) .GE. XMAX  .OR.
     &              XLOCA( IDX ) .LE. XORIG .OR.
     &              YLOCA( IDX ) .GE. YMAX  .OR.
     &              YLOCA( IDX ) .LE. YORIG      ) ENDX( I ) = 0
            ENDIF

   22   CONTINUE

C.........  Create pointer IEND for sources with same FIP, PLANT, STACK
        I = 0
        DOWHILE( I .LE. NELPT-1 )
            I = I + 1
            IF( ENDX( I ) .NE. 0 ) THEN
                IDX = INDXIN( ENDX( I ) ) 

                J   = I + 1
                IGO = .TRUE.
                DOWHILE( IGO .AND. J .LE. NELPT )
                    IF( ENDX( J ) .NE. 0 ) THEN
                        JDX = INDXIN( ENDX( J ) )

                        IF( IFIP  ( IDX ) .EQ. IFIP  ( JDX ) .AND.
     &                      IPLANT( IDX ) .EQ. IPLANT( JDX ) .AND.
     &                      ISTACK( IDX ) .EQ. ISTACK( JDX ) ) THEN

                            IEND( I ) = IEND( I ) + 1
                            IEND( J ) = 0
                            J = J + 1
                        ELSE 
                            I   = J - 1
                            IGO = .FALSE. 
                        ENDIF
                    ELSE 
                        J = J + 1
                    ENDIF
                ENDDO
            ENDIF
        ENDDO

C.........  Prompt for name of output file

        ODEV = PROMPTFFILE(
     &        'Enter logical name for output file ',
     &        .FALSE., .TRUE., 'UAMPTS', 'SMK2PTS' )

C.........  Write file for use as input to PTSRCE
C.........  Control packet:

        WRITE( ODEV, 93010 ) PNAME( 1 )
        WRITE( ODEV, 93000 ) 'PTSOURCE  '
        WRITE( ODEV, 93000 ) FNOTE

C... .....  Count number of sources which will be output:

        NPOUT = 0
        DO 23 I = 1, NELPT
            IF( ENDX( I ) .NE. 0 .AND. IEND( I ) .NE. 0 )
     &          NPOUT = NPOUT + 1
   23   CONTINUE

C.........  Insert correct number of sources into 1st line of ESTRING
        BUFFER = ESTRING( 1 )
        WRITE( ESTRING( 1 ), 94000 ) BUFFER( 1:20 ), NPOUT,
     &                               BUFFER( 31:ISL )

        DO 25 I = 1, MXS
            BUFFER = ESTRING( I )
            M = TRIMLEN( BUFFER )
            WRITE( ODEV, 93000 ) ( BUFFER( 1:M ) )
   25   CONTINUE
 
        DO 33 I = 1, NSPECS
            WRITE( ODEV, 93020, IOSTAT=IOS ) SPCNAM( I )
   33   CONTINUE

        WRITE( ODEV, 93030, IOSTAT=IOS ) INT( IBGDAT ), INT( BEGTIM ),
     &                                   INT( IENDAT ), INT( ENDTIM )
        WRITE( ODEV, 93000 ) 'END' 

C.........  Region packet:
 
        WRITE( ODEV, 93010 ) PNAME( 2 )
        WRITE( ODEV, 93040 ) UTMX, UTMY, IZONE
        WRITE( ODEV, 93040 ) XORIG, YORIG

        IF( DXCELL .LT. 10. ) THEN    ! Print out for Lat-Lon
            WRITE( ODEV, 93045 ) DXCELL, DYCELL
        ELSE
            WRITE( ODEV, 93040 ) DXCELL, DYCELL
        ENDIF

        WRITE( ODEV, 93050 ) NXCELL, NYCELL, NZCELL
        WRITE( ODEV, 93060 ) NZLOWR, NZUPPR, HTSUR, 
     &                        HTLOWR, HTUPPR 
        WRITE( ODEV, 93000 ) 'END' 

C.........  Point Sources packet:

        WRITE( ODEV, 93010 ) PNAME( 3 )

C.............  Loop through elevated sources

        ICNT = 0 
        DO 44 I = 1, NELPT

            IF( ENDX( I ) .NE. 0 .AND. IEND( I ) .NE. 0 ) THEN
                ICNT = ICNT + 1

C.................  Set pointer
                IDX = INDXIN( ENDX( I ) ) 

C.................  Adjust units

                VH = STKVE( IDX ) * 3600.     ! m/s -> m/hr

C.................  Write records - format for 1st write depends on 
C.................  Lat-Lon or not

                IF( ABS( XLOCA( IDX ) ) .LT. 360. ) THEN

                    WRITE( ODEV, 93065, IOSTAT=IOS )
     &               ICNT, 'STD       ', XLOCA( IDX ), YLOCA( IDX ),
     &               IPLANT( IDX ), ISTACK( IDX ), IFIP( IDX )

                ELSE
                    WRITE( ODEV, 93070, IOSTAT=IOS )
     &               ICNT, 'STD       ', XLOCA( IDX ), YLOCA( IDX ),
     &               IPLANT( IDX ), ISTACK( IDX ), IFIP( IDX )
                ENDIF

                WRITE( ODEV, 93080, IOSTAT=IOS )
     &               STKHT( IDX ), STKDM( IDX ), STKTK( IDX ), VH

            ENDIF
   44   CONTINUE

        WRITE( ODEV, 93000 ) 'END' 

C.........  Initialize time and date for time step after JTIME and JDATE

        JDATE   = BDATE
        JTIME   = BTIME
        JDATEP1 = JDATE
        JTIMEP1 = JTIME
        CALL NEXTIME( JDATEP1, JTIMEP1, TSTEP )

C.........  Loop through time periods

        DO 55 T = 1, NSTEPS

C.............  Read emissions input file data

            IF( .NOT. READ3( UNAME, ALLVAR3, ALLAYS3, 
     &                   JDATE, JTIME, NELVIN ) ) THEN
                CALL M3EXIT( 'SMK2PTS', JDATE, JTIME,
     &                       'Could not read time step from ' 
     &                       // UNAME, 2 )
            ENDIF

C.............  Ensure that the number of sources has stayed the same
            
            IF( T .GT. 2 .AND. NELVIN .NE. NELVPREV ) THEN
                WRITE( MESG,94010 ) 'Number of elevated sources '//
     &                              'changed.  Cannot continue! '
                CALL M3EXIT( 'SMK2PTS', JDATE, JTIME, MESG, 2 )
            ENDIF
            NELVPREV = NELVIN

C.............  Time Interval packet:
C.............  Set required time parameters

            IBT = JTIME  /100
            IET = JTIMEP1/100
            IBD = JDATE   - 1900000
            IED = JDATEP1 - 1900000 

            WRITE( ODEV, 93010 ) PNAME( 4 )
            WRITE( ODEV, 93050 ) IBD, IBT, IED, IET

C.............  Method packet: 
C.............  Provide same output as PSTPNT

            WRITE( ODEV, 93010 ) PNAME( 5 )
            WRITE( ODEV, 93000 ) 
     &            'STD       ALL       EMVALUES  0.        50000.'
            WRITE( ODEV, 93000 ) 'END' 

C.............  Vertical Method packet:
C.............  Provide same output as PSTPNT with variable plume-rise method

            WRITE( ODEV, 93010 ) PNAME( 6 )
            WRITE( ODEV, 93000 ) 
     &            'STD       ALL       ' // VTYPE //
     &            ' 0.       10000.'
            WRITE( ODEV, 93000 ) 'END' 

C.............  Emissions Values packet:

            WRITE( ODEV, 93010 ) PNAME( 7 )
            WRITE( ODEV, 93000 ) 'ALL       ALL            0.000'

            ICNT = 0
            DO 53 I = 1, NELPT

C.................  Ensure that some species have non-zero emissions
                IF( ENDX( I ) .NE. 0 .AND. IEND( I ) .NE. 0 ) THEN

                    ICNT = ICNT + 1

                    DO 52 ISPC = 1, NSPECS  !  UAM species list

C.........................  Ensure that variable is in input file
                        IF( ISPTR( ISPC ) .NE. 0 ) THEN

C.............................  Sum up emissions in case of duplicate source
                            ESUM = 0.
                            DO 50 J = I, I + IEND( I ) - 1
                                IDX  = ENDX( J )
                                ESUM = ESUM + EMIS( IDX, ISPTR( ISPC ) )
   50                       CONTINUE

C.............................  Ensure emissions are not zero
                            IF( ESUM .NE. 0. ) THEN

C.................................  Write emissions in correct format
                                WRITE( ODEV, 93090 ) 
     &                                 ICNT, SPCNAM( ISPC ), ESUM
                            ENDIF
                        ENDIF
   52               CONTINUE
                ENDIF
   53       CONTINUE

            WRITE( ODEV, 93000 ) 'END' 
            WRITE( ODEV, 93000 ) 'ENDTIME'

C.............  Get next time and date

            JDATE = JDATEP1
            JTIME = JTIMEP1
            CALL NEXTIME( JDATEP1, JTIMEP1, TSTEP) 

   55   CONTINUE

      CALL M3EXIT( 'SMK2PTS', 0, 0,
     &             'Normal completion of PROGRAM SMK2PTS', 0 )


C******************    ERROR MESSAGES     ******************************


C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx
 
93000   FORMAT( A )
93010   FORMAT( A60 )
93020   FORMAT( A10 )
93030   FORMAT( 4I10 )
93040   FORMAT(	2F10.0, I10 )
93045   FORMAT( 2F10.7 )
93050   FORMAT( 4I10 )
93060   FORMAT( 2I10, 3F10.0 )
93065   FORMAT( I10, A10, F10.5, F10.5, I10.5, I10, I10.5 )
93070   FORMAT( I10, A10, F10.0, F10.0, I10.5, I10, I10.5 )
93080   FORMAT( F10.1, F10.2, F10.1, F10.0 )
93090   FORMAT( I10, A10, F10.3 )

C...........   Internal buffering formats.............94xxx
 
94000   FORMAT( A, I10, A )
94010   FORMAT( A )

      END

