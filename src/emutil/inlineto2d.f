
        PROGRAM INLINETO2D

C***********************************************************************
C  program body starts at line 
C
C  DESCRIPTION:
C    Program INLINETO2D reads STACK_GROUPS and INLN I/O API files and
C    converts them into a 2-D emission file for QA
C
C  PRECONDITIONS REQUIRED:  
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C    Original by G. Pouliot 2/25/2010
C    Revised 11/16/2011 to use GRIDDESC information for output file 
C    and ignore grid info from input files
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id: inlineto2d.f,v 1.19 2007/07/11 19:30:43 bbaek Exp $
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
C***********************************************************************
 
C.........  MODULES for public variables
C.........  This module contains the global variables for the 3-d grid
        USE MODSOURC, ONLY: XLOCA, YLOCA

        USE MODGRID, ONLY: GRDNM, COORD, GDTYP, P_ALP, P_BET, P_GAM,
     &                     XCENT, YCENT, XORIG, YORIG, XCELL, YCELL,
     &                     NCOLS, NROWS, NGRID

        USE MODGRDLIB

        IMPLICIT NONE

C...........   INCLUDES:
        INCLUDE 'EMCNST3.EXT'
        INCLUDE 'PARMS3.EXT'
        INCLUDE 'IODECL3.EXT'
        INCLUDE 'FDESC3.EXT'
      
C...........   EXTERNAL FUNCTIONS
        CHARACTER(2)  CRLF
        LOGICAL       ENVYN
        INTEGER       GETFLINE
        LOGICAL       GETYN       
        INTEGER       INDEX1
        INTEGER       LBLANK
        INTEGER       PROMPTFFILE
        CHARACTER(16) PROMPTMFILE
        INTEGER       SEC2TIME
        INTEGER       SECSDIFF
        LOGICAL       DSCM3GRD

        EXTERNAL CRLF, ENVYN, GETFLINE, GETYN, INDEX1, LBLANK,
     &           PROMPTFFILE, PROMPTMFILE, SEC2TIME, SECSDIFF, DSCM3GRD

C.........  LOCAL PARAMETERS and their descriptions:

        CHARACTER(50), PARAMETER :: 
     &  CVSW = '$Name: SMOKEv4.8.1_Jan2021$' ! CVS release tag

C...........   LOCAL VARIABLES and their descriptions:
C...........   Emissions arrays
        REAL,      ALLOCATABLE :: EIN ( : ,: )
        REAL,      ALLOCATABLE :: EOUT( : ,:, :)
        REAL,      ALLOCATABLE :: LAT(:), LON(:)
        INTEGER,   ALLOCATABLE :: LMAJOR(:)
        INTEGER,   ALLOCATABLE :: ROW(:), COL(:), ISTACK(:)
        INTEGER,   ALLOCATABLE :: NEW_COL(:), NEW_ROW(:)

C...........   Input file descriptors
        INTEGER       DURATA, DURATB ! no. time steps
        INTEGER       NCOLSA, NCOLSB! no. columns file a, file b, output file
        INTEGER       NROWSA, NROWSB  ! no. rows
        INTEGER       NVARSA, NVARSB ! no. variables
        INTEGER       SDATEA, SDATEB   ! start date
        INTEGER       STIMEA, STIMEB   ! start time
        INTEGER       NLAYSA, NLAYSB ! number of layers in the file
        INTEGER       NCOLS_N, NROWS_N

        REAL, ALLOCATABLE :: EMIS_SUM(:,:)
        CHARACTER(16), ALLOCATABLE :: VNAMEA( :),  VNAMEB(:) ! variable names
        CHARACTER(16), ALLOCATABLE :: VUNITA( :),  VUNITB(:) ! variable units
        CHARACTER(80), ALLOCATABLE :: VDESCA( :),  VDESCB(:) ! var descrip
        INTEGER      , ALLOCATABLE :: VTYPEA( :),  VTYPEB(:) ! var type

C...........   Other local variables 
        INTEGER       C, F, J, K, L, L1, L2, NL, V, T, S ! pointers and counters
        INTEGER       ISTART, IEND, ICT,V_R, V_I
        INTEGER       J_LOOP
        INTEGER       DUMMY                      ! dummy value for use 
        INTEGER       EDATE                      ! ending julian date
        INTEGER       ETIME                      ! ending time HHMMSS
        INTEGER       NSRC                       ! no of sources 
        INTEGER       IOS                        ! i/o status
        INTEGER       IREC                       ! line number count
        INTEGER       JDATE                      ! iterative julian date
        INTEGER       JTIME                      ! iterative time HHMMSS
        INTEGER       LB                         ! leading blanks counter
        INTEGER       LE                         ! location of end of string
        INTEGER       MXNFIL_1                   ! max no. of stack groups files
        INTEGER       MXNFIL_2                   ! max no. of emission files
        INTEGER       MXNFIL                     ! max no. of both
        INTEGER       NFILE                      ! no. of 2-d input files
        INTEGER       NSTEPS                     ! no. of output time steps
        INTEGER       NVOUT(2)                   ! no. of output variables
        INTEGER       NVOUT_I(2)                 ! no. of integer type output variables        
        INTEGER       NVOUT_R(2)                 ! no. of real type output variables        
        INTEGER       RDATE                      ! reference date
        INTEGER       SAVLAYS                    ! number of layers
        INTEGER       SDATE                      ! starting julian date
        INTEGER       SECS                       ! tmp seconds
        INTEGER       SECSMAX                    ! seconds maximum
        INTEGER       SECSMIN                    ! seconds minimum
        INTEGER       STIME                      ! starting time HHMMSS
        INTEGER       STEPS                      ! tmp number of steps
        INTEGER       TIMET                      ! tmp time from seconds
        INTEGER       TSTEP                      ! time step
        INTEGER       VLB                        ! VGLVS3D lower bound 
        INTEGER       TMAX                       ! max time counter
        INTEGER       LDEV, S_CT

        CHARACTER(16)  ONAME                     ! ONAME
        CHARACTER(16)  FDESC                     ! tmp file description
        CHARACTER(16)  NAM                       ! tmp file name
        CHARACTER(16)  VNM                       ! tmp variable name
        CHARACTER(256) LINE                      ! input buffer
        CHARACTER(256) MESG                      ! message field
        CHARACTER(15)  RPTCOL                    ! single column in report line
        CHARACTER(300) RPTLINE                   ! line of report file

        LOGICAL    :: EFLAG   = .FALSE.   ! error flag
        LOGICAL    :: FIRST3D = .TRUE.    ! true: first 3-d file not yet input
        LOGICAL    :: LFLAG   = .FALSE.   ! true iff 3-d file input
        LOGICAL    :: TFLAG   = .FALSE.   ! true: grid didn't match
        LOGICAL    :: LATLON_IN_FILE 

        CHARACTER(80) :: GDESC        ! grid description
        CHARACTER(IOVLEN3) COORD3D    ! coordinate system name
        CHARACTER(IOVLEN3) COORUN3D   ! coordinate system projection units
        CHARACTER(80) :: GDNAM

        CHARACTER(16) :: PROGNAME = 'INLINETO2D' ! program name

C***********************************************************************
C   begin body of program INLINETO2D
 
        LDEV = INIT3()
 
C.........  Write out copyright, version, web address, header info, and prompt
C           to continue running the program.
        CALL INITEM( LDEV, CVSW, PROGNAME )

        LATLON_IN_FILE = .TRUE.

        IF ( .NOT. OPEN3( 'INLN', FSREAD3, PROGNAME )) THEN
            MESG = 'Could not open file INLN ' 
            CALL M3MSG2( MESG )
            MESG = 'Ending program "INLINETO2D".'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF   

        IF ( .NOT. DESC3( 'INLN' ) ) THEN
            MESG = 'Could not get description of file INLN ' 
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ELSE
            NROWSA = NROWS3D
            NCOLSA = NCOLS3D
            NLAYSA = NLAYS3D
            NVARSA = NVARS3D
            SDATEA = SDATE3D
            STIMEA = STIME3D
            DURATA = MXREC3D
            TSTEP = TSTEP3D

            ALLOCATE (VNAMEA(  NVARSA), STAT=IOS)
            ALLOCATE (VUNITA(  NVARSA), STAT=IOS)
            ALLOCATE (VDESCA(  NVARSA), STAT=IOS)
            ALLOCATE (VTYPEA(  NVARSA), STAT=IOS)

            DO V = 1, NVARSA
                VNAMEA( V ) = VNAME3D( V )
                VUNITA( V ) = UNITS3D( V )
                VDESCA( V ) = VDESC3D( V )
                VTYPEA( V ) = VTYPE3D( V )
            END DO
        END IF

        NSRC = NROWSA
        ALLOCATE (EIN( NSRC, NVARSA), STAT=IOS)            
        CALL CHECKMEM( IOS, 'EIN', PROGNAME )
        ALLOCATE (EMIS_SUM(NVARSA,3))
            
        IF ( .NOT. OPEN3( 'STACK_GROUPS', FSREAD3, PROGNAME )) THEN
            MESG = 'Could not open file INLN ' 
            CALL M3MSG2( MESG )
            MESG = 'Ending program "INLINETO2D".'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF   
 
        IF ( .NOT. DESC3( 'STACK_GROUPS' ) ) THEN
            MESG = 'Could not get description of file STACK_GROUPS ' 
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ELSE
            NROWSB = NROWS3D
            NCOLSB = NCOLS3D
            NLAYSB = NLAYS3D
            NVARSB = NVARS3D
            SDATEB = SDATE3D
            STIMEB = STIME3D
            DURATB = MXREC3D

            ALLOCATE (VNAMEB(  NVARSB), STAT=IOS)
            ALLOCATE (VUNITB(  NVARSB), STAT=IOS)
            ALLOCATE (VDESCB(  NVARSB), STAT=IOS)
            ALLOCATE (VTYPEB(  NVARSB), STAT=IOS)

            DO V = 1, NVARSB
                VNAMEB( V ) = VNAME3D( V )
                VUNITB( V ) = UNITS3D( V )
                VDESCB( V ) = VDESC3D( V )
                VTYPEB( V ) = VTYPE3D( V )
            END DO

        END IF

        ALLOCATE (COL(NSRC))
        ALLOCATE (ROW(NSRC))
        ALLOCATE (XLOCA(NSRC))            
        ALLOCATE (YLOCA(NSRC))
        ALLOCATE (ISTACK(NSRC))
        ALLOCATE (LAT(NSRC))
        ALLOCATE (LON(NSRC))
        ALLOCATE (LMAJOR(NSRC))
        ALLOCATE (NEW_COL(NSRC))
        ALLOCATE (NEW_ROW(NSRC))

        IF( .NOT. READ3( 'STACK_GROUPS', 'COL', 1, SDATEB, STIMEB, COL)) THEN
           MESG = 'Could not read  COL from STACK_GROUPS'
           CALL M3EXIT( PROGNAME, 0, 0,MESG, 2 )
        ENDIF

        IF( .NOT. READ3( 'STACK_GROUPS', 'ROW', 1, SDATEB, STIMEB, ROW)) THEN
           MESG = 'Could not read  ROW from STACK_GROUPS'
           CALL M3EXIT( PROGNAME, 0, 0,MESG, 2 )
        ENDIF

        IF( .NOT. READ3( 'STACK_GROUPS', 'XLOCA', 1,SDATEB, STIMEB, XLOCA)) THEN
           MESG = 'Could not read  XLOCA from STACK_GROUPS'
           CALL M3EXIT( PROGNAME, 0, 0,MESG, 2 )
        ENDIF

        IF( .NOT. READ3( 'STACK_GROUPS', 'YLOCA', 1, SDATEB, STIMEB, YLOCA)) THEN
           MESG = 'Could not read  YLOCA from STACK_GROUPS'
           CALL M3EXIT( PROGNAME, 0, 0,MESG, 2 )
        ENDIF

        IF (LATLON_IN_FILE) THEN
            IF( .NOT. READ3( 'STACK_GROUPS', 'LATITUDE', 1, SDATEB, STIMEB, LAT)) THEN
               MESG = 'Could not read  LAT from STACK_GROUPS'
               CALL M3EXIT( PROGNAME, 0, 0,MESG, 2 )
            ENDIF

            IF( .NOT. READ3( 'STACK_GROUPS', 'LONGITUDE', 1, SDATEB, STIMEB, LON)) THEN
               MESG = 'Could not read  LON from STACK_GROUPS'
               CALL M3EXIT( PROGNAME, 0, 0,MESG, 2 )
            ENDIF
        ENDIF ! latlon in file check
                        
        IF( .NOT. READ3( 'STACK_GROUPS', 'ISTACK', 1, SDATEB, STIMEB, ISTACK)) THEN
            MESG = 'Could not read  ISTACK from STACK_GROUPS'
            CALL M3EXIT( PROGNAME, 0, 0,MESG, 2 )
        ENDIF

        IF( .NOT. READ3( 'STACK_GROUPS', 'LMAJOR', 1, SDATEB, STIMEB, LMAJOR)) THEN
            MESG = 'Could not read  LMAJOR from STACK_GROUPS'
            CALL M3EXIT( PROGNAME, 0, 0,MESG, 2 )
        ENDIF
                        
C***** GET GRID INFORMATION FOR OUTPUT FILE AND OPEN OUTPUT FILE
        IF( .NOT. DSCM3GRD( GDNAM3D, GDESC, COORD3D, GDTYP3D, COORUN3D,
     &                P_ALP3D, P_BET3D, P_GAM3D, XCENT3D, YCENT3D,
     &                XORIG3D, YORIG3D, XCELL3D, YCELL3D,
     &                NCOLS3D, NROWS3D, NTHIK3D ) ) THEN

            MESG = 'Could not get Models-3 grid description.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ELSE
            GDTYP = GDTYP3D
            P_ALP = P_ALP3D
            P_BET = P_BET3D
            P_GAM = P_GAM3D
            XCENT = XCENT3D
            YCENT = YCENT3D
            NCOLS = NCOLS3D
            NROWS = NROWS3D
            NGRID = NCOLS * NROWS

        END IF    

        IF (LATLON_IN_FILE) THEN
            XLOCA(1:NSRC) = LON(1:NSRC)
            YLOCA(1:NSRC) = LAT(1:NSRC)

            GDNAM = GDNAM3D

            CALL CONVRTXY( NSRC, GDTYP, GDNAM, P_ALP, P_BET, P_GAM, 
     &                     XCENT, YCENT, XLOCA, YLOCA )

            DO S = 1, NSRC
                IF( .NOT. INGRID( XLOCA(S), YLOCA(S), NCOLS_N, NROWS_N,
     &                            COL(S), ROW(S) ) ) THEN
                    CYCLE  ! To end of loop
                ENDIF  
            ENDDO

        ENDIF    ! latlon in file check
             
        ALLOCATE (EOUT(NCOLS, NROWS, NVARSA))
        EOUT(1:NCOLS,1:NROWS,1:NVARSA) = 0.0
        
C.........  Set up layer structure for output file
        SDATE3D = SDATEA
        STIME3D = STIMEA
        NVARS3D = NVARSA
        GDTYP3D = GDTYP
        P_ALP3D = P_ALP            
        P_BET3D = P_BET
        P_GAM3D = P_GAM
        XCENT3D = XCENT
        YCENT3D = YCENT
        NCOLS3D = NCOLS
        NROWS3D = NROWS       

C.........  Set up layer structure for output file
        NLAYS3D = 1
        VGLVS3D = 0.     ! initialize array
                DO V = 1, NVARSA
                    VNAME3D( V ) = VNAMEA( V )
                    UNITS3D( V ) = VUNITA( V )
                    VDESC3D( V ) = VDESCA( V )
                    VTYPE3D( V ) = VTYPEA( V )
                END DO

        ONAME = PROMPTMFILE( 
     &          'Enter logical name for OUTPUT file',
     &          FSUNKN3, 'OUTFILE', PROGNAME )        
            
        JDATE = SDATEA
        JTIME = STIMEA
        TMAX = DURATA
        
        EMIS_SUM(1:NVARSA,1:3) = 0.0

        DO T = 1, TMAX

            EOUT(1:NCOLS,1:NROWS,1:NVARSA) = 0.0

            DO V = 1, NVARSA
               IF( .NOT. READ3( 'INLN', VNAMEA(V), 1, JDATE,  
     &                           JTIME, EIN(1,V))) THEN
                   MESG = 'Could not read "' // VNAMEA(V) //
     &                    '" from file INLN "' 
                   CALL M3EXIT( PROGNAME, RDATE, JTIME, MESG, 2 )
               ENDIF
               
               S_CT = 0
               DO S = 1, NSRC
                 IF (COL(S) .lt. 1) CYCLE
                 IF (ROW(S) .lt. 1) CYCLE
                 IF (COL(S) .gt. ncols) CYCLE
                 IF (ROW(S) .gt. nrows) CYCLE
                 S_CT = S_CT + 1
                 EOUT(COL(S),ROW(S),V) = EOUT(COL(S),ROW(S),V) + EIN(S,V)
                 EMIS_SUM(V,3) = EMIS_SUM(V,3) + EIN(S,V)
               END DO

               EMIS_SUM(V,1) = EMIS_SUM(V,1) + SUM(EIN(1:NSRC,V))
               EMIS_SUM(V,2) = EMIS_SUM(V,2) + SUM(EOUT(1:NCOLS,1:NROWS,V))

               IF( .NOT. WRITE3( ONAME, VNAMEA(V), JDATE,  
     &                           JTIME, EOUT(1,1,V))) THEN
                    MESG = 'Could not write "' // VNAMEA(V) //
     &                     '" to file OUTFILE "' 
                    CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
               ENDIF
               
            END DO   ! loop through variables
            
            CALL NEXTIME( JDATE, JTIME, TSTEP )     

        END DO       ! loop through timesteps

        CALL M3EXIT( PROGNAME, 0, 0, ' ', 0)
    
C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )
 
C...........   Formatted file I/O formats............ 93xxx


93000   FORMAT(  A )

93010   FORMAT( A15 )

93020   FORMAT( I15 )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I7, :, 1X ) )

94020   FORMAT( A, :, I3, :, 1X, 10 ( A, :, F8.5, :, 1X ) )

        END PROGRAM INLINETO2D
