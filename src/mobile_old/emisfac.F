
        PROGRAM EMISFAC

C***********************************************************************
C  program body starts at line 246
C
C  DESCRIPTION:
C      This program drives the MOBILE5 program by supplying an ambient 
C      temperature and an input file name from which mobile5 derives 
C      all of its input information.  Each call to MOBILE is a single 
C      temperature solution for all input scenarios.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Prototype 11/95 by MRH
C
C  FUTURE NEEDS:
C      (1) Need to create better analysis of IPSIR combo emission factors to 
C          ensure that no loops and better yet, no commons using commons.
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

        INCLUDE 'CHDGEN3.EXT'    !  Mobile include paramaters
        INCLUDE 'MBDIMS3.EXT'    !  Mobile include paramaters
        INCLUDE 'PARMS3.EXT'     !  I/O API parameters
        INCLUDE 'IODECL3.EXT'    !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'     !  I/O API file description data structures.
        INCLUDE 'NETCDF.EXT'     !  NetCDF declarations

C...........   EXTERNAL FUNCTIONS and their descriptions:

        CHARACTER*2  CRLF 
        INTEGER      ENVINT
        REAL         ENVREAL
        LOGICAL      ENVYN
        INTEGER      FIND1
        INTEGER      GETEFILE
        INTEGER      GETNUM
        REAL         GETREAL
        LOGICAL      GETYN 
        INTEGER      INDEX1
        INTEGER      LOCAT1
        INTEGER      PROMPTFFILE
        CHARACTER*16 PROMPTMFILE
        INTEGER      SEC2TIME
        INTEGER      STR2INT
        INTEGER      TRIMLEN

        EXTERNAL     CRLF, ENVINT, ENVYN, FIND1, GETEFILE, GETNUM, 
     &               GETREAL, GETYN, INDEX1, LOCAT1, PROMPTFFILE, 
     &               PROMPTMFILE, SEC2TIME, STR2INT, TRIMLEN

C...........   EXTERNAL BLOCK DATA for mobile:

        EXTERNAL      BD01, BD02, BD03, BD04, BD05, BD06, BD10,
     &                BD11, BD12, BD13, BD14, BD15, BD16, BD17,
     &                BD18, BD19, BD20, BD21, BD22, BD23, BD24,
     &                BD25, BD26, BD27, BD28, BD29, BD30, BD31,
     &                BD32, BD33, BD34, BD35, BD36, BD37, BD38

C...........   LOCAL PARAMETERS and their descriptions:

        CHARACTER*5   BLANK5
        PARAMETER(    BLANK5 = ' ' )

C...........   LOCAL VARIABLES and their descriptions:

C...........   Simulation-specific MOBILE5 parmeter list variables
        INTEGER IPSIR ( MXPSI )      ! Sorted   scenario-specific PSI list
        INTEGER MPNUM                ! Actual count of scenario-specific PSIs
        INTEGER P1, P2, P3, P4( 24 ) ! Dummies

C...........   MOBILE5 parameter x-ref variables
        INTEGER MRNUM                     ! Actual no. of normal x-ref entries
        INTEGER INDEXA( MXPREF )          ! Index for sorting x-ref scheme index
        INTEGER IPSI  ( MXPREF )          ! Sorted   x-ref scheme index
        INTEGER IPSIA ( MXPREF )          ! Unsorted x-ref scheme index
        INTEGER ISTYPE( MXPREF, 2     )   ! Type-specifier (1), and pointer (2)
        INTEGER NIDX  ( MXPREF, NNFAC )   ! Normal x-ref entries
        INTEGER CIDX  ( MXPREF, NCFAC )   ! Control Normal x-ref entries
        INTEGER M5PSI ( MXM5PSI )         ! Temp input array for MOBILE5 packets

C...........   Variables for storing PSI descriptors:

        REAL    PCNTPRO( MXNPT, MXPREF )
        INTEGER IPLIST ( MXNPT, MXPREF )
        INTEGER ISLINE ( MXPREF )         ! Line number of MOBILE5 packet for PSIs
        INTEGER ISROOT ( MXPREF )         ! Root PSI for multi-scenario MOBILE5 pckt
        INTEGER ISCNT  ( MXPREF )         ! Count of PSIs multi-scenario
        INTEGER MPTYPE

C...........   Master Non-Diurnal EFs list variables
        REAL    EXH_CO ( MAXTMP, NVTYPE )  ! Non-Diurn Emis fac: CO exhaust
        REAL    EXH_NOX( MAXTMP, NVTYPE )  ! Non-Diurn Emis fac: NOX exhaust
        REAL    EXH_VOC( MAXTMP, NVTYPE )  ! Non-Diurn Emis fac: VOC exhaust
        REAL    EVP_VOC( MAXTMP, NVTYPE )  ! Non-Diurn Emis fac: VOC evaporative
        REAL    RNL_VOC( MAXTMP, NVTYPE )  ! Non-Diurn Emis fac: VOC running
        REAL    RST_VOC( MAXTMP, NVTYPE )  ! Non-Diurn Emis fac: VOC resting
        REAL    RFL_VOC( MAXTMP, NVTYPE )  ! Non-Diurn Emis fac: VOC refueling

        COMMON / EFS1 / EXH_CO , EXH_NOX, EXH_VOC, EVP_VOC, RNL_VOC, 
     &                  RST_VOC, RFL_VOC

C...........   Master Diurnal EFs list variables

        REAL    WDL_VOC( NMINMAX, NVTYPE )  ! Diurnal Emis fac: VOC weighted
        REAL    DNL_VOC( NMINMAX, NVTYPE )  ! Diurnal Emis fac: VOC diurnal
        REAL    HOT_VOC( NMINMAX, NVTYPE )  ! Diurnal Emis fac: VOC hot soak
        REAL    CRC_VOC( NMINMAX, NVTYPE )  ! Diurnal Emis fac: VOC crank case

        COMMON / EFS2 / WDL_VOC, DNL_VOC, HOT_VOC, CRC_VOC

C...........   Variables for targeted update of EFs
        INTEGER NUPDAT                      ! No. of PSIs targeted for update
        INTEGER UPDATE ( MXPSI )            ! Sorted   targeted PSI list
        INTEGER UPDATA ( MXPSI )            ! Unsorted targted PSI list
        INTEGER UPIDXA ( MXPSI )            ! Targeted PSI list index

C...........   Variables to carry over from one MOBILE5 call to the next

        INTEGER      IREAD
        REAL         CR12HC
        REAL         CR12CO
        CHARACTER*40 TNAME
        COMMON /IMPAR7/ TNAME,IREAD
        COMMON /IM12HC/ CR12HC(19,20,5,2)
        COMMON /IM12CO/ CR12CO(19,20,5,2)

C...........   Variables for storing temperature index:

        INTEGER   NTREF                     ! No of temperature index entries
        INTEGER   MMTC   (         MXPSI )  ! Counter for current temp combo
        INTEGER   MMTREF (         MXPSI )  ! EF reference IDs with temp index
        INTEGER   NTIDX  (         MXPSI )  ! number of indices per EF ref
        INTEGER   TIDXPSI( MXTMMI, MXPSI )  ! temperature indices

C...........   Unit numbers and logical names

        INTEGER       EDEV
        INTEGER       LDEV
        INTEGER       MODEV
        INTEGER       NLDEV
        INTEGER       NRDEV
        INTEGER       SDEV    ! Optional input for targeted update on EF files
        CHARACTER*16  DNAME   ! Input  diurnal master EF file
        CHARACTER*16  ENAME   ! Input  non-diu master EF file
        CHARACTER*16  NNAME   ! Output diurnal master EF file
        CHARACTER*16  ONAME   ! Output non-diu master EF file

C...........   Other local variables

        REAL    EFACT ( MXNDET, NVTYPE )  ! Temporary non-diurnal EFs 
        REAL    DFACT ( MXDET , NVTYPE )  ! Temporary diurnal EFs
        REAL    EFSAVE( NVTYPE, MXM5PSI, SAVTMP, MXEMIS ) 
                                          ! Saves EFs from multi-scen MOB5 run
        REAL    TEMPIN( MAXTMP )          ! Input temperatures in deg F

        REAL    BTEMP                     ! Beginning temp for non-diurnal
        REAL    ETEMP                     ! Ending temp for non-diurnal
        REAL    MAXTL                     ! intermediate temperature index value
        REAL    TF                        ! temporary temperature
        REAL    TF_MAX                    ! temporary maximum temperature
        REAL    TF_MIN                    ! temporary minimum temperature

        INTEGER CNTPSI                    ! Counter of PSIs on MOBILE5 packets of NPREF
        INTEGER DEPCNT                    ! Temporary count for combo emis facs
        INTEGER DEPOINT                   ! Temporary pointer to combo-only lists
        INTEGER EMTIDX( MXEMIS )          ! Index to EMTLST for actual EFs
        INTEGER I, II, IM, IPM            ! Counters and pointers
        INTEGER J, JJ, K, KK              ! Counters and pointers
        INTEGER L, L1, L2, N, T           ! Counters and pointers
        INTEGER IERR                      ! MOBILE5 error value
        INTEGER IOS                       ! I/O status
        INTEGER IPARM                     ! temporary PSI value 
        INTEGER IPTIM                     ! IPARM converted to time for file ops
        INTEGER IREC                      ! Line number when reading ASCII files
        INTEGER ITEMP                     ! Indicator of tmpr status for MOBILE5 
        INTEGER JYEAR                     ! year to calculate emission factors
        INTEGER MA, MB                    ! indicator for MOBILE5A or MOBILE5b packets
        INTEGER MIMICHK                   ! for checking infile MINT_MIN
        INTEGER MIMACHK                   ! for checking infile MINT_MAX
        INTEGER MAMICHK                   ! for checking infile MAXT_MIN
        INTEGER MAMACHK                   ! for checking infile MAXT_MAX
        INTEGER NTEMP                     ! Actual no. of temps for eval of EFs
        INTEGER P5( MXTMMI )              ! Dummy for reading EF-based TMMI vals
        INTEGER PSI                       ! Temporary PSI
        INTEGER ROOT       
        INTEGER TMMI                      ! temporary min/max temp index
        INTEGER TMMIDV                    ! temporary index divided by MAXTL
        INTEGER VOLORG                    ! Index of Mobile5 hydrocarbon type

        LOGICAL DFLAG                     ! Calculate diurnal factors when true
        LOGICAL EFLAG                     ! Error exists when true
        LOGICAL FIRSTMPR                  ! First tmptr of current PSI when true
        LOGICAL MOB5AON                   ! MOBILE5A packet has been found
        LOGICAL MOB5BON                   ! MOBILE5B packet has been found
        LOGICAL PROMPTF                   ! iff PROMPTFLAG E.V. true or not defined
        LOGICAL RFLAG                     ! Diurnal facs _were_ calc'd when true
        LOGICAL SFLAG( MXPSI )            ! Update diurnal combos when true
        LOGICAL TFLAG                     ! Non-diu facs _were_ calc'd when true
        LOGICAL VFLAG( MXPSI )            ! Update non-diu combos when true
        LOGICAL UFLAG                     ! Update EFs when true

        CHARACTER*3   VOLDEFLT            !  Default Mobile hydrocarbon
        CHARACTER*3   VOLNAM              !  name of Mobile5 hydrocarbon
        CHARACTER*8   M5VRSION            !  Run time MOBILE version
        CHARACTER*30  BUFFER              !  input buffer
        CHARACTER*256 LINE                !  line input buffer
        CHARACTER*256 MESG                !  message buffer

C...........   Commons linked with MOBILE subroutine
        COMMON / MULTISC / EFSAVE, M5VRSION 
        COMMON / TOVRIDE / TF, TF_MIN, TF_MAX

C***********************************************************************
C   begin body of program EMISFAC

        LDEV = INIT3()

        CALL INITEM( LDEV )
 
        WRITE( *,92000 )
     &  ' ',
     &  'Program EMISFAC drives the MOBILE5 program by supplying',
     &  'a range of ambient temperatures and a scenario-specific',
     &  'MOBILE5 parameter file MPLIST. Using multiple calls',
     &  'to MOBILE5, EMISFAC creates a master emission factors',
     &  'table for the specified temperatures and input parameter',
     &  'conbinations.',
     &  ' ',
     &  'You will need to enter the logical names for the input and',
     &  'output files (and to have set them prior to program launch,',
     &  'using "setenv <logicalname> <pathname>").',
     &  'Optional checking that all species are reported for each ',
     &  'source may be turned on via "setenv RAW_SRC_CHECK Y".',
     &  ' ',
     &  'You may use END_OF-FILE (control-D) to quit the program',
     &  'during logical-name entry. Default responses are given in',
     &  'brackets [LIKE THIS] and can be accepted by hitting the',
     &  '<RETURN> key.',
     &  ' '

        IF ( .NOT. GETYN( 'Continue with program?', .TRUE. ) ) THEN
            CALL M3EXIT( 'EMISFAC', 0, 0, 'Ending program', 2 )
        END IF

C.......   Get environment variable logical variables
        PROMPTF = ENVYN ( 'PROMPTFLAG', 'Prompt for inputs or not',
     &                    .TRUE., IOS )

C.........  Turn off verbose NetCDF to elimiate error messages when using
C.........  CHECK3 to check for available time steps
        CALL NCPOPT( 0 )

C.........  Get MOBILE5 input file (no prompt)
        MODEV = GETEFILE( 'MO5AIN', .FALSE., .TRUE., 'EMISFAC' )

C.........  Prompt for diurnal master list EFs input file for modeling year
        ENAME = PROMPTMFILE(
     &          'Enter logical name for NON-DIURNAL EF input '//
     &          'file or "NONE"', FSRDWR3, 'MEFSND_IN', 'EMISFAC' )

C.........  Prompt for diurnal master list EFs input file for modeling year
        DNAME = PROMPTMFILE(
     &          'Enter logical name for DIURNAL EF input ' //
     &          'file or "NONE"', FSRDWR3, 'MEFSD_IN', 'EMISFAC' )

        SDEV = -9
        IF( DNAME .NE. 'NONE' .OR. ENAME .NE. 'NONE' ) THEN

            SDEV = PROMPTFFILE(
     &          'Enter logical name for TARGETED EF UPDATE ' //
     &          'file or "NONE"', .TRUE., .TRUE., 'MEFUPD', 'EMISFAC' )

        ENDIF

C.........  Prompt for EF-based file of temp combination indices
        EDEV  = PROMPTFFILE(
     &          'Enter logical name for output EF-REF/TEMP file',
     &          .TRUE., .TRUE., 'MEFTEMP', 'EMISFAC' )

C.........  Prompt for MPLIST (scenario specific parameter scheme list)
        NLDEV = PROMPTFFILE(
     &          'Enter logical name for PARAMETER LIST input file',
     &          .TRUE., .TRUE., 'MPLIST', 'EMISFAC' )
        
C.........  Prompt for MPREF (general MOBILE5 parameter x-ref file)
        NRDEV = PROMPTFFILE( 
     &          'Enter logical name for MPREF input file',
     &          .TRUE., .TRUE., 'MPREF', 'EMISFAC' )

C.........  If starting with existing non-diurnal emission factors file,
C.........  set parameters that would otherwise get from prompt
        IF( ENAME .NE. 'NONE' ) THEN

C.............  Get file description
            IF( .NOT. DESC3( ENAME ) ) THEN

                WRITE( MESG, 94010 ) 
     &                 'Could not get description of ' // ENAME
                CALL M3EXIT( 'EMISFAC', 0, 0, MESG, 2 )

            ENDIF

C.........  Set start and ending temperature depending on input file
C.............  Set output file name to same as input
            ONAME = ENAME

C.............  Set modeling year
            JYEAR = SDATE3D / 1000

            WRITE( MESG,94010 ) 
     &          'NOTE: Emission factors are for year', JYEAR
            CALL M3MSG2( MESG )

C.............  Check that file delta temperature is same as compiled value
            IF( XCELL3D .NE. TINTV ) THEN
                WRITE( MESG, 94030 ) 
     &                 'Delta temperature for file ' // 
     &                 ENAME( 1:TRIMLEN( ENAME ) ) // CRLF()// BLANK5//
     &                 '  not equal to compiled value of', TINTV
                CALL M3EXIT( 'EMISFAC', 0, 0, MESG, 2 )

            ENDIF

            NTEMP = STR2INT( FDESC3D( 5 ) )

            BTEMP = - SNGL( XORIG3D )
            ETEMP = BTEMP + TINTV * FLOAT( NTEMP - 1 ) 

            DO 22 I = 1, NTEMP
                TEMPIN( I ) = BTEMP + ( I - 1 ) * TINTV
22          CONTINUE

C.............  Set the volatile pollutant used in this file

            L = TRIMLEN( FDESC3D( 7 ) )
            J = INDEX1 ( FDESC3D( 7 )( 1:L ), NUMVOL, VOLTYP )

            IF( J .GT. 0 ) THEN
                VOLORG = J
                CALL M3MSG2( 
     &               'NOTE: Emissions factors are for organic pollutant'
     &               // ' "' // VOLTYP(J)( 1:TRIMLEN( VOLTYP(J) ) )
     &               // '"' )
            ELSE

                WRITE( MESG,94010 ) 
     &             'ERROR: Volatile pollutant type "' //
     &             FDESC3D( 7 )( 1:L ) // '" in file "' //
     &             ENAME( 1:TRIMLEN( ENAME ) ) // 
     &             '" is not a valid volatile pollutant type.'

                CALL M3EXIT( 'EMISFAC', 0, 0, MESG, 2 )
                  
            ENDIF

C.........  If ENAME does _not_ exist, build the parameters from scratch
        ELSE 

C.......   Get default characteristic for output file:
C.......   If we're going to prompt, then set the defaults as best we can
C.......      otherwise, use environment variables to set defaults

            EFLAG = .FALSE.
            IF( PROMPTF ) THEN
                VOLDEFLT = VOLTYP( 3 )
                JYEAR    = 1988
                BTEMP    = 50.
                ETEMP    = 110.

            ELSE

                CALL ENVSTR( 'MB_HC_TYPE', 'Volatile pollutant type', 
     &                       VOLTYP( 3 ) , VOLDEFLT, IOS )
                IF( IOS .NE. 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94020 ) 'MB_HC_TYPE'
                    CALL M3MSG2( MESG )
                ENDIF

                JYEAR = ENVINT( 'EF_YEAR', 
     &                          'Emission factors year', 1988, IOS )
                IF( IOS .NE. 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94020 ) 'EF_YEAR'
                    CALL M3MSG2( MESG )
                ENDIF

                BTEMP = ENVREAL( 'EF_MIN_TEMP', 'Minimum emission ' //
     &                           'factors temperature', 50., IOS )
                IF( IOS .NE. 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94020 ) 'EF_MIN_TEMP'
                    CALL M3MSG2( MESG )
                ENDIF

                ETEMP = ENVREAL( 'EF_MAX_TEMP', 'Minimum emission ' //
     &                           'factors temperature', 100., IOS )
                IF( IOS .NE. 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94020 ) 'EF_MAX_TEMP'
                    CALL M3MSG2( MESG )

                ENDIF
     
                IF( EFLAG ) THEN
                    CALL M3EXIT( 'EMISFAC', 0, 0,
     &                   'Bad environment variable setting(s)', 2 )
                ENDIF

            ENDIF

C.............  Write out list of legal hydrocarbon species
            MESG = VOLTYP( 1 )
            DO 33 I = 2, NUMVOL
                WRITE( MESG, 94010 ) MESG( 1:TRIMLEN( MESG ) ) // 
     &                           CRLF() // BLANK5 // VOLTYP( I )
33          CONTINUE 

            CALL M3MSG2( MESG )

C.........  Prompt for type of volatile pollutant to use (numbers match
C.........  MOBILE5 numbers)

            DO 40 I = 5, 1, -1
                CALL GETSTR( 'Enter type of hydrocarbon emissions', 
     &                       VOLDEFLT, VOLNAM )
 
                VOLORG = INDEX1( VOLNAM, NUMVOL, VOLTYP )

                IF( VOLORG .GT. 0 ) GO TO 41  ! Valid reponse found

                WRITE( MESG,94010 ) 'Invalid response, get', I-1, 
     &                              'more tries'
                CALL M3MSG2( MESG )
 
40          CONTINUE

            WRITE( MESG,94010 ) 'Bad hydrocarbon emissions type'
            CALL M3EXIT( 'EMISFAC', 0, 0, MESG, 2 )

41          CONTINUE

C............  Find VOLNAM in EMTLST for setting output variable names

            EMTIDX( 1 ) = 1   ! For CO
            EMTIDX( 2 ) = 2   ! For NOX
            K = 2

            DO 44 I = 3, MXEALL

                L = INDEX( EMTLST( I ), '_' )
                IF( EMTLST( I )(L+1:7 ) .EQ. VOLNAM ) THEN
                    K = K + 1
                    EMTIDX( K ) = I
                ENDIF

44          CONTINUE

C.............  Prompt for modeling year
            JYEAR = GETNUM( 0, 9999, JYEAR,
     &             'Enter modeling year for emission factors (YYYY)' )

            CALL M3MSG2( 
     &       'Enter temperatures for non-diurnal emission factors:' )

C.............  Prompt for start temperature
            BTEMP = GETREAL( MINT_MIN, MINT_MAX, BTEMP,
     &                   'Enter minimum temperature [deg F]' )

C.............   Prompt for end temperature
            ETEMP = GETREAL( BTEMP, MAXT_MAX, ETEMP, 
     &                   'Enter maximum temperature [deg F]' )

C.............  Set temperatures
            I = 0
50          CONTINUE

                TF = BTEMP + TINTV * I

                IF( TF .GT. ETEMP ) GO TO 51

                I = I + 1
                TEMPIN( I ) = TF

            GO TO 50

51          CONTINUE
            NTEMP = I

C.............  Ensure no overflow of temperatures
            IF( NTEMP .GT. MAXTMP ) THEN

                WRITE( MESG,94010 )
     &                 'Number of temperatures requested:', NTEMP,
     &                 CRLF() // BLANK5 //
     &                 'Max dimensioned number (MAXTMP):', MAXTMP,
     &                 CRLF() // BLANK5 //
     &                 'Max number of temperatures exceeded.'
                CALL M3EXIT( 'EMISFAC', 0, 0, MESG, 2 )

            ENDIF

C.............  Set up non-diurnal master list file description

            FTYPE3D = GRDDED3
            P_ALP3D = 0.D0
            P_BET3D = 0.D0
            P_GAM3D = 0.D0
            XCENT3D = 0.D0
            YCENT3D = 0.D0
            XORIG3D = -DBLE( BTEMP )  ! Negative starting temp (- for PAVE)
            YORIG3D = 26.D0           ! Set to prevent PAVE crash
            XCELL3D = DBLE( TINTV )   ! Temperature increment
            YCELL3D = 1.D0            ! Set to prevent PAVE crash
            SDATE3D = JYEAR * 1000 + 1
            STIME3D = 1
            TSTEP3D = 1
            NVARS3D = MXNDET
            NCOLS3D = MAXTMP
            NROWS3D = NVTYPE 
            NLAYS3D = 1
            NTHIK3D = 1
            GDTYP3D = LATGRD3           ! Set to prevent PAVE crash
            VGTYP3D = IMISS3
            VGTOP3D = BADVAL3
            GDNAM3D = ' ' 

            FDESC3D( 1 ) = 'Non-diurnal Master Emissions Factors list'
            WRITE  ( MESG, 94070 ) 'for temperatures from', BTEMP, 
     &                             'to', ETEMP
            FDESC3D( 2 ) = MESG

            WRITE  ( MESG, 94070 ) 'on increments of', TINTV
            FDESC3D( 3 ) = MESG

            WRITE  ( MESG, 94070 ) 'Total number of temperatures='
            FDESC3D( 4 ) = MESG

            WRITE  ( MESG, 94080 ) NTEMP
            FDESC3D( 5 ) = MESG

            FDESC3D( 6 ) = ' '

            FDESC3D( 7 ) = VOLTYP( VOLORG )

            DO 55 J = 8, MXDESC3
                FDESC3D( J ) = ' '
55          CONTINUE 

            K = 0
            DO 61 J = 1, MXNDET
                K = K + 1
                VNAME3D( J ) = EMTLST ( EMTIDX( K ) )
                UNITS3D( J ) = EMTUNIT( EMTIDX( K ) )
                VDESC3D( J ) = EMTDESC( EMTIDX( K ) )
                VTYPE3D( J ) = M3REAL
61          CONTINUE
 
            ONAME = PROMPTMFILE(
     &        'Enter logical name for output NON-DIURNAL MASTER EF ' //
     &        'output file', FSNEW3, 'MEFSND', 'EMISFAC' )

        ENDIF  ! File exists already or not

C.........  Set output diurnal name if input file exists
        IF( DNAME .NE. 'NONE' ) THEN  ! Master list exists

            IF( .NOT. DESC3( DNAME ) ) THEN

                WRITE( MESG, 94010 ) 
     &                 'Could not get description of ' // DNAME
                CALL M3EXIT( 'EMISFAC', 0, 0, MESG, 2 )

            ENDIF

C.............  Set output file name to same as input
            NNAME  = DNAME

C.............  Check modeling year to ensure same as ENAME
            I = SDATE3D / 1000
            IF( I .NE. JYEAR ) THEN

                WRITE( MESG,94010 ) 'The year of input file ' //
     &              DNAME( 1:TRIMLEN( DNAME ) ) // ' is', I,
     &              CRLF() // BLANK5 //
     &              '  but the year of input file ' //
     &              ENAME( 1:TRIMLEN( ENAME ) ) // ' is', JYEAR
                CALL M3MSG2( MESG )
                CALL M3EXIT( 'EMISFAC', 0, 0, 
     &                       'Incompatible input files', 2 )

            ENDIF 

C.............  Check volatile organic in this file is same as non-diurnal
            L = TRIMLEN( FDESC3D( 7 ) )
            J = INDEX1 ( FDESC3D( 7 )( 1:L ), NUMVOL, VOLTYP )

            IF( J .NE. VOLORG ) THEN

                WRITE( MESG,94010 ) 
     &             'ERROR: Volatile pollutant type "' //
     &             FDESC3D( 7 )( 1:L ) // '" in file "' //
     &             DNAME( 1:TRIMLEN( DNAME ) ) // 
     &             '" is different from "' // VOLTYP( VOLORG ) //
     &             '" in file "' // ENAME( 1:TRIMLEN( ENAME ) ) //
     &             '".'

                CALL M3EXIT( 'EMISFAC', 0, 0, MESG, 2 )
                  
            ENDIF

C.............  Check that file min/max temperatures is same as compiled value
            MIMICHK = STR2INT( FDESC3D( 2 )( 12:21 ) )
            MIMACHK = STR2INT( FDESC3D( 3 )( 12:21 ) )
            MAMICHK = STR2INT( FDESC3D( 4 )( 12:21 ) )
            MAMACHK = STR2INT( FDESC3D( 5 )( 12:21 ) )

            IF( MIMICHK .NE. MINT_MIN .OR.
     &          MIMACHK .NE. MINT_MAX .OR.
     &          MAMICHK .NE. MAXT_MIN .OR.
     &          MAMACHK .NE. MAXT_MAX      ) THEN

                WRITE( MESG,94010 ) 
     &           '         Program  File' // CRLF() // BLANK5 //
     &           'MINT_MIN', MINT_MIN, ' ', MIMICHK, CRLF() // BLANK5 //
     &           'MINT_MAX', MINT_MAX, ' ', MIMACHK, CRLF() // BLANK5 //
     &           'MAXT_MIN', MAXT_MIN, ' ', MAMICHK, CRLF() // BLANK5 //
     &           'MAXT_MAX', MAXT_MAX, ' ', MAMACHK

                CALL M3MSG2( MESG )

                WRITE( MESG,94010 ) 
     &             'ERROR: File temperature dimensions in file "' //
     &             DNAME( 1:TRIMLEN( DNAME ) ) // 
     &             '" are different from program compiled values.'

                CALL M3EXIT( 'EMISFAC', 0, 0, MESG, 2 )

            ENDIF

        ELSE

C.............  Set up diurnal master List file description

            FTYPE3D = GRDDED3
            P_ALP3D = 0.D0
            P_BET3D = 0.D0
            P_GAM3D = 0.D0
            XCENT3D = 0.D0
            YCENT3D = 0.D0
            XORIG3D = -DBLE( MINT_MIN ) ! Negative starting temperature (- for PAVE)
            YORIG3D = 26.D0             ! Set for PAVE
            XCELL3D = DBLE( TINTV )     ! Temperature increment
            YCELL3D = 1.0               ! Set for PAVE
            SDATE3D = JYEAR * 1000 + 1
            STIME3D = 1
            TSTEP3D = 1
            NVARS3D = MXDET
            NCOLS3D = NMINMAX
            NROWS3D = NVTYPE 
            NLAYS3D = 1
            NTHIK3D = 1
            GDTYP3D = LATGRD3           ! Set for PAVE
            VGTYP3D = IMISS3
            VGTOP3D = BADVAL3
            GDNAM3D = ' ' 

            FDESC3D( 1 ) = 'Diurnal Master Emissions Factors list'
            WRITE( FDESC3D( 2 ),94030 ) 'MINT_MIN = ', MINT_MIN, 'deg F'
            WRITE( FDESC3D( 3 ),94030 ) 'MINT_MAX = ', MINT_MAX, 'deg F'
            WRITE( FDESC3D( 4 ),94030 ) 'MAXT_MIN = ', MAXT_MIN, 'deg F'
            WRITE( FDESC3D( 5 ),94030 ) 'MAXT_MAX = ', MAXT_MAX, 'deg F'
            WRITE( FDESC3D( 6 ),94030 ) 'INCREMENT= ', TINTV, 'deg F'

            FDESC3D( 7 ) = VOLTYP( VOLORG )

            DO 66 J = 8, MXDESC3
                FDESC3D( J ) = ' '
66          CONTINUE 

            K = MXNDET
            DO 77 J = 1, MXDET
                K = K + 1
                VNAME3D( J ) = EMTLST ( EMTIDX( K ) )
                UNITS3D( J ) = EMTUNIT( EMTIDX( K ) )
                VDESC3D( J ) = EMTDESC( EMTIDX( K ) )
                VTYPE3D( J ) = M3REAL
77          CONTINUE
 
            NNAME = PROMPTMFILE(
     &         'Enter logical name for output DIURNAL MASTER EF ' //
     &         'output file', FSNEW3, 'MEFSD', 'EMISFAC' )

        ENDIF

C...................................................................
C.........  Read in and process MPREF  .............................
C...................................................................

        I       = 0
        IREC    = 0
        J       = 0
        JJ      = 0 
        MA      = 0
        MB      = 0
        EFLAG   = .FALSE.
        MOB5AON = .FALSE.
        MOB5BON = .FALSE.
        M5VRSION= ' '
88      CONTINUE   ! Start loop for reading MPREF 

C.............  Read line as character string, which can be the
C.............  first line of 2-line records, or can have a MOBILE5* 
C.............  packet in it.
            READ( NRDEV, 93100, END=99, IOSTAT=IOS ) LINE

            IREC = IREC + 1

            IF( IOS .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 )
     &            'Error ', IOS, 'reading MPREF file at line', IREC
                CALL M3MESG( MESG )

                GO TO 88    ! To head of read loop
            ENDIF

C.............  Check for MOBILE5 packets
            MA = INDEX( LINE, 'MOBILE5A' )
            MB = INDEX( LINE, 'MOBILE5B' )
            L2 = TRIMLEN( LINE )

C.............  Process MPREF entry for MOBILE5 input type
            IF( MA .GT. 0 .OR. MB .GT. 0 ) THEN

C.................  Mobile5a-specific input info
                IF( MA .GT. 0 ) THEN
                    MOB5AON  = .TRUE.
                    M5VRSION = 'MOBILE5a'
                    L1 = MA + 8

                    IF( MOB5BON ) THEN
                        EFLAG = .TRUE.
                        WRITE( MESG,94010 )
     &                    'Line', IREC, 'of MPREF: MOBILE5A ' //
     &                    'packet found, but MOBILE5B already used'
                        CALL M3MESG( MESG )
                    ENDIF

C.................  Mobile5b-specific input info
                ELSE
                    MOB5BON = .TRUE.
                    M5VRSION = 'MOBILE5b'
                    L1 = MB + 8

                    IF( MOB5AON ) THEN
                        EFLAG = .TRUE.
                        WRITE( MESG,94010 )
     &                    'Line', IREC, 'of MPREF: MOBILE5B ' //
     &                    'packet found, but MOBILE5A already used'
                        CALL M3MESG( MESG )
                    ENDIF
                ENDIF

C.................  Set flag for input record type
                MPTYPE = -1

C.................  Store count of PSIs on this line
                READ( LINE( L1:L2 ), *, IOSTAT=IOS ) CNTPSI

                IF( IOS .NE. 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 )
     &                'Bad MOBILE5* packet in MPREF file at line', IREC
                    CALL M3MESG( MESG )
                    GO TO 88    ! To head of read loop

                ELSEIF( CNTPSI .LE. 0 .OR. CNTPSI .GT. MXM5PSI ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG, 94010 )
     &                'Max no. PSIs per MOBILE5 packet (MXM5PSI)=',
     &                MXM5PSI, CRLF() // BLANK5 //
     &                '  exceeded at line', IREC, ', file MPREF'
                    CALL M3MESG( MESG )
                    GO TO 88    ! To head of read loop

                ELSE
                    READ( LINE( L1:L2 ), *, IOSTAT=IOS ) 
     &                    CNTPSI, ( M5PSI( K ), K=1, CNTPSI )

                    DO 91 K = 2, CNTPSI
                        IF( M5PSI( K ) .NE. M5PSI( K-1 ) + 1 ) THEN
                            EFLAG = .TRUE.
                            WRITE( MESG, 94010 )
     &                        'PSIs for MOBILE5 packet are not ' //
     &                        'sequential at line', IREC, ', file MPREF'
                            CALL M3MESG( MESG )
                            GO TO 88    ! To head of read loop
                        ENDIF
91                  CONTINUE

                ENDIF

C.............  If we're still reading pure and combo formats
            ELSEIF( .NOT. MOB5AON .AND. .NOT. MOB5BON ) THEN
                CNTPSI = 1        ! Count of new PSIs is only one

                READ( LINE, *, IOSTAT=IOS ) P1, MPTYPE

                IF( IOS .NE. 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 )
     &                'Error ', IOS, 'reading MPREF file at line', IREC
                    CALL M3MESG( MESG )
                ENDIF

C.............  Otherwise, go to head of read loop because we're skipping
C.............  over the meat of the Mobile5* input files
            ELSE
                GO TO 88 

            ENDIF

C.............  If no overflow for main PSI storage array
            IF( I + CNTPSI .LE. MXPREF ) THEN

C................  For MOBILE5* input type, collect all of the PSI records
                IF( MPTYPE .EQ. -1 ) THEN

                    DO 95 K = 1, CNTPSI

                        I             = I + 1
                        IPSIA ( I )   = M5PSI( K )
                        INDEXA( I )   = I
                        ISTYPE( I,1 ) = MPTYPE
                        ISTYPE( I,2 ) = K
                        ISLINE( I )   = IREC
                        ISROOT( I )   = M5PSI( 1 )
                        ISCNT ( I )   = CNTPSI
    
95                  CONTINUE

C................  For non-combo (normal) PSIs, increment and store properties
                ELSEIF( MPTYPE .EQ. 0 ) THEN
                    I = I + CNTPSI

                    IPSIA ( I )   = P1
                    INDEXA( I )   = I
                    ISTYPE( I,1 ) = MPTYPE
                    ISCNT ( I )   = 1

                    J = J + 1 
                    ISTYPE( I,2 ) = J

                    READ( NRDEV, * ) ( NIDX( J,K ), K=1,NNFAC ),
     &                               ( CIDX( J,K ), K=1,NCFAC )

                    IREC = IREC + 1

C................  For combo PSIs, increment and store properties
                ELSEIF( MPTYPE .GT. 0 ) THEN

                    I = I + CNTPSI

                    IPSIA ( I )   = P1
                    INDEXA( I )   = I
                    ISTYPE( I,1 ) = MPTYPE
                    ISCNT ( I )   = 1

C.....................  Make sure that combination EF doesn't use more than
C                       the legal number of other PSIs
                    IF( MPTYPE .GT. MXNPT ) THEN

                        EFLAG = .TRUE.

                        WRITE( MESG, 94010 ) 
     &                   'Problem at line', IREC, ', file MPREF:' //
     &                   CRLF() // BLANK5 // '  ' //
     &                   'Attmptd # PSIs used in combo        :',MPTYPE,
     &                   CRLF() // BLANK5 // '  ' //
     &                   'Max     # PSIs used in combo (MXNPT):',MXNPT
                        CALL M3MESG( MESG )

C.....................  No overflow, set combo EF properties...
                    ELSE
                        JJ = JJ + 1
                        ISTYPE( I,2 ) = JJ

                        READ( NRDEV, *, IOSTAT=IOS ) 
     &                                   ( PCNTPRO( K,JJ ), 
     &                                     IPLIST ( K,JJ ), K=1,MPTYPE )

                        IREC = IREC + 1

                        IF( IOS .NE. 0 ) THEN
                            EFLAG = .TRUE.
                            WRITE( MESG,94010 )
     &                        'Error ', IOS,
     &                        'reading MPREF file at line', IREC
                            CALL M3MESG( MESG )
 
                        ENDIF

                    ENDIF

                ELSE

                    EFLAG = .TRUE.
                    WRITE( MESG, 94010 ) 
     &                 'PSI type flag illegal value:', MPTYPE,
     &                 CRLF() // BLANK5 //
     &                 '  (profile', P1, ', line', IREC, ', file MPREF)'
                    CALL M3MESG( MESG )

                ENDIF

            ENDIF

        GO TO 88    ! End loop for reading MPREF
99      CONTINUE

        MRNUM  = I

        IF( MRNUM .GT. MXPREF ) THEN

            EFLAG = .TRUE.
            WRITE( MESG, 94010 )
     &       'Records count for master PSIs        :', MRNUM ,
     &       CRLF() // BLANK5 //
     &       'Max dimensioned record-count (MXPREF):', MXPREF,
     &       CRLF() // BLANK5 //
     &       'Max record-count exceeded.'
            CALL M3MSG2( MESG )

        ENDIF

        IF( EFLAG ) THEN
            CALL M3EXIT( 'EMISFAC', 0, 0,
     &                   'Error in MPREF file', 2 )
        END IF

C.........  Sort master PSI list
        CALL SORTI1( MRNUM, INDEXA, IPSIA )

C.........  Create sorted list for searching purposes only
        DO 101 I = 1, MRNUM
           IPSI( I ) = IPSIA( INDEXA( I ) )
101     CONTINUE

C...................................................................
C.........  Read targeted EF update file if it exists ............
C...................................................................

        NUPDAT = 0
        J      = 0
        IF( SDEV .GT. 0 ) THEN
            IREC  = 0
            EFLAG = .FALSE.

111         CONTINUE

                READ( SDEV, 93100, END=117, IOSTAT=IOS ) BUFFER

                IREC = IREC + 1
                L = INDEX( BUFFER, '-' ) 
 
                IF( IOS .NE. 0 ) THEN

                    EFLAG = .TRUE.
                    WRITE( MESG,94010 )
     &                  'Error ', IOS,
     &                  'reading TARGETED EF UPDATE file at line', IREC
                    CALL M3MESG( MESG )
                    GO TO 111           ! To head of read loop

                ELSEIF( L .EQ. 0 ) THEN
                    K  = STR2INT( BUFFER( 1:TRIMLEN( BUFFER ) ) )
                    KK = K 

                ELSE
                    K  = STR2INT( BUFFER( 1  :L-1 ) )
                    KK = STR2INT( BUFFER( L+1:TRIMLEN( BUFFER ) ) )

                ENDIF

                DO 113 JJ = K, KK

                    J = J + 1

                    IF( J .LE. MXPSI ) THEN
                        UPIDXA( J ) = J
                        UPDATA( J ) = JJ
                    ENDIF

113             CONTINUE

            GO TO 111 ! To head of read loop

117         CONTINUE

            NUPDAT = MIN( J, MXPSI )

C.............  Sort list that we have so far
            CALL SORTI1( NUPDAT, UPIDXA, UPDATA )

            DO 121 I = 1, NUPDAT
                UPDATE( I ) = UPDATA( UPIDXA( I ) )
121         CONTINUE

        ENDIF         ! End update input file

C.........  Include PSIs that are first in multi-scenario input
C.........  if PSIs that are NOT first are in update list
        DO 169 I = 1, MRNUM

            MPTYPE = ISTYPE( INDEXA( I ), 1 )
            PSI    = IPSIA ( INDEXA( I )    )
            ROOT   = ISROOT( INDEXA( I )    )

            IF( MPTYPE .EQ. -1 .AND. PSI .GT.  1 ) THEN

C.................  Check if UPDATE list exists and PSI is not in list
                IF( NUPDAT .GT. 0 .AND.
     &              FIND1( PSI , NUPDAT, UPDATE ) .GT. 0 .AND.
     &              FIND1( ROOT, NUPDAT, UPDATE ) .LE. 0       ) THEN

                    J = J + 1

                    IF( J .LE. MXPSI ) THEN
                        UPIDXA( J ) = J
                        UPDATA( J ) = ROOT
                    ENDIF

                ENDIF
            ENDIF

169     CONTINUE

        NUPDAT = J

        IF( NUPDAT .GT. MXPSI ) THEN

            EFLAG = .TRUE.
            WRITE( MESG, 94010 )
     &             'Records count for PSI updates       :', NUPDAT,
     &             CRLF() // BLANK5 //
     &             'Max dimensioned record-count (MXPSI):', MXPSI,
     &             CRLF() // BLANK5 //
     &             'Max record-count exceeded.'
            CALL M3MSG2( MESG )
        ENDIF

        IF( EFLAG ) THEN
            CALL M3EXIT( 'EMISFAC', 0, 0,
     &                   'Error creating update list', 2 )

        ENDIF

        CALL SORTI1( NUPDAT, UPIDXA, UPDATA )

        DO 171 I = 1, NUPDAT
            UPDATE( I ) = UPDATA( UPIDXA( I ) )
171     CONTINUE

C...................................................................
C.........  Read in and process MPLIST file             ............
C...................................................................

C.........  Create SORTED list of UNIQUE parameter scheme indices (PSIs)
C.........  Do this sorting and uniqueness no to avoid a _huge_ array

        I     = 0 
        IREC  = 0 
        EFLAG = .FALSE.
177     CONTINUE   ! Start loop for reading MPLIST

            READ( NLDEV, *, END=188 ) P1, P2, P3, ( P4( JJ ), JJ=1,24 )

            IREC = IREC + 1

C.............  Loop through all entries on line (for all hours)
            DO 181 JJ = 1, 24
 
                IM = LOCAT1( P4( JJ ), I, IPSIR )

C.................  If PSI not found, and not first entry...
                IF( IM .GT. 0 .AND. I .NE. 0 ) THEN
 
                    I = I + 1

                    IF( I .GT. MXPSI ) THEN        ! Overflow will occur
                        WRITE( MESG, 94010 )
     &                     'Maximum number of unqiue PSIs', MXPSI,
     &                     CRLF() // BLANK5 //
     &                     '  exceeded at line', IREC, ', file MPLIST'
                        CALL M3MESG( MESG )
 
                    ELSE

C......................... Shift list by 1, starting from bottom of list
                        DO 179 K = I, IM+1, -1
                            IPSIR( K ) = IPSIR( K-1 )
179                     CONTINUE

C......................... Inert new PSI into sorted list 
                        IPSIR( IM ) = P4( JJ )
                    ENDIF

C.................  First entry...
                ELSEIF( I .EQ. 0 ) THEN
                    I = I + 1
                    IPSIR( I ) = P4( JJ )
 
                ENDIF
 
181         CONTINUE    ! End loop on hours
 
        GO TO 177       ! End loop for reading MPLIST
188     CONTINUE

C.........  Set (current) count of scenario-specific PSIs
        MPNUM = I

        IF( MPNUM .EQ. 0 ) THEN

            CALL M3EXIT( 'EMISFAC', 0, 0, 
     &                   'Empty MPLIST file not allowed', 1 )

        ENDIF

C...................................................................
C.........  Ensure that all PSI on which others depend are in IPSIR
C...................................................................

        DO 199 N = 1, MRNUM

            MPTYPE = ISTYPE( N,1 )

C.............  For combination emission factors...
            IF( MPTYPE .GT. 0 ) THEN

C.................  Set pointer to combo-only lists
                DEPOINT = ISTYPE( N,2 )

C.................  Loop through all factors depended on by Ith PSI
                DO 197 J = 1,MPTYPE

C.....................  Get Jth PSI from list for Ith PSI
                    P1 = IPLIST( J, DEPOINT )

C.....................  Locate Jth PSI in sorted scenario-specific PSI list
                    IM = LOCAT1( P1, MPNUM, IPSIR )

                    IF( IM .GT. 0 ) THEN      ! Jth PSI not in scen-spec list

                        MPNUM = MPNUM + 1

                        IF( MPNUM .GT. MXPSI ) THEN
                            WRITE( MESG, 94010 )
     &                         'Maximum number of unqiue PSIs', MXPSI,
     &                         CRLF() // BLANK5 //
     &                         '  exceeded during list update for ' //
     &                         'combo EFs'
                            CALL M3MESG( MESG )

                        ELSE

C............................. Shift list by 1, starting from bottom of list
                            DO 195 K = MPNUM, IM+1, -1
                                IPSIR( K ) = IPSIR( K-1 )
195                         CONTINUE
    
C............................. Inert new PSI into sorted list 
                            IPSIR( IM ) = P1
                        ENDIF

                    ENDIF

197             CONTINUE

            ENDIF

199     CONTINUE
        
        IF( MPNUM .GT. MXPSI ) THEN

            EFLAG = .TRUE.
            WRITE( MESG, 92010 )
     &             'Records count for unique PSIs       :', MPNUM ,
     &             'Max dimensioned record-count (MXPSI):', MXPSI,
     &             'Max record-count exceeded.'
        ENDIF

        IF( EFLAG ) THEN
            CALL M3EXIT( 'EMISFAC', 0, 0, 
     &                   'Error in MPLIST file', 1 )

        ENDIF

C...................................................................
C.........  Read EF-based file of temp combo indices ...............
C...................................................................

        I     = 0
        IREC  = 0
        EFLAG = .FALSE.
211     CONTINUE   ! Start loop for reading MEFTEMP

            READ( EDEV, *, END=222, IOSTAT=IOS ) 
     &            P1, P2, ( P5( JJ ), JJ = 1, P2 )

            IREC = IREC + 1
            I    = I + 1

            IF( IOS .NE. 0 ) THEN

                EFLAG = .TRUE.
                WRITE( MESG,94010 )
     &              'Error ', IOS,
     &              'reading TEMPERATURE INDEX file at line', I
                CALL M3MESG( MESG )

            ELSEIF( I .GT. MXPSI ) THEN
                WRITE( MESG, 94010 )
     &             'Maximum number PSI-based temperature indices', 
     &             MXPSI, CRLF() // BLANK5 //
     &             '  exceeded at line', IREC, ', file MEFTEMP'
                CALL M3MESG( MESG )

            ELSE    

                MMTREF( I ) = P1
                NTIDX ( I ) = P2
                DO 215 J = 1, P2
                    TIDXPSI( J,I ) = P5( J )
215             CONTINUE

            ENDIF

        GO TO 211

222     CONTINUE
        NTREF = I 

        IF( MPNUM .GT. MXPSI ) THEN

            EFLAG = .TRUE.
            WRITE( MESG, 94010 )
     &       'Records count of temperature indices:', NTREF,
     &       CRLF() // BLANK5 //
     &       'Max dimensioned record-count (MXPSI):', MXPSI,
     &       CRLF() // BLANK5 //
     &       'Max record-count exceeded.'
            CALL M3MSG2( MESG )

        ENDIF

C.........  Check if the number of entries in the temperature-index file
C.........  is the same as the number in MPLIST.  If not - give warning.
        IF( MPNUM .NE. NTREF ) THEN

            WRITE( MESG,94010 ) 
     &                 'The number of entries in MEFTEMP is', NTREF,
     &                 CRLF() // BLANK5 //
     &                 'The number of entries in MPLIST  is', MPNUM
            CALL M3MSG2( MESG )

            WRITE( MESG,94010 ) 
     &        '*****************************************************'//
     &         CRLF() // BLANK5 //
     &        'WARNING: The parameter scheme indices in MEFTEMP are'//
     &         CRLF() // BLANK5 //
     &        '         inconsistent with those in the MPLIST file'//
     &         CRLF() // BLANK5 //
     &        '*****************************************************'
            CALL M3MSG2( MESG )

        ENDIF

        IF ( EFLAG ) THEN
            CALL M3EXIT( 'EMISFAC', 0, 0,
     &                   'Error in MEFTEMP file', 2 )
        END IF

C.........  Update temperature indices lists so pure PSIs include
C.........  indices for the combo PSIs that use them

        EFLAG = .FALSE.
        DO 233 I = 1, NTREF

C.............  Find PSI for temperature index in main PSI list
            K = FIND1 ( MMTREF( I ), MRNUM, IPSI )

            IF( K .LE. 0 ) THEN
                WRITE( MESG, 94010 ) 
     &            'Parameter scheme index', MMTREF( I ),
     &            'from MEFTEMP not found in MPREF'
                CALL M3EXIT('EMISFAC', 0, 0, MESG, 2 )
            ENDIF

            DEPCNT  = ISTYPE( INDEXA( K ),1 )
            DEPOINT = ISTYPE( INDEXA( K ),2 )

            IF( K .LE. 0 ) THEN  ! Error
                EFLAG = .TRUE.
                WRITE( MESG,94010 )
     &               'Error: Temperature index PSI', MMTREF( I ),
     &               'not found in current MPREF file.'
                CALL M3MESG( MESG )

C.............  If PSI is a combination format
            ELSEIF( DEPCNT .GT. 0 ) THEN

C.................  Loop through all PSIs used by PSI( I )
                DO 232 J = 1, DEPCNT

C.....................  Find each PSI used in temperature index
                    L = FIND1( IPLIST( J,DEPOINT ), NTREF, MMTREF )

C.....................  Loop through temp indices for PSI( I ) 
                    DO 231 T = 1, NTIDX( I )

C.........................  Find temp index location in list for PSI(L)
                        KK = LOCAT1( TIDXPSI( T,I ), NTIDX( L ), 
     &                               TIDXPSI( 1,L )              )

C.........................  Initialize temperature index list if needed
                        IF( NTIDX( L ) .EQ. 0 ) THEN

                            NTIDX( L ) = 1
                            TIDXPSI( 1,L ) = TIDXPSI( T,I )

C.........................  Insert temperature index in list for I
                        ELSEIF( KK .GT. 0 ) THEN

C.............................  Check for overflow...
                            IF( NTIDX( L ) .EQ. MXTMMI ) THEN
                                WRITE( MESG, 94010 )
     &                             'Maximum dimension MXTMMI=', MXTMMI,
     &                             'was exceeded'
                                CALL M3EXIT( 'EMISFAC', 0, 0, MESG, 2 )
 
C.............................  If no overflow
                            ELSE
 
C.................................  Add one to count of TMMIs for PSI(L)
                                NTIDX( L ) = NTIDX( L ) + 1
 
C.................................  Shift stored TMMIs down one for PSl(L)
                                DO 229 JJ = NTIDX( L ), KK+1, -1
                                    TIDXPSI( JJ,L ) = TIDXPSI( JJ-1,L )
229                             CONTINUE
 
C.................................  Insert new TMMI for PSI(L) from PSI(I)
                                TIDXPSI( KK,L ) = TIDXPSI( T,I )
 
                            ENDIF  ! End overflow check
                        ENDIF

231                 CONTINUE
232             CONTINUE

            ENDIF 

233     CONTINUE   !  ... Phew!

        IF ( EFLAG ) THEN
            CALL M3EXIT( 'EMISFAC', 0, 0,
     &                   'Bad TEMPERATURE INDEX file', 2 )
        END IF

C.........  Initialize constants for loops
        MAXTL = INT( ( MAXT_MAX - MAXT_MIN ) / TINTV ) + 1

C...................................................................
C.........  Create new pure EFs where needed
C...................................................................

        DO 333 IPM = 1, MPNUM

C.............  Get indice from MPLIST
            IPARM = IPSIR( IPM )
            IPTIM = SEC2TIME( IPARM )

            IF( IPARM .GT. 24*3600 ) THEN
                WRITE( MESG, 94010 ) 'Value of PSIs cannot exceed 86400'
                CALL M3EXIT('EMISFAC', 0, 0, MESG, 2 )
            ENDIF

C.............  Search for MPLIST indice in MPREF index
            J = FIND1( IPARM, MRNUM, IPSI )

            IF( J .LE. 0 ) THEN
                WRITE( MESG, 94010 ) 
     &            'Parameter scheme index', IPARM,
     &            'from MPLIST not found in MPREF'
                CALL M3EXIT('EMISFAC', 0, 0, MESG, 2 )
            ENDIF

C.............  Set count of PSI on which current PSI depends
            DEPCNT  = ISTYPE( INDEXA( J ),1 )
            DEPOINT = ISTYPE( INDEXA( J ),2 )
            CNTPSI  = ISCNT ( INDEXA( J )   )

C.............  Recognize only indice representing pure EF calculations
C.............  or MOBILE5* input files
            IF( DEPCNT .EQ. 0 .OR. DEPCNT .EQ. -1 ) THEN

C.................  Initialize flags to say EFs were not calculated
                TFLAG = .FALSE.  ! Non-diurnal
                RFLAG = .FALSE.  ! Diurnal

C.................  Search for MPLIST indice in PSI list used in temperature index
                K = FIND1( IPARM, NTREF, MMTREF )

                IF( K .LE. 0 ) THEN
                    WRITE( MESG, 94010 ) 
     &                'Parameter scheme index', IPARM,
     &                'from MPLIST not found in MEFTEMP'
                    CALL M3EXIT('EMISFAC', 0, 0, MESG, 2 )
                ENDIF
 
C.................  Initialize MIN/MAX temperature indice counter
                MMTC( K ) = 1

C.................  Write message to screen if non-MOBILE5 packet type
                IF( DEPCNT .EQ. 0 ) THEN
                    WRITE( MESG, 94010 ) 'Processing PSI', IPARM, '...'
                    CALL M3MSG2( MESG )

                ELSEIF( DEPCNT .EQ. -1 .AND. DEPOINT .EQ. 1 ) THEN
                    WRITE( MESG, 94010 ) 
     &                   'Processing PSI', IPARM, 'through', 
     &                   IPARM + ISCNT( INDEXA ( J ) ) - 1, '...'
                    CALL M3MSG2( MESG )

                ENDIF

C.................  Set update flag using search of MPLIST indice in update list
                UFLAG = .FALSE.
                IF ( NUPDAT .GT. 0 .AND.
     &               FIND1( IPARM, NUPDAT, UPDATE ) .GT. 0 ) 
     &               UFLAG = .TRUE.

C.................  Check to see any diurnals for this profile indice exist
                IF( .NOT. CHECK3
     $              ( DNAME, ALLVAR3, SDATE3D, IPTIM ) ) THEN

C.....................  Init diurnal arrays b/c blanks will exist on write3
                    DO 255 I = 1, NVTYPE
                        DO 253 JJ = 1, NMINMAX
                            WDL_VOC( JJ, I ) = BADVAL3
                            DNL_VOC( JJ, I ) = BADVAL3
                            HOT_VOC( JJ, I ) = BADVAL3
                            CRC_VOC( JJ, I ) = BADVAL3
253                     CONTINUE
255                 CONTINUE

C.....................  Diurnal emission factors were "calculated"
                    RFLAG = .TRUE.

C.................  If exists, initialize by reading existing values
                ELSE

                    IF( .NOT. READ3( DNAME, VNAME3D( 1 ), ALLAYS3,
     &                                SDATE3D, IPTIM, WDL_VOC   ) ) THEN
                        WRITE( MESG,94010 ) 
     &                         'Could not read record', IPARM, 
     &                         'of variable WDL_VOC from ' // DNAME
                        CALL M3EXIT( 'EMISFAC', 0, 0, MESG, 2 )
                    ENDIF

                    IF( .NOT. READ3( DNAME, VNAME3D( 2 ), ALLAYS3,
     &                                SDATE3D, IPTIM, DNL_VOC   ) ) THEN
                        WRITE( MESG,94010 ) 
     &                         'Could not read record', IPARM, 
     &                         'of variable DNL_VOC from ' // DNAME
                        CALL M3EXIT( 'EMISFAC', 0, 0, MESG, 2 )
                    ENDIF

                    IF( .NOT. READ3( DNAME, VNAME3D( 3 ), ALLAYS3,
     &                                SDATE3D, IPTIM, HOT_VOC   ) ) THEN
                        WRITE( MESG,94010 ) 
     &                         'Could not read record', IPARM, 
     &                         'of variable HOT_VOC from ' // DNAME
                        CALL M3EXIT( 'EMISFAC', 0, 0, MESG, 2 )
                    ENDIF

                    IF( .NOT. READ3( DNAME, VNAME3D( 4 ), ALLAYS3,
     &                                SDATE3D, IPTIM, CRC_VOC   ) ) THEN
                        WRITE( MESG,94010 ) 
     &                         'Could not read record', IPARM, 
     &                         'of variable CRC_VOC from ' // DNAME
                        CALL M3EXIT( 'EMISFAC', 0, 0, MESG, 2 )
                    ENDIF
                ENDIF

C.................  Specify that more has changed than just temperature
C.................  This should not be set inside non-diurnal and diurnal
C.................  sections individually, because ITEMP = 1 will reset TCNT
C.................  in MOBILE subroutine.
                ITEMP = 1

C.................  Check for Non-diurnal EFs for profile indice
C.................  Also check if an update has been specifically requested
                IF( .NOT. CHECK3 ( ENAME, ALLVAR3, SDATE3D, IPTIM ) .OR.
     &              UFLAG ) THEN

                    FIRSTMPR = .TRUE. ! Specify first temperature for MOINBLD

C.....................  Loop over temperature range for non-diurnal factors
                    DO 266 T = 1, NTEMP

                        TMMI   = TIDXPSI( MMTC( K ),K )

C.........................  Get temperature index for PSI(K) and count(K)
C.........................  Check if count has not yet reached NTIDX (the number
C.........................  for this PSI), and diurnal emissions are not zero
                        IF( MMTC( K ) .LE. NTIDX( K ) .AND. 
     &                      TMMI      .NE. 0                ) THEN

                            TMMIDV = INT( TMMI / MAXTL ) 
                            TF_MIN = MINT_MIN + TINTV * TMMIDV
                            TF_MAX = MAXT_MIN + TINTV *
     &                                       ( TMMI - MAXTL*TMMIDV - 1)
                        ELSE
                            TF_MIN = 0.
                            TF_MAX = 0.

                        ENDIF

C..........................  Initialize flag to _not_ calculate diurnal factors
                        DFLAG = .FALSE.

C.........................  Increase counter for this PSI if ambient temp is 
C.........................  evenly between min/max, otherwise set min/max 
                        IF( TEMPIN( T ) .EQ. 
     &                      TF_MIN + ( TF_MAX - TF_MIN ) / 2. ) THEN

                            MMTC( K ) = MMTC( K ) + 1  ! Set for next time

C.............................  Diurnal emission factors to be calculated 
                            DFLAG = .TRUE. 

C.........................  Set min/max so 2*TINTV less/greater than ambient
                        ELSE

C.............................  Temporarily round current temp to nearest in 
C.............................  allowed DIURNAL set
                            TF = MINT_MIN + TINTV *
     &                           NINT( ( TEMPIN(T) - MINT_MIN )/TINTV ) 

C.............................  Set new min and max based on rounded temp
                            TF_MIN = MAX( TF - 2.*TINTV, MINT_MIN )
                            TF_MIN = MIN( TF_MIN, MINT_MAX )

                            TF_MAX = MIN( TF + 2.*TINTV, MAXT_MAX )

C............................. Set new temperature index for new min/max
                            TMMI = MAXTL *
     &                           NINT( ( TF_MIN-MINT_MIN ) / TINTV ) +
     &                           NINT( ( TF_MAX-MAXT_MIN ) / TINTV ) + 1

                        ENDIF

C.........................  Set ambient temperature
                        TF = TEMPIN( T )

C.........................  Initialize error count
                        IERR = 0 

C.......................  For pure EFs from tables, create input for MOBILE5.
C.......................  Also, turn off first temperature after first call. 
C.......................  Finally, create emission factors.
                        IF( DEPCNT .EQ. 0 ) THEN

                            REWIND( MODEV )                  
                            CALL MOINBLD( MODEV,IPARM ,DEPOINT, NIDX,
     &                                    CIDX ,VOLORG,JYEAR, FIRSTMPR )

                            FIRSTMPR = .FALSE.

                            REWIND( MODEV )                  
                            CALL MOBILE( IERR , ITEMP , MODEV,
     &                                   CNTPSI, EFACT, DFACT  )

C.......................  For MOBILE5* input when the PSI is the first in a
C.......................  multi-scenario run, set input file to MPREF file, and
C.......................  skip the appropriate number of lines.  Finally, create
C.......................  emission factors and store extras.
                        ELSEIF( DEPCNT  .EQ. -1 .AND. 
     &                          DEPOINT .EQ.  1       ) THEN
                            
                            REWIND( NRDEV )                  
                            CALL SKIPL ( NRDEV, ISLINE( INDEXA( J ) ) )

                            CALL MOBILE( IERR , ITEMP , NRDEV,
     &                                   CNTPSI, EFACT, DFACT  )

C.......................  For MOBILE5* input when the PSI is NOT the first in a
C.......................  multi-scenario run, retrieve EFs from storage
                        ELSEIF( DEPCNT  .EQ. -1 .AND. 
     &                          DEPOINT .GT.  1       ) THEN

                            DO 259 I = 1, NVTYPE
                                EFACT( 1,I )= EFSAVE( I,DEPOINT,T,1  )
                                EFACT( 2,I )= EFSAVE( I,DEPOINT,T,2  )
                                EFACT( 3,I )= EFSAVE( I,DEPOINT,T,3  )
                                EFACT( 4,I )= EFSAVE( I,DEPOINT,T,4  )
                                EFACT( 5,I )= EFSAVE( I,DEPOINT,T,5  )
                                EFACT( 6,I )= EFSAVE( I,DEPOINT,T,6  )
                                EFACT( 7,I )= EFSAVE( I,DEPOINT,T,7  )
                                DFACT( 1,I )= EFSAVE( I,DEPOINT,T,8  )
                                DFACT( 2,I )= EFSAVE( I,DEPOINT,T,9  )
                                DFACT( 3,I )= EFSAVE( I,DEPOINT,T,10 )
                                DFACT( 4,I )= EFSAVE( I,DEPOINT,T,11 )
259                         CONTINUE

                        ENDIF

C.........................  Loop through all vehicle types
                        DO 261 I = 1, NVTYPE

C.............................  Put EFACT into non-diurnal arrays
                            EXH_CO ( T, I ) = EFACT( 1, I )
                            EXH_NOX( T, I ) = EFACT( 2, I )
                            EXH_VOC( T, I ) = EFACT( 3, I )
                            EVP_VOC( T, I ) = EFACT( 4, I )
                            RNL_VOC( T, I ) = EFACT( 5, I )
                            RST_VOC( T, I ) = EFACT( 6, I )
                            RFL_VOC( T, I ) = EFACT( 7, I )

C.............................  Put DFACT into diurnal arrays and flip switch 
C.............................  that diurnal emission factors were calculated
                            IF( DFLAG ) THEN
                                WDL_VOC( TMMI, I ) = DFACT( 1, I )
                                DNL_VOC( TMMI, I ) = DFACT( 2, I )
                                HOT_VOC( TMMI, I ) = DFACT( 3, I )
                                CRC_VOC( TMMI, I ) = DFACT( 4, I )

                                RFLAG = .TRUE. 
                            ENDIF

261                     CONTINUE  ! Vehicle types

C.........................  Non-diurnal emission factors were calculated 
                        TFLAG = .TRUE. 

266                 CONTINUE  ! T, temperature for non-diurnal factors

C.....................  Add non-diurnal EFs to master list 
                    IF( .NOT. WRITE3( ONAME, ALLVAR3, SDATE3D,
     &                                IPTIM, EXH_CO            )) THEN

                        WRITE( MESG, 94010 ) 
     &                    'Could not write non-diurnal EF set', IPARM,
     &                    'to ' // ONAME( 1:TRIMLEN( ONAME ) )
                        CALL M3EXIT('EMISFAC', 0, 0, MESG, 2 )

                    ENDIF  ! Write to non-diurnal master list

                ELSE

                    WRITE( MESG, 94010 )
     &                'Non-diurnal emission factors for set', IPARM,
     &                CRLF() // BLANK5 //
     &                '  already exist in ' // ENAME( 1:TRIMLEN(ENAME) )
                    CALL M3MESG( MESG )

                ENDIF      ! Non-Diurnal EFs have been created or not

                VFLAG( K ) = .TRUE.  ! Initialize - default: update non-diu EFs
                SFLAG( K ) = .TRUE.  ! Initialize - default: update diurnal EFs
                FIRSTMPR   = .TRUE.  ! Specific first temperature for MOINBLD

C.................  Create remaining diurnal factors for IPARM by
C.................  looping through remaining temperatures
                DO 288 T = MMTC( K ), NTIDX( K )

                    DFLAG = .FALSE.  ! Update for each new T combo

C.....................  Get temperature index for PSI(K) and count(K)
                    TMMI   = TIDXPSI( T,K )

C.....................  Check if diurnal emissions for this PSI are zero
                    IF( TMMI .NE. 0 ) THEN
                        TMMIDV = INT( TMMI / MAXTL ) 

C.........................  Specify min and max temperatures
                        TF_MIN = MINT_MIN + TINTV * TMMIDV

                        TF_MAX = MAXT_MIN + TINTV *
     &                               ( TMMI - MAXTL*TMMIDV - 1)

C.........................  Check if EFs for this TMMI is missing or if
C.........................  targeted EF update is required
                        IF( WDL_VOC( TMMI,1 ) .LT. 0 .OR. UFLAG ) 
     &                      DFLAG = .TRUE.

C.........................  If diurnal EF has not yet been calculated for this
C.........................  min/max temperature combination
                        IF( DFLAG ) THEN

                          TF = TF_MIN + ( TF_MAX - TF_MIN ) / 2.0

C...........................  Initialize error count
                          IERR = 0 

C.........................  For pure EFs from tables, create input for MOBILE5.
C.........................  Also, turn off first temperature after first call. 
C.........................  Finally, create emission factors.
                          IF( DEPCNT .EQ. 0 ) THEN

                            REWIND( MODEV )                  
                            CALL MOINBLD( MODEV,IPARM ,DEPOINT, NIDX,
     &                                    CIDX ,VOLORG,JYEAR, FIRSTMPR )

                            FIRSTMPR = .FALSE.

                            REWIND( MODEV )                  
                            CALL MOBILE( IERR , ITEMP , MODEV,
     &                                   CNTPSI, EFACT, DFACT  )

C.........................  For MOBILE5* input when the PSI is the first in a
C.........................  multi-scenario run, set input file to MPREF file, and
C.........................  skip the appropriate number of lines.  Finally, create
C.........................  emission factors and store extras.
                          ELSEIF( DEPCNT  .EQ. -1 .AND. 
     &                            DEPOINT .EQ.  1       ) THEN

                            REWIND( NRDEV )                  
                            CALL SKIPL ( NRDEV, ISLINE( INDEXA( J ) ) )

                            CALL MOBILE( IERR , ITEMP , NRDEV,
     &                                   CNTPSI, EFACT, DFACT  )

C.........................  For MOBILE5* input when the PSI is NOT the first in a
C.........................  multi-scenario run, retrieve EFs from storage
                          ELSEIF( DEPCNT  .EQ. -1 .AND. 
     &                            DEPOINT .GT.  1       ) THEN

                            DO 271 I = 1, NVTYPE
                                DFACT( 1,I )= EFSAVE( I,DEPOINT,T,8  )
                                DFACT( 2,I )= EFSAVE( I,DEPOINT,T,9  )
                                DFACT( 3,I )= EFSAVE( I,DEPOINT,T,10 )
                                DFACT( 4,I )= EFSAVE( I,DEPOINT,T,11 )
271                         CONTINUE

                          ENDIF

C...........................  Put DFACT into diurnal arrays
                          DO 277 I = 1, NVTYPE

                            WDL_VOC( TMMI, I ) = DFACT( 1, I )
                            DNL_VOC( TMMI, I ) = DFACT( 2, I )
                            HOT_VOC( TMMI, I ) = DFACT( 3, I )
                            CRC_VOC( TMMI, I ) = DFACT( 4, I )

277                       CONTINUE
 
C...........................  Diurnal emission factors were calculated 
                          RFLAG = .TRUE. 

                        ENDIF  ! If need to call MOBILE or not
                    ENDIF      ! If valid min/max temp combo or not
 
288             CONTINUE       ! Loop through temperature indices for IPARM

C.................  If factors were not calculated, then don't output them
                IF( .NOT. TFLAG ) VFLAG( K ) = .FALSE.  ! Non-diurnal
                IF( .NOT. RFLAG ) SFLAG( K ) = .FALSE.  ! Diurnal

C.................  Write out diurnal emission factors for IPARM
C.................  NOTE: Won't output for UPDATE if there are NO diurnal
C.................        temperature variations.
                IF( SFLAG( K ) ) THEN
                    IF( .NOT. WRITE3( NNAME, ALLVAR3, SDATE3D,
     &                                IPTIM, WDL_VOC )         ) THEN

                        WRITE( MESG, 94010 ) 
     &                    'Could not write diurnal EF set', IPARM,
     &                    'to ' // NNAME( 1:TRIMLEN( NNAME ) )
                        CALL M3EXIT('EMISFAC', 0, 0, MESG, 2 )

                    ENDIF  ! Write to diurnal master list
                ELSEIF( DNAME .NE. 'NONE' ) THEN
                    WRITE( MESG, 94010 )
     &                'Diurnal     emission factors for set', IPARM,
     &                CRLF() // BLANK5 //
     &                '  already exist in ' // DNAME( 1:TRIMLEN(DNAME) )
                    CALL M3MESG( MESG )

                ELSE
                    WRITE( MESG, 94010 )
     &                'Error: Check MEFTEMP file for zero TMMI ' //
     &                'values for EFs set', IPARM
                    CALL M3EXIT('EMISFAC', 0, 0, MESG, 2 )

                ENDIF      ! update of diurnal file is needed or not

            ENDIF          ! EFs are pure or combo

333     CONTINUE           ! IPM

C...................................................................
C.........  Create new mixed EFs where needed
C...................................................................

        DO 401 IPM = 1, MPNUM

C.............  Get indice from MPLIST
            IPARM = IPSIR( IPM )
            IPTIM = SEC2TIME( IPARM )

C.............  Search for MPLIST indice in MPREF index
            J      = FIND1 ( IPARM, MRNUM, IPSI )
            DEPCNT = ISTYPE( INDEXA( J ),1 )
            DEPOINT= ISTYPE( INDEXA( J ),2 )

C.............  Recognize only indice representing combo EF calculations
            IF( DEPCNT .GT. 0 ) THEN

C.................  Write message to screen
                WRITE( MESG, 94010 ) 'Processing PSI', IPARM, '...'
                CALL M3MSG2( MESG )

C.................  Set update flag using search of MPLIST indice in update list
                UFLAG = .FALSE.
                IF ( NUPDAT .GT. 0 .AND.
     &               FIND1( IPARM, NUPDAT, UPDATE ) .GT. 0 ) 
     &               UFLAG = .TRUE.

C.................  Update diurnals or non-diurnal EFs if any related 
C.................  pure EFs were updated
                DO 381 I = 1, DEPCNT         ! Loop has early exit if needed

                    II    = FIND1( IPLIST( I,DEPOINT ), NTREF, MMTREF )
                    TFLAG = VFLAG( II )  ! Non-diurnal
                    RFLAG = SFLAG( II )  ! Diurnal

C.....................  If either non-diurnal or diurnal are required,
C.....................  update both in efcombo, because running it is fast
                    IF( UFLAG .OR. TFLAG .OR. RFLAG ) THEN

C.........................  Create combo profile through common blocks 
C.........................  EFS1 and EFS2
                        CALL EFCOMBO( ONAME, NNAME, NTEMP, 
     &                                DEPCNT, PCNTPRO( 1, DEPOINT ),
     &                                IPLIST ( 1, DEPOINT ), SDATE3D )

                        GO TO 383  ! Exit out of loop

                    ENDIF

381             CONTINUE
383             CONTINUE

C.................  Store emission factors in master lists
                IF( TFLAG ) THEN
                    IF( .NOT. WRITE3( ONAME, ALLVAR3, SDATE3D,
     &                                IPTIM, EXH_CO )         ) THEN

                        WRITE( MESG, 94010 ) 
     &                        'Could not write non-diurnal EF set',
     &                        IPARM, 'to ' // ONAME( 1:TRIMLEN(ONAME) )
                        CALL M3EXIT('EMISFAC', 0, 0, MESG, 2 )

                    ENDIF  ! Write to non-diurnal master list
                ELSE

                    WRITE( MESG, 94010 )
     &                'Non-diurnal emission factors for set', IPARM,
     &                CRLF() // BLANK5 //
     &                '  already exist in ' // ENAME( 1:TRIMLEN(ENAME) )
                    CALL M3MESG( MESG )

                ENDIF      ! Missing non-diurnal for this IPARM or not

                IF( RFLAG ) THEN
                    IF( .NOT. WRITE3( NNAME, ALLVAR3, SDATE3D,
     &                                IPTIM, WDL_VOC )        ) THEN

                        WRITE( MESG, 94010 )  
     &                    'Could not write diurnal EF set',
     &                    IPARM, 'to ' // NNAME( 1:TRIMLEN(NNAME) )
                        CALL M3EXIT('EMISFAC', 0, 0, MESG, 2 )

                    ENDIF  ! Write to diurnal master list

                ELSE

                    WRITE( MESG, 94010 )
     &                'Diurnal     emission factors for set', IPARM,
     &                CRLF() // BLANK5 //
     &                '  already exist in ' // DNAME( 1:TRIMLEN(DNAME) )
                    CALL M3MESG( MESG )

                ENDIF      ! Missing diurnal for this IPARM or not

            ENDIF          ! EFs are combo or pure

401     CONTINUE           ! End loop on parameter schemes

        MESG = 'Normal completion of EMISFAC'
        CALL M3EXIT( 'EMISFAC', 0, 0, MESG, 0 )
 
C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )
92010   FORMAT( 5X, A, I5 / 20( 5X, A / ) )
92020   FORMAT( 20( 5X, A, 1X, I5 ) )

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( 12( 3( 9(F6.3, 1X)/ ) // ) )
93100   FORMAT( A )
93200   FORMAT( 3I6, 2X, 5E12.4 )
93210   FORMAT( 18X, 2X, 5E12.4 )
c temp
93999   FORMAT( F6.2,1X, 8(F6.2,1X) )

C...........   Internal buffering formats............ 94xxx
 
94010   FORMAT( 20 ( A, :, I5, :, 2X ) )
94020   FORMAT( 'Env var PROMPTFLAG is N, but ', A, ' is bad' )
94030   FORMAT( A, F10.3, 1X, A )
94060   FORMAT( 8( F6.2, ', ' ) )
94070   FORMAT( 5( A, 1X, F6.1, 1X, A, 1X, F6.1, 1X ) )
94080   FORMAT( I8 )

        END

