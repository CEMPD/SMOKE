C Version "@(#)$Id$ $Source$ $Date$ 

        PROGRAM EMISFAC

C***********************************************************************
C  program body starts at line
C
C  DESCRIPTION:
C      This program drives the MOBILE5A program by supplying an ambient 
C      temperature and an input file name from which mobile5a derives 
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

        IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'CHDIMS3.EXT'    !  Mobile include paramaters
        INCLUDE 'MBDIMS3.EXT'    !  Mobile include paramaters
        INCLUDE 'PARMS3.EXT'     !  I/O API parameters
        INCLUDE 'IODECL3.EXT'    !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'     !  I/O API file description data structures.
        INCLUDE 'NETCDF.EXT'     !  NetCDF declarations

C...........   EXTERNAL FUNCTIONS and their descriptions:

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

        EXTERNAL     FIND1, GETEFILE, GETNUM, GETREAL, GETYN, 
     &               INDEX1, LOCAT1, PROMPTFFILE, PROMPTMFILE, 
     &               SEC2TIME, STR2INT, TRIMLEN

C...........   EXTERNAL BLOCK DATA for mobile:

        EXTERNAL      BD01, BD02, BD03, BD04, BD05, BD06, BD10,
     &                BD11, BD12, BD13, BD14, BD15, BD16, BD17,
     &                BD18, BD19, BD20, BD21, BD22, BD23, BD24,
     &                BD25, BD26, BD27, BD28, BD29, BD30, BD31,
     &                BD32, BD33, BD34, BD35, BD36, BD37, BD38

C...........   LOCAL VARIABLES and their descriptions:

C...........   Simulation-specific MOBILE5a parmeter list variables
        INTEGER IPSIDX( MXPSI )      ! Sorting index for IPSIRA
        INTEGER IPSIR ( MXPSI )      ! Sorted   scenario-specific PSI list
        INTEGER MPNUM                ! Actual count of scenario-specific PSIs
        INTEGER P1, P2, P3, P4( 24 ) ! Dummies

C...........   MOBILE5a parameter x-ref variables
        INTEGER MRNUM                     ! Actual no. of normal x-ref entries
        INTEGER INDEXA( MXPREF )          ! Index for sorting x-ref scheme index
        INTEGER IPSI  ( MXPREF )          ! Sorted   x-ref scheme index
        INTEGER IPSIA ( MXPREF )          ! Unsorted x-ref scheme index
        INTEGER ISTYPE( MXPREF, 2     )   ! Type-specifier (1), and pointer (2)
        INTEGER NIDX  ( MXPREF, NNFAC )   ! Normal x-ref entries
        INTEGER CIDX  ( MXPREF, NCFAC )   ! Control Normal x-ref entries

C...........   Master Non-Diurnal EFs list variables
        REAL    EXH_CO ( MAXTMP, NVTYPE )  ! Non-Diurn Emis fac: CO exhaust
        REAL    EXH_NOX( MAXTMP, NVTYPE )  ! Non-Diurn Emis fac: NOX exhaust
        REAL    EXH_VOC( MAXTMP, NVTYPE )  ! Non-Diurn Emis fac: VOC exhaust
        REAL    EVP_VOC( MAXTMP, NVTYPE )  ! Non-Diurn Emis fac: VOC evaporative
        REAL    RNL_VOC( MAXTMP, NVTYPE )  ! Non-Diurn Emis fac: VOC running
        REAL    RST_VOC( MAXTMP, NVTYPE )  ! Non-Diurn Emis fac: VOC resting

        COMMON / EFS1 / EXH_CO , EXH_NOX, EXH_VOC, EVP_VOC, RNL_VOC, 
     &                  RST_VOC

C...........   Master Diurnal EFs list variables

        REAL    WDL_VOC( NMINMAX, NVTYPE )  ! Diurnal Emis fac: VOC weighted
        REAL    DNL_VOC( NMINMAX, NVTYPE )  ! Diurnal Emis fac: VOC diurnal
        REAL    HOT_VOC( NMINMAX, NVTYPE )  ! Diurnal Emis fac: VOC hot soak
        REAL    CRC_VOC( NMINMAX, NVTYPE )  ! Diurnal Emis fac: VOC crank case

        COMMON / EFS2 / WDL_VOC, DNL_VOC, HOT_VOC, CRC_VOC

C...........   Variables for targeted update of EFs
        INTEGER NUPDAT                      ! No. of PSIs targeted for update
        INTEGER UPDATE ( MXPSI )            ! Targeted PSI list

C...........   Variables to carry over from one MOBILE5 call to the next

        INTEGER      IREAD
        REAL         CR12HC
        REAL         CR12CO
        CHARACTER*40 TNAME
        COMMON /IMPAR7/ TNAME,IREAD
        COMMON /IM12HC/ CR12HC(19,20,5,2)
        COMMON /IM12CO/ CR12CO(19,20,5,2)

C...........   Variables for storing PSI descriptors:

        REAL    PCNTPRO( MXNPT, MXPREF )
        INTEGER IPLIST ( MXNPT, MXPREF )
        INTEGER MPTYPE, MCOMBO, MPURE

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
        REAL    TEMPIN( MAXTMP )          ! Input temperatures in deg F

        REAL    BTEMP                     ! Beginning temp for non-diurnal
        REAL    ETEMP                     ! Ending temp for non-diurnal
        REAL    MAXTL                     ! intermediate temperature index value
        REAL    TF                        ! temporary temperature
        REAL    TF_MAX                    ! temporary maximum temperature
        REAL    TF_MIN                    ! temporary minimum temperature

        INTEGER DEPCNT                    ! Temporary count for combo emis facs
        INTEGER DEPOINT                   ! Temporary pointer to combo-only lists
        INTEGER EMTIDX( MXEMIS )          ! Index to EMTLST for actual EFs
        INTEGER I, II, IM, IPM            ! Counters and pointers
        INTEGER J, JJ, K, KK, L, N, T     ! Counters and pointers
        INTEGER IERR                      ! MOBILE5 error value
        INTEGER IOS                       ! I/O status
        INTEGER IPARM                     ! temporary PSI value 
        INTEGER IPTIM                     ! IPARM converted to time for file ops
        INTEGER IREC                      ! Line number when reading ASCII files
        INTEGER ITEMP                     ! Indicator of tmpr status for MOBILE5 
        INTEGER JYEAR                     ! year to calculate emission factors
        INTEGER NTEMP                     ! Actual no. of temps for eval of EFs
        INTEGER P5( MXTMMI )              ! Dummy for reading EF-based TMMI vals
        INTEGER TMMI                      ! temporary min/max temp index
        INTEGER TMMIDV                    ! temporary index divided by MAXTL
        INTEGER VOLORG                    ! Index of Mobile5 hydrocarbon type

        LOGICAL DFLAG                     ! Calculate diurnal factors when true
        LOGICAL EFLAG                     ! Error exists when true
        LOGICAL RFLAG                     ! Diurnal facs _were_ calc'd when true
        LOGICAL SFLAG( MXPSI )            ! Update diurnal combos when true
        LOGICAL TFLAG                     ! Non-diu facs _were_ calc'd when true
        LOGICAL VFLAG( MXPSI )            ! Update non-diu combos when true
        LOGICAL UFLAG                     ! Update EFs when true

        CHARACTER*4   VOLNAM              !  name of Mobile5 hydrocarbon
        CHARACTER*30  BUFFER              !  input buffer
        CHARACTER*256 MESG                !  message buffer

C***********************************************************************
C   begin body of program EMISFAC

        LDEV = INIT3()
 
        WRITE( *,92000 )
     &  ' ',
     &  'Program EMISFAC drives the MOBILE5a program by supplying',
     &  'a range of ambient temperatures and a scenario-specific',
     &  'MOBILE5 parameter file MPLIST. Using multiple calls',
     &  'to MOBILE5a, EMISFAC creates a master emission factors',
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
            WRITE( *,92000 ) ' ', 'Ending program EMISFAC', ' '
            CALL EXIT( 0 )
        END IF

C.........  Turn off verbose NetCDF to elimiate error messages when using
C.........  CHECK3 to check for available time steps
        CALL NCPOPT( 0 )

C.........  Get MOBILE5a input file (no prompt)
        MODEV = GETEFILE( 'MO5AIN', .FALSE., .TRUE., 'EMISFAC' )

C.........  Prompt for diurnal master list EFs input file for modeling year
        ENAME = PROMPTMFILE(
     &          'Enter logical name for NON-DIURNAL EF input '//
     &          'file or "NONE"', FSRDWR3, 'MEFSND', 'EMISFAC' )

C.........  Prompt for diurnal master list EFs input file for modeling year
        DNAME = PROMPTMFILE(
     &          'Enter logical name for DIURNAL EF input ' //
     &          'file or "NONE"', FSRDWR3, 'MEFSD', 'EMISFAC' )

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
        
C.........  Prompt for MPREF (general MOBILE5a parameter x-ref file)
        NRDEV = PROMPTFFILE( 
     &          'Enter logical name for MPREF input file',
     &          .TRUE., .TRUE., 'MPREF', 'EMISFAC' )

C.........  Write out list of legal hydrocarbon species
        WRITE( *,92000 ) ' '
        DO 22 I = 1, NUMVOL
            CALL M3MESG( '     ' // VOLTYP( I ) )
   22   CONTINUE 

C.........  Prompt for type of organic compound to use (numbers match
C.........  MOBILE5 numbers

        DO 23 I = 5, 1, -1
            CALL GETSTR( 'Enter type of hydrocarbon emissions', 
     &                   VOLTYP( 3 ), VOLNAM )
 
            VOLORG = INDEX1( VOLNAM, NUMVOL, VOLTYP )

            IF( VOLORG .GT. 0 ) GO TO 25  ! Valid reponse found

            WRITE( MESG,94010 ) 'Invalid response, get', I-1, 
     &                          'more tries'
            CALL M3MESG( MESG )
 
   23   CONTINUE

        WRITE( MESG,94010 ) 'Bad hydrocarbon emissions type'
        CALL M3EXIT( 'EMISFAC', 0, 0, MESG, 2 )

   25   CONTINUE

C.........  Find VOLNAM in EMTLST for setting output variable names

        EMTIDX( 1 ) = 1   ! For CO
        EMTIDX( 2 ) = 2   ! For NOX
        K = 2

        DO 26 I = 3, MXEALL

            L = INDEX( EMTLST( I ), '_' )
            IF( EMTLST( I )(L+1:7 ) .EQ. VOLNAM ) THEN
                K = K + 1
                EMTIDX( K ) = I
            ENDIF

   26   CONTINUE

C.........  Set start and ending temperature depending on input file
C.........  If it exists, use the parameters from the file.
        IF( ENAME .NE. 'NONE' ) THEN

C.............  Get file description
            IF( .NOT. DESC3( ENAME ) ) THEN

                WRITE( MESG, 94010 ) 
     &                 'Could not get description of ' // ENAME
                CALL M3EXIT( 'EMISFAC', 0, 0, MESG, 2 )

            ENDIF

C.............  Set output file name to same as input
            ONAME = ENAME

C.............  Set modeling year
            JYEAR = SDATE3D / 1000

            WRITE( MESG,94010 ) 
     &          'Emission factors are being calculated for', JYEAR
            CALL M3MESG( MESG )

C.............  Check that file delta temperature is same as compiled value
            IF( XCELL3D .NE. TINTV ) THEN
                WRITE( MESG, 94030 ) 
     &                 'Delta temperature for file ' // 
     &                 ENAME( 1:TRIMLEN( ENAME ) ) // 
     &                 'not equal to compiled value of', TINTV
                CALL M3EXIT( 'EMISFAC', 0, 0, MESG, 2 )

            ENDIF

            NTEMP = STR2INT( FDESC3D( 5 ) )

            BTEMP = XORIG3D
            ETEMP = BTEMP + TINTV * FLOAT( NTEMP - 1 ) 

            DO 28 I = 1, NTEMP
                TEMPIN( I ) = BTEMP + ( I - 1 ) * TINTV
   28       CONTINUE

C.........  If it does not exist, build the parameters from scratch
        ELSE 

C.............  Prompt for modeling year
            JYEAR = 1988
            JYEAR = GETNUM( 0, 9999, JYEAR,
     &             'Enter modeling year for emission factors (YYYY)' )

C.............  Prompt for start temperature
            BTEMP = GETREAL( 1.,110., 1.,
     &                   'Enter minimum temperature [deg F]' )

C.............   Prompt for end temperature
            ETEMP = GETREAL( BTEMP,110., 110., 
     &                   'Enter maximum temperature [deg F]' )

C.............  Set temperatures
            I = 0
   31       CONTINUE

                TF = BTEMP + TINTV * I

                IF( TF .GT. ETEMP ) GO TO 32

                I = I + 1
                TEMPIN( I ) = TF

            GO TO 31

   32       CONTINUE
            NTEMP = I

C.............  Set up non-diurnal master list file description

            FTYPE3D = GRDDED3
            P_ALP3D = DBLE( IMISS3 )
            P_BET3D = DBLE( BADVAL3 )
            P_GAM3D = DBLE( BADVAL3 )
            XCENT3D = 0.0D0
            YCENT3D = 0.0D0
            XORIG3D = DBLE( BTEMP )     ! Starting temperature
            YORIG3D = DBLE( BADVAL3 )
            XCELL3D = DBLE( TINTV )     ! Temperature increment
            YCELL3D = DBLE( BADVAL3 )
            SDATE3D = JYEAR * 1000 + 1
            STIME3D = 1
            TSTEP3D = 1
            NVARS3D = MXNDET
            NCOLS3D = MAXTMP
            NROWS3D = NVTYPE 
            NLAYS3D = 1
            NTHIK3D = 1
            GDTYP3D = IMISS3
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

            K = 0
            DO 33 J = 1, MXNDET
                K = K + 1
                VNAME3D( J ) = EMTLST ( EMTIDX( K ) )
                UNITS3D( J ) = EMTUNIT( EMTIDX( K ) )
                VDESC3D( J ) = EMTDESC( EMTIDX( K ) )
                VTYPE3D( J ) = M3REAL
   33       CONTINUE
 
            ONAME = PROMPTMFILE(
     &        'Enter logical name for output NON-DIURNAL MASTER EF ' //
     &        'output file', FSNEW3, 'MEFSND', 'EMISFAC' )

        ENDIF

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

                WRITE( *,92020 ) 'The year of input file ' //
     &              DNAME( 1:TRIMLEN( DNAME ) ) // ' is', I,
     &              'but the year of input file ' //
     &              ENAME( 1:TRIMLEN( ENAME ) ) // ' is', JYEAR
                CALL M3EXIT( 'EMISFAC', 0, 0, 
     &                       'Incompatible input files', 2 )

            ENDIF 

        ELSE

C.............  Set up diurnal master List file description

            FTYPE3D = GRDDED3
            P_ALP3D = DBLE( IMISS3 )
            P_BET3D = DBLE( BADVAL3 )
            P_GAM3D = DBLE( BADVAL3 )
            XCENT3D = 0.0D0
            YCENT3D = 0.0D0
            XORIG3D = DBLE( MINT_MIN )     ! Starting temperature
            YORIG3D = DBLE( BADVAL3 )
            XCELL3D = DBLE( TINTV )     ! Temperature increment
            YCELL3D = DBLE( BADVAL3 )
            SDATE3D = JYEAR * 1000 + 1
            STIME3D = 1
            TSTEP3D = 1
            NVARS3D = MXDET
            NCOLS3D = NMINMAX
            NROWS3D = NVTYPE 
            NLAYS3D = 1
            NTHIK3D = 1
            GDTYP3D = IMISS3
            VGTYP3D = IMISS3
            VGTOP3D = BADVAL3
            GDNAM3D = ' ' 

            FDESC3D( 1 ) = 'Diurnal Master Emissions Factors list'
            WRITE( FDESC3D( 2 ),94030 ) 'MINT_MIN = ', MINT_MIN, 'deg F'
            WRITE( FDESC3D( 3 ),94030 ) 'MINT_MAX = ', MINT_MAX, 'deg F'
            WRITE( FDESC3D( 4 ),94030 ) 'MAXT_MIN = ', MAXT_MIN, 'deg F'
            WRITE( FDESC3D( 5 ),94030 ) 'MAXT_MAX = ', MAXT_MAX, 'deg F'
            WRITE( FDESC3D( 6 ),94030 ) 'INCREMENT= ', TINTV, 'deg F'

            K = MXNDET
            DO 34 J = 1, MXDET
                K = K + 1
                VNAME3D( J ) = EMTLST ( EMTIDX( K ) )
                UNITS3D( J ) = EMTUNIT( EMTIDX( K ) )
                VDESC3D( J ) = EMTDESC( EMTIDX( K ) )
                VTYPE3D( J ) = M3REAL
   34       CONTINUE
 
            NNAME = PROMPTMFILE(
     &         'Enter logical name for output DIURNAL MASTER EF ' //
     &         'output file', FSNEW3, 'MEFSD', 'EMISFAC' )

        ENDIF
C...................................................................
C.........  Read targeted EF update file if it exists ............
C...................................................................

        NUPDAT = 0
        IF( SDEV .GT. 0 ) THEN
            IREC  = 0
            J     = 0
            EFLAG = .FALSE.

   35       CONTINUE

                READ( SDEV, 93100, END=38, IOSTAT=IOS ) BUFFER

                IREC = IREC + 1
                L = INDEX( BUFFER, '-' ) 
 
                IF( IOS .NE. 0 ) THEN

                    EFLAG = .TRUE.
                    WRITE( MESG,94010 )
     &                  'Error ', IOS,
     &                  'reading TARGETED EF UPDATE file at line', IREC
                    CALL M3MESG( MESG )
                    GO TO 35            ! To head of read loop

                ELSEIF( L .EQ. 0 ) THEN
                    K  = STR2INT( BUFFER( 1:TRIMLEN( BUFFER ) ) )
                    KK = K 

                ELSE
                    K  = STR2INT( BUFFER( 1  :L-1 ) )
                    KK = STR2INT( BUFFER( L+1:TRIMLEN( BUFFER ) ) )

                ENDIF

                DO 37 JJ = K, KK

                    J = J + 1

                    IF( J .GT. MXPSI ) THEN
                        WRITE( MESG, 94010 )
     &                     'Maximum number of EF updates', MXPSI, 
     &                     'exceeded at line', IREC, ', file MEFUPD'
                        CALL M3MESG( MESG )

                    ELSE

                        UPDATE( J ) = JJ

                        IF( UPDATE( J ) .LE. UPDATE( J-1 ) ) THEN

                            EFLAG = .TRUE.
                            WRITE( MESG,94010 )
     &                          'Error: Entries in MEFUPD are not ' //
     &                          'uniquely sorted in ascending order'
                            CALL M3MESG( MESG )

                        ENDIF

                    ENDIF

   37           CONTINUE

            GO TO 35  ! To head of read loop

   38       CONTINUE
            NUPDAT = J

            IF( MPNUM .GT. MXPSI ) THEN

                EFLAG = .TRUE.
                WRITE( *, 92000 ) ' '
                WRITE( *, 92010 )
     &                 'Records count for PSI updates       :', NUPDAT,
     &                 'Max dimensioned record-count (MXPSI):', MXPSI,
     &                 'Max record-count exceeded.'
            ENDIF

            IF( EFLAG ) THEN
                CALL M3EXIT( 'EMISFAC', 0, 0,
     &                       'Error in MEFUPD file', 2 )

            ENDIF

        ENDIF

C...................................................................
C.........  Read in and process MPREF  .............................
C...................................................................
        I     = 0
        IREC  = 0
        J     = 0
        JJ    = 0 
        EFLAG = .FALSE.
   55   CONTINUE   ! Start loop for reading MPREF 

C.............  Read first line of 2-line records
            READ( NRDEV, *, END=61, IOSTAT=IOS ) P1, MPTYPE

            I    = I + 1
            IREC = IREC + 1

            IF( IOS .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 )
     &            'Error ', IOS, 'reading MPREF file at line', IREC
                CALL M3MESG( MESG )

            ELSEIF( I .GT. MXPREF ) THEN
                WRITE( MESG, 94010 )
     &                'Maximum number of master PSIs', MXPREF,
     &                 'exceeded at line', IREC, ', file MPREF'
                CALL M3MESG( MESG )

C.............  If no overflow and no I/O errors...
            ELSE

                IPSIA ( I )   = P1
                INDEXA( I )   = I
                ISTYPE( I,1 ) = MPTYPE

C................  For non-combo (normal) PSIs, increment and store properties
                IF( MPTYPE .EQ. 0 ) THEN
                    J = J + 1 
                    ISTYPE( I,2 ) = J

                    READ( NRDEV, * ) ( NIDX( J,K ), K=1,NNFAC ),
     &                               ( CIDX( J,K ), K=1,NCFAC )

                    IREC = IREC + 1

C................  For combo PSIs, increment and store properties
                ELSEIF( MPTYPE .GT. 0 ) THEN

C.....................  Make sure that combination EF doesn't use more than
C                       the legal number of other PSIs
                    IF( MPTYPE .GT. MXNPT ) THEN

                        EFLAG = .TRUE.

                        WRITE( MESG, 94010 ) 
     &                         'Problem at line', IREC, ', file MPREF:'
                        CALL M3MESG( MESG )

                        WRITE( *, 92010 ) 
     &                   'Attmptd # PSIs used in combo        :',MPTYPE,
     &                   'Max     # PSIs used in combo (MXNPT):',MXNPT

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
     &                 '(profile', P1, ', line', IREC, ', file MPREF)'
                    CALL M3MESG( MESG )

                ENDIF

            ENDIF

        GO TO 55    ! End loop for reading MPREF
   61   CONTINUE

        MRNUM  = I
        MPURE  = J  ! Note: J  <= I
        MCOMBO = JJ ! Note: JJ <= I

        IF( MRNUM .GT. MXPREF ) THEN

            EFLAG = .TRUE.
            WRITE( *, 92000 ) ' '
            WRITE( *, 92010 )
     &             'Records count for master PSIs        :', MRNUM ,
     &             'Max dimensioned record-count (MXPREF):', MXPREF,
     &             'Max record-count exceeded.'
        ENDIF

        IF( EFLAG ) THEN
            CALL M3EXIT( 'EMISFAC', 0, 0,
     &                   'Error in MPREF file', 2 )
        END IF

C.........  Sort master PSI list
        CALL SORTI1( MRNUM, INDEXA, IPSIA )

        DO 63 I = 1, MRNUM
           IPSI( I ) = IPSIA( INDEXA( I ) )
   63   CONTINUE

C...................................................................
C.........  Read in and process MPLIST file             ............
C...................................................................

C.........  Create SORTED list of UNIQUE parameter scheme indices (PSIs)
C.........  Do this sorting and uniqueness no to avoid a _huge_ array

        I     = 0 
        IREC  = 0 
        EFLAG = .FALSE.
   44   CONTINUE   ! Start loop for reading MPLIST

            READ( NLDEV, *, END=48 ) P1, P2, P3, ( P4( JJ ), JJ=1,24 )

            IREC = IREC + 1

C.............  Loop through all entries on line (for all hours)
            DO 47 JJ = 1, 24
 
                IM = LOCAT1( P4( JJ ), I, IPSIR )

C.................  If PSI not found, and not first entry...
                IF( IM .GT. 0 .AND. I .NE. 0 ) THEN
 
                    I = I + 1

                    IF( I .GT. MXPSI ) THEN        ! Overflow will occur
                        WRITE( MESG, 94010 )
     &                     'Maximum number of unqiue PSIs', MXPSI,
     &                     'exceeded at line', IREC, ', file MPLIST'
                        CALL M3MESG( MESG )
 
                    ELSE

C......................... Shift list by 1, starting from bottom of list
                        DO 46 K = I, IM+1, -1
                            IPSIR( K ) = IPSIR( K-1 )
   46                   CONTINUE

C......................... Inert new PSI into sorted list 
                        IPSIR( IM ) = P4( JJ )
                    ENDIF

C.................  First entry...
                ELSEIF( I .EQ. 0 ) THEN
                    I = I + 1
                    IPSIR( I ) = P4( JJ )
 
                ENDIF
 
   47       CONTINUE    ! End loop on hours
 
        GO TO 44        ! End loop for reading MPLIST
   48   CONTINUE

C.........  Set (current) count of scenario-specific PSIs
        MPNUM = I

        IF( MPNUM .EQ. 0 ) THEN

            CALL M3EXIT( 'EMISFAC', 0, 0, 
     &                   'Empty MPLIST file not allowed', 1 )

        ENDIF

C.........  Ensure that all PSI on which others depend are in IPSIR
        DO 69 N = 1, MRNUM

            MPTYPE = ISTYPE( N,1 )

C.............  For combination emission factors...
            IF( MPTYPE .GT. 0 ) THEN

C.................  Set pointer to combo-only lists
                DEPOINT = ISTYPE( N,2 )

C.................  Loop through all factors depended on by Ith PSI
                DO 68 J = 1,MPTYPE

C.....................  Get Jth PSI from list for Ith PSI
                    P1 = IPLIST( J, DEPOINT )

C.....................  Locate Jth PSI in sorted scenario-specific PSI list
                    IM = LOCAT1( P1, MPNUM, IPSIR )

                    IF( IM .GT. 0 ) THEN      ! Jth PSI not in scen-spec list

                        MPNUM = MPNUM + 1

                        IF( MPNUM .GT. MXPSI ) THEN
                            WRITE( MESG, 94010 )
     &                         'Maximum number of unqiue PSIs', MXPSI,
     &                         'exceeded during list update for ' //
     &                         'combo EFs'
                            CALL M3MESG( MESG )

                        ELSE

C............................. Shift list by 1, starting from bottom of list
                            DO 66 K = MPNUM, IM+1, -1
                                IPSIR( K ) = IPSIR( K-1 )
   66                       CONTINUE
    
C............................. Inert new PSI into sorted list 
                            IPSIR( IM ) = P1
                        ENDIF

                    ENDIF

   68           CONTINUE

            ENDIF

   69   CONTINUE
        
        IF( MPNUM .GT. MXPSI ) THEN

            EFLAG = .TRUE.
            WRITE( *, 92000 ) ' '
            WRITE( *, 92010 )
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
   39   CONTINUE   ! Start loop for reading MEFTEMP

            READ( EDEV, *, END=41, IOSTAT=IOS ) 
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
     &             MXPSI, 'exceeded at line', IREC, ', file MEFTEMP'
                CALL M3MESG( MESG )

            ELSE    

                MMTREF( I ) = P1
                NTIDX ( I ) = P2
                DO 40 J = 1, P2
                    TIDXPSI( J,I ) = P5( J )
   40           CONTINUE

            ENDIF

        GO TO 39

   41   CONTINUE
        NTREF = I 

        IF( MPNUM .GT. MXPSI ) THEN

            EFLAG = .TRUE.
            WRITE( *, 92000 ) ' '
            WRITE( *, 92010 )
     &             'Records count of temperature indices:', NTREF ,
     &             'Max dimensioned record-count (MXPSI):', MXPSI,
     &             'Max record-count exceeded.'
        ENDIF

C.........  Check if the number of entries in the temperature-index file
C.........  is the same as the number in MPLIST.  If not - give warning.
        IF( MPNUM .NE. NTREF ) THEN

            WRITE( *,92000 ) ' ' 
            WRITE( *,92010 ) 
     &                 'The number of entries in MEFTEMP is', NTREF 
            WRITE( *,92010 ) 
     &                 'The number of entries in MPLIST  is', MPNUM
            WRITE( *,92000 ) 
     &        ' ', 
     &        '*****************************************************',
     &        'WARNING: The parameter scheme indices in MEFTEMP are',
     &        '         inconsistent with those in the MPLIST file',
     &        '*****************************************************',
     &        ' '
        ENDIF

        IF ( EFLAG ) THEN
            CALL M3EXIT( 'EMISFAC', 0, 0,
     &                   'Error in MEFTEMP file', 2 )
        END IF

C.........  Update temperature indices lists so pure PSIs include
C.........  indices for the combo PSIs that use them

        EFLAG = .FALSE.
        DO 75 I = 1, NTREF

C.............  Find PSI for temperature index in main PSI list
            K = FIND1 ( MMTREF( I ), MRNUM, IPSI )

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
                DO 73 J = 1, DEPCNT

C.....................  Find each PSI used in temperature index
                    L = FIND1( IPLIST( J,DEPOINT ), NTREF, MMTREF )

C.....................  Loop through temp indices for PSI( I ) 
                    DO 72 T = 1, NTIDX( I )

C.........................  Find temp index location in list for PSI(L)
                        KK = LOCAT1( TIDXPSI( T,I ), NTIDX( L ), 
     &                               TIDXPSI( 1,L )              )

C.........................  Insert temperature index in list for I
                        IF( KK .GT. 0 ) THEN

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
                                DO 71 JJ = NTIDX( L ), KK+1, -1
                                    TIDXPSI( JJ,L ) = TIDXPSI( JJ-1,L )
   71                           CONTINUE
 
C.................................  Insert new TMMI for PSI(L) from PSI(I)
                                TIDXPSI( KK,L ) = TIDXPSI( T,I )
 
                            ENDIF  ! End overflow check
                        ENDIF

   72               CONTINUE
   73           CONTINUE

            ENDIF 

   75   CONTINUE   !  ... Phew!

        IF ( EFLAG ) THEN
            CALL M3EXIT( 'EMISFAC', 0, 0,
     &                   'Bad TEMPERATURE INDEX file', 2 )
        END IF

C.........  Initialize constants for loops
        MAXTL = INT( ( MAXT_MAX - MAXT_MIN ) / TINTV ) + 1

C.........  Create new pure EFs where needed
        DO 201 IPM = 1, MPNUM

C.............  Get indice from MPLIST
            IPARM = IPSIR( IPSIDX( IPM ) )
            IPTIM = SEC2TIME( IPARM )

            IF( IPARM .GT. 24*3600 ) THEN
                WRITE( MESG, 94010 ) 
     &            'Number of parameter combinations cannot exceed 86400'
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

C.............  Search for MPLIST indice in PSI list used in temperature index
            K = FIND1( IPARM, NTREF, MMTREF )

            IF( K .LE. 0 ) THEN
                WRITE( MESG, 94010 ) 
     &            'Parameter scheme index', IPARM,
     &            'from MPLIST not found in MEFTEMP'
                CALL M3EXIT('EMISFAC', 0, 0, MESG, 2 )
            ENDIF
 
C.............  Write message to screen
            WRITE( *,92000 ) ' '
            WRITE( *,92020 ) 'Processing PSI', IPARM, '...'

C.............  Initialize MIN/MAX temperature indice counter
            MMTC( K ) = 1

C.............  Initialize flags to say EFs were not calculated
            TFLAG = .FALSE.  ! Non-diurnal
            RFLAG = .FALSE.  ! Diurnal

C.............  Set count of PSI on which current PSI depends
            DEPCNT  = ISTYPE( INDEXA( J ),1 )
            DEPOINT = ISTYPE( INDEXA( J ),2 )

C.............  Recognize only indice representing pure EF calculations
            IF( DEPCNT .EQ. 0 ) THEN

C.................  Set update flag using search of MPLIST indice in update list
                UFLAG = .FALSE.
                IF ( NUPDAT .GT. 0 .AND.
     &               FIND1( IPARM, NUPDAT, UPDATE ) .GT. 0 ) 
     &               UFLAG = .TRUE.

C.................  Check to see any diurnals for this profile indice exist
                IF( .NOT. CHECK3
     $              ( DNAME, ALLVAR3, SDATE3D, IPTIM ) ) THEN

C.....................  Init diurnal arrays b/c blanks will exist on write3
                    DO 81 I = 1, NVTYPE
                        DO 80 JJ = 1, NMINMAX
                            WDL_VOC( JJ, I ) = BADVAL3
                            DNL_VOC( JJ, I ) = BADVAL3
                            HOT_VOC( JJ, I ) = BADVAL3
                            CRC_VOC( JJ, I ) = BADVAL3
   80                   CONTINUE
   81               CONTINUE

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

C.................  Check for Non-diurnal EFs for profile indice
C.................  Also check if an update has been specifically requested
                IF( .NOT. CHECK3 ( ENAME, ALLVAR3, SDATE3D, IPTIM ) .OR.
     &              UFLAG ) THEN

C.....................  Loop over temperature range for non-diurnal factors
                    ITEMP = 1 ! Specify that more has changed than just temp
                    DO 101 T = 1, NTEMP

                        TMMI   = TIDXPSI( MMTC( K ),K )

C.........................  Get temperature index for PSI(K) and count(K)
C.............................  Check if count has not yet reached NTIDX and
C.............................  diurnal emissions are not zero
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

C.........................  Increase counter if using min/max, otherwise
C.........................  reset min/max current temp is not between them
                        IF( TEMPIN( T ) .EQ. 
     &                      TF_MIN + ( TF_MAX - TF_MIN ) / 2. ) THEN

                            MMTC( K ) = MMTC( K ) + 1  ! Set for next time

C.............................  Diurnal emission factors to be calculated 
                            DFLAG = .TRUE. 

                        ELSE

C.............................  Round current temp to nearest in diurnal set
                            TF = MINT_MIN + TINTV *
     &                           NINT( ( TEMPIN(T) - MINT_MIN )/TINTV ) 

C.............................  Set new min and max based on rounded temp
                            TF_MIN = MAX( TF - TINTV, MINT_MIN )
                            TF_MIN = MIN( TF_MIN, MINT_MAX )

                            TF_MAX = MIN( TF + TINTV, MAXT_MAX )

C............................. Set new temperature index for new min/max
                            TMMI = MAXTL *
     &                           NINT( ( TF_MIN-MINT_MIN ) / TINTV ) +
     &                           NINT( ( TF_MAX-MAXT_MIN ) / TINTV ) + 1

                        ENDIF

                        TF = TEMPIN( T )

C.........................  Create input file for MOBILE5a
                        REWIND( MODEV )                  
                        CALL MOINBLD( MODEV, IPARM , DEPOINT, NIDX,
     &                                CIDX , VOLORG, JYEAR, 
     &                                TF   , TF_MIN, TF_MAX        )

C.........................  Initialize error count
                        IERR = 0 

C.........................  Create EFs for current NPSI, CPSI, temperature
                        REWIND( MODEV )                  
                        CALL MOBILE ( IERR, ITEMP, MODEV, EFACT, DFACT )

                        DO 88 I = 1, NVTYPE

C.............................  Put EFACT into non-diurnal arrays
                            EXH_CO ( T, I ) = EFACT( 1, I )
                            EXH_NOX( T, I ) = EFACT( 2, I )
                            EXH_VOC( T, I ) = EFACT( 3, I )
                            EVP_VOC( T, I ) = EFACT( 4, I )
                            RNL_VOC( T, I ) = EFACT( 5, I )
                            RST_VOC( T, I ) = EFACT( 6, I )

C.............................  Put DFACT into diurnal arrays
                            IF( DFLAG ) THEN
                                WDL_VOC( TMMI, I ) = DFACT( 1, I )
                                DNL_VOC( TMMI, I ) = DFACT( 2, I )
                                HOT_VOC( TMMI, I ) = DFACT( 3, I )
                                CRC_VOC( TMMI, I ) = DFACT( 4, I )

C.................................  Diurnal emission factors were calculated 
                                RFLAG = .TRUE. 

                            ENDIF

   88                   CONTINUE

C.........................  Non-diurnal emission factors were calculated 
                        TFLAG = .TRUE. 

  101               CONTINUE  ! T, temperature for non-diurnal factors

C.....................  Add non-diurnal EFs to master list 
                    IF( .NOT. WRITE3
     &                  ( ONAME, ALLVAR3, SDATE3D, IPTIM, EXH_CO) ) THEN

                        WRITE( MESG, 94010 ) 
     &                    'Could not write non-diurnal EF set', IPARM,
     &                    'to ' // ONAME( 1:TRIMLEN( ONAME ) )
                        CALL M3EXIT('EMISFAC', 0, 0, MESG, 2 )

                    ENDIF  ! Write to non-diurnal master list

                ELSE

                    WRITE( MESG, 94010 )
     &                'Non-diurnal emission factors for set', IPARM,
     &                'already exist in ' // ENAME( 1:TRIMLEN(ENAME) )
                    CALL M3MESG( MESG )

                ENDIF      ! Non-Diurnal EFs have been created or not

                ITEMP = 1            ! More has changed than just temperature

                VFLAG( K ) = .TRUE.  ! Initialize - default: update non-diu EFs
                SFLAG( K ) = .TRUE.  ! Initialize - default: update diurnal EFs

C.................  Create remaining diurnal factors for IPARM by
C.................  looping through remaining temperatures
                DO 151 T = MMTC( K ), NTIDX( K )

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

C.............................  Create input file for MOBILE5a
                            REWIND( MODEV )                  
                            CALL MOINBLD( MODEV, IPARM , DEPOINT, NIDX,
     &                                    CIDX , VOLORG, JYEAR, 
     &                                    TF   , TF_MIN, TF_MAX        )

C.............................  Initialize error count
                            IERR = 0 

C.............................  Create EFs current temperature combo
                            REWIND( MODEV )                  
                            CALL MOBILE ( IERR, ITEMP, MODEV, 
     &                                    EFACT, DFACT        )

C.............................  Put DFACT into diurnal arrays
                            DO 121 I = 1, NVTYPE

                                WDL_VOC( TMMI, I ) = DFACT( 1, I )
                                DNL_VOC( TMMI, I ) = DFACT( 2, I )
                                HOT_VOC( TMMI, I ) = DFACT( 3, I )
                                CRC_VOC( TMMI, I ) = DFACT( 4, I )

  121                       CONTINUE
 
C.............................  Diurnal emission factors were calculated 
                            RFLAG = .TRUE. 

                        ENDIF  ! If need to call MOBILE or not
                    ENDIF      ! If valid min/max temp combo or not
 
  151           CONTINUE       ! Loop through temperature indices for IPARM

C.................  If factors were not calculated, then don't output them
                IF( .NOT. TFLAG ) VFLAG( K ) = .FALSE.  ! Non-diurnal
                IF( .NOT. RFLAG ) SFLAG( K ) = .FALSE.  ! Diurnal

C.................  Write out diurnal emission factors for IPARM
C.....................  Add non-diurnal EFs to master list 
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
     &                'already exist in ' // DNAME( 1:TRIMLEN(DNAME) )
                    CALL M3MESG( MESG )

                ELSE
                    WRITE( MESG, 94010 )
     &                'Error: Check MEFTEMP file for zero TMMI ' //
     &                'values for EFs set', IPARM
                    CALL M3EXIT('EMISFAC', 0, 0, MESG, 2 )

                ENDIF      ! update of diurnal file is needed or not

            ENDIF          ! EFs are pure or combo

  201   CONTINUE           ! IPM

C.........  Create new mixed EFs where needed
        DO 401 IPM = 1, MPNUM

C.............  Get indice from MPLIST
            IPARM = IPSIR( IPSIDX( IPM ) )
            IPTIM = SEC2TIME( IPARM )

C.............  Search for MPLIST indice in MPREF index
            J      = FIND1 ( IPARM, MRNUM, IPSI )
            DEPCNT = ISTYPE( INDEXA( J ),1 )
            DEPOINT= ISTYPE( INDEXA( J ),2 )

C.............  Recognize only indice representing combo EF calculations
            IF( DEPCNT .GT. 0 ) THEN

C.................  Write message to screen
                WRITE( *,92000 ) ' '
                WRITE( *,92020 ) 'Processing PSI', IPARM, '...'

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

  381           CONTINUE
  383           CONTINUE

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
     &                'already exist in ' // ENAME( 1:TRIMLEN(ENAME) )
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
     &                'already exist in ' // DNAME( 1:TRIMLEN(DNAME) )
                    CALL M3MESG( MESG )

                ENDIF      ! Missing diurnal for this IPARM or not

            ENDIF          ! EFs are combo or pure

  401   CONTINUE           ! End loop on parameter schemes

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
 
94010   FORMAT( 10 ( A, :, I5, :, 2X ) )
94030   FORMAT( A, F10.3, 1X, A )
94060   FORMAT( 8( F6.2, ', ' ) )
94070   FORMAT( 5( A, 1X, F6.1, 1X, A, 1X, F6.1, 1X ) )
94080   FORMAT( I8 )

        END

