
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
C***********************************************************************

        IMPLICIT NONE

C...........   INCLUDES:

c       INCLUDE 'EFS.EXT'        !  Emissions factors include paramaters
        INCLUDE 'CHDIMS3.EXT'
        INCLUDE 'MBDIMS3.EXT'    !  Mobile include paramaters
        INCLUDE 'PARMS3.EXT'     !  I/O API parameters
        INCLUDE 'IODECL3.EXT'    !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'     !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS and their descriptions:

        INTEGER      FIND1
        INTEGER      GETNUM
        REAL         GETREAL
        LOGICAL      GETYN 
        INTEGER      LOCAT1
        INTEGER      PROMPTFFILE
        CHARACTER*16 PROMPTMFILE
        INTEGER      TRIMLEN

        EXTERNAL     FIND1, GETNUM, GETREAL, GETYN, LOCAT1, 
     &               PROMPTFFILE, PROMPTMFILE, TRIMLEN

C...........   LOCAL VARIABLES and their descriptions:

C...........   Simulation-specific MOBILE5a parmeter list variables
        INTEGER IPSIR ( NLIST*24 )       ! Parameter scheme index for 24 hours
        INTEGER MPNUM                    ! Actual no. of normal parm list recs
        INTEGER P1, P2, P3, P4( 24 )     ! Dummies

C...........   MOBILE5a parameter x-ref variables
        INTEGER MRNUM                    ! Actual no. of normal x-ref entries
        INTEGER IPSI  ( NPREF )          ! Normal x-ref scheme index
        INTEGER NIDX  ( NPREF, NNFAC )   ! Normal x-ref entries
        INTEGER CIDX  ( NPREF, NCFAC )   ! Control Normal x-ref entries
        INTEGER P1

        COMMON / MPARMS / MRNUM, IPSI, NIDX, CIDX

C...........   Master EFs list variables
        INTEGER NMSTR                      ! Actual no. of master EF records
        REAL    EXH_CO ( MAXTMP, NVTYPE )  ! Emis factor: CO exhaust
        REAL    EXH_NOX( MAXTMP, NVTYPE )  ! Emis factor: NOX exhaust
        REAL    EXH_VOC( MAXTMP, NVTYPE )  ! Emis factor: VOC exhaust
        REAL    EVP_VOC( MAXTMP, NVTYPE )  ! Emis factor: VOC evaporative
        REAL    RNL_VOC( MAXTMP, NVTYPE )  ! Emis factor: VOC running
        REAL    RST_VOC( MAXTMP, NVTYPE )  ! Emis factor: VOC resting
        REAL    RF1_VOC( MAXTMP, NVTYPE )  ! Emis factor: VOC refuel /mile
        REAL    RF2_VOC( MAXTMP, NVTYPE )  ! Emis factor: VOC refuel /gallon
        REAL    HOT_VOC( MAXTMP, NVTYPE )  ! Emis factor: VOC hot soak
        REAL    CRC_VOC( MAXTMP, NVTYPE )  ! Emis factor: VOC crank case
        REAL    RHR_VOC( MAXTMP, NVTYPE )  ! Emis factor: VOC resting /hour

        COMMON / EFS / EXH_CO , EXH_NOX, EXH_VOC, EVP_VOC, RNL_VOC, 
     &                 RST_VOC, RF1_VOC, RF2_VOC, HOT_VOC, CRC_VOC,
     &                 RHR_VOC

C...........   Variables for saving for MOBILE5
        INTEGER      IREAD
        REAL         CR12HC
        REAL         CR12CO
        CHARACTER*40 TNAME
        COMMON /IMPAR7/ TNAME,IREAD
        COMMON /IM12HC/ CR12HC(19,20,5,2)
        COMMON /IM12CO/ CR12CO(19,20,5,2)

C...........   Local variables
        REAL    BTEMP
        REAL    ETEMP
        REAL    TEMP
        REAL    TINTV
        REAL    EFACT ( MXEMIS, NVTYPE )  ! EFs for each temperature
        REAL    TEMPIN( MAXTMP )          ! Input temperatures in deg F
        INTEGER ITEMP                     ! Indicator for MOBILE5 for only temp changing
        INTEGER NTEMP                     ! Actual no. of temps for eval of EFs
        INTEGER IPSIDX( NPREF )           ! Sorting index for IPSIR

C...........   Other local variables

        REAL    PCNTPRO( MXNPT, NPREF )
        INTEGER IPLIST ( MXNPT, NPREF )
        INTEGER ISTYPE ( NPREF, 2     )
        INTEGER MPTYPE, J, JJ, MCOMBO, MPURE

        INTEGER T
        INTEGER I, IERR, IM, IPARM, IPM, 
     &          J, JJ, JDATE, JTIME, JYEAR, K

        INTEGER LDEV, MODEV, NLDEV, NRDEV

        CHARACTER*16  MNAME, ONAME
        CHARACTER*256 MESG

C***********************************************************************
C   begin body of program EMISFAC

        LDEV = INIT3()
 
        WRITE( *,92000 )
     &  ' ',
     &  'Program EMISFAC drives the MOBILE5a program by supplying',
     &  'a range of ambient temperatures and a scenario-specific',
     &  'MOBILE5 parameter file MPLIST. Using multiple calls',
     &  'to MOBILE5a, EMISFAC creates a master emissions factors',
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

C.........  Prompt for MOBILE5a input file
        MODEV = PROMPTFFILE(
     &          'Enter logical name for temporary MOBILE5a input file',
     &          .FALSE., .TRUE., 'MO5AIN', 'EMISFAC' )

C.........  Prompt for MEFS (master list EFs) input file for modeling year
        MNAME = PROMPTMFILE(
     &          'Enter logical name for MASTER EF file or "NONE"',
     &          FSUNKN3, 'MEFS', 'EMISFAC' )

C.........  Prompt for MPLIST (scenario specific parameter scheme list)
        NLDEV = PROMPTFFILE(
     &          'Enter logical name for PARAMETER LIST input file',
     &          .TRUE., .TRUE., 'MPLIST', 'EMISFAC' )
        
C.........  Prompt for MPREF (general MOBILE5a parameter x-ref file)
        NRDEV = PROMPTFFILE( 
     &          'Enter logical name for MPREF input file',
     &          .TRUE., .TRUE., 'MPREF', 'EMISFAC' )

C.........  Prompt for modeling year
        JYEAR = 1988
        JYEAR = GETNUM( 0, 9999, JYEAR,
     &         'Enter modeling year for emissions factors (YYYY)' )

C.........  Prompt for start temperature
        BTEMP = GETREAL( 1.,99., 1.,
     &                   'Enter minimum temperature [deg F]' )

C.........  Prompt for end temperature
        ETEMP = GETREAL( BTEMP,99., 99., 
     &                   'Enter maximum temperature [deg F]' )

C.........  Prompt for temperature increment
        TINTV = GETREAL( 0.0001, 50., 2., 
     &                   'Enter temperature increment [deg F]' )

C.........  Set temperatures ..................................
        I = 0 
        TEMP = BTEMP
        DOWHILE( TEMP .LE. ETEMP )

            I = I + 1
            TEMPIN( I ) = TEMP
            TEMP = BTEMP + TINTV * I 

        ENDDO
        NTEMP = I

C.........  Get number in Master list if list exists
        IF( MNAME .NE. 'NONE' ) THEN  ! Master list exists

C.............  Get file description
            IF( .NOT. DESC3( MNAME ) ) THEN

                WRITE( MESG, 94010 ) 
     &                 'Could not get description of Master list'
                CALL M3EXIT( 'EMISFAC', 0, 0, MESG, 2 )

            ENDIF

            NMSTR = MXREC3D
            ONAME = MNAME

cmh         JDATE = SDATE3D
cmh         JTIME = 0
cmh         CALL NEXTIME( JDATE, JTIME, MXREC3D )
 
        ELSE

C.............  Set up Master List file description

            FTYPE3D = GRDDED3
            P_ALP3D = DBLE( IMISS3 )
            P_BET3D = DBLE( AMISS3 )
            P_GAM3D = DBLE( AMISS3 )
            XCENT3D = 0.0D0
            YCENT3D = 0.0D0
            XORIG3D = DBLE( BTEMP )     ! Starting temperature
            YORIG3D = DBLE( AMISS3 )
            XCELL3D = DBLE( TINTV )     ! Temperature increment
            YCELL3D = DBLE( AMISS3 )
            SDATE3D = JYEAR * 1000 + 1
            STIME3D = 1
            TSTEP3D = 1
            NVARS3D = MXEMIS
            NCOLS3D = MAXTMP
            NROWS3D = NVTYPE 
            NLAYS3D = 1
            NTHIK3D = 1
            GDTYP3D = IMISS3
            VGTYP3D = IMISS3
            VGTOP3D = AMISS3
            GDNAM3D = '' 

            FDESC3D( 1 ) = 'Emissions Factors Master list for'
            WRITE  ( MESG, 94070 ) 'temperatures from', BTEMP, 
     &                             'to', ETEMP
            FDESC3D( 2 ) = MESG

            WRITE  ( MESG, 94070 ) 'on increments of', TINTV
            FDESC3D( 3 ) = MESG

            WRITE  ( MESG, 94070 ) 'Total number of temperatures='
            FDESC3D( 4 ) = MESG

            WRITE  ( MESG, 94080 ) NTEMP
            FDESC3D( 5 ) = MESG

            DO 33 J = 1, MXEMIS
                VNAME3D( J ) = EMTLST ( J )
                UNITS3D( J ) = EMTUNIT( J )
                VDESC3D( J ) = EMTDESC( J )
                VTYPE3D( J ) = M3REAL
   33       CONTINUE
 
            ONAME = PROMPTMFILE(
     &              'Enter logical name for MASTER EF output file',
     &              FSUNKN3, 'MEFS', 'EMISFAC' )

cmh         JDATE = SDATE3D
cmh         JTIME = 0
            NMSTR = 0

        ENDIF

        I   = 0 
C.........  Read in MPLIST .......................................
   44   CONTINUE   ! Start loop for reading MPLIST

            READ( NLDEV, *, END=48 ) P1, P2, P3, ( P4( JJ ), JJ=1,24 )

C.............  Create sorted list of unique parameter scheme indices
            DO 47 JJ = 1, 24

                IM = LOCAT1( P4( JJ ), I, IPSIR )

                IF( IM .GT. 0 ) THEN
                    I = I + 1

                    IF( I .GT. NLIST ) THEN
                        WRITE( MESG, 94010 ) 
     &                     'Maximum dimension NLIST=', NLIST, 
     &                     'was exceeded while reading MPLIST'
                        CALL M3EXIT( 'EMISFAC', 0, 0, MESG, 1 )

                    ELSE
                        DO 46 K = I, IM+1, -1
                            IPSIR( K ) = IPSIR( K-1 )
   46                   CONTINUE

                        IPSIR( IM ) = P4( JJ )
                    ENDIF
                ELSEIF( I. EQ. 0 ) THEN

                    I = I + 1
                    IPSIR( I ) = P4( JJ )

                ENDIF

   47       CONTINUE

        GO TO 44   ! End loop for reading MPLIST
   48   MPNUM = I

C.........  Read in MPREF ..........................................
        I  = 0
        J  = 0
        JJ = 0 
   55   CONTINUE   ! Start loop for reading MPREF 

            READ( NRDEV, *, END=61 ) P1, MPTYPE
            I = I + 1

            IF( I .GT. NPREF ) THEN
                WRITE( MESG, 94010 ) 
     &                 'Maximum dimension NPREF=', NPREF, 
     &                 'was exceeded while reading MPREF'
                CALL M3EXIT( 'EMISFAC', 0, 0, MESG, 1 )
            ELSE

                IPSI  ( I ) = P1
                ISTYPE( I,1 ) = MPTYPE

                IF( MPTYPE .EQ. 0 ) THEN
                    J = J + 1 
                    ISTYPE( I,2 ) = J

                    READ( NRDEV, * ) ( NIDX( J,K ), K=1,NNFAC ),
     &                               ( CIDX( J,K ), K=1,NCFAC )

                ELSEIF( MPTYPE .GT. 0 ) THEN

                    IF( MPTYPE .GT. MXNPT ) THEN

                        WRITE( MESG, 94010 ) 
     &                     'Maximum dimension MXNPT=', MXNPT, 
     &                     'was exceeded while reading profile', P1,
     &                     'in file MPREF'
                        CALL M3EXIT( 'EMISFAC', 0, 0, MESG, 1 )

                    ELSE
                        JJ = JJ + 1
                        ISTYPE( I,2 ) = JJ

                        READ( NRDEV, * ) ( PCNTPRO( K,JJ ), 
     &                                     IPLIST ( K,JJ ), K=1,MPTYPE )
                    ENDIF

                ELSE

                    WRITE( MESG, 94010 ) 
     &                 'Profile type flag illegal value:', MPTYPE,
     &                 'for MPREF profile', P1
                    CALL M3EXIT( 'EMISFAC', 0, 0, MESG, 1 )

                ENDIF

            ENDIF

        GO TO 55   ! End loop for reading MPREF
   61   MRNUM  = I
        MPURE  = J
        MCOMBO = JJ

C.........  Ensure that all PSI on which others depend are in IPSIR

        DO 69 I = 1, MRNUM

            MPTYPE = ISTYPE( I,1 )
            IF( MPTYPE .GT. 0 ) THEN

                DO 68 J = 1,MPTYPE

                    P1 = IPLIST( J, ISTYPE( I,2 ) )
                    IM = LOCAT1( P1, MPNUM, IPSIR )

                    IF( IM .GT. 0 ) THEN
                        MPNUM = MPNUM + 1

                        IF( MPNUM .GT. NLIST ) THEN
                            WRITE( MESG, 94010 ) 
     &                         'Maximum dimension NLIST=', NLIST, 
     &                         'was exceeded while reading MPLIST'
                            CALL M3EXIT( 'EMISFAC', 0, 0, MESG, 2 )

                        ELSE
                            DO 66 K = MPNUM, IM+1, -1
                                IPSIR( K ) = IPSIR( K-1 )
   66                       CONTINUE
    
                            IPSIR( IM ) = P1
                        ENDIF

                    ELSEIF( IM. EQ. 0 ) THEN

                        MPNUM = MPNUM + 1
                        IPSIR( MPNUM ) =  P1

                    ENDIF

   68           CONTINUE

            ENDIF

   69   CONTINUE
        
C.........  Initialize sorting index
        DO 70 I = 1, MPNUM
            IPSIDX( I ) = I
   70   CONTINUE

C.........  Sort scenario specific IPSI array
        CALL SORTI1( MPNUM, IPSIDX, IPSIR )

C.........  Create new pure EFs where needed
        DO 201 IPM = 1, MPNUM

C.............  Get indice from MPLIST
            IPARM = IPSIR( IPSIDX( IPM ) )

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

C.............  Recognize only indice representing pure EF calculations
            IF( ISTYPE( J,1 ) .EQ. 0 ) THEN

C.................  Check to see if already have EFs for profile indice
                IF( .NOT. CHECK3
     $              ( ONAME, ALLVAR3, SDATE3D, IPARM ) ) THEN

C.....................  Loop over temperature range
                    ITEMP = 1 ! Specify that more has changed than just temp
                    DO 101 T = 1, NTEMP

C.........................  Create input file for MOBILE5a
                        REWIND( MODEV )                  
                        CALL MOINBLD( MODEV, IPARM,
     &                                JYEAR, TEMPIN( T ) )

C.........................  Create EFs for current NPSI, CPSI, temperature
                        REWIND( MODEV )                  
                        CALL MOBILE ( IERR, ITEMP, MODEV, EFACT )

C.....................  Put EFACT into Emissions factors
                        DO 77 J = 1, NVTYPE

                            EXH_CO ( T, J ) = EFACT( 1 , J )
                            EXH_NOX( T, J ) = EFACT( 2 , J )
                            EXH_VOC( T, J ) = EFACT( 3 , J )
                            EVP_VOC( T, J ) = EFACT( 4 , J )
                            RNL_VOC( T, J ) = EFACT( 5 , J )
                            RST_VOC( T, J ) = EFACT( 6 , J )
                            RF1_VOC( T, J ) = EFACT( 7 , J )
                            RF2_VOC( T, J ) = EFACT( 8 , J )
                            HOT_VOC( T, J ) = EFACT( 9 , J )
                            CRC_VOC( T, J ) = EFACT( 10, J )
                            RHR_VOC( T, J ) = EFACT( 11, J )

   77                   CONTINUE

  101               CONTINUE  ! T, temperature

C.....................  Add EFs to master list 
                    IF( .NOT. WRITE3
     &                  ( ONAME, ALLVAR3, SDATE3D, IPARM, EXH_CO) ) THEN

                        WRITE( MESG, 94010 ) 
     &                    'Could not write EF set', IPARM,
     &                    'to ' // TNAME( 1:TRIMLEN( TNAME ) )
                        CALL M3EXIT('EMISFAC', 0, 0, MESG, 2 )

                    ENDIF  ! Write to master list

                ELSE

                    WRITE( MESG, 94010 )
     &                'Emissions factors for indice', IPARM,
     &                'already exists in ' 
     &                // ONAME( 1:TRIMLEN( ONAME ) )
                    CALL M3MESG( MESG )

                ENDIF      ! EFs have been created or not

            ENDIF          ! EFs are pure or combo

  201   CONTINUE           ! IPM

C.........  Create new mixed EFs where needed
        DO 401 IPM = 1, MPNUM

C.............  Get indice from MPLIST
            IPARM = IPSIR( IPSIDX( IPM ) )

C.............  Search for MPLIST indice in MPREF index
            J  = FIND1( IPARM, MRNUM, IPSI )
            JJ = ISTYPE( J,2 )

C.............  Recognize only indice representing combo EF calculations
            IF( ISTYPE( J,1 ) .GT. 0 ) THEN

C.................  Check to see if already have EFs for profile indice
                IF( .NOT. CHECK3
     &              ( ONAME, ALLVAR3, SDATE3D, IPARM ) ) THEN

C.....................  Create combo profile through common block EFS
                    CALL EFCOMBO( ONAME, NTEMP, ISTYPE( J,1 ), 
     &                            PCNTPRO( 1, JJ ), 
     &                            IPLIST ( 1, JJ ), SDATE3D )

C.....................  Store emission factor in main list

                    IF( .NOT. WRITE3
     &                  ( ONAME, ALLVAR3, SDATE3D, IPARM, EXH_CO) ) THEN

                        WRITE( MESG, 94010 ) 
     &                    'Could not write EF set', IPARM,
     &                    'to ' // TNAME( 1:TRIMLEN( TNAME ) )
                        CALL M3EXIT('EMISFAC', 0, 0, MESG, 2 )

                    ENDIF  ! Write to master list

                ELSE

                    WRITE( MESG, 94010 )
     &                'Emissions factors for indice', IPARM,
     &                'already exists in ' 
     &                // ONAME( 1:TRIMLEN( ONAME ) )
                    CALL M3MESG( MESG )

                ENDIF      ! EFs have been created or not

            ENDIF          ! EFs are combo or pure

  401   CONTINUE           ! IPM

        MESG = 'Normal completion of EMISFAC'
        CALL M3EXIT( 'EMISFAC', 0, 0, MESG, 0 )
 
C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )
92010   FORMAT( 5X, A, I5 / 20( 5X, A / ) )

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( 12( 3( 9(F6.3, 1X)/ ) // ) )
93200   FORMAT( 3I6, 2X, 5E12.4 )
93210   FORMAT( 18X, 2X, 5E12.4 )
c temp
93999   FORMAT( F6.2,1X, 8(F6.2,1X) )

C...........   Internal buffering formats............ 94xxx
 
94010   FORMAT( 10 ( A, :, I5, :, 2X ) )
94060   FORMAT( 8( F6.2, ', ' ) )
94070   FORMAT( 5( A, 1X, F6.1, 1X, A, 1X, F6.1, 1X ) )
94080   FORMAT( I8 )

        END

