C ****************************************************************************
C *
C * Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C *                System
C * File: @(#)$Id$
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
C ****************************************************************************/

        PROGRAM SPCPMAT

C***********************************************************************
C  program body starts at line 242
C
C  DESCRIPTION:
C       Construct point source chemical speciation matrix from data 
C       contained in EPS-style chem-split and chem-xref files.
C       Uses fixed NOtoNOX  factor of 0.0197  moles  NO/gm NOX (from EPS)
C       Uses fixed NO2toNOX factor of 0.00217 moles NO2/gm NOX  "
C
C  PRECONDITIONS REQUIRED:
C       Input data for chem-xref and sorted chem-split files;
C       M3IO point-sources file from program GRDAMAT
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       Models-3 I/O
C       FIND1, FIND2, FIND3, GETEFILE, GETYN, TRIMLEN, RDSPROF
C
C  REVISION  HISTORY:
C       Prototype  3/95 by CJC.
C       Version   11/95 by CJC sorts the chem-xref input file
C
C***********************************************************************

      IMPLICIT NONE

C...........   INCLUDES:
        
        INCLUDE 'PTDIMS3.EXT'   !  point-source dimensioning parameters
        INCLUDE 'CHDIMS3.EXT'   !  emis chem parms (inventory + model)
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.


C...........   EXTERNAL FUNCTIONS and their descriptions:
        
        INTEGER         FIND1, FIND2, FIND3, FIND4
        LOGICAL         GETYN
        INTEGER         GETMENU 
        INTEGER         LBLANK
        INTEGER         PROMPTFFILE
        CHARACTER*16    PROMPTMFILE
        INTEGER         RDSPROF
        INTEGER         STR2INT
        REAL            STR2REAL
        INTEGER         TRIMLEN
        
        EXTERNAL        FIND1, FIND2, FIND3, FIND4, GETYN, GETMENU, 
     &                  LBLANK, PROMPTFFILE, PROMPTMFILE, RDSPROF,
     &                  STR2INT, STR2REAL, TRIMLEN


C...........   PARAMETERS and their descriptions:

        REAL        TON2GM
        PARAMETER ( TON2GM = 907184.74 )
        
C...........   LOCAL VARIABLES and their descriptions:

        INTEGER     FIP ! state/county ID
        INTEGER     SCC ! (Source Category Code, of course)
        INTEGER     PID ! plant ID
        INTEGER     SID ! stack ID
        INTEGER     XID ! stack ID
        INTEGER     LID ! last ID read
                                 
C...........   Speciation profiles table
        
        REAL        VREACT
        REAL        ETHADJ
        REAL        FRMADJ
        REAL        ALLADJ
        REAL        VOLADJ   ! volatile conversion factor (e.g. ROG to TOG)
       
        INTEGER     NPROFS                      !  actual number of profiles
        INTEGER     SXREF( MXSPRO )             !  profile indexes (sorted)
        REAL        SCOEF( MXSPRO, NMPOL )      ! profile coefficients


C.......   Actual-SCC  table

        INTEGER      NSCCA
        INTEGER      ISCCA( NPSRC )  !  actutally-occurring SCC table


C...........   Speciation cross-reference table:
C.......   Before sorting:

        INTEGER     NREFA
        INTEGER     INDEXA( NPSRC )
        INTEGER     FIPSFA( NPSRC )
        INTEGER     SCCSFA( NPSRC )
        INTEGER     PLANTA( NPSRC )
        INTEGER     STACKA( NPSRC )
        INTEGER     SPXRFA( NPSRC )
        REAL        VREACA( NPSRC )
        REAL        ETHADA( NPSRC )
        REAL        FRMADA( NPSRC )
        REAL        ALLADA( NPSRC )

C.......   After sorting -- Five parts:   first depends on scc,fip,plant,stack

        INTEGER     NREF1
        INTEGER     FIPSF1( NSREF )
        INTEGER     SCCSF1( NSREF )
        INTEGER     PLANT1( NSREF )
        INTEGER     STACK1( NSREF )
        INTEGER     SPXRF1( NSREF )
        REAL        VREAC1( NSREF )
        REAL        ETHAD1( NSREF )
        REAL        FRMAD1( NSREF )
        REAL        ALLAD1( NSREF )

C.......   Second part depends on scc,fip,plant

        INTEGER     NREF2
        INTEGER     FIPSF2( NSREF )
        INTEGER     SCCSF2( NSREF )
        INTEGER     PLANT2( NSREF )
        INTEGER     SPXRF2( NSREF )
        REAL        VREAC2( NSREF )
        REAL        ETHAD2( NSREF )
        REAL        FRMAD2( NSREF )
        REAL        ALLAD2( NSREF )

C.......   Third part depends on scc,fip

        INTEGER     NREF3
        INTEGER     FIPSF3( NSREF )
        INTEGER     SCCSF3( NSREF )
        INTEGER     SPXRF3( NSREF )
        REAL        VREAC3( NSREF )
        REAL        ETHAD3( NSREF )
        REAL        FRMAD3( NSREF )
        REAL        ALLAD3( NSREF )

C.......   Fourth part depends on scc, state fip

        INTEGER     NREF4
        INTEGER     FIPSF4( NSREF )
        INTEGER     SCCSF4( NSREF )
        INTEGER     SPXRF4( NSREF )
        REAL        VREAC4( NSREF )
        REAL        ETHAD4( NSREF )
        REAL        FRMAD4( NSREF )
        REAL        ALLAD4( NSREF )

C.......   Fifth part depends only on scc

        INTEGER     NREF5
        INTEGER     SCCSF5( NSREF )
        INTEGER     SPXRF5( NSREF )
        REAL        VREAC5( NSREF )
        REAL        ETHAD5( NSREF )
        REAL        FRMAD5( NSREF )
        REAL        ALLAD5( NSREF )


C...........   Volatile Pollutant Type Conversion Table 
C........ After sorting, two parts

C........ First part depends on scc, state fip

        INTEGER     NCNV1
        INTEGER     FIPSV1( MXSCNV )  !  conversion factor STATE ID
        INTEGER     SCCSV1( MXSCNV )  !  conversion factor SCC
        REAL        CONVV1( MXSCNV )  !  conversion factor

C........ Second part depends on scc only

        INTEGER     NCNV2
        INTEGER     SCCSV2( MXSCNV )  !  conversion factor SCC
        REAL        CONVV2( MXSCNV )  !  conversion factor

C.......   Speciation matrix:
        
        REAL        SPCFACS( NPSRC, NMPOL )    !  speciation coefficients

C...........   Area Sources Table
        
        INTEGER      IFIP( NPSRC )  !  source FIPS (county) ID
        INTEGER      ISCC( NPSRC )  !  source SCC
        INTEGER      IPLT( NPSRC )  !  source plant ID
        INTEGER      ISTK( NPSRC )  !  source stack ID

C...........   Other local variables
        
        INTEGER         S, F, I, J, V       !  loop counters.
        INTEGER         CHOICE              !  Choice from menu call
        INTEGER         CDEV                !  for actual-SCC file
        INTEGER         KDEV                !  volatile form conversion file
        INTEGER         SDEV                !  for surrogate coeff file
        INTEGER         XDEV                !  for surrogate xref  file
        INTEGER         LDEV                !  log-device
        
        INTEGER         IOS                 !  I/O status
        INTEGER         IREC                !  input line (record) number
        INTEGER         I1, I2, I3, I4, I5  !  for counting XREFs
        INTEGER         LS,EC

        REAL            FAC

        CHARACTER*5     INVPRO  !  volatile inventory pollutant in spec profile
        CHARACTER*5     POLV    !  volatile inventory pollutant in inventory
        CHARACTER*16    PNAME   !  logical name for point-source      input file
        CHARACTER*16    MNAME   !  logical name for spec matrix output file
        CHARACTER*256   MESG    !  message buffer for M3MESG(), M3EXIT()
        CHARACTER*80    LINE    !  input line from XREF file
        
        LOGICAL         EFLAG   !  error flag
        LOGICAL         KFLAG   !  if volatile format conversion factor or not


C***********************************************************************
C   begin body of program SPCPMAT
        
        LDEV = INIT3()
        
        WRITE( *,92000 ) 
     &  ' ',
     &  'Program SPCPMAT to take the netCDF POINT SOURCES file, the',
     &  'sorted ACTUAL SCC and SPECIATION PROFILES files, and the',
     &  'SPECIATION CROSS-REFERENCE file, and produce the ',
     &  'SPECIATION MATRIX file.',
     &  ' ',
     &  'You will need to enter the logical names for the input and',
     &  'output files (and to have set them prior to program launch,',
     &  'using "setenv <logicalname> <pathname>").  Input files must',
     &  'have been sorted as indicated, prior to program execution.',
     &  ' ',
     &  'You may use END_OF-FILE (control-D) to quit the program',
     &  'during logical-name entry. ',
     &  ' '
        
        IF ( .NOT. GETYN( 'Continue with program?', .TRUE. ) ) THEN
            CALL M3EXIT( 'SPCPMAT', 0, 0, 'Ending program SPCPMAT', 2 )
        END IF


C.......   Get file name; open input point sources file

        PNAME = PROMPTMFILE( 
     &          'Enter logical name for POINT SOURCE input file >> ',
     &          FSREAD3, 'PNTS', 'SPCPMAT' )

        CDEV = PROMPTFFILE( 
     &           'Enter logical name for ACTUAL SCC file >> ',
     &           .TRUE., .TRUE., 'PSCC', 'SPCPMAT' )

        XDEV = PROMPTFFILE( 
     &           'Enter logical name for SPECIATION XREF file >> ',
     &           .TRUE., .TRUE., 'PSREF', 'SPCPMAT' )

        SDEV = PROMPTFFILE( 
     &           'Enter logical name for SPECIATION PROFILES file >> ',
     &           .TRUE., .TRUE., 'PSPRO', 'SPCPMAT' )

C...........   Read in emissions indexing data IFIP, ASC7, ASC3.
        
        IF ( .NOT. READ3( PNAME, 'IFIP', ALLAYS3, 0, 0, IFIP ) ) THEN
            CALL M3EXIT( 'SPCPMAT', 0, 0, 
     &                   'Error reading IFIP from file' // PNAME, 2 )
        END IF

        IF ( .NOT. READ3( PNAME, 'ISCC', ALLAYS3, 0, 0, ISCC ) ) THEN
            CALL M3EXIT( 'SPCPMAT', 0, 0, 
     &                   'Error reading ISCC from file' // PNAME, 2 )
        END IF

        IF ( .NOT. READ3( PNAME, 'IPLANT', ALLAYS3, 0, 0, IPLT ) ) THEN
            CALL M3EXIT( 'SPCPMAT', 0, 0, 
     &                   'Error reading IPLANT from file' // PNAME, 2 )
        END IF

        IF ( .NOT. READ3( PNAME, 'ISTACK', ALLAYS3, 0, 0, ISTK ) ) THEN
            CALL M3EXIT( 'SPCPMAT', 0, 0, 
     &                   'Error reading ISTACK from file' // PNAME, 2 )
        END IF


C.......   Read the ACTUAL-SCC FILE

        WRITE( *,92000 ) ' ', 'Reading ACTUAL SCC file...', ' '

        IREC  =  0
        I     =  0
        EFLAG = .FALSE.
        LID   = -1

22      CONTINUE        !  head of the CDEV-read loop

            READ( CDEV, *, END=33, IOSTAT=IOS ) SCC

            IREC = IREC + 1

            IF ( IOS .NE. 0 ) THEN

                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 
     &              'I/O error', IOS, 
     &              'reading ACTUAL SCC file at line', IREC
                CALL M3MESG( MESG )

            ELSE IF ( I .GE. NPSRC ) THEN 

                EFLAG = .TRUE.
                WRITE( MESG,94010 )                         
     &              'Max SCC table size', NPSRC, 
     &              'exceeded at line', IREC
                CALL M3MESG( MESG )

            ELSE IF ( LID .GE. SCC ) THEN 

                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 
     &              'SCC table out of order at line', IREC
                CALL M3MESG( MESG )

            ELSE

                I   = I + 1
                LID = SCC
                ISCCA( I ) = SCC

            END IF

            GO TO  22   !  to head of loop

33      CONTINUE        !  end of the CDEV-read loop

        IF ( EFLAG ) THEN
            MESG = 'Error reading ACTUAL SCC file'
            CALL M3EXIT( 'SPCPMAT', 0, 0, MESG, 2 )

        END IF

        NSCCA = I


C.......   Read the SPECIATION PROFILES FILE
C...............   ASSUMPTION:  sorted by profile ID

        WRITE( *,92000 ) 
     &       ' ', 'Reading SPECIATION PROFILES file...', ' '

        NPROFS = RDSPROF( SDEV , NMPOL, EMNAM, MXSPRO, 
     &                    SXREF, SCOEF, INVPRO )

        IF ( NPROFS .LE. 0 ) THEN
            CALL M3EXIT( 'SPCMMAT', 0, 0,
     &                   'Error reading SPECIATION PROFILE file.', 2 )
        END IF

C.......   Compare inventory pollutant in the speciation profile with
C.......   volatile pollutant in inventory.  If not the same, prompt for
C.......   a conversion file.

C.......   Prompt for which inventory pollutant is volatile
        
        CHOICE = GETMENU( NIPOL, 1, 
     &                    'Which inventory pollutant is volatile?',
     &                    EINAM )

        POLV = EINAM( CHOICE )

        IF( INVPRO .NE. POLV ) THEN

            MESG = 'Enter logical name for ' //
     &              POLV ( 1:TRIMLEN( POLV ) ) // ' TO ' //
     &              INVPRO ( 1:TRIMLEN( INVPRO ) ) // 
     &             ' CONVERSION file >> '
            KDEV = PROMPTFFILE( 
     &             MESG, .TRUE., .TRUE., 'PSCNV', 'SPCPMAT' )

            KFLAG = .TRUE.
             
        ENDIF

C.......   If needed, read the volatile pollutant conversion file

        IF( KFLAG ) THEN

            CALL RDCONV( KDEV, FIPSV1, SCCSV1, CONVV1, NCNV1,
     &                   SCCSV2, CONVV2, NCNV2 )

        ENDIF

C.......   Read the speciation cross-reference file

        WRITE( *,92000 ) 
     &      ' ', 'Reading SPECIATION XREF file...', ' '

        IREC  = 0
        I     = 0
        EFLAG = .FALSE.

77      CONTINUE        !  head of the XDEV-read loop

            IREC = IREC + 1
            READ( XDEV, 93000, END=88, IOSTAT=IOS ) LINE

            IF ( IOS .NE. 0 ) THEN

                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 
     &              'I/O error ', IOS, 
     &              'reading SPECIATION XREF file at line', IREC
                CALL M3MESG( MESG )

            ENDIF

            LS = LBLANK ( LINE( 17:26 ) )   ! Number of leading spaces
            EC = TRIMLEN( LINE( 17:26 ) )   ! Size without ending blanks

C............. Convert character SCC field to integer SCC number
            IF( LS .EQ. 0 .AND. EC .EQ. 10 ) THEN  ! Area x-ref entry
                GO TO 77

            ELSEIF( LINE(17:17) .GT. '9' ) THEN    ! Motor vehicle x-ref entry
                GO TO 77

            ELSE
                SCC = STR2INT( LINE( 17:26 ) )

            ENDIF


            IF( SCC .GT. 0 ) THEN
                J = FIND1( SCC, NSCCA, ISCCA )  ! Specific SCC

            ELSE
                J = 1                           ! Default SCC

            ENDIF

            IF ( J .LE. 0 ) THEN  ! SCC code not found in SCC LIST file
                GO TO  77

            ELSE IF ( I .GE. NPSRC ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 )
     &                 'Max XREF table size', NPSRC,
     &                 'exceeded at line', IREC
                CALL M3MESG( MESG )

            ELSE
                I = I + 1
                INDEXA( I ) = I
                FIPSFA( I ) = STR2INT ( LINE( 1:5 ) )
                SCCSFA( I ) = SCC
                PLANTA( I ) = STR2INT ( LINE( 7:11 ) )
                STACKA( I ) = STR2INT ( LINE( 13:15 ) )
                SPXRFA( I ) = STR2INT ( LINE( 28:32 ) )
                VREACA( I ) = STR2REAL( LINE( 34:38 ) )
                ETHADA( I ) = STR2REAL( LINE( 55:59 ) )
                FRMADA( I ) = STR2REAL( LINE( 61:65 ) )
                ALLADA( I ) = STR2REAL( LINE( 67:71 ) )

            END IF

        GO TO  77   !  to head of the SDEV-read loop

88      CONTINUE        !  end of the SDEV-read loop

        IF ( EFLAG ) THEN
            MESG = 'Error reading SPECIATION XREF file'
            CALL M3EXIT( 'SPCPMAT', 0, 0, MESG, 2 )
        END IF

        NREFA = I

C.........  Sort and process the speciation cross-reference file

        CALL SORTI4( NREFA, INDEXA, FIPSFA, SCCSFA, PLANTA, STACKA )

        WRITE( *,92000 ) 
     &      ' ', 'Processing SPECIATION XREF file...', ' '

        I1    = 0
        I2    = 0
        I3    = 0
        I4    = 0
        I5    = 0
        EFLAG = .FALSE.

        DO  99  I = 1, NREFA

            J = INDEXA( I )
            FIP    = FIPSFA( J )
            PID    = PLANTA( J )
            SID    = STACKA( J )
            SCC    = SCCSFA( J )

            IF ( FIP .LE. 0  ) THEN     !  scc only

                I5 = I5 + 1
                IF ( I5 .LE. NSREF ) THEN
                    SPXRF5( I5 ) = SPXRFA( J )
                    SCCSF5( I5 ) = SCC
                    VREAC5( I5 ) = VREACA( J )
                    ETHAD5( I5 ) = ETHADA( J )
                    FRMAD5( I5 ) = FRMADA( J )
                    ALLAD5( I5 ) = ALLADA( J )
                END IF

            ELSE IF ( MOD( FIP,1000 ) .LE. 0  ) THEN

                I4 = I4 + 1
                IF ( I4 .LE. NSREF ) THEN
                    SPXRF4( I4 ) = SPXRFA( J )
                    FIPSF4( I4 ) = FIP / 1000   !  state code part only
                    SCCSF4( I4 ) = SCC
                    VREAC4( I4 ) = VREACA( J )
                    ETHAD4( I4 ) = ETHADA( J )
                    FRMAD4( I4 ) = FRMADA( J )
                    ALLAD4( I4 ) = ALLADA( J )
                END IF

            ELSE IF ( PID .LE. 0  ) THEN

                I3 = I3 + 1
                IF ( I3 .LE. NSREF ) THEN
                    SPXRF3( I3 ) = SPXRFA( J ) 
                    FIPSF3( I3 ) = FIP
                    SCCSF3( I3 ) = SCC
                    VREAC3( I3 ) = VREACA( J )
                    ETHAD3( I3 ) = ETHADA( J )
                    FRMAD3( I3 ) = FRMADA( J )
                    ALLAD3( I3 ) = ALLADA( J )
                END IF 

            ELSE IF ( SID .LE. 0  ) THEN

                I2 = I2 + 1
                IF ( I2 .LE. NSREF ) THEN
                    SPXRF2( I2 ) = SPXRFA( J ) 
                    FIPSF2( I2 ) = FIP
                    PLANT2( I2 ) = PID
                    SCCSF2( I2 ) = SCC
                    VREAC2( I2 ) = VREACA( J )
                    ETHAD2( I2 ) = ETHADA( J )
                    FRMAD2( I2 ) = FRMADA( J )
                    ALLAD2( I2 ) = ALLADA( J )
                END IF

            ELSE

                I1 = I1 + 1
                IF ( I1 .LE. NSREF ) THEN
                    SPXRF1( I1 ) = SPXRFA( J ) 
                    FIPSF1( I1 ) = FIP
                    PLANT1( I1 ) = PID
                    STACK1( I1 ) = SID
                    SCCSF1( I1 ) = SCC
                    VREAC1( I1 ) = VREACA( J )
                    ETHAD1( I1 ) = ETHADA( J )
                    FRMAD1( I1 ) = FRMADA( J )
                    ALLAD1( I1 ) = ALLADA( J )
                END IF

            END IF

99      CONTINUE        !  end of the xref-processing loop


        IF ( I1 .GT. NSREF ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 )  
     &          'Max XREF count ', NSREF, 
     &          'exceeded in XREF file:  count=', I1
            CALL M3MESG( MESG )
        ENDIF
        IF ( I2 .GT. NSREF ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 )  
     &          'Max XREF count ', NSREF, 
     &          'exceeded in XREF file:  count=', I2
            CALL M3MESG( MESG )
        END IF
        IF ( I3 .GT. NSREF ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 )  
     &          'Max XREF count ', NSREF, 
     &          'exceeded in XREF file:  count=', I3
            CALL M3MESG( MESG )
        END IF
        IF ( I4 .GT. NSREF ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 )  
     &          'Max XREF count ', NSREF, 
     &          'exceeded in XREF file:  count=', I4
            CALL M3MESG( MESG )
        END IF
        IF ( I5 .GT. NSREF ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 )  
     &          'Max XREF count ', NSREF, 
     &          'exceeded in XREF file:  count=', I5
            CALL M3MESG( MESG )
        END IF

        IF ( EFLAG ) THEN
            MESG = 'Error reading SPECIATION XREF file'
            CALL M3EXIT( 'SPCPMAT', 0, 0, MESG, 2 )
        END IF

        NREF1 = I1
        NREF2 = I2
        NREF3 = I3
        NREF4 = I4
        NREF5 = I5

C.......   Read in description of point-source file
        
        IF ( .NOT. DESC3( PNAME ) ) THEN
            MESG = 'Could not get description of file ' // PNAME
            CALL M3EXIT( 'SPCPMAT', 0, 0, MESG, 2 )

        END IF
        
C.......   Open speciation matrix file:  reuse much of point file description
        
        NVARS3D = NMPOL
        FTYPE3D = GRDDED3
        DO  11  V = 1, NMPOL

            I = INDEX( VNAME3D( V ), '_to_' )
            VNAME3D( V ) = SFNAM( V )
            UNITS3D( V ) = 'moles/gm'
            VDESC3D( V ) = 'Conversion factor, ' // 
     &              SFNAM( V )( 1:I-1 ) // ' to ' // 
     &              SFNAM( V )( I+4:16 )
            VTYPE3D( V ) = M3REAL
11      CONTINUE

        FDESC3D( 1 ) = 
     &     'Matrix of speciation factors, by source, for converting '
        FDESC3D( 2 ) = 
     &     'inventory species into model species for the CB-IV chemical'
        FDESC3D( 3 ) = 
     &    'mechanism.'
        DO  45  V = 4, MXDESC3
            FDESC3D( V ) = ' '
45      CONTINUE

        MNAME = PROMPTMFILE( 
     &    'Enter logical name for SPECIATION MATRIX output file >> ',
     &          FSUNKN3, 'PSMAT', 'SPCPMAT' )


C.......   Compute the speciation coefficient matrix
C...........   ASSUMPTION:  everything sorted by FIP,SCC,PLANT,STACK 
C...........   or by XREF ID.

C...........   Finding the correct profile is a
C...........       Five stage process:
C...........       by      FIP,SCC,plant,stack;
C...........       else by FIP,SCC,plant;
C...........       else by FIP,SCC;
C...........       else by state,SCC;
C...........       else by SCC only.

        EFLAG = .FALSE.
        DO  144  S = 1, NPSRC
                
            FIP = IFIP( S )
            SCC = ISCC( S )
            PID = IPLT( S )
            SID = ISTK( S )

            F = FIND4( FIP, SCC, PID, SID, 
     &                 NREF1, FIPSF1, SCCSF1, PLANT1, STACK1 )

            IF ( F .GT. 0 ) THEN

                XID    = SPXRF1( F )
                ETHADJ = ETHAD1( F )
                FRMADJ = FRMAD1( F )
                ALLADJ = ALLAD1( F )
                VREACT = VREAC1( F )

            ELSE        !  no stack,plant,fip,scc match

                F = FIND3( FIP, SCC, PID,
     &                     NREF2, FIPSF2, SCCSF2, PLANT2 )
                   
                IF ( F .GT. 0 ) THEN
                   
                    XID    = SPXRF2( F )
                    ETHADJ = ETHAD2( F )
                    FRMADJ = FRMAD2( F )
                    ALLADJ = ALLAD2( F )
                    VREACT = VREAC2( F )
                   
                ELSE            !  no plant,fip,scc match
                     
                    F = FIND2( FIP, SCC, 
     &                         NREF3, FIPSF3, SCCSF3 )
                       
                    IF ( F .GT. 0 ) THEN
                       
                        XID    = SPXRF3( F )
                        ETHADJ = ETHAD3( F )
                        FRMADJ = FRMAD3( F )
                        ALLADJ = ALLAD3( F )
                        VREACT = VREAC3( F )
                       
                    ELSE        !  no state+county,scc match
                         
                        F = FIND2( (FIP/1000)*1000, SCC, 
     &                             NREF4, FIPSF4, SCCSF4 )
                           
                        IF ( F .GT. 0 ) THEN
                           
                            XID    = SPXRF4( F )
                            ETHADJ = ETHAD4( F )
                            FRMADJ = FRMAD4( F )
                            ALLADJ = ALLAD4( F )
                            VREACT = VREAC4( F )
                           
                        ELSE    !  no state,scc match
                             
                            F = FIND1( SCC, NREF5, SCCSF5 )
                               
                            IF ( F .GT. 0 ) THEN
                               
                                XID    = SPXRF5( F )
                                ETHADJ = ETHAD5( F )
                                FRMADJ = FRMAD5( F )
                                ALLADJ = ALLAD5( F )
                                VREACT = VREAC5( F )
                               
                            ELSE        !  no scc match:  WARNING
                                
                                WRITE( MESG,94010 )
     &                              'No XREF for FIP', FIP,
     &                              'SCC', SCC, 'plant', PID,
     &                              'stack', SID
                                CALL M3MESG( MESG )

                                WRITE( MESG,94010 ) 
     &                              '     Using default profile 0000'
                                CALL M3MESG( MESG )

                                XID    = 0
                                ETHADJ = 1.0
                                FRMADJ = 0.0
                                ALLADJ = 1.0
                                VREACT = 1.0
                                 
                            END IF      !  if find4:fip,scc,pid,scc succeeded or not
                               
                        END IF      !  if find4:fip,scc,pid,scc succeeded or not
                           
                    END IF      !  if find4:fip,scc,pid,scc succeeded or not
                       
                END IF      !  if find4:fip,scc,pid,scc succeeded or not
                   
            END IF      !  if find4:fip,scc,pid,scc succeeded or not

C...........   Finding the the correct volatile conversion factor is a
C...........       Two stage process:
C...........       by      FIP(state),SCC
C...........       by      SCC only

            VOLADJ = 1.0     ! Initialize volatile conversion factor
    
            IF( KFLAG ) THEN

                F = FIND2( (FIP/1000)*1000, SCC, NCNV1, 
     &                     FIPSV1, SCCSV1 )

                IF ( F .GT. 0 ) THEN

                    VOLADJ = CONVV1( F )

                ELSE                      ! no state,scc match

                    F = FIND1( SCC, NCNV2, SCCSV2 )

                    IF ( F .GT. 0 ) THEN

                        VOLADJ = CONVV2( F )

                    ELSE                  !  no scc match:  WARNING

                        WRITE( MESG,94010 )
     &                   'Default ' // POLV( 1:TRIMLEN( POLV ) ) // 
     &                   ' TO ' // INVPRO ( 1:TRIMLEN( INVPRO ) ) // 
     &                  ' of 1.0 used for FIP= ', FIP, 'SCC=', SCC
                        CALL M3MESG( MESG )

                        CALL M3MESG( MESG )

                    ENDIF  ! scc only match or not
                ENDIF      ! state and scc match or not
            ENDIF          ! KFLAG is true

C...........   Now use this cross-reference index to compute this source's
C...........   entries in the speciation matrix

            I   = FIND1( XID, NPROFS, SXREF )

            IF ( I .LE. 0 ) THEN   !  error:  XREF index not found
           
                WRITE( MESG,94010 ) 
     &              'Profile', XID, 
     &              'not found for FIP', FIP,
     &              'SCC', SCC, 'plant', PID, 'stack', SID
                CALL M3MESG( MESG )

                WRITE( MESG,94010 ) 
     &              '     Using default profile 0000'
                CALL M3MESG( MESG )

                I = FIND1( 0, NPROFS, SXREF )

                IF( I .LE. 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 
     &                  'ERROR: Default profile 0000 not found'
                    CALL M3MESG( MESG )
                    GO TO 144            ! to end of do-loop

                END IF  !  if i<0
            END IF      !  if i<0

            FAC = TON2GM * ALLADJ
            SPCFACS( S, 1 ) = FAC * COFAC     ! 1/molwt for CO
            SPCFACS( S, 2 ) = FAC * NOFAC     ! 90% of 1/molwt for NO2
            SPCFACS( S, 3 ) = FAC * NO2FAC    ! 10% of 1/molwt for NO2

            !  if TCH supported:  use FAC = FAC * VREACT here
            DO  143  V = 4, NMPOL
                SPCFACS( S,V ) = FAC * VOLADJ * SCOEF( I,V )
143         CONTINUE

            SPCFACS( S, 8 ) = FRMADJ * SPCFACS( S, 8 )
            SPCFACS( S,10 ) = ETHADJ * SPCFACS( S,10 )
           
144     CONTINUE        ! end loop on point sources, computing speciation matrix

        IF ( EFLAG ) THEN
            MESG = 'Error computing SPECIATION MATRIX'
            CALL M3EXIT( 'SPCPMAT', 0, 0, MESG, 2 )
        END IF
        
C.......   Write out the speciation matrix:
        
        WRITE( *,92000 ) ' ', 'Writing out SPECIATION MATRIX...', ' '

        IF ( .NOT. WRITE3( MNAME, 'ALL', 0, 0, SPCFACS ) ) THEN
            MESG = 'Error writing SPECIATION MATRIX file'
            CALL M3EXIT( 'SPCPMAT', 0, 0, MESG, 2 )
        END IF



999   CONTINUE          !  spcamat completed:  exit program
     
      MESG = 'Normal completion of SPCPMAT' 
      CALL M3EXIT( 'SPCPMAT', 0, 0, MESG, 0 )


C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )

92010   FORMAT ( 5X , A, :, I10 )


C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

93010   FORMAT( A16 )

93030   FORMAT( I5, I6, I4, 1X, A10, 1X, I5, F6.3, 15X, 3 F6.3 ) ! for ASREF


C...........   Internal buffering formats............ 94xxx

94010   FORMAT ( 10 ( A, :, I10, :, 2X ) )

        END

