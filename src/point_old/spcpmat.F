
C Version "@(#)$Id$ $Source$ $Date$ 

        PROGRAM SPCPMAT

C***********************************************************************
C  program body starts at line 182
C
C  DESCRIPTION:
C       Construct point source chemical speciation matrix from data 
C       contained in EPS-style chem-split and chem-xref files.
C       Uses fixed NOtoNOX  factor of 0.0197  moles  NO/gm NOX (from EPS)
C       Uses fixed NO2toNOX factor of 0.00217 moles NO2/gm NOX  "
C
C  PRECONDITIONS REQUIRED:
C       Sorted, cut-down input data for chem-split and chem-xref files;
C       M3IO point-sources file from program GRDAMAT
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       Models-3 I/O
C       FIND1, FIND2, FIND3, GETEFILE, GETYN, TRIMLEN
C
C  REVISION  HISTORY:
C       Prototype  3/95 by CJC.
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
        INTEGER         PROMPTFFILE
        CHARACTER*16    PROMPTMFILE
        INTEGER         STR2INT
        REAL            STR2REAL
        INTEGER         TRIMLEN
        
        EXTERNAL        FIND1, FIND2, FIND3, FIND4, GETYN, 
     &                  PROMPTFFILE, PROMPTMFILE, 
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
                                 
C...........   Speciation profiles table
        
        REAL        FACS( NMPOL )    !  input record:  speciation coefficients
        REAL        VREACT
        REAL        ETHADJ
        REAL        FRMADJ
        REAL        ALLADJ
        
        INTEGER     NPROFS                      !  actual number of profiles
        INTEGER     SXREF( NSPRP )              !  profile indexes (sorted)
        REAL        SCOEF( NSPRP, NMPOL )       ! profile coefficients


C...........   Speciation cross-reference table:  
C.......   Five parts:   first part depends on scc,fip,plant,stack

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


C.......   Speciation matrix:
        
        REAL        SPCFACS( NPSRC, NMPOL )    !  speciation coefficients


C...........   Area Sources Table
        
        INTEGER      IFIP( NPSRC )  !  source FIPS (county) ID
        INTEGER      ISCC( NPSRC )  !  source SCC
        INTEGER      IPLT( NPSRC )  !  source plant ID
        INTEGER      ISTK( NPSRC )  !  source stack ID


C...........   Other local variables
        
        INTEGER         S, F, I, V          !  loop counters.
        INTEGER         SDEV                !  for surrogate coeff file
        INTEGER         XDEV                !  for surrogate xref  file
        INTEGER         LDEV                !  log-device
        
        INTEGER         IOS                 !  I/O status
        INTEGER         IREC                !  input line (record) number
        INTEGER         I1, I2, I3, I4, I5  !  for counting XREFs

        REAL            FAC

        CHARACTER*16    PNAME   !  logical name for point-source      input file
        CHARACTER*16    MNAME   !  logical name for spec matrix output file
        CHARACTER*120   MESG    !  message buffer for M3ERR()
        CHARACTER*80    LINE    !  input line from XREF file


C***********************************************************************
C   begin body of program SPCPMAT
        
        LDEV = INIT3()
        
        WRITE( *,92000 ) 
     &  ' ',
     &  'Program SPCPMAT to take the netCDF POINT SOURCES file, the',
     &  'sorted SPECIATION PROFILES file, and the sorted SPECIATION',
     &  'CROSS-REFERENCE file, and produce the SPECIATION MATRIX ',
     &  'file.',
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

        XDEV = PROMPTFFILE( 
     &           'Enter logical name for SPECIATION XREF file >> ',
     &           .TRUE., .TRUE., 'PSREF', 'SPCPMAT' )

        SDEV = PROMPTFFILE( 
     &           'Enter logical name for SPECIATION PROFILES file >> ',
     &           .TRUE., .TRUE., 'PSPRO', 'SPCPMAT' )

C.......   Read in description of point-source file
        
        IF ( .NOT. DESC3( PNAME ) ) THEN
            CALL M3ERR( 'SPCPMAT', 0, 0, 
     &                  'Could not get description of file ' // PNAME,
     &                  .TRUE. )
        END IF
        
        
C.......   Open speciation matrix file:  reuse much of point file description
        
        NVARS3D = NMPOL
        FTYPE3D = GRDDED3
        DO  44  V = 1, NMPOL

            I = INDEX( VNAME3D( V ), '_to_' )
            VNAME3D( V ) = SFNAM( V )
            UNITS3D( V ) = 'moles/gm'
            VDESC3D( V ) = 'Conversion factor, ' // 
     &              SFNAM( V )( 1:I-1 ) // ' to ' // 
     &              SFNAM( V )( I+4:16 )
            VTYPE3D( V ) = M3REAL
44      CONTINUE

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


C...........   Read in emissions indexing data IFIP, ASC7, ASC3.
        
        IF ( .NOT. READ3( PNAME, 'IFIP', ALLAYS3, 0, 0, IFIP ) ) THEN
            CALL M3ERR( 'SPCPMAT', 0, 0, 
     &                  'Error reading IFIP from file' // PNAME,
     &                  .TRUE. )
        END IF

        IF ( .NOT. READ3( PNAME, 'ISCC', ALLAYS3, 0, 0, ISCC ) ) THEN
            CALL M3ERR( 'SPCPMAT', 0, 0, 
     &                  'Error reading ISCC from file' // PNAME,
     &                  .TRUE. )
        END IF

        IF ( .NOT. READ3( PNAME, 'IPLANT', ALLAYS3, 0, 0, IPLT ) ) THEN
            CALL M3ERR( 'SPCPMAT', 0, 0, 
     &                  'Error reading IPLANT from file' // PNAME,
     &                  .TRUE. )
        END IF

        IF ( .NOT. READ3( PNAME, 'ISTACK', ALLAYS3, 0, 0, ISTK ) ) THEN
            CALL M3ERR( 'SPCPMAT', 0, 0, 
     &                  'Error reading ISTACK from file' // PNAME,
     &                  .TRUE. )
        END IF


C.......   Read the SPECIATION PROFILES FILE
C...............   ASSUMPTION:  sorted by FIP, ASC
        
        WRITE( *,92000 ) 
     &       ' ', 'Reading SPECIATION PROFILES file...', ' '

        IREC  =  0
        
101     CONTINUE        !  head of the SDEV-read loop
        
            READ( SDEV, *, END=122, IOSTAT=IOS ) 
     &          SID, ( FACS( V ), V = 4, NMPOL )

            IREC = IREC + 1
            IF ( IOS .NE. 0 ) THEN
                WRITE( *,92010 ) 
     &              'Unit # ', SDEV, 
     &              'IOSTAT ', IOS, 
     &              'Line # ', IREC, ' '
                CALL M3ERR( 'SPCPMAT', 0, 0, 
     &              'I/O error reading SPECIATION PROFILE file.',
     &              .TRUE. )
            ELSE IF ( IREC .GT. NSPRP ) THEN 
                WRITE( *,92010 ) 
     &              'Record number  ', IREC,
     &              'Dimensioned max', NSPRP, ' '
                CALL M3ERR( 'SPCPMAT', 0, 0, 
     &              'Max speciation profile table size NSPRP exceeded',
     &              .TRUE. )
            END IF
            
            SXREF( IREC ) = SID
            DO  111  V = 4, NMPOL
                SCOEF( IREC,V ) = FACS( V )
111         CONTINUE

            GO TO  101  !  to head of loop
      
122     CONTINUE        !  end of the SDEV-read loop

        NPROFS = IREC



C.......   Read the speciation cross-reference file
C...............   ASSUMPTION:  sorted by FIP and ASC (i.e., <PID,ID3> pair)

        WRITE( *,92000 ) 
     &      ' ', 'Reading SPECIATION XREF file...', ' '

        IREC = 0
        I1   = 0
        I2   = 0
        I3   = 0
        I4   = 0
        I5   = 0

133     CONTINUE        !  head of the XDEV-read loop
        
            READ( XDEV, 93000, END=135, IOSTAT=IOS ) LINE
            IREC = IREC + 1

            IF ( IOS .NE. 0 ) THEN
            
                WRITE( MESG,94050 ) 
     &              'Error ', IOS, 
     &              'reading SPECIATION XREF file at line', IREC
                CALL M3EXIT( 'SPCPMAT', 0, 0, MESG, 2 )

            END IF
            
            FIP    = STR2INT ( LINE(  1: 5 ) )
            PID    = STR2INT ( LINE(  7:11 ) )
            SID    = STR2INT ( LINE( 13:15 ) )
            SCC    = STR2INT ( LINE( 17:26 ) )
            XID    = STR2INT ( LINE( 28:32 ) )
            VREACT = STR2REAL( LINE( 34:38 ) )
            ETHADJ = STR2REAL( LINE( 55:59 ) )
            FRMADJ = STR2REAL( LINE( 61:65 ) )
            ALLADJ = STR2REAL( LINE( 67:71 ) )
            
            IF ( ETHADJ .LT. AMISS3 ) ETHADJ = 1.0
            IF ( FRMADJ .LT. AMISS3 ) FRMADJ = 1.0
            IF ( ALLADJ .LT. AMISS3 ) ALLADJ = 1.0
            
            IF ( FIP .EQ. 0  ) THEN     !  scc only

                I5 = I5 + 1
                IF ( I5 .LE. NSREF ) THEN
                    SPXRF5( I5 ) = XID 
                    SCCSF5( I5 ) = SCC
                    VREAC5( I5 ) = VREACT
                    ETHAD5( I5 ) = ETHADJ
                    FRMAD5( I5 ) = FRMADJ
                    ALLAD5( I5 ) = ALLADJ
                END IF

            ELSE IF ( MOD( FIP,1000 ) .EQ. 0  ) THEN
                
                I4 = I4 + 1
                IF ( I4 .LE. NSREF ) THEN
                    SPXRF4( I4 ) = XID 
                    FIPSF4( I4 ) = FIP / 1000   !  state code part only
                    SCCSF4( I4 ) = SCC
                    VREAC4( I4 ) = VREACT
                    ETHAD4( I4 ) = ETHADJ
                    FRMAD4( I4 ) = FRMADJ
                    ALLAD4( I4 ) = ALLADJ
                END IF

            ELSE IF ( PID .LE. 0  ) THEN

                I3 = I3 + 1
                IF ( I3 .LE. NSREF ) THEN
                    SPXRF3( I3 ) = XID 
                    FIPSF3( I3 ) = FIP
                    SCCSF3( I3 ) = SCC
                    VREAC3( I3 ) = VREACT
                    ETHAD3( I3 ) = ETHADJ
                    FRMAD3( I3 ) = FRMADJ
                    ALLAD3( I3 ) = ALLADJ
                END IF

            ELSE IF ( SID .LE. 0  ) THEN

                I2 = I2 + 1
                IF ( I2 .LE. NSREF ) THEN
                    SPXRF2( I2 ) = XID 
                    FIPSF2( I2 ) = FIP
                    PLANT2( I2 ) = PID
                    SCCSF2( I2 ) = SCC
                    VREAC2( I2 ) = VREACT
                    ETHAD2( I2 ) = ETHADJ
                    FRMAD2( I2 ) = FRMADJ
                    ALLAD2( I2 ) = ALLADJ
                END IF

            ELSE

                I1 = I1 + 1
                IF ( I1 .LE. NSREF ) THEN
                    SPXRF1( I1 ) = XID 
                    FIPSF1( I1 ) = FIP
                    PLANT1( I1 ) = PID
                    STACK1( I1 ) = SID
                    SCCSF1( I1 ) = SCC
                    VREAC1( I1 ) = VREACT
                    ETHAD1( I1 ) = ETHADJ
                    FRMAD1( I1 ) = FRMADJ
                    ALLAD1( I1 ) = ALLADJ
                END IF

            END IF

            GO TO  133	!  to head of the SDEV-read loop

135     CONTINUE        !  end of the SDEV-read loop


        IF      ( I1 .GT. NSREF ) THEN
            WRITE( MESG,94050 )  
     &          'Max XREF count ', NSREF, 
     &          'exceeded in XREF file:  count=', I1
            CALL M3EXIT( 'SPCPMAT', 0, 0, MESG, 2 )
        ELSE IF ( I2 .GT. NSREF ) THEN
            WRITE( MESG,94050 )  
     &          'Max XREF count ', NSREF, 
     &          'exceeded in XREF file:  count=', I2
            CALL M3EXIT( 'SPCPMAT', 0, 0, MESG, 2 )
        ELSE IF ( I3 .GT. NSREF ) THEN
            WRITE( MESG,94050 )  
     &          'Max XREF count ', NSREF, 
     &          'exceeded in XREF file:  count=', I3
            CALL M3EXIT( 'SPCPMAT', 0, 0, MESG, 2 )
        ELSE IF ( I4 .GT. NSREF ) THEN
            WRITE( MESG,94050 )  
     &          'Max XREF count ', NSREF, 
     &          'exceeded in XREF file:  count=', I4
            CALL M3EXIT( 'SPCPMAT', 0, 0, MESG, 2 )
        ELSE IF ( I5 .GT. NSREF ) THEN
            WRITE( MESG,94050 )  
     &          'Max XREF count ', NSREF, 
     &          'exceeded in XREF file:  count=', I5
            CALL M3EXIT( 'SPCPMAT', 0, 0, MESG, 2 )
        END IF
        
        NREF1 = I1
        NREF2 = I2
        NREF3 = I3
        NREF4 = I4
        NREF5 = I5
        
        
C.......   Compute the speciation coefficient matrix

        DO  144  S = 1, NPSRC
                
            FIP = IFIP( S )
            SCC = ISCC( S )
            PID = IPLT( S )
            SID = ISTK( S )

C...........   Find the relevant cross-reference table entry:            
C...........   ASSUMPTION:  everything sorted by FIP,SCC,PLANT,STACK 
C...........   or by XREF ID.
C...........   Five stage process:
C...........       by      FIP,SCC,plant,stack;
C...........       else by FIP,SCC,plant;
C...........       else by FIP,SCC;
C...........       else by state,SCC;
C...........       else by SCC only.

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
                         
                        F = FIND2( FIP/1000, SCC, 
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
                               
                            ELSE        !  no scc match:  ERROR
                                
                                WRITE( MESG,94050 )
     &                              'No XREF for FIP', FIP,
     &                              'SCC', SCC,
     &                              'plant', PID,
     &                              'stack', SID
                                CALL M3EXIT( 'SPCPMAT', 0,0, MESG, 2 )
                                 
                            END IF      !  if find4:fip,scc,pid,scc succeeded or not
                               
                        END IF      !  if find4:fip,scc,pid,scc succeeded or not
                           
                    END IF      !  if find4:fip,scc,pid,scc succeeded or not
                       
                END IF      !  if find4:fip,scc,pid,scc succeeded or not
                   
            END IF      !  if find4:fip,scc,pid,scc succeeded or not


C...........   Now use this cross-reference index to compute this source's
C...........   entries in the speciation matrix

            I   = FIND1( XID, NPROFS, SXREF )

            IF ( I .GT. 0 ) THEN
           
                FAC = TON2GM * ALLADJ
                SPCFACS( S, 1 ) = FAC * COFAC     ! 1/molwt for CO
                SPCFACS( S, 2 ) = FAC * NOFAC     ! 90% of 1/molwt for NO2
                SPCFACS( S, 3 ) = FAC * NO2FAC    ! 10% of 1/molwt for NO2

                !  if TCH supported:  use FAC = FAC * VREACT here
                DO  143  V = 4, NMPOL
                    SPCFACS( S,V ) = FAC * SCOEF( I,V )
143             CONTINUE

                SPCFACS( S, 8 ) = FRMADJ * SPCFACS( S, 8 )
                SPCFACS( S,10 ) = ETHADJ * SPCFACS( S,10 )
           
            ELSE    !  error:  XREF index not found
           
                WRITE( MESG,94050 ) 
     &              'XREF index ', XID, 
     &              'not found for FIP', FIP,
     &              'SCC', SCC,
     &              'plant', PID,
     &              'stack', SID
                CALL M3EXIT( 'SPCPMAT', 0,0, MESG, 2 )
           
            END IF  !  if i>0, or not

144     CONTINUE        ! end loop on point sources, computing speciation matrix

        
C.......   Write out the speciation matrix:
        
        WRITE( *,92000 ) ' ', 'Writing out SPECIATION MATRIX...', ' '

        IF ( .NOT. WRITE3( MNAME, 'ALL', 0, 0, SPCFACS ) ) THEN
            CALL M3EXIT( 'SPCPMAT', 0, 0, 
     &                   'Error writing SPECIATION MATRIX file.', 2 )
        END IF



999   CONTINUE          !  spcamat completed:  exit program
      
      CALL M3EXIT( 'SPCPMAT', 0, 0, 'Normal completion.', 0 )


C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )

92010   FORMAT ( 5X , A, :, I10 )


C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

93010   FORMAT( A16 )


C...........   Internal buffering formats............ 94xxx

94050   FORMAT( 5( A, :, I9, :, 2X ) )


        END

