C copied by: mhouyoux
C origin: emspoint.F 4.3

         PROGRAM PCPOINT

C***********************************************************************
C  program body starts at line 332
C
C  DESCRIPTION:
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C
C****************************************************************************/
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
C***************************************************************************

        IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   PARAMETERS

        CHARACTER*50    SCCSW
        PARAMETER     ( SCCSW = '@(#)$Id$' )

C...........   EXTERNAL FUNCTIONS and their descriptions:
        
        CHARACTER*2             CRLF
        LOGICAL                 ENVYN
        LOGICAL                 GETYN
        INTEGER                 GETISIZE
        INTEGER                 GETFLINE
        INTEGER                 GETFMTPT
        INTEGER                 GETTZONE
        INTEGER                 INDEX1
        CHARACTER*5             ORDERNAM
        INTEGER                 PROMPTFFILE
        INTEGER                 STR2INT
        INTEGER                 TRIMLEN

        EXTERNAL CRLF, ENVYN, GETYN, GETISIZE, GETFLINE, GETFMTPT,
     &           GETTZONE, INDEX1, ORDERNAM, PROMPTFFILE, , STR2INT,
     &           TRIMLEN

C...........   LOCAL VARIABLES and their descriptions:
C...............   Time zone tables:  FIP-independent; state-only; state-county

        INTEGER                TZONE0
        INTEGER                NZS         !  no of state-specific time zones
        INTEGER                NZF         !  no of FIP-specific time zones
        INTEGER, ALLOCATABLE:: TZONST( : ) !  state-specific zones
        INTEGER, ALLOCATABLE:: TFIPST( : ) !  state FIPS codes (2 digit)
        INTEGER, ALLOCATABLE:: TZONEF( : ) !  fip-specific zones
        INTEGER, ALLOCATABLE:: TFIPEF( : ) !  state/county FIPS codes (5 digit)

C.........  Sorted list of point sources for SMOKE inventory file
        INTEGER                NPSRC         !  actual source count
        REAL   , ALLOCATABLE:: POLALL( :,: ) !  emissions-specific variables

        CHARACTER(LEN=BLRLEN3), ALLOCATABLE:: CBLRID( : ) !  boiler ID
        CHARACTER(LEN=DSCLEN3), ALLOCATABLE:: CPDESC( : ) !  plant description

        CHARACTER(LEN=SRCLEN3), ALLOCATABLE:: CSOURC( : ) !  concatonated source
        CHARACTER(LEN=SRCLEN3)                LSRCCHR     !  previous CSOURC
        CHARACTER(LEN=SRCLEN3)                TSRCCHR     !  tmporary CSOURC

C...........  IDA Output buffers
C...........  NOTE: Each pollutant requires 52 columns
        CHARACTER(LEN=248), ALLOCATABLE:: CHRBUF( : )   ! IDA characteristics
        CHARACTER(LEN=52 ), ALLOCATABLE:: POLBUF( :,: ) ! IDA pol-specific

C...........  Inventory pollutants actually in the inventory
        INTEGER                               NIPOL      ! Actual no of inv pols
        INTEGER               , ALLOCATABLE:: EIIDX( : ) ! pos in full inven arr
        CHARACTER(LEN=IOVLEN3), ALLOCATABLE:: EINAM( : ) ! Name of actual pols
        CHARACTER(LEN=IOVLEN3)                CPOL       ! Tmp pollutant code

C...........  Variable names in matrices per inventory pollutant
        CHARACTER(LEN=IOVLEN3), ALLOCATABLE:: CPVNAM( :,: )

C...........  Matrix variables not matching any inventory pollutants
        CHARACTER(LEN=IOVLEN3) CVODDLST( MXVARS3D ) 

C...........   File units and logical/physical names

        INTEGER         DDEV    !  Unit number for output IDA formatted file
        INTEGER         LDEV    !  log-device
        INTEGER         SDEV    !  for ASCII input inventory file

        CHARACTER(LEN=NAMLEN3) CNAMEA( MXCMAT )! unsorted control/proj matrices
        CHARACTER(LEN=NAMLEN3) CNAME ( MXCMAT )! sorted   control/proj matrices
        CHARACTER(LEN=NAMLEN3) ENAME !  emissions input inventory logical name
        CHARACTER(LEN=NAMLEN3) MNAME !  temporary control/proj matrix name
        CHARACTER(LEN=NAMLEN3) ONAME !  emissions output inventory logical name

C...........   Other local variables

        REAL            EMISI   !  Temporary inverse emission value
        REAL            EMISN   !  Temporary new emission value
        REAL            EMISO   !  Temporary old emission value
                                
        INTEGER         S, I, J, K, L, LK, LS, V !  counters and indices
        INTEGER         L1, L2           !  counters and indices
        INTEGER         FIP     ! Temporary FIPS code

        INTEGER         INVFMT  !  Inventory format code
        INTEGER         IOS     !  I/O status
        INTEGER         NFIPLIN !  Number of lines in ZDEV
        INTEGER         NLINE   !  Number of lines for a file
        INTEGER         TZONE   !  tmp time zone

        LOGICAL         CFLAG   !  velocity recalc: TRUE iff VELOC_RECALC = Y
        LOGICAL         DFLAG   !  input verification:  TRUE iff ERROR
        LOGICAL      :: EFLAG = .FALSE.  !  TRUE iff ERROR
        LOGICAL         SFLAG   !  input verification:  report missing species
        LOGICAL         TFLAG   !  TRUE if temporal x-ref output
        LOGICAL         WFLAG   !  input verification:  convert bad lat-lons

        CHARACTER*5     TPOLPOS !  Temporary pollutant position
        CHARACTER*300   MESG    !  text for M3EXIT()

        CHARACTER*16 :: PROGNAME= 'PCPOINT' !  program name

C***********************************************************************
C   begin body of program PCPOINT

        LDEV = INIT3()

C.........  Write out copywrite, version, web address, header info, and prompt
C           to continue running the program.
        CALL INITEM( LDEV, SCCSW, PROGNAME )

C.........  Get environment variables that control the program
        PROMPTF = ENVYN ( 'PROMPTFLAG', 'Prompt for inputs or not',
     &                    .TRUE., IOS )

C.........  Set default settings, getting from environment if PROMPTF is false
        IF( PROMPTF ) THEN
            NCMAT = 1

        ELSE
            NCMAT = ENVINT( 'NUM_PT_CTLMAT', 
     &                      'Number of control and projection matrices',
     &                      1, IOS )

            IF( IOS .NE. 0 ) THEN
                WRITE( MESG,94010 ) 'WARNING: Environment variable ' //
     &                 'NUM_PT_CTLMAT is not set or is invalid.' //
     &                 CRLF() // BLANK16 // 'Default value of', NCMAT,
     &                 'will be used.'
                CALL M3MSG2( MESG )
            ENDIF

        ENDIF

C.........  Prompt for and open I/O API inventory input file
        ENAME= PROMPTMFILE( 
     &       'Enter logical name for the POINT I/O API INVENTORY file',
     &       FSREAD3, 'PNTS', PROGNAME )

C.........  Prompt for and open ASCII inventory input file
        SDEV= PROMPTFFILE( 
     &      'Enter logical name for the POINT ASCII INVENTORY file',
     &      .TRUE., .TRUE., 'PSRC', PROGNAME )

C.........  Read NetCDF header and store info
        IF( .NOT. DESC3( ENAME ) ) THEN
            MESG = 'Could not read description for "' //
     &             ENAME( 1:TRIMLEN( ENAME ) ) // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
 
C.........  Store number of sources, number of non-pollutant variables, number
C           of variables per pollutant, and pollutant names
        ELSE
            NPSRC  = NROWS3D
            NPVARS = NVARS3D
            NNONPV = GETIFDSC( FDESC3, '/NON POLLUTANT/' )
            NVPERP = GETIFDSC( FDESC3, '/PER POLLUTANT/' )

C.............  Allocate memory for inventory pollutant names
            NIPOL = ( NPVARS - NNONPV ) / NVPERP
            ALLOCATE( EINAM( NIPOL ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EINAM', PROGNAME )

C.............  Store pollutant names
            I = 0 
            DO V = NNONPV + 1, NPVARS, NVPERP
                I = I + 1
                EINAM( I ) = VNAME3D( V )
            ENDDO
 
        END IF

C.........  Quick check number of sources in SDEV file
        I = GETFLINE( SDEV, 'Point ASCII inventory file')
        I = I - 1 ! Adjust for header
        IF( I .LE. 0 ) THEN
            MESG = 'ERROR: Empty ASCII inventory file!'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        ELSEIF( I .NE. NPSRC ) THEN
            WRITE( MESG,94010 ) 'ERROR: Mismatched inventory files!' //
     &             CRLF() // BLANK16 // 
     &             'Number of sources in NetCDF file:', NPSRC,
     &             CRLF() // BLANK16 // 
     &             'Number of sources in ASCII  file:', I
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        ENDIF

C.........  Get number of control and/or projection matrices
        MESG = 'Enter number of control and/or projection matrices'
        NCMAT = GETNUM( 0, 10000, NCMAT, MESG )

C.........  Allocate memory based on number of control/projection matrices
        ALLOCATE( CNAMEA( NCMAT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CNAMEA', PROGNAME )
        ALLOCATE( CTYPEA( NCMAT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CTYPEA', PROGNAME )
        ALLOCATE( CINDXA( NCMAT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CINDXA', PROGNAME )
        ALLOCATE( CCNTRA( NCMAT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CCNTRA', PROGNAME )
        ALLOCATE( CNAME( NCMAT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CNAME', PROGNAME )
        ALLOCATE( CTYPE( NCMAT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CTYPE', PROGNAME )
        ALLOCATE( NCPVARS( NCMAT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'NCPVARS', PROGNAME )
        ALLOCATE( CPVNAMS( MXVARS3, NCMAT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CPVNAMS', PROGNAME )

C.........  Loop through potential control and projection matrices
        DO I = 1, NCMAT 

C.............  Generate default name
            WRITE( MNAME, '(A5,I2.2)' ) 'PCMAT', I

C.............  Prompt for name 
            CBUF5 = ORDERNAM( I )
            L = TRIMLEN( CBUF5 )
            WRITE( MESG, 94010 ) 'Enter logical name for ' //
     &             CBUF5( 1:L ) // ' CONTROL/PROJECTION MATRIX file'

            MNAME = PROMPTMFILE( MESG, FSREAD3, MNAME, PROGNAME )

            IF( .NOT. DESC3( MNAME ) ) THEN
                EFLAG = .TRUE.
                MESG = 'Could not read description for "' //
     &                 MNAME( 1:TRIMLEN( MNAME ) ) // '"'
                CALL M3MSG2( MESG )
            ENDIF

            CNAMEA( I ) = MNAME  ! Store name in unsorted list
            CTYPEA( I ) = GETIFDSC( FDESC3, '/CTYPE/' )
            CINDXA( I ) = I
            CCNTRA( I ) = I

C.............  Compare number of sources to NPSRC 
            BUFF1 = CBUF5( 1:L ) // ' control/projection file'
            BUFF2 = 'Point inventory'
            CALL CHKISIZ( MNAME, BUFF1, BUFF2, NPSRC, IOS )
            IF( IOS .NE. 0 ) EFLAG = .TRUE.

C.............  Initialize array
            DO J = 1, NIPOL
                CPVNAM( J,I ) = EMCMISS3
            ENDDO

C.............  Store number of variables in current matrix
            NCPVAR( I ) = NVARS3D

C.............  Interpret variable names and compare to pollutant list.  Keep
C               track of those that are not in the inventory file.
            N = 0
            M = 0
            DO V = 1, NVARS3D
 
                CPVNAMS( V,I ) = VNAME3D( V )

                K = CHKCPVAR( VNAME3D( V ), EINAM, NIPOL ) 

                IF( K .EQ. 0 ) THEN   ! Count matrices with "all"
                    M = M + 1

                ELSEIF( K .EQ. -1 ) THEN   ! Store names don't match with EINAM
                    N = N + 1
                    CVODDLST( N ) = VNAME3D( V )

                ELSEIF( K .EQ. -2 )    ! Error
                    EFLAG = .TRUE.

                ELSEIF( K .NE. -3 )    ! Bad value, -3 lets us skip variable
                    WRITE( MESG,94010 ) 'Invalid value "', K,
     &                     '" returned from subroutine CHKCPVAR'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

                ENDIF

            ENDDO

C.............  Report pollutant-specific factors that do not match pollutants
C               in inventory. Do inside matrix loop so can report by matrix.
            IF( N .GT. 0 ) THEN
                WRITE( MESG,94010 ) 'In matrix', I, 'some factors ' //
     &               'do not apply to any pollutants in the inventory:'
                CALL M3MESG( MESG )

                DO L = 1, N
                    MESG = BLANK5 // CVODDLST( L )
                    CALL M3MESG( MESG )
                ENDDO
            ENDIF

        ENDDO  ! End loop on matrices, I

        NCMAT = I - 1  ! Store final number of matrices
        NCALL = M

C.........  Abort if error occurred
        IF( EFLAG ) THEN
            MESG = 'ERROR: Problem with control or projection matrices'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C.........  Give warning if no control or projection matrices entered
        ELSEIF( NCMAT .EQ. 0 ) THEN
            MESG = 'WARNING: No control or projection matrices. '
            CALL M3MSG2( MESG )

        ENDIF
        
C.........  Sort control matrices in order of precedence, and sort sorted list.
C           The sort must make sure that when two matrices of the same type
C           are present, the input order is maintained.
        CALL SORTI2( NCMAT, CINDXA, CTYPEA, CCNTRA )

        DO I = 1, NCMAT
            J = CINDXA( I )
            CNAME( I ) = CNAMEA( J )
            CTYPE( I ) = CTYPEA( J )
        ENDDO

C.........  Allocate memory for SMOKE inventory arrays
        ALLOCATE( CHARVAL( NPSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CHARVAL', PROGNAME )  
        ALLOCATE( POLALL( NPSRC, NVPERP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'POLALL', PROGNAME )
        ALLOCATE( CBLRID( NPSRC ), STAT=IOS )

C.........  Allocate memory for IDA buffer arrays:
        ALLOCATE( CHRBUF( NPSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CHRBUF', PROGNAME )  
        ALLOCATE( POLBUF( NPSRC,NIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'POLBUF', PROGNAME )  

C.........  Allocate memory for control arrays
        ALLOCATE( CFACALL( NPSRC, NCALL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CFACALL', PROGNAME )  
        ALLOCATE( CFAC( NPSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CFAC', PROGNAME )  
        ALLOCATE( IDXALL( NCMAT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IDXALL', PROGNAME )  

C.........  Re-read header for input inventory file
        IF( .NOT. DESC3( ENAME ) ) THEN
            MESG = 'Could not read description for "' //
     &             ENAME( 1:TRIMLEN( ENAME ) ) // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ENDIF

C.........  Prompt file output SMOKE file (or NONE)
        MESG = 'Enter logical name for the SMOKE I/O API OUTPUT ' //
               'file or "NONE"'
        ONAME = PROMPTMFILE( MESG, FSUNKN3, 'PNTSOUT', PROGNAME )

C.........  Prompt file output IDA file (or NONE)
        MESG = 'Enter logical name for the IDA-FORMATTED OUTPUT ' //
               'file or "NONE"'

        DDEV = PROMPTFFILE( MESG, .FALSE., .TRUE., 'PTIDA', PROGNAME )

C.........  Read in control matrix variables that are "all".  First, store
C           position of data in storage array, then read.
        J = 0
        DO I = 1, NCMAT

            K = INDEX1( 'all', NCPVAR(I), CPVNAM(1,I) )
            IF( K .GT. 0 ) THEN
                J = J + 1
                IDXALL( I ) = J
                IF( .NOT. READ3( CNAME( I ), 'all', ALLAYS3, 
     &                           0, 0, CFACALL( 1,J ) ) ) THEN
                    MESG = 'ERROR: Could not read variable "all" ' //
     &                     'from file.'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                ENDIF

            ELSE
                IDXALL( I ) = 0
            ENDIF

        ENDDO ! End loop on matrices

C.........  For SMOKE output, copy all non-pollutant variables from I/O API
        IF( ONAME .NE. 'NONE' ) THEN

C.............  Loop through I/O API non-emissions pollutants and write to 
C               output file
            DO V = 1, NNONPV

                VARBUF = VNAME3D( V )
                CALL READWR3( ENAME, ONAME, VARBUF, ALLAYS3,
     &                        SDATE3D, STIME3D, VTYPE3D(V), NPSRC, IOS )
 
                IF( IOS .NE. 0 ) EFLAG = .TRUE.

            ENDDO ! End of loop on non-pollutant variables

        END IF
                
C.........  Abort if error occurred
        IF( EFLAG ) THEN
            MESG = 'ERROR: Problem with read or write of I/O API files.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ENDIF

C.........  For SMOKE or IDA output, process emission values and write
C           pollutant-based variables
        IF( ONAME .NE. 'NONE' .OR. DDEV .GT. 0 ) THEN

C.............  Loop through inventory pollutants
            DO J = 1, NIPOL

                VARBUF = EINAM( J )                     

C.................  Read inventory pollutant emissions from I/O API file
                IS = NNONPV + ( J-1 ) * NVPERP + 1
                CALL RDINVPOL( ENAME, NPSRC, NVPERP, VNAME3D( IS ), 
     &                         POLALL, IOS ) 

                IF( IOS .GT. 0 ) EFLAG = .TRUE.

C.................  Loop through matrices, if any
                DO I = 1, NCMAT

C.....................  Search for pollutant name in list for this matrix
                    K = INDEX1( VARBUF, NCPVAR(I), CPVNAM(1,I) )

C.....................  Read in pollutant-specific array 
                    IF( K .GT. 0 ) THEN

                        IF( .NOT. READ3( CNAME( I ), VARBUF, 
     &                                   ALLAYS3, 0, 0, CFAC ) ) THEN
                            MESG = 'ERROR: Could not read ' //
     &                             'variable ' //
     &                             VARBUF( 1:TRIMLEN( VARBUF ) ) //
     &                             ' from file.'
                            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                        ENDIF

                    ENDIF

C.....................  Set default factors
                    ALLFAC = 1.
                    PSFAC  = 1.

C.....................  Change operation based on type of matrix, then
C                       loop through sources, applying factors as needed
                    IF( CTYPE(I) .EQ. CTYPPROJ .AND. K .GT. 0 ) THEN! Projection

                        DO S = 1, NPSRC
                            POLLALL(S,1) = POLLALL(S,1) * CFAC(S)
                        ENDDO

                    ELSEIF( CTYPE(I) .EQ. CTYPREAC ) THEN ! Reactivity 

                        DO S = 1, NPSRC
                            POLLALL(S,1) = POLLALL(S,1) * CFAC(S)
                        ENDDO

                    ELSEIF( CTYPE(I) .EQ. CTYPADD .AND. K .GT. 0 ) THEN ! Add

                        DO S = 1, NPSRC
                            POLLALL(S,1) = POLLALL(S,1) + CFAC(S)
                        ENDDO
                             
                    ELSEIF( CTYPE(I) .EQ. CTYPMULT ) THEN ! Multiply (standard)

                        DO S = 1, NPSRC
                            J = INDXALL( I )
                            IF( J .GT. 0 ) ALLFAC = CFACALL( S,J )
                            IF( K .GT. 0 )  PSFAC = CFAC( S )

                            POLLALL(S,1) = POLLALL(S,1)*ALLFAC*CFAC
                        ENDDO

                    ELSEIF( CTYPE(I) .EQ. CTYPOVER .AND. K .GT. 0 ) THEN!Overwrt

                        DO S = 1, NPSRC
                            POLLALL(S,1) = CFAC(S)
                        ENDDO

                    ENDIF
                        
                ENDDO  ! End loop on control/projection matrices

C.................  Write out pollutant-based variables to SMOKE file
                IF( ONAME .NE. 'NONE' ) THEN
                    CALL WRINVPOL( ONAME, NPSRC, NVPERP, VNAME3D( IS ),
     &                             POLALL, IOS )

                ENDIF

                IF( IOS .GT. 0 ) EFLAG = .TRUE.

C.................  Write out pollutant-based variables to buffer for IDA
C.................  The buffer approach is being used so that all of the 
C                   pollutants do not have to be stored in memory at the same
C                   time.
                IF( DDEV .GT. 0 ) THEN
                    
                    CALL WRIDAPOL( 'POINT', NPSRC, POLBUF( 1,J ), 
     &                             NVPERP, POLALL, IOS )
                ENDIF

                IF( IOS .GT. 0 ) EFLAG = .TRUE.

            ENDDO  ! End loop on inventory pollutants
            
        ENDIF     ! End of section for both SMOKE and IDA output

        IF( EFLAG ) THEN
            MESG = 'ERROR: Could not read and write all pollutant-' //
     &             'specific variables'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ENDIF

C.........  For IDA output
        IF( DDEV .GT. 0 ) THEN

            CALL WRIDAOUT( 'POINT', ENAME, SDEV, NPSRC, NIPOL, 
     &                     POLBUF, IOS )

            IF( IOS .GT. 0 ) EFLAG = .TRUE.

        ENDIF
      
        IF( EFLAG ) THEN
            MESG = 'ERROR: Could not read and write all source ' //
     &             'characteristics for IDA output.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ENDIF

C.........  End program successfully
        MESG = ' '
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 0 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )

92010   FORMAT( 5X, A, :, I10 )


C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

93500   FORMAT( A248, <NIPOL>(A52) )


C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94020   FORMAT( A, 1X, I5.5, 1X, A, 1X, I8.8, 1X,
     &          A, I6, X, A, I6, X, A, :, I6 )

94040   FORMAT( A, I2.2 )

94060   FORMAT( 10( A, :, E10.3, :, 1X ) )

94080   FORMAT( '************  ', A, I7, ' ,  ' , A, I12 )
 
        END

