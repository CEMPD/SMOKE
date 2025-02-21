
        SUBROUTINE BLDREPIDX( SLNAME, SSNAME )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      The BLDREPIDX sets the data-reading and data-aggregating indices
C      for all reports and for the Smkreport program in general.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created 7/2000 by M Houyoux
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
C***********************************************************************

C...........   MODULES for public variables
C.........  This module contains Smkreport-specific settings
        USE MODREPRT, ONLY: TFLAG, PRRPTFLG, NREPORT, SLFLAG, SSFLAG,
     &                      MXOUTDAT, MXINDAT, OUTDNAM, RPT_, ALLRPT,
     &                      INDNAM, SPCPOL, NSPCPOL,
     &                      SPFLAG, SKFLAG

C.........  This module contains report arrays for each output bin
        USE MODREPBN, ONLY: NSVARS, LV1, LV2, LV3, TODOUT, ETPNAM,
     &                      DATNAM, INVIDX, TPRIDX, INVTOPRJ, INVTOCMU,
     &                      TOSOUT, SPCNAM, ETPSPCNAM, PRCSPCNAM,
     &                      SUMETPNAM, SUMPOLNAM, SUMSPCNAM, SPCOUT, 
     &                      SPCTOINV, SPCTOTPR, SPCIDX, TPACTIDX, 
     &                      SLUNIT, SSUNIT, NMSPC, EMNAM, NSPCIN

C.........  This module contains the temporal profile tables
        USE MODTMPRL, ONLY: NTPDAT, TPNAME, TPDESC

C.........  This module contains the control packet data and control matrices
        USE MODCNTRL, ONLY: NVPROJ, PNAMPROJ, NVCMULT, PNAMMULT

C.........  This module contains the information about the source category
        USE MODINFO, ONLY: NIPPA, NIACT, ACTVTY, EANAM

C.........  This module is required for the FileSetAPI
        USE MODFILESET

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'SETDECL.EXT'   !  FileSetAPI function declarations

C...........   EXTERNAL FUNCTIONS and their descriptions
        CHARACTER(2)  CRLF
        INTEGER       INDEX1

        EXTERNAL   CRLF, INDEX1

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN) :: SLNAME     ! mole-based spec matrix name
        CHARACTER(*), INTENT (IN) :: SSNAME     ! mass-based spec matrix name

C...........   Local paramaters
        CHARACTER(7), PARAMETER :: RPRTPRE = 'ProjRep'   ! Check prjctn prefix
        CHARACTER(8), PARAMETER :: DIFFPRE = 'ProjDiff'  ! Check prjctn diff

C...........   Sorting index for creating unique species list
        INTEGER  :: SRTIDX( NSVARS )

C...........   Other local variables
        INTEGER          E, I, I1, I2, I3, J, K, L, L2, N, V ! counters and indices

        INTEGER          IDX     !  temporary index
        INTEGER          IOS     !  i/o status
        INTEGER          LS      !  length of speciation name joiner
        INTEGER          LT      !  length of emission type name joiner
        INTEGER          MXDATALL!  maximum input data for memory allocation
        INTEGER       :: NDATA=0  !  tmp number of data variables per report
        INTEGER       :: NDATALL=0 !  no. of all input data
        INTEGER          NPOL    !  no. of pollutants

        LOGICAL       :: SKIP    = .FALSE. !  true: skipping selecting output species 
        LOGICAL       :: ANYOUT  = .FALSE. !  true: data select will be output
        LOGICAL       :: EFLAG   = .FALSE. !  true: error found
        LOGICAL       :: SFLAG   = .FALSE. !  true: speciation

        CHARACTER(300)   MESG              !  message buffer

        CHARACTER(IOVLEN3), ALLOCATABLE :: POLNAM( : )     !  tmp pol name

        CHARACTER(LV1)      ABUF       !  tmp activity
        CHARACTER(LV2)      EBUF       !  tmp emission type
        CHARACTER(LV1)      PBUF       !  previous species
        CHARACTER(LV1)      SBUF       !  tmp species
        CHARACTER(LV3)      VBUF       !  tmp speciation variable name

        CHARACTER(16) :: PROGNAME = 'BLDREPIDX' ! program name

C***********************************************************************
C   begin body of subroutine BLDREPIDX

C.........  Set the maximum number of output variables, depending on whether 
C           hourly data are in use for any reports
C.........  Also, increase further if any projection or control checks
C           are being run
        MXDATALL = NIPPA
        IF( TFLAG ) MXDATALL = NIPPA + NTPDAT
        IF( PRRPTFLG ) MXDATALL = MXDATALL * 3

C.........  Allocate memory for the data indexing and labeling arrays...

C.........  Pollutant/activity/emission type arrays
        ALLOCATE( TODOUT( MXDATALL, NREPORT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TODOUT', PROGNAME )
        ALLOCATE( ETPNAM( MXDATALL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ETPNAM', PROGNAME )
        ALLOCATE( DATNAM( MXDATALL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'DATNAM', PROGNAME )
        ALLOCATE( INVIDX( MXDATALL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INVIDX', PROGNAME )
        ALLOCATE( TPRIDX( MXDATALL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TPRIDX', PROGNAME )
        ALLOCATE( INVTOPRJ( MXDATALL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INVTOPRJ', PROGNAME )
        ALLOCATE( INVTOCMU( MXDATALL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INVTOCMU', PROGNAME )

c NOTE: will have to adjust MXDATALL to be larger when a CHECK instruction
C   N: is included

        TODOUT%ETP   = 0        ! array
        TODOUT%DAT   = 0        ! array
        TODOUT%AGG   = 0        ! array
        TODOUT%SPCYN = .FALSE.  ! array
        TODOUT%PRYN  = .FALSE.  ! array
        TODOUT%CUYN  = .FALSE.  ! array
        TODOUT%CRYN  = .FALSE.  ! array
        ETPNAM       = ' '      ! array
        DATNAM       = ' '      ! array
        INVIDX       = 0        ! array
        TPRIDX       = 0        ! array
        INVTOPRJ     = 0        ! array
        INVTOCMU     = 0        ! array

C.........  Speciation variable arrays
        ALLOCATE( TOSOUT( NSVARS, NREPORT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TOSOUT', PROGNAME )
        ALLOCATE( SPCNAM( NSVARS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SPCNAM', PROGNAME )
        ALLOCATE( POLNAM( NSVARS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'POLNAM', PROGNAME )
        ALLOCATE( ETPSPCNAM( NSVARS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ETPSPCNAM', PROGNAME )
        ALLOCATE( PRCSPCNAM( NSVARS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'PRCSPCNAM', PROGNAME )
        ALLOCATE( SUMETPNAM( NSVARS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SUMETPNAM', PROGNAME )
        ALLOCATE( SUMPOLNAM( NSVARS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SUMPOLNAM', PROGNAME )
        ALLOCATE( SUMSPCNAM( NSVARS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SUMSPCNAM', PROGNAME )
        ALLOCATE( SPCOUT( NSVARS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SPCOUT', PROGNAME )
        ALLOCATE( SPCTOINV( NSVARS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SPCTOINV', PROGNAME )
        ALLOCATE( SPCTOTPR( NSVARS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SPCTOTPR', PROGNAME )
        ALLOCATE( SPCIDX( NSVARS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SPCIDX', PROGNAME )
        ALLOCATE( SLUNIT( NSVARS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SLUNIT', PROGNAME )
        ALLOCATE( SSUNIT( NSVARS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SSUNIT', PROGNAME )
        ALLOCATE( SKFLAG( NSVARS,NREPORT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SKFLAG', PROGNAME )
 
        TOSOUT%SPC    = 0       ! array
        TOSOUT%ETPSPC = 0       ! array
        TOSOUT%PRCSPC = 0       ! array
        TOSOUT%SUMETP = 0       ! array
        TOSOUT%SUMPOL = 0       ! array
        TOSOUT%SUMSPC = 0       ! array
        TOSOUT%AGG    = 0       ! array
        SPCNAM        = ' '     ! array
        POLNAM        = ' '     ! array
        ETPSPCNAM     = ' '     ! array
        PRCSPCNAM     = ' '     ! array
        SUMETPNAM     = ' '     ! array
        SUMPOLNAM     = ' '     ! array
        SUMSPCNAM     = ' '     ! array
        SPCOUT        = .FALSE. ! array
        SPCTOINV      = 0       ! array
        SPCTOTPR      = 0       ! array
        SPCIDX        = 0       ! array
        SKFLAG        = .FALSE. ! array

C.........  Set the length of the emission type joiner
        LT = LEN_TRIM( ETJOIN )
        LS = LEN_TRIM( SPJOIN )
C.........  Populate the emission type list and the pollutant list from the
C           reporting bins module
C.........  Will add variables when a CHECK instruction is included
c NOTE: Will have to index differently so that the CHECK variables can be
C    N: added in during the loop.
        N = 0
        DO V = 1, NIPPA

            EBUF = EANAM( V )

            K = INDEX( EBUF, ETJOIN )   ! Look for emission type joiner

            N = N + 1
            IF( K .GT. 0 ) THEN         ! Store emission type and pol from it
                L2 = LEN_TRIM( EBUF )
                ETPNAM( N ) = EBUF

C..............  See if pollutant projection is being checked and add
C                for all pollutants (since that's how it's implemented
C                in Cntlmat)
                IF( PRRPTFLG ) THEN
                    N = N + 1
                    ETPNAM( N ) = RPRTPRE // EBUF

                    N = N + 1
                    ETPNAM( N ) = DIFFPRE // EBUF
                ENDIF

            ELSE                        ! Store pollutant only
                DATNAM( N ) = EBUF
                
C..............  See if pollutant projection is being checked and add
C                for all pollutants (since that's how it's implemented
C                in Cntlmat)
                IF( PRRPTFLG ) THEN
                    N = N + 1
                    DATNAM( N ) = RPRTPRE // EBUF
                    N = N + 1
                    DATNAM( N ) = DIFFPRE // EBUF
                ENDIF

            END IF

C.............  Store index from "all" list to inventory file variables
            INVIDX( N ) = V

        END DO

        NDATALL = N

C.........  If temporal allocated used during program...
        IF( TFLAG ) THEN

C.............  Allocate memory for activity index
            ALLOCATE( TPACTIDX( NTPDAT ), STAT=IOS )
            CALL CHECKMEM( IOS, 'TPACTIDX', PROGNAME )
            TPACTIDX = 0

C.............  Loop through hourly emissions data variables
            N = NIPPA            
            DO V = 1, NTPDAT

                N = N + 1
                EBUF = TPNAME( V )

                K  = INDEX( EBUF, ETJOIN )   ! Look for emission type joiner

C.................  Store emission type and/or pollutant
                IF( K .GT. 0 ) THEN
                    L2 = LEN_TRIM( EBUF )
                    ETPNAM( N ) = EBUF
                    DATNAM( N ) = EBUF( K+LT:L2 )

                ELSE
                    DATNAM( N ) = EBUF

                END IF

C.................  Store index from "all" list to hourly file variables
                TPRIDX( N ) = V

C.................  Determine if any variables from hourly file are associated 
C                   with activities in the inventory file.
                ABUF = ' '
                L2   = LEN_TRIM( TPDESC( V ) )

                K = INDEX( TPDESC( V ), 'from' )
                IF( K .GT. 0 ) ABUF = ADJUSTL( TPDESC( V )( K+4:L2 ) )

                J = INDEX1( ABUF, NIACT, ACTVTY )
                IF( J .GT. 0 ) THEN
                    TPACTIDX( V ) = J
                END IF

            END DO

            NDATALL = N

        END IF

C.........  If speciation is used during program...
        SFLAG = ( SLFLAG .OR. SSFLAG )

        IF( SFLAG ) THEN

C.............  Initialize sorting index for species names
            SRTIDX = 0     ! array

C.............  Get header of mole speciation matrix
            IF( SLFLAG ) THEN

                IF ( .NOT. DESCSET( SLNAME, ALLFILES ) ) THEN

                    MESG = 'Could not get description of file "' //
     &                     SLNAME( 1:LEN_TRIM( SLNAME ) ) // '"'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

                ENDIF

C.................  Store units
                CALL STORE_VUNITSET( 1, 1, NVARSET, SLUNIT )

            END IF

C.............  Get header of mass speciation matrix
            IF( SSFLAG ) THEN

                IF ( .NOT. DESCSET( SSNAME, ALLFILES ) ) THEN

                    MESG = 'Could not get description of file "' //
     &                     SSNAME( 1:LEN_TRIM( SSNAME ) ) // '"'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

                ENDIF

C.................  Store units
                CALL STORE_VUNITSET( 1, 1, NVARSET, SSUNIT )

            END IF

C.............  Populate the speciation variable lists
C.............  Populate the speciation units lists
C.............  Count the number of unique species
            NMSPC = 0
            DO V = 1, NSVARS

                VBUF = VDESCSET( V )

                J  = INDEX( VBUF, ETJOIN )   ! Look for emission type joiner
                K  = INDEX( VBUF, SPJOIN )   ! Look for speciation joiner
                L2 = LEN_TRIM( VBUF )

                IF( K .LE. 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 'ERROR: Speciation joiner "'//
     &                     SPJOIN( 1:LS ) // '" is not found for ' //
     &                     'speciation variable', V
                    CALL M3MSG2( MESG )
                    CYCLE 

                END IF

                SBUF = VBUF( K+LS:L2 )
                EBUF = VBUF( 1:K-1 )

C.................  Find pollutant or emission type in list and store index
                I1 = INDEX1( EBUF, NDATALL, DATNAM  ) ! Look in data
                I2 = INDEX1( EBUF, NDATALL, ETPNAM )  ! Look in emission types
                IF( TFLAG ) THEN
                    I3 = INDEX1( EBUF, NTPDAT, TPNAME )   ! Look in hourly file
                ELSE
                    I3 = 0
                END IF

C.................  For species based on inventory pollutant, set index
                IF( I1 .GT. 0 ) THEN

                    SPCTOINV( V ) = MAX( I1,0 )

C.................  For species based on activity, set index and write note
                ELSE IF( I2 .GT. 0 .AND. I3 .GT. 0 ) THEN
                    IF ( TPACTIDX( I3 ) .GT. 0 ) THEN

                        SPCTOTPR( V ) = I3 + NIPPA

                        I2 = INDEX1( SBUF, V, SPCNAM )
                        IF( I2 .LE. 0 ) THEN
                            L2 = LEN_TRIM( SBUF )
                            MESG = 'NOTE: Species "' // SBUF(1:L2) // 
     &                             '" created based on activity data.'
                            CALL M3MESG( MESG )
                        END IF
                    END IF

C.................  For emission type not based on activity
                ELSE IF( I2 .GT. 0 ) THEN
                    SPCTOINV( V ) = I2

C.................  Otherwise, skip species 
                ELSE
                    MESG = 'WARNING: Speciation variable "' //
     &                     VBUF( 1:L2 ) // '" does not match any ' //
     &                     CRLF()// BLANK10// 'data in emissions' // 
     &                     'input file. Species will be skipped.'
                    CALL M3MSG2( MESG )
                    CYCLE

                END IF

C.................  Update index for species-to-temporal pollutant
                IF( I3 .GT. 0 ) THEN
                    SPCTOTPR( V ) = I3 + NIPPA
                ENDIF

C.................  Data variable is emission type
                IF( J .GT. 0 ) THEN

                    SPCNAM   ( V )= SBUF
                    ETPSPCNAM( V )= EBUF
                    PRCSPCNAM( V )= 'S-'// VBUF( 1:J-1 )// ETJOIN// SBUF
                    SUMETPNAM( V )= 'S-'// EBUF
                    SUMPOLNAM( V )= 'I-'// VBUF( J+LT:K-1 )
                    SUMSPCNAM( V )= 'S-'// TRIM( SBUF )

C.................  No emission type
                ELSE

                    SPCNAM   ( V ) = SBUF
                    SUMPOLNAM( V ) = 'I-' // EBUF
                    SUMSPCNAM( V ) = 'S-' // TRIM( SBUF )

                END IF

                SRTIDX( V ) = V
                POLNAM( V ) = EBUF

C.................  Count unique species.  Look for this species in previous 
C                   list and add one to count of not found.
                K = INDEX1( SBUF, V-1, SPCNAM )
                IF( K .LE. 0 ) NMSPC = NMSPC + 1

            END DO  ! End of loop for speciation variables

C.............  Allocate memory for unique species list
            ALLOCATE( EMNAM( NMSPC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EMNAM', PROGNAME )
            EMNAM = ' '
        
C.............  Create sorted list of unique species names
            CALL SORTIC( NSVARS, SRTIDX, SPCNAM )

            PBUF = ' '
            K = 0
            DO V = 1, NSVARS
                J = SRTIDX( V )
 
C.................  Skip any zeroes in the sorted index due to unused species
                IF( J == 0 ) CYCLE
 
                SBUF = SPCNAM( J )

                IF( SBUF .NE. PBUF ) THEN
                    K = K + 1
                    EMNAM( K ) = SBUF
                END IF

                PBUF = SBUF

            END DO  ! End creating list of unique species

C.............  If species name in SPCNAM is the same as a pollutant name, 
C               then set it to blank
            DO V = 1, NSVARS
                K = INDEX1( SPCNAM( V ), NIPPA, EANAM )
                IF( K .GT. 0 )  THEN
                    SPCNAM( V ) = ' '
                END IF
            END DO  ! End creating list of unique species

        END IF      ! End if speciation

C.........  Determine if any reports did not have a DATA instruction
        I = 0  ! added by GAP 1/17/07
        DO V = 1, NREPORT
            IF( ALLRPT( V )%NUMDATA < 0 ) THEN
                I = ALLRPT( V )%NUMDATA
                EXIT            
            END IF
        END DO
!        I = MINVAL( ALLRPT%NUMDATA )

C.........  If a report did not have a DATA instruction, then maximum output
C           variables depends on whether we have speciation or not
        MXOUTDAT = MXINDAT
        IF( I .LT. 0 ) THEN
            MXOUTDAT = MAX( MXOUTDAT, NDATALL )
            IF( SFLAG ) MXOUTDAT = MAX( MXOUTDAT, NDATALL + NMSPC )
        END IF

C.........  Allocate memory of reading and output data names
C.........  UNC-IE: Dec 2023: Increase allocation of OUTDNAM
c       ALLOCATE( OUTDNAM( MXOUTDAT, NREPORT ), STAT=IOS )
        ALLOCATE( OUTDNAM(max(MXOUTDAT,NSVARS+1), NREPORT ),
     &            STAT=IOS )
        CALL CHECKMEM( IOS, 'OUTDNAM', PROGNAME )
        OUTDNAM = ' '

C.........  Go through reports, and assign indices to from data to output 
C           columns
        DO N = 1, NREPORT

            RPT_ = ALLRPT( N )
            NDATA = RPT_%NUMDATA

C.............  If no data values specified, then output all available 
C               pol/act/emission types and species.  This is the default
C               behavior of the reports.
            IF( NDATA .LT. 0 ) THEN

C.................  Set the default output data depending on hourly inputs or not
                IF( RPT_%USEHOUR ) THEN
                    J = NTPDAT
                    NPOL = NTPDAT
                    OUTDNAM( 1:NTPDAT,N ) = TPNAME( 1:NTPDAT )
                ELSE
                    J = NIPPA
                    NPOL = NIPPA
                    OUTDNAM( 1:NIPPA,N ) = EANAM( 1:NIPPA )
                END IF

C.................  If species are used for this report, set them for output
C.................  If species is same as pollutant, then add prefix
                IF( RPT_%USESSMAT .OR.
     &              RPT_%USESLMAT      ) THEN

                    DO V = 1, NMSPC
                        J = J + 1
                        SBUF = EMNAM( V )

                        K = INDEX1( SBUF, NDATALL, DATNAM )

                        IF( K .GT. 0 ) THEN
                            IF (SPFLAG) THEN
                              OUTDNAM( J,N ) = SBUF
                            ELSE
                              OUTDNAM( J,N ) = 'S-' // SBUF
                            END IF
                            I2 = INDEX1( SBUF, NPOL, OUTDNAM( 1:NPOL,N )) 
                            IF ( I2 .GT. 0) THEN
                                IF (SPFLAG) SKFLAG(I2,N) = .TRUE. ! UNC-IE H.Tran: Mark this inventory species for skipping
                                OUTDNAM( I2,N ) =  'I-' // SBUF
                            END IF
                        ELSE
                            OUTDNAM( J,N ) = SBUF
                        END IF
                    END DO

                END IF

                ALLRPT( N )%NUMDATA = J
                NDATA = J

C.............  Otherwise, set output data values as input data values
            ELSE

                IF( RPT_%BYSPC ) THEN   ! BY SPCCODE poll+ add associated species....

                    NDATA = 1      ! add default firt pol name
                    DO V = 1, NSVARS
                        IF( INDNAM( 1,N ) == POLNAM( V ) ) NDATA = NDATA + 1
                    END DO
C                    NSPCPOL = NDATA
                    ALLRPT( N )%NUMDATA = NDATA

C                    DEALLOCATE( OUTDNAM, SPCPOL )
C   UNC-IE Dec 2023: Comment out the three lines below this comment lines 
C                    deallocate OUTDNAM could lead to array out of bound issue
C                    when there are multiple configurations of reported species in REPCONFIG
c                   DEALLOCATE( OUTDNAM )
c                   ALLOCATE( OUTDNAM( NDATA, NREPORT ), STAT=IOS )
c                   CALL CHECKMEM( IOS, 'OUTDNAM', PROGNAME )
C                    ALLOCATE( SPCPOL( NDATA ), STAT=IOS )
C                    CALL CHECKMEM( IOS, 'SPCPOL', PROGNAME )

                    L = 1
                    OUTDNAM( 1,N )     = INDNAM( 1,N )
                    ALLRPT( N )%SPCPOL = INDNAM( 1,N )
                    DO V = 1, NSVARS
                        IF( INDNAM( 1,N ) == POLNAM( V ) ) THEN
                            L = L + 1
                            OUTDNAM( L,N ) = SPCNAM( V ) 
                        END IF
                    END DO
C                    SPCPOL = OUTDNAM( :,N )

                ELSE

                    OUTDNAM( 1:NDATA,N ) = INDNAM( 1:NDATA,N )

                END IF

            END IF

C.............  Loop through requested data for this report
            DO I = 1, NDATA

                SKIP   = .FALSE.
                ANYOUT = .FALSE.

C.................  Loop through names of emission types, emis, or act
                DO E = 1, NDATALL

C.....................  Skip "all" list entry if index is inventory and
C                       report is temporal or if index is temporal and
C                       report is inventory
                    IF(      RPT_%USEHOUR .AND. INVIDX(E) .GT. 0 ) CYCLE
                    IF(.NOT. RPT_%USEHOUR .AND. TPRIDX(E) .GT. 0 ) CYCLE

C.....................  To emission type column
C.....................  Only when emission type is not based on activity
                    IF( OUTDNAM( I,N ) .EQ. ETPNAM( E ) ) THEN
                        TODOUT( E,N )%ETP = I
                        TODOUT( E,N )%AGG = 1
                        ANYOUT = .TRUE.
                        SKIP   = .TRUE.
                    END IF

C.....................  To pollutant/activity column
                    IF( OUTDNAM( I,N ) .EQ. DATNAM( E ) .OR.  
     &                      OUTDNAM( I,N ) .EQ. ('I-'//DATNAM( E ) )) THEN
                        TODOUT( E,N )%DAT = I
                        TODOUT( E,N )%AGG = 1
                        ANYOUT = .TRUE.
                        SKIP   = .TRUE.
                    END IF

C..................  If projection matrix applies to this report
                    IF( RPT_%USEPRMAT ) THEN
                        IF( PNAMPROJ( 1 ) .EQ. 'pfac' ) THEN
                            TODOUT( E,N )%PRYN = .TRUE.
                            INVTOPRJ( E ) = 1
                        ELSE
                            IDX = INDEX1( DATNAM(E), NVPROJ, PNAMPROJ )
                            IF( IDX .GT. 0 ) THEN
                                TODOUT( E,N )%PRYN = .TRUE.
                                INVTOPRJ( E ) = IDX
                            END IF
                        END IF
                    END IF

C..................  If mult. control matrix applies to this report
                    IF( RPT_%USECUMAT ) THEN
                        IF( PNAMMULT( 1 ) .EQ. 'all' ) THEN
                            TODOUT( E,N )%CUYN = .TRUE.
                            INVTOCMU( E ) = 1
                        ELSE
                            IDX = INDEX1( DATNAM(E), NVCMULT, PNAMMULT )
                            IF( IDX .GT. 0 ) THEN
                                TODOUT( E,N )%CUYN = .TRUE.
                                INVTOCMU( E ) = IDX
                            END IF
                        END IF
                    END IF

                END DO  ! End loop over emission types, emis, or act

                IF( SKIP ) CYCLE  ! skip once output species is chosen

C.................  If speciation for current report
                IF( RPT_%USESLMAT .OR.
     &              RPT_%USESSMAT      ) THEN

C.....................  Loop through names of speciation variables
                    DO V = 1, NSVARS

C.........................  Get species to inventory data index
                        E = SPCTOINV( V )

C.....................  Set inventory variable as having speciation
                       IF( E .GT. 0 )
     &                      TODOUT( E,N )%SPCYN = ( RPT_%USESLMAT .OR.
     &                                              RPT_%USESSMAT    )

C.........................  Get species to hourly data index
                        E = SPCTOTPR( V )

C......................  Set hourly variable as having speciation
                        IF( E .GT. 0 ) 
     &                      TODOUT( E,N )%SPCYN= ( RPT_%USESLMAT .OR.
     &                                             RPT_%USESSMAT    )

C.........................  Skip species that don't match inventory data
C.........................  To species column
                        IF( OUTDNAM( I,N ) .EQ. SPCNAM( V ) ) THEN
                            TOSOUT( V,N )%SPC = I
                            TOSOUT( V,N )%AGG = 1
                            IF( .NOT. SPCOUT( V ) ) NSPCIN = NSPCIN + 1
                            SPCOUT( V ) = .TRUE.
                            ANYOUT = .TRUE.
                        END IF

C.........................  To emission-type/species column
                        IF( OUTDNAM( I,N ) .EQ. ETPSPCNAM( V ) ) THEN 
                            TOSOUT( V,N )%ETPSPC = I
                            TOSOUT( V,N )%AGG = 1
                            IF( .NOT. SPCOUT( V ) ) NSPCIN = NSPCIN + 1
                            SPCOUT( V ) = .TRUE.
                            ANYOUT = .TRUE.
                        END IF

C.........................  To process/species column
                        IF( OUTDNAM( I,N ) .EQ. PRCSPCNAM( V ) ) THEN 
                            TOSOUT( V,N )%PRCSPC = I
                            TOSOUT( V,N )%AGG = 1
                            IF( .NOT. SPCOUT( V ) ) NSPCIN = NSPCIN + 1
                            SPCOUT( V ) = .TRUE.
                            ANYOUT = .TRUE.

C.........................  To post-speciation summed emission type column
                        ELSEIF( OUTDNAM( I,N ) .EQ. SUMETPNAM(V) ) THEN 
                            TOSOUT( V,N )%SUMETP = I
                            TOSOUT( V,N )%AGG = 1
                            IF( .NOT. SPCOUT( V ) ) NSPCIN = NSPCIN + 1
                            SPCOUT( V ) = .TRUE.
                            ANYOUT = .TRUE.
                        END IF

C.........................  To post-speciation summed pollutant column
C.........................  Also for species with same name as pollutants
                        IF( OUTDNAM( I,N ) .EQ. SUMPOLNAM( V ) ) THEN
                            TOSOUT( V,N )%SUMPOL = I
                            TOSOUT( V,N )%AGG = 1
                            IF( .NOT. SPCOUT( V ) ) NSPCIN = NSPCIN + 1
                            SPCOUT( V ) = .TRUE.
                            ANYOUT = .TRUE.

C.........................  To post-speciation summed species column
                        ELSEIF( OUTDNAM( I,N ) .EQ. SUMSPCNAM(V) ) THEN 
                            TOSOUT( V,N )%SUMSPC = I
                            TOSOUT( V,N )%AGG = 1
                            IF( .NOT. SPCOUT( V ) ) NSPCIN = NSPCIN + 1
                            SPCOUT( V ) = .TRUE.
                            ANYOUT = .TRUE.
                        END IF

                    END DO   ! End loop over speciation variables

                END IF

C.................  Give warning if no matches
                IF( .NOT. ANYOUT ) THEN
                    L = LEN_TRIM( OUTDNAM( I,N ) )
                    WRITE( MESG,94010 ) 'WARNING: Skipping requested '//
     &                     'output data "'// OUTDNAM( I,N )( 1:L ) // 
     &                     '" for report', N, CRLF() // BLANK10 //
     &                     'because no match found with inputs.'
                    CALL M3MSG2( MESG )
                END IF

            END DO       ! End loop over selected data

        END DO           ! End loop over reports

C.........  Set global-to-input index for speciation factors
        K = 0
        DO V = 1, NSVARS

            IF( SPCOUT( V ) ) THEN
                K = K + 1
                SPCIDX( V ) = K
            END IF

        END DO

C.........  If there was any error, exit 
        IF( EFLAG ) THEN
             MESG = 'Problem setting up input data to output columns.'
             CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I10, :, 1X ) )

C******************  INTERNAL SUBPROGRAMS  *****************************
 
        CONTAINS
 
C.............  This subprogram stores I/O API NetCDF variable units from
C               a file set into a local array based on indices in subprogram call.
            SUBROUTINE STORE_VUNITSET( ISTART, INCRMT, NUNIT, UNITS )

C.............  Subprogram arguments
            INTEGER      ISTART        ! starting position in VUNITSET of names
            INTEGER      INCRMT        ! increment of VUNITSET for names
            INTEGER      NUNIT         ! number of units
            CHARACTER(*) UNITS( NUNIT )! stored variable units

C.............  Local variables
            INTEGER  I, J, L

C----------------------------------------------------------------------

            UNITS = ' '

            J = ISTART
            DO I = 1, NUNIT

                L = LEN_TRIM( VUNITSET( J ) )
                UNITS( I ) = VUNITSET( J )( 1:L )
                J = J + INCRMT

            END DO
 
            END SUBROUTINE STORE_VUNITSET

        END SUBROUTINE BLDREPIDX

