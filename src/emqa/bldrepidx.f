
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
C COPYRIGHT (C) 2002, MCNC Environmental Modeling Center
C All Rights Reserved
C  
C See file COPYRIGHT for conditions of use.
C  
C Environmental Modeling Center
C MCNC
C P.O. Box 12889
C Research Triangle Park, NC  27709-2889
C  
C smoke@emc.mcnc.org
C  
C Pathname: $Source$
C Last updated: $Date$ 
C  
C***********************************************************************

C...........   MODULES for public variables
C.........  This module contains Smkreport-specific settings
        USE MODREPRT

C.........  This module contains report arrays for each output bin
        USE MODREPBN

C.........  This module contains the temporal profile tables
        USE MODTMPRL

C.........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS and their descriptions
        CHARACTER*2   CRLF
        INTEGER       INDEX1

        EXTERNAL   CRLF, INDEX1

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN) :: SLNAME     ! mole-based spec matrix name
        CHARACTER(*), INTENT (IN) :: SSNAME     ! mass-based spec matrix name

C...........   Sorting index for creating unique species list
        INTEGER  :: SRTIDX( NSVARS )

C...........   Other local variables
        INTEGER          E, I, I1, I2, I3, J, K, L, L2, N, V ! counters and indices

        INTEGER          IOS     !  i/o status
        INTEGER          LS      !  length of speciation name joiner
        INTEGER          LT      !  length of emission type name joiner
        INTEGER          MXDATALL!  maximum input data for memory allocation
        INTEGER          NDATA   !  tmp number of data variables per report
        INTEGER          NDATALL !  no. of all input data

        LOGICAL       :: ANYOUT  = .FALSE. !  true: data select will be output
        LOGICAL       :: EFLAG   = .FALSE. !  true: error found
        LOGICAL       :: SFLAG   = .FALSE. !  true: speciation

        CHARACTER*300    MESG              !  message buffer

        CHARACTER(LEN=LV1)      ABUF       !  tmp activity
        CHARACTER(LEN=LV2)      EBUF       !  tmp emission type
        CHARACTER(LEN=LV1)      PBUF       !  previous species
        CHARACTER(LEN=LV1)      SBUF       !  tmp species
        CHARACTER(LEN=LV3)      VBUF       !  tmp speciation variable name

        CHARACTER*16 :: PROGNAME = 'BLDREPIDX' ! program name

C***********************************************************************
C   begin body of subroutine BLDREPIDX

C.........  Set the maximum number of outputs variables, depending on whether 
C           hourly data are in use for any reports
        MXDATALL = NIPPA
        IF( TFLAG ) MXDATALL = NIPPA + NTPDAT

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

        TODOUT%ETP   = 0        ! array
        TODOUT%DAT   = 0        ! array
        TODOUT%AGG   = 0        ! array
        TODOUT%SPCYN = .FALSE.  ! array
        ETPNAM       = ' '      ! array
        DATNAM       = ' '      ! array
        INVIDX       = 0        ! array
        TPRIDX       = 0        ! array

C.........  Speciation variable arrays
        ALLOCATE( TOSOUT( NSVARS, NREPORT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TOSOUT', PROGNAME )
        ALLOCATE( SPCNAM( NSVARS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SPCNAM', PROGNAME )
        ALLOCATE( ETPSPCNAM( NSVARS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ETPSPCNAM', PROGNAME )
        ALLOCATE( PRCSPCNAM( NSVARS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'PRCSPCNAM', PROGNAME )
        ALLOCATE( SUMETPNAM( NSVARS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SUMETPNAM', PROGNAME )
        ALLOCATE( SUMPOLNAM( NSVARS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SUMPOLNAM', PROGNAME )
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
 
        TOSOUT%SPC    = 0       ! array
        TOSOUT%ETPSPC = 0       ! array
        TOSOUT%PRCSPC = 0       ! array
        TOSOUT%SUMETP = 0       ! array
        TOSOUT%SUMPOL = 0       ! array
        TOSOUT%AGG    = 0       ! array
        SPCNAM        = ' '     ! array
        ETPSPCNAM     = ' '     ! array
        PRCSPCNAM     = ' '     ! array
        SUMETPNAM     = ' '     ! array
        SUMPOLNAM     = ' '     ! array
        SPCOUT        = .FALSE. ! array
        SPCTOINV      = 0       ! array
        SPCTOTPR      = 0       ! array
        SPCIDX        = 0       ! array

C.........  Set the length of the emission type joiner
        LT = LEN_TRIM( ETJOIN )
        LS = LEN_TRIM( SPJOIN )

C.........  Populate the emission type list and the pollutant list from the
C           reporting bins module
        DO V = 1, NIPPA

            EBUF = EANAM( V )

            K = INDEX( EBUF, ETJOIN )   ! Look for emission type joiner

            IF( K .GT. 0 ) THEN         ! Store emission type and pol from it
        	L2 = LEN_TRIM( EBUF )
                ETPNAM( V ) = EBUF
                DATNAM( V ) = EBUF( K+LT:L2 )

            ELSE                        ! Store pollutant only
                DATNAM( V ) = EBUF
                
            END IF

C.............  Store index from "all" list to inventory file variables
            INVIDX( V ) = V

        END DO

        NDATALL = NIPPA

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

                IF ( .NOT. DESC3( SLNAME ) ) THEN

        	    MESG = 'Could not get description of file "' //
     &                     SLNAME( 1:LEN_TRIM( SLNAME ) ) // '"'
        	    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

                ENDIF

C.................  Store units
                CALL STORE_VUNITS( 1, 1, NVARS3D, SLUNIT )

            END IF

C.............  Get header of mass speciation matrix
            IF( SSFLAG ) THEN

                IF ( .NOT. DESC3( SSNAME ) ) THEN

        	    MESG = 'Could not get description of file "' //
     &                     SSNAME( 1:LEN_TRIM( SSNAME ) ) // '"'
        	    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

                ENDIF

C.................  Store units
                CALL STORE_VUNITS( 1, 1, NVARS3D, SSUNIT )

            END IF

C.............  Populate the speciation variable lists
C.............  Populate the speciation units lists
C.............  Count the number of unique species
            NMSPC = 0
            DO V = 1, NSVARS

        	VBUF = VDESC3D( V )

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
                I3 = INDEX1( EBUF, NTPDAT, TPNAME )   ! Look in hourly file

C.................  For species based on inventory pollutant, set index
                IF( I1 .GT. 0 ) THEN

                    SPCTOINV( V ) = MAX( I1,0 )

C.................  For species based on activity, set index and write note
                ELSE IF( I2 .GT. 0 .AND. TPACTIDX( I3 ) .GT. 0 ) THEN

                    SPCTOTPR( V ) = I3 + NIPPA

                    I2 = INDEX1( SBUF, V, SPCNAM )
                    IF( I2 .LE. 0 ) THEN
                        L2 = LEN_TRIM( SBUF )
                        MESG = 'NOTE: Species "' // SBUF(1:L2) // 
     &                         '" created based on activity data.'
                        CALL M3MESG( MESG )
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
                    ETPSPCNAM( V )= VBUF
                    PRCSPCNAM( V )= 'S-'// VBUF( 1:J-1 )// ETJOIN// SBUF
                    SUMETPNAM( V )= 'S-'// EBUF
                    SUMPOLNAM( V )= 'S-'// VBUF( J+LT:K-1 )

C.................  No emission type
                ELSE

                    SPCNAM   ( V ) = SBUF
                    SUMPOLNAM( V ) = 'S-' // EBUF

                END IF

                SRTIDX( V ) = V

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
        	IF( K .GT. 0 ) SPCNAM( V ) = ' '
            END DO  ! End creating list of unique species

        END IF      ! End if speciation

C.........  Determine if any reports did not have a DATA instruction
        I = MINVAL( ALLRPT%NUMDATA )

C.........  If a report did not have a DATA instruction, then maximum output
C           variables depends on whether we have speciation or not
        MXOUTDAT = MXINDAT
        IF( I .LT. 0 ) THEN
            MXOUTDAT = MAX( MXOUTDAT, NDATALL )
            IF( SFLAG ) MXOUTDAT = MAX( MXOUTDAT, NDATALL + NMSPC )
        END IF

C.........  Allocate memory of reading and output data names
        ALLOCATE( OUTDNAM( MXOUTDAT, NREPORT ), STAT=IOS )
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
                    OUTDNAM( 1:NTPDAT,N ) = TPNAME( 1:NTPDAT )
                ELSE
                    J = NIPPA
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
                	    OUTDNAM( J,N ) = 'S-' // SBUF
                        ELSE
                	    OUTDNAM( J,N ) = SBUF
                        END IF
                    END DO

                END IF

                ALLRPT( N )%NUMDATA = J
                NDATA = J

C.............  Otherwise, set output data values as input data values
            ELSE

                OUTDNAM( 1:NDATA,N ) = INDNAM( 1:NDATA,N )

            END IF

C.............  Loop through requested data for this report
            DO I = 1, NDATA

                ANYOUT = .FALSE.

C.................  Loop through names of emission types
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
                    END IF

C.....................  To pollutant/activity column
                    IF( OUTDNAM( I,N ) .EQ. DATNAM( E ) ) THEN
                        TODOUT( E,N )%DAT = I
                        TODOUT( E,N )%AGG = 1
                        ANYOUT = .TRUE.
                    END IF

                END DO

C.................  If speciation for current report
                IF( RPT_%USESLMAT .OR.
     &              RPT_%USESSMAT      ) THEN

C.....................  Loop through names of speciation variables
                    DO V = 1, NSVARS

C.........................  Get species to inventory data index
                        E = SPCTOINV( V )

C..........................  Set inventory variable as having speciation
                        IF( E .GT. 0 )
     &                      TODOUT( E,N )%SPCYN= ( RPT_%USESLMAT .OR.
     &                                             RPT_%USESSMAT    )

C.........................  Get species to hourly data index
                        E = SPCTOTPR( V )

C..........................  Set hourly variable as having speciation
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
                	END IF

C.........................  To post-speciation summed emission type column
                	IF( OUTDNAM( I,N ) .EQ. SUMETPNAM( V ) ) THEN 
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
 
C.............  This subprogram stores I/O API NetCDF variable units into
C               a local array based on indices in subprogram call.
            SUBROUTINE STORE_VUNITS( ISTART, INCRMT, NUNIT, UNITS )

C.............  Subprogram arguments
            INTEGER      ISTART        ! starting position in UNITS3D of names
            INTEGER      INCRMT        ! increment of UNITS3D for names
            INTEGER      NUNIT         ! number of units
            CHARACTER(*) UNITS( NUNIT )! stored variable units

C.............  Local variables
            INTEGER  I, J, L

C----------------------------------------------------------------------

            UNITS = ' '

            J = ISTART
            DO I = 1, NUNIT

                L = LEN_TRIM( UNITS3D( J ) )
                UNITS( I ) = UNITS3D( J )( 1:L )
                J = J + INCRMT

            END DO
 
            END SUBROUTINE STORE_VUNITS

        END SUBROUTINE BLDREPIDX

