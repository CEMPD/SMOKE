
        SUBROUTINE ASGNBINS( RCNT )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C       The ASGNBINS routine is responsible for assigning a bin number to each
C       output record.  A bin is a group of records for which the emissions,
C       activity, or emission-type data will be summed to provide a single
C       record in a report.
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
C COPYRIGHT (C) 2001, MCNC--North Carolina Supercomputing Center
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
C***********************************************************************

C...........   MODULES for public variables
C...........   This module is the inventory arrays
        USE MODSOURC

C.........  This module contains the lists of unique source characteristics
        USE MODLISTS

C.........  This module contains Smkreport-specific settings
        USE MODREPRT

C.........  This module contains report arrays for each output bin
        USE MODREPBN

C.........  This module contains arrays for plume-in-grid and major sources
        USE MODELEV

C.........  This module contains the arrays for state and county summaries
        USE MODSTCY

C.........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........  EXTERNAL FUNCTIONS and their descriptions:
        INTEGER    INDEX1
        INTEGER    FIND1
        INTEGER    FINDC

        EXTERNAL   INDEX1, FIND1, FINDC

C...........   SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN) :: RCNT    ! current report number

C...........   Local parameters
        INTEGER, PARAMETER :: BUFLEN = 101 + SCCLEN3 + SPNLEN3

C...........   Sorting arrays
        INTEGER              , ALLOCATABLE :: SORTIDX( : )

        CHARACTER(LEN=BUFLEN), ALLOCATABLE :: SORTBUF( : )

C...........   Local variables
        INTEGER         B, C, F, I, J, K, LB, S

        INTEGER         COL               ! tmp column number
        INTEGER         DIUID             ! tmp diurnal profile number
        INTEGER         FIP               ! tmp country/state/county
        INTEGER         IOS               ! i/o status
        INTEGER         MONID             ! tmp monthly profile number
        INTEGER         NDATA             ! no. output data columns for current
        INTEGER         RCL               ! tmp road class code
        INTEGER         ROW               ! tmp row number
        INTEGER         SRCID             ! tmp source ID
        INTEGER         SRGID1            ! tmp primary surrogate ID
        INTEGER         SRGID2            ! tmp fallback surrogate ID
        INTEGER         WEKID             ! tmp weekly profile number

        CHARACTER*1            ESTAT      ! tmp elevated status
        CHARACTER*60           FMTBUF     ! format buffer
        CHARACTER*300          MESG       ! message buffer

        CHARACTER(LEN=BUFLEN)  BUFFER     ! sorting info buffer
        CHARACTER(LEN=BUFLEN)  LBUF       ! previous sorting info buffer
        CHARACTER(LEN=SCCLEN3) SCC        ! tmp SCC
        CHARACTER(LEN=SPNLEN3) SPCID      ! tmp speciation profile

        CHARACTER*16 :: PROGNAME = 'ASGNBINS' ! program name

C***********************************************************************
C   begin body of subroutine ASGNBINS

C.........  Set report-specific local settings
        NDATA   = ALLRPT( RCNT )%NUMDATA
        RPT_    = ALLRPT( RCNT )
        LREGION = ( RPT_%BYCNTY .OR. RPT_%BYSTAT .OR. RPT_%BYCNRY )

C.........  Allocate (and deallocate) memory for sorting arrays
        IF( ALLOCATED( SORTIDX ) ) DEALLOCATE( SORTIDX, SORTBUF )

        ALLOCATE( SORTIDX( NOUTREC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SORTIDX', PROGNAME )
        ALLOCATE( SORTBUF( NOUTREC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SORTBUF', PROGNAME )

C.........  Build format statement for writing the sorting buffer
C           (building it in case SCC width changes in the future)
        WRITE( FMTBUF, '(A,I2.2,A,I2.2,A)' ) 
     &         '(4I8, A', SCCLEN3, ',5I8,A', SPNLEN3, ',I8,A)'

C.........  Initialize local variables for building sorting array for this 
C           report
        COL    = 0
        ROW    = 0
        SRCID  = 0
        FIP    = 0
        RCL    = 0
        SRGID1 = 0
        SRGID2 = 0
        MONID  = 0
        WEKID  = 0
        DIUID  = 0
        SPCID  = ' '
        SCC    = ' '
        ESTAT  = ' '

C.........  Create a sorting array for all output records
        DO I = 1, NOUTREC

C.............  If BY CELL, insert X-cell and then Y-cell 
            IF( RPT_%BYCELL ) THEN
                C = OUTCELL( I )
                ROW = C / NCOLS          ! note: integer math
                IF( MOD( C, NCOLS ) .GT. 0. ) ROW = ROW + 1
                COL = C - ( ROW-1 ) * NCOLS
            END IF

C.............  If BY SOURCE, then nothing else needed
            IF( RPT_%BYSRC ) THEN
                SRCID = OUTSRC( I )
                FIP   = IFIP( SRCID )
                SCC   = CSCC( SRCID )

            ELSE

C.................  If BY COUNTY insert region code (state and country not 
C                   needed)
                IF( RPT_%BYCNTY ) THEN

                    FIP = IFIP( OUTSRC( I ) )

C.................  If BY STATE, insert region code with trailing zeros 
C                   (country not needed)
                ELSE IF( RPT_%BYSTAT ) THEN

                    FIP = IFIP( OUTSRC( I ) )
                    FIP = ( FIP / 1000 ) * 1000    ! integer math

C.................  If BY COUNTRY, insert region code with trailing zeros
                ELSE IF( RPT_%BYCNRY ) THEN

                    FIP = IFIP( OUTSRC( I ) )
                    FIP = ( FIP / 100000 ) * 100000    ! integer math

                END IF  ! End by county, state, or country

C.................  If BY SCC10, insert full SCC, and nothing else needed
                IF( RPT_%BYSCC ) THEN

                    SCC = CSCC( OUTSRC( I ) )

C.................  If BY ROADCLASS, insert roadclass code        
                ELSE IF( RPT_%BYRCL ) THEN

                    RCL = IRCLAS( OUTSRC( I ) )

                END IF  ! End by source or by roadclass


            END IF      ! End by source or not

C.................  If by surrogates IDs, insert them depending on resolution
            IF( RPT_%BYSRG ) THEN

                SRGID2 = SRGID( OUTSRC( I ), 2 )
                IF( RPT_%SRGRES .EQ. 1 ) THEN
                    SRGID1 = SRGID( OUTSRC( I ), 1 )
                END IF

            END IF  ! End by surrogate IDs

C.................  If by monthly temporal profile, insert them
            IF( RPT_%BYMON ) MONID = IMON( OUTSRC( I ) )

C.................  If by weekly temporal profile, insert them
            IF( RPT_%BYWEK ) WEKID = IWEK( OUTSRC( I ) )

C.................  If by diurnal temporal profile, insert them
            IF( RPT_%BYDIU ) DIUID = IDIU( OUTSRC( I ) )

C.................  If by speciation profile, insert them based on the pollutant
C                   specified
            IF( RPT_%BYSPC ) THEN
                J = INDEX1( RPT_%SPCPOL, NSPCPOL, SPCPOL )

C.................  Unless coding error, RPT_%SPCPOL should be found
                IF ( J .LE. 0 ) THEN

                    MESG = 'INTERNAL ERROR: Pollutant "'// RPT_%SPCPOL//
     &                     'not found in list created from REPCONFIG.'
                    CALL M3MSG2( MESG )
                    CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

                ELSE
                    SPCID = SPPROF( OUTSRC(I),J )                    
                END IF
            END IF

C.................  If BY ELEVSTAT, insert elevated status code
            IF( RPT_%BYELEV ) THEN

                IF( LPING( OUTSRC( I ) ) ) THEN       ! PinG
                    ESTAT = 'P'
                ELSE IF( LMAJOR( OUTSRC( I ) ) ) THEN ! Elevated
                    ESTAT = 'E'
                ELSE                                      ! Low-level
                    ESTAT = 'L'
                END IF

            END IF  ! End by elevated status

C.............  Store sorting information for current record
            WRITE( BUFFER,FMTBUF ) COL, ROW, SRCID, FIP, SCC, 
     &                             SRGID1, SRGID2, MONID, WEKID, DIUID,
     &                             SPCID, RCL, ESTAT

            SORTIDX( I ) = I
            SORTBUF( I ) = BUFFER

        END DO   ! End loop I over output records

C.........  Sort sorting array
        CALL SORTIC( NOUTREC, SORTIDX, SORTBUF )

C.........  Assign bins to output records based on sorting array
        LBUF = ' '
        B = 0
        DO I = 1, NOUTREC

            J = SORTIDX( I )
            IF( SORTBUF( J ) .NE. LBUF ) THEN
                B = B + 1
                LBUF = SORTBUF( J )
            END IF

            OUTBIN( J ) = B

        END DO

        NOUTBINS = B

C.........  If memory is allocated for bin arrays, then deallocate
        IF( ALLOCATED( BINCOIDX  ) ) DEALLOCATE( BINCOIDX )
        IF( ALLOCATED( BINSTIDX  ) ) DEALLOCATE( BINSTIDX )
        IF( ALLOCATED( BINCYIDX  ) ) DEALLOCATE( BINCYIDX )
        IF( ALLOCATED( BINREGN   ) ) DEALLOCATE( BINREGN )
        IF( ALLOCATED( BINSMKID  ) ) DEALLOCATE( BINSMKID )
        IF( ALLOCATED( BINSCC    ) ) DEALLOCATE( BINSCC )
        IF( ALLOCATED( BINSRGID1 ) ) DEALLOCATE( BINSRGID1 )
        IF( ALLOCATED( BINSRGID2 ) ) DEALLOCATE( BINSRGID2 )
        IF( ALLOCATED( BINSNMIDX ) ) DEALLOCATE( BINSNMIDX )
        IF( ALLOCATED( BINRCL    ) ) DEALLOCATE( BINRCL )
        IF( ALLOCATED( BINX      ) ) DEALLOCATE( BINX )
        IF( ALLOCATED( BINY      ) ) DEALLOCATE( BINY )
        IF( ALLOCATED( BINELEV   ) ) DEALLOCATE( BINELEV )
        IF( ALLOCATED( BINDATA   ) ) DEALLOCATE( BINDATA )

C.........  Allocate memory for bins
        IF( RPT_%BYCONAM ) ALLOCATE( BINCOIDX ( NOUTBINS ), STAT=IOS )
        IF( RPT_%BYSTNAM ) ALLOCATE( BINSTIDX ( NOUTBINS ), STAT=IOS )
        IF( RPT_%BYCYNAM ) ALLOCATE( BINCYIDX ( NOUTBINS ), STAT=IOS )
        IF( LREGION      ) ALLOCATE( BINREGN  ( NOUTBINS ), STAT=IOS )
        IF( RPT_%BYSRC   ) ALLOCATE( BINSMKID ( NOUTBINS ), STAT=IOS )
        IF( RPT_%BYSCC   ) ALLOCATE( BINSCC   ( NOUTBINS ), STAT=IOS )
        IF( RPT_%SCCNAM  ) ALLOCATE( BINSNMIDX( NOUTBINS ), STAT=IOS )
        IF( RPT_%SRGRES .EQ. 1 ) 
     &                     ALLOCATE( BINSRGID1( NOUTBINS ), STAT=IOS )
        IF( RPT_%SRGRES .GE. 1 ) 
     &                     ALLOCATE( BINSRGID2( NOUTBINS ), STAT=IOS )
        IF( RPT_%BYMON   ) ALLOCATE( BINMONID ( NOUTBINS ), STAT=IOS )
        IF( RPT_%BYWEK   ) ALLOCATE( BINWEKID ( NOUTBINS ), STAT=IOS )
        IF( RPT_%BYDIU   ) ALLOCATE( BINDIUID ( NOUTBINS ), STAT=IOS )
        IF( RPT_%BYSPC   ) ALLOCATE( BINSPCID ( NOUTBINS ), STAT=IOS )
        IF( RPT_%BYRCL   ) ALLOCATE( BINRCL   ( NOUTBINS ), STAT=IOS )
        IF( RPT_%BYCELL  ) ALLOCATE( BINX     ( NOUTBINS ), STAT=IOS )
        IF( RPT_%BYCELL  ) ALLOCATE( BINY     ( NOUTBINS ), STAT=IOS )
        IF( RPT_%BYELEV  ) ALLOCATE( BINELEV  ( NOUTBINS ), STAT=IOS )

        ALLOCATE( BINDATA( NOUTBINS, NDATA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'BINDATA', PROGNAME )
        IF( NDATA .GT. 0 ) BINDATA = 0  ! array

C.........  Populate the bin characteristic arrays (not the data array)

        LB = 0
        DO I = 1, NOUTREC

            J = SORTIDX( I )
            B = OUTBIN( J )
            BUFFER = SORTBUF( J )

            IF( B .NE. LB ) THEN

                READ( BUFFER,FMTBUF ) 
     &                COL, ROW, SRCID, FIP, SCC, SRGID1, SRGID2, 
     &                MONID, WEKID, DIUID, SPCID, RCL, ESTAT

C.................  Store region code
                IF( LREGION ) BINREGN( B ) = FIP

C.................  Store country name index
                IF( RPT_%BYCONAM ) THEN
                    F = FIP / 100000
                    K = FIND1( F, NCOUNTRY, CTRYCOD )
                    BINCOIDX( B ) = K
                END IF

C.................  Store state name index
                IF( RPT_%BYSTNAM ) THEN
                    F = ( FIP / 1000 ) * 1000        ! In case by-county also
                    K = FIND1( F, NSTATE, STATCOD )
                    BINSTIDX( B ) = K
                END IF

C.................  Store county name index
                IF( RPT_%BYCYNAM ) THEN
                    K = FIND1( FIP, NCOUNTY, CNTYCOD )
                    BINCYIDX( B ) = K
                END IF

C.................  Store road class
                IF( RPT_%BYRCL ) BINRCL( B ) = RCL

C.................  Store SMOKE ID
                IF( RPT_%BYSRC ) BINSMKID( B ) = SRCID

C.................  Store SCC
                IF( RPT_%BYSCC ) BINSCC( B ) = SCC

C.................  Store SCC name index
                IF( RPT_%SCCNAM ) THEN
                    K = FINDC( SCC, NINVSCC, INVSCC )
                    BINSNMIDX( B ) = K
                END IF

C.................  Store surrogate codes
                IF( RPT_%SRGRES .EQ. 1 ) BINSRGID1( B ) = SRGID1
                IF( RPT_%SRGRES .GE. 1 ) BINSRGID2( B ) = SRGID2

C.................  Store temporal profiles
                IF( RPT_%BYMON ) BINMONID( B ) = MONID
                IF( RPT_%BYWEK ) BINWEKID( B ) = WEKID
                IF( RPT_%BYDIU ) BINDIUID( B ) = DIUID

C.................  Store speciation profiles
                IF( RPT_%BYSPC ) BINSPCID( B ) = SPCID

C.................  Store x-cell and y-cell
                IF( RPT_%BYCELL ) BINX( B ) = COL
                IF( RPT_%BYCELL ) BINY( B ) = ROW

C.................  Store Elevated status
                IF( RPT_%BYELEV ) BINELEV( B ) = ESTAT

            END IF

        END DO   ! End loop I over output records

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I10, :, 1X ) )

        END SUBROUTINE ASGNBINS
 
