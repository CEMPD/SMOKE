
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
C...........   This module is the inventory arrays
        USE MODSOURC, ONLY: CSOURC, IFIP, CSCC, IRCLAS, SRGID, IMON,
     &                      IWEK, IDIU, SPPROF, ISIC, CMACT, CNAICS,
     &                      CSRCTYP, CORIS

C.........  This module contains the lists of unique source characteristics
        USE MODLISTS, ONLY: NINVSCC, INVSCC, NINVSIC, INVSIC, NINVMACT,
     &                      INVMACT, NINVNAICS, INVNAICS

C.........  This module contains Smkreport-specific settings
        USE MODREPRT, ONLY: RPT_, LREGION, AFLAG, ALLRPT, NSPCPOL,
     &                      SPCPOL, STKX, STKY, LOC_BEGP, LOC_ENDP

C.........  This module contains report arrays for each output bin
        USE MODREPBN, ONLY: NOUTREC, NOUTBINS, BINBAD, BINCOIDX,
     &                      BINSTIDX, BINCYIDX, BINREGN, BINSMKID,
     &                      BINSCC, BINSRGID1, BINSRGID2, BINSNMIDX,
     &                      BINRCL, BINMONID, BINWEKID, BINDIUID,
     &                      BINSPCID, BINPLANT, BINX, BINY, BINELEV,
     &                      BINPOPDIV, BINDATA, OUTBIN, OUTCELL,OUTSRC,
     &                      BINSIC, BINSICIDX, BINMACT, BINMACIDX,
     &                      BINNAICS, BINNAIIDX, BINSRCTYP, BINORIS,
     &                      BINORSIDX

C.........  This module contains the global variables for the 3-d grid
        USE MODGRID, ONLY: NCOLS

C.........  This module contains arrays for plume-in-grid and major sources
        USE MODELEV, ONLY: LPING, LMAJOR

C.........  This module contains the arrays for state and county summaries
        USE MODSTCY, ONLY: NCOUNTRY, CTRYCOD, NSTATE, STATCOD, NCOUNTY,
     &                     CNTYCOD, CTRYPOPL, STATPOPL, CNTYPOPL,
     &                     NORIS, ORISLST

C.........  This module contains the information about the source category
        USE MODINFO, ONLY: CATEGORY
        
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
        INTEGER, PARAMETER :: BUFLEN = 101 + SCCLEN3 + SICLEN3 + SPNLEN3
     &                                     + MACLEN3 + NAILEN3 + STPLEN3
     &                                     + ORSLEN3
        INTEGER, PARAMETER :: PTSCCLEV( NSCCLV3 ) =
     &                        ( / 1, 3, 6, 8 / )
        INTEGER, PARAMETER :: ARSCCLEV( NSCCLV3 ) =
     &                        ( / 2, 4, 7, 10 / )
     
C...........   Sorting arrays
        INTEGER          , ALLOCATABLE :: SORTIDX( : )

        CHARACTER(BUFLEN), ALLOCATABLE :: SORTBUF( : )

C...........   Local variables
        INTEGER         B, C, F, I, J, K, L, LB, S

        INTEGER         COL               ! tmp column number
        INTEGER         DIUID             ! tmp diurnal profile number
        INTEGER         FIP               ! tmp country/state/county
        INTEGER         IOS               ! i/o status
        INTEGER         MONID             ! tmp monthly profile number
        INTEGER         NDATA             ! no. output data columns for current
        INTEGER         PREVFIP           ! previous FIPs code
        INTEGER         PREVSRCID         ! previous source ID
        INTEGER         RCL               ! tmp road class code
        INTEGER         ROW               ! tmp row number
        INTEGER         SIC               ! tmp SIC
        INTEGER         SRCID             ! tmp source ID
        INTEGER         SRGID1            ! tmp primary surrogate ID
        INTEGER         SRGID2            ! tmp fallback surrogate ID
        INTEGER         WEKID             ! tmp weekly profile number

        CHARACTER              ESTAT      ! tmp elevated status
        CHARACTER(60)          FMTBUF     ! format buffer
        CHARACTER(300)         MESG       ! message buffer

        CHARACTER(5)       SCCTYPE    ! tmp determination of SCC type
        CHARACTER(BUFLEN)  BUFFER     ! sorting info buffer
        CHARACTER(BUFLEN)  LBUF       ! previous sorting info buffer
        CHARACTER(SCCLEN3) SCC        ! tmp SCC
        CHARACTER(MACLEN3) MACT       ! tmp MACT
        CHARACTER(NAILEN3) NAICS      ! tmp NAICS
        CHARACTER(ORSLEN3) ORIS       ! tmp ORIS
        CHARACTER(STPLEN3) SRCTYP     ! tmp SRCTYP
        CHARACTER(SPNLEN3) SPCID      ! tmp speciation profile
        CHARACTER(PLTLEN3) PLANT      ! tmp plant ID
        CHARACTER(PLTLEN3) PREVPLT    ! previous plant ID

        CHARACTER(16) :: PROGNAME = 'ASGNBINS' ! program name

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
      WRITE( FMTBUF,'(A,I2.2,A,I1,A,I2.2,A,I2.2,A,I1,A,I1,A,I1,
     &                A,I1,A)') 
     &    '(4I8,A',SCCLEN3,',I',SICLEN3,',5I8,A', SPNLEN3,',A',
     &    PLTLEN3,',A',ORSLEN3,',I8,A,A', MACLEN3,',A', NAILEN3,',A', 
     &    STPLEN3, ')'

C.........  Initialize local variables for building sorting array for this 
C           report
        COL    = 0
        ROW    = 0
        SRCID  = 0
        FIP    = 0
        RCL    = 0
        SIC    = 0
        SRGID1 = 0
        SRGID2 = 0
        MONID  = 0
        WEKID  = 0
        DIUID  = 0
        SPCID  = ' '
        PLANT  = ' '
        SCC    = ' '
        MACT   = ' '
        NAICS  = ' '
        ORIS   = ' '
        SRCTYP = ' '
        ESTAT  = ' '

C.........  Create a sorting array for all output records
        DO I = 1, NOUTREC

C.............  If BY CELL, insert X-cell and then Y-cell 
            IF( RPT_%BYCELL ) THEN
                IF( .NOT. AFLAG ) THEN
                    C = OUTCELL( I )
                    ROW = 1 + INT( (C-1) / NCOLS )
                    COL = C - NCOLS * ( ROW-1 )
                ELSE
                    ROW = STKY( I )
                    COL = STKX( I )
                END IF
            END IF

C.............  If BY SOURCE, then nothing else needed
            IF( RPT_%BYSRC ) THEN
                SRCID = OUTSRC( I )
                FIP   = IFIP( SRCID )
                IF( .NOT. AFLAG ) THEN
                    SCC   = CSCC( SRCID )
                    IF( RPT_%BYSIC ) SIC = ISIC( SRCID )
                END IF

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

C.................  If BY SCC, insert SCC based on SCCRES (resolution) set in input file
                IF( RPT_%BYSCC ) THEN

                    SCC = CSCC( OUTSRC( I ) )

C.....................  Get rid of leading zeros
                    IF( SCC( 1:2 ) == '00' ) SCC = SCC(3:SCCLEN3) // '  '

                    IF( RPT_%SCCRES < 4 ) THEN

                        L = LEN( TRIM( SCC ) )
                        SELECT CASE( L )

C.........................  If 10 characters long, then nonpoint SCC
                        CASE( 10 )
                            SCCTYPE = 'AREA'

C.........................  If 8 characters long, then point SCC
                        CASE( 8 )
                            SCCTYPE = 'POINT'

C.........................  If SCC length less than full SCC, then assume based on category
                        CASE DEFAULT
                            SCCTYPE = CATEGORY

                        END SELECT

C.........................  Set SCC to use for bins based on level from REPCONFIG
                        IF( SCCTYPE == 'POINT' ) THEN
                            SCC = SCC( 1:PTSCCLEV( RPT_%SCCRES ) )

                        ELSE IF( SCCTYPE == 'AREA' .OR.
     &                           SCCTYPE == 'MOBILE'    ) THEN
                            SCC = SCC( 1:ARSCCLEV( RPT_%SCCRES ) )

                        END IF
                    END IF
                END IF  ! end if by SCC
            END IF  ! end if by source

C.............  If BY ROADCLASS, insert roadclass code        

            IF( RPT_%BYRCL ) THEN

                RCL = IRCLAS( OUTSRC( I ) )

            END IF  ! End by source or by roadclass

C.................  If BY SIC, insert full SIC
            IF( RPT_%BYSIC ) SIC = ISIC( OUTSRC( I ) )
 
C.................  If BY MACT, insert full MACT
            IF( RPT_%BYMACT ) THEN
                IF( .NOT. ASSOCIATED( CMACT ) ) THEN
                    MESG = 'ERROR: BY MACT is requested, but ' //
     &                    'MACT is not present in ASCII inventory file'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                ELSE
                    MACT = CMACT( OUTSRC( I ) )
                END IF
            END IF

C.................  If BY NAICS, insert full NAICS
            IF( RPT_%BYNAICS ) THEN
                IF( .NOT. ASSOCIATED( CNAICS ) ) THEN
                    MESG = 'ERROR: BY NAICS is requested, but ' //
     &                    'NAICS is not present in ASCII inventory file'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                ELSE
                   NAICS = CNAICS( OUTSRC( I ) )
                END IF
            END IF

C.................  If BY ORIS, insert full ORIS
            IF( RPT_%BYORIS ) THEN
                IF( .NOT. ALLOCATED( CORIS ) ) THEN
                    MESG = 'ERROR: BY ORIS is requested, but ' //
     &                    'ORIS is not present in ASCII inventory file'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                ELSE
                   ORIS = CORIS( OUTSRC( I ) )
                END IF
            END IF

C.................  If BY SRCTYP, insert full Source type
            IF( RPT_%BYSRCTYP ) THEN
                IF( .NOT. ASSOCIATED( CSRCTYP ) ) THEN
                    MESG = 'ERROR: BY SRCTYP is requested, but ' //
     &                    'source type code is not present in ASCII ' //
     &                    'inventory file'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                ELSE
                    SRCTYP = CSRCTYP( OUTSRC( I ) )
                END IF
            END IF

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

C.................  If BY PLANT, get plant ID and set same source
C                   number until plant changes
            IF( RPT_%BYPLANT ) THEN
                S = OUTSRC( I )
                PLANT = CSOURC( S ) (LOC_BEGP(2):LOC_ENDP(2))

C...............  If this is the same plant, then set the old source
C                 ID so that the bins will still be "by plant"
                IF ( IFIP( S ) .EQ. PREVFIP .AND. 
     &               PLANT     .EQ. PREVPLT       ) THEN
                    SRCID = PREVSRCID

C...............  If this is a different plant, the reset SRCID to
C                 be the first source for the current plant
                ELSE
                    SRCID = S
                    PREVFIP   = IFIP( S )
                    PREVPLT   = PLANT
                    PREVSRCID = S
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
            WRITE( BUFFER,FMTBUF ) COL, ROW, SRCID, FIP, SCC, SIC,
     &                             SRGID1, SRGID2, MONID, WEKID, DIUID,
     &                             SPCID, PLANT, ORIS, RCL, ESTAT, MACT,
     &                             NAICS, SRCTYP

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
        IF( ALLOCATED( BINBAD    ) ) DEALLOCATE( BINBAD )
        IF( ALLOCATED( BINCOIDX  ) ) DEALLOCATE( BINCOIDX )
        IF( ALLOCATED( BINSTIDX  ) ) DEALLOCATE( BINSTIDX )
        IF( ALLOCATED( BINCYIDX  ) ) DEALLOCATE( BINCYIDX )
        IF( ALLOCATED( BINREGN   ) ) DEALLOCATE( BINREGN )
        IF( ALLOCATED( BINSMKID  ) ) DEALLOCATE( BINSMKID )
        IF( ALLOCATED( BINSCC    ) ) DEALLOCATE( BINSCC )
        IF( ALLOCATED( BINSIC    ) ) DEALLOCATE( BINSIC )
        IF( ALLOCATED( BINMACT   ) ) DEALLOCATE( BINMACT )
        IF( ALLOCATED( BINMACIDX ) ) DEALLOCATE( BINMACIDX )
        IF( ALLOCATED( BINNAICS  ) ) DEALLOCATE( BINNAICS )
        IF( ALLOCATED( BINNAIIDX ) ) DEALLOCATE( BINNAIIDX )
        IF( ALLOCATED( BINORIS   ) ) DEALLOCATE( BINORIS )
        IF( ALLOCATED( BINORSIDX ) ) DEALLOCATE( BINORSIDX )
        IF( ALLOCATED( BINSRCTYP ) ) DEALLOCATE( BINSRCTYP)
        IF( ALLOCATED( BINSRGID1 ) ) DEALLOCATE( BINSRGID1 )
        IF( ALLOCATED( BINSRGID2 ) ) DEALLOCATE( BINSRGID2 )
        IF( ALLOCATED( BINSNMIDX ) ) DEALLOCATE( BINSNMIDX )
        IF( ALLOCATED( BINRCL    ) ) DEALLOCATE( BINRCL )
        IF( ALLOCATED( BINMONID  ) ) DEALLOCATE( BINMONID )
        IF( ALLOCATED( BINWEKID  ) ) DEALLOCATE( BINWEKID )
        IF( ALLOCATED( BINDIUID  ) ) DEALLOCATE( BINDIUID )
        IF( ALLOCATED( BINSPCID  ) ) DEALLOCATE( BINSPCID )
        IF( ALLOCATED( BINPLANT  ) ) DEALLOCATE( BINPLANT )
        IF( ALLOCATED( BINX      ) ) DEALLOCATE( BINX )
        IF( ALLOCATED( BINY      ) ) DEALLOCATE( BINY )
        IF( ALLOCATED( BINELEV   ) ) DEALLOCATE( BINELEV )
        IF( ALLOCATED( BINPOPDIV ) ) DEALLOCATE( BINPOPDIV )
        IF( ALLOCATED( BINDATA   ) ) DEALLOCATE( BINDATA )

C.........  Allocate memory for bins
        ALLOCATE( BINBAD( NOUTBINS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'BINBAD', PROGNAME )
        BINBAD = 0    ! array

        IF( RPT_%BYCONAM ) THEN
            ALLOCATE( BINCOIDX ( NOUTBINS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BINCOIDX', PROGNAME )
        ENDIF
        IF( RPT_%BYSTNAM ) THEN
            ALLOCATE( BINSTIDX ( NOUTBINS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BINSTIDX', PROGNAME )
        ENDIF
        IF( RPT_%BYCYNAM ) THEN
            ALLOCATE( BINCYIDX ( NOUTBINS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BINCYIDX', PROGNAME )
        ENDIF
        IF( LREGION      ) THEN
            ALLOCATE( BINREGN  ( NOUTBINS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BINREGN', PROGNAME )
        ENDIF
        IF( RPT_%BYSRC .OR. RPT_%BYPLANT ) THEN
            ALLOCATE( BINSMKID ( NOUTBINS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BINSMKID', PROGNAME )
        ENDIF
        IF( RPT_%BYSCC   ) THEN
            ALLOCATE( BINSCC   ( NOUTBINS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BINSCC', PROGNAME )
        ENDIF
        IF( RPT_%SCCNAM  ) THEN
            ALLOCATE( BINSNMIDX( NOUTBINS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BINSNMIDX', PROGNAME )
        ENDIF
        IF( RPT_%BYSIC   ) THEN
            ALLOCATE( BINSIC   ( NOUTBINS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BINSIC', PROGNAME )
        ENDIF
        IF( RPT_%SICNAM   ) THEN
            ALLOCATE( BINSICIDX( NOUTBINS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BINSICIDX', PROGNAME )
        ENDIF
                
        IF( RPT_%BYMACT   ) THEN
            ALLOCATE( BINMACT   ( NOUTBINS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BINMACT', PROGNAME )
        ENDIF
        
        IF( RPT_%MACTNAM   ) THEN
            ALLOCATE( BINMACIDX( NOUTBINS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BINMACIDX', PROGNAME )
        ENDIF
        
        IF( RPT_%BYNAICS   ) THEN
            ALLOCATE( BINNAICS   ( NOUTBINS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BINNAICS', PROGNAME )
        ENDIF

        IF( RPT_%NAICSNAM   ) THEN
            ALLOCATE( BINNAIIDX( NOUTBINS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BINNAIIDX', PROGNAME )
        ENDIF

        IF( RPT_%BYORIS   ) THEN
            ALLOCATE( BINORIS   ( NOUTBINS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BINORIS', PROGNAME )
        ENDIF

        IF( RPT_%ORISNAM   ) THEN
            ALLOCATE( BINORSIDX( NOUTBINS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BINORSIDX', PROGNAME )
        ENDIF

        IF( RPT_%BYSRCTYP   ) THEN
            ALLOCATE( BINSRCTYP   ( NOUTBINS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BINSRCTYP', PROGNAME )
        ENDIF
        
        IF( RPT_%SRGRES .EQ. 1 ) THEN
            ALLOCATE( BINSRGID1( NOUTBINS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BINSRGID1', PROGNAME )
        ENDIF
        IF( RPT_%SRGRES .GE. 1 ) THEN
            ALLOCATE( BINSRGID2( NOUTBINS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BINSRGID2', PROGNAME )
        ENDIF
        IF( RPT_%BYMON   ) THEN
            ALLOCATE( BINMONID ( NOUTBINS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BINMONID', PROGNAME )
        ENDIF
        IF( RPT_%BYWEK   ) THEN
            ALLOCATE( BINWEKID ( NOUTBINS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BINWEKID', PROGNAME )
        ENDIF
        IF( RPT_%BYDIU   ) THEN
            ALLOCATE( BINDIUID ( NOUTBINS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BINDIUID', PROGNAME )
        ENDIF
        IF( RPT_%BYSPC   ) THEN
            ALLOCATE( BINSPCID ( NOUTBINS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BINSPCID', PROGNAME )
        ENDIF
        IF( RPT_%BYPLANT ) THEN
            ALLOCATE( BINPLANT ( NOUTBINS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BINPLANT', PROGNAME )
        ENDIF
        IF( RPT_%BYRCL   ) THEN
            ALLOCATE( BINRCL   ( NOUTBINS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BINRCL', PROGNAME )
        ENDIF
        IF( RPT_%BYCELL  ) THEN
            ALLOCATE( BINX     ( NOUTBINS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BINX', PROGNAME )
        ENDIF
        IF( RPT_%BYCELL  ) THEN
            ALLOCATE( BINY     ( NOUTBINS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BINY', PROGNAME )
        ENDIF
        IF( RPT_%BYELEV  ) THEN
            ALLOCATE( BINELEV  ( NOUTBINS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BINELEV', PROGNAME )
        ENDIF

        ALLOCATE( BINPOPDIV( NOUTBINS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'BINPOPDIV', PROGNAME )
        BINPOPDIV = 1.    ! array

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
     &                COL, ROW, SRCID, FIP, SCC, SIC, SRGID1, SRGID2, 
     &                MONID, WEKID, DIUID, SPCID, PLANT, ORIS, RCL, 
     &                ESTAT, MACT, NAICS, SRCTYP

C.................  Store region code
                IF( LREGION ) BINREGN( B ) = FIP

C.................  Store country name index. Note for population that some
C                   form of "by region is required"
                IF( RPT_%BYCONAM ) THEN
                    F = ( FIP / 100000 ) * 100000
                    K = FIND1( F, NCOUNTRY, CTRYCOD )
                    BINCOIDX( B ) = K

C.....................  If using population normalization, initialize with 
C                       country population
                    IF ( RPT_%NORMPOP ) THEN
                        IF( CTRYPOPL( K ) . GT. 0. ) THEN
                            BINPOPDIV( B ) = 1. / CTRYPOPL(K)

C.........................  If population data unavailable, then flags bins
C                           that will not be able to have normalization by pop.
                        ELSE
                            BINBAD   ( B ) = 100  ! cpde for bad pop
                            BINPOPDIV( B ) = 1.
                        END IF
                    END IF
                
                END IF

C.................  Store state name index. Note for population that some
C                   form of "by region is required"
                IF( RPT_%BYSTNAM ) THEN
                    F = ( FIP / 1000 ) * 1000        ! In case by-county also
                    K = FIND1( F, NSTATE, STATCOD )
                    BINSTIDX( B ) = K

C.....................  If using population normalization, reset with 
C                       state population
                    IF ( RPT_%NORMPOP ) THEN
                        IF( STATPOPL( K ) . GT. 0. ) THEN
                            BINPOPDIV( B )= 1. / STATPOPL(K)
C.........................  If population data unavailable, then flags bins
C                           that will not be able to have normalization by pop.
                        ELSE
                            BINBAD   ( B ) = 100  ! code for bad pop
                            BINPOPDIV( B ) = 1.
                        END IF
                    END IF

                END IF

C.................  Store county name index. Note for population that some
C                   form of "by region is required"
                IF( RPT_%BYCYNAM ) THEN
                    K = FIND1( FIP, NCOUNTY, CNTYCOD )
                    BINCYIDX( B ) = K

C.....................  If using population normalization, reset with 
C                       county population
                    IF ( RPT_%NORMPOP ) THEN
                        IF( CNTYPOPL( K ) . GT. 0. ) THEN
                            BINPOPDIV( B )= 1. / CNTYPOPL(K)
C.........................  If population data unavailable, then flags bins
C                           that will not be able to have normalization by pop.
                        ELSE
                            BINBAD   ( B ) = 100  ! code for bad pop
                            BINPOPDIV( B ) = 1.
                        END IF
                    END IF

                END IF

C.................  Store road class
                IF( RPT_%BYRCL ) BINRCL( B ) = RCL

C.................  Store SMOKE ID
                IF( RPT_%BYSRC ) BINSMKID( B ) = SRCID

C.................  Store SCC
                IF( RPT_%BYSCC ) BINSCC( B ) = SCC

C.................  Store SCC name index (for full name, regardless of SCC truncation.
C                   Note: have confirmed that using the OUTSRC(J) index with CSCC maps
C                   to SCC from BUFFER properly
                IF( RPT_%SCCNAM ) THEN
                    K = FINDC( CSCC( OUTSRC(J) ), NINVSCC, INVSCC )
                    BINSNMIDX( B ) = K
                END IF

C.................  Store SIC
                IF( RPT_%BYSIC ) BINSIC( B ) = SIC

C.................  Store SIC name index
                IF( RPT_%SICNAM ) THEN
                    K = FIND1( SIC, NINVSIC, INVSIC )
                    BINSICIDX( B ) = K
                END IF

C.................  Store MACT
                IF( RPT_%BYMACT ) BINMACT( B ) = MACT

C.................  Store MACT name index
                IF( RPT_%MACTNAM ) THEN
                    K = FINDC( MACT, NINVMACT, INVMACT )
                    BINMACIDX( B ) = K
                END IF

C.................  Store NAICS
                IF( RPT_%BYNAICS ) BINNAICS( B ) = NAICS

C.................  Store NAICS name index
                IF( RPT_%NAICSNAM ) THEN
                    K = FINDC( NAICS, NINVNAICS, INVNAICS )
                    BINNAIIDX( B ) = K
                END IF

C.................  Store ORIS
                IF( RPT_%BYORIS ) BINORIS( B ) = ORIS

C.................  Store ORIS name index
                IF( RPT_%ORISNAM ) THEN
                    K = FINDC( ORIS, NORIS, ORISLST )
                    BINORSIDX( B ) = K
                END IF

C.................  Store SRCTYP
                IF( RPT_%BYSRCTYP ) BINSRCTYP( B ) = SRCTYP

C.................  Store surrogate codes
                IF( RPT_%SRGRES .EQ. 1 ) BINSRGID1( B ) = SRGID1
                IF( RPT_%SRGRES .GE. 1 ) BINSRGID2( B ) = SRGID2

C.................  Store temporal profiles
                IF( RPT_%BYMON ) BINMONID( B ) = MONID
                IF( RPT_%BYWEK ) BINWEKID( B ) = WEKID
                IF( RPT_%BYDIU ) BINDIUID( B ) = DIUID

C.................  Store speciation profiles
                IF( RPT_%BYSPC ) BINSPCID( B ) = SPCID

C.................  Store plant ID code
                IF( RPT_%BYPLANT ) THEN
                    BINPLANT( B ) = PLANT
                    BINSMKID( B ) = SRCID   ! Needed for plant names
                END IF

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
 
