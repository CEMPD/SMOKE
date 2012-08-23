
        SUBROUTINE WRPDEMIS( DAYFLAG, JDATE, JTIME, TIDX, NPDSRC, NVAR,
     &                       NVASP, FNAME, PFLAG, EAIDX, SPIDX, 
     &                       LASTSTEP, PDIDX, PDDATA, EFLAG )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C      Subroutines: I/O API subroutine
C
C  REVISION  HISTORY:
C      Created 12/99 by M. Houyoux
C
C**************************************************************************
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
C***************************************************************************

C.........  MODULES for public variables
C...........   This module is the inventory arrays
        USE MODSOURC, ONLY: CSOURC

C.........  This module contains the information about the source category
        USE MODINFO, ONLY: NSRC, NCHARS, EANAM, NCOMP, VAR_FORMULA,
     &                     CHKPLUS, CHKMINUS, 
     &                     VIN_A, VIN_B, VNAME

C.........  This module contains data for day- and hour-specific data
        USE MODDAYHR, ONLY: PDTOTL, NPDPT, IDXSRC, SPDIDA, CODEA,
     &                      EMISVA, DYTOTA, CIDXA

C.........  This module contains the lists of unique inventory information
        USE MODLISTS, ONLY: FIREFLAG, NUNIQCAS, UNIQCAS

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C.........  EXTERNAL FUNCTIONS
        CHARACTER(2) CRLF
        INTEGER      ENVINT
        LOGICAL      ENVYN

        EXTERNAL     CRLF, ENVINT, ENVYN


C.........  SUBROUTINE ARGUMENTS
        LOGICAL     , INTENT  (IN) :: DAYFLAG              ! true: day-, false: hour-spec
        INTEGER     , INTENT  (IN) :: JDATE                ! Julian date
        INTEGER     , INTENT  (IN) :: JTIME                ! time HHMMSS
        INTEGER     , INTENT  (IN) :: TIDX                 ! time index
        INTEGER     , INTENT  (IN) :: NPDSRC               ! no. part-day srcs
        INTEGER     , INTENT  (IN) :: NVAR                 ! no. pol/act vars
        INTEGER     , INTENT  (IN) :: NVASP                ! no. pol/act/special
        CHARACTER(*), INTENT  (IN) :: FNAME                ! output file name
        LOGICAL     , INTENT  (IN) :: PFLAG                ! true: gen profiles
        INTEGER     , INTENT  (IN) :: EAIDX( NVAR )        ! pol/act index
        INTEGER     , INTENT  (IN) :: SPIDX( MXSPDAT )     ! special var index
        LOGICAL     , INTENT  (IN) :: LASTSTEP             ! true: last timestep
        INTEGER     , INTENT (OUT) :: PDIDX ( NPDSRC )     ! sparse src index
        REAL        , INTENT (OUT) :: PDDATA( NPDSRC,NVASP)! sparse data storage
        LOGICAL     , INTENT (OUT) :: EFLAG                ! true: error found

C...........   Local allocatable arrays
        INTEGER, ALLOCATABLE, SAVE :: EAIDX2( : )    ! reverse index for EAIDX
        INTEGER, ALLOCATABLE, SAVE :: IDNNOTE( : )   ! flag to note for Inventory Data names renamed to same SMOKE data name
        LOGICAL, ALLOCATABLE, SAVE :: NOMISS( :,: )

C...........   Local arrays
        INTEGER, ALLOCATABLE, SAVE :: SPIDX2( : )
        INTEGER, ALLOCATABLE       :: SIDX( : )      ! start index of a source
        INTEGER, ALLOCATABLE       :: EIDX( : )      ! end index of a source

C...........   LOCAL PARAMETERS
        CHARACTER(16), PARAMETER :: FORMEVNM = 'SMKINVEN_FORMULA'

C...........   Other local variables
        INTEGER          I, J, K, L2, LN, LS, N, S, V, V2
        INTEGER          II, JJ, VV, F, CK

        INTEGER          IOS                  ! i/o status
        INTEGER, SAVE :: NWARN = 0            ! warning count
        INTEGER, SAVE :: NWARN1= 0            ! warning count
        INTEGER, SAVE :: MXEA                 ! maximum pol/var # in EAIDX
        INTEGER, SAVE :: MXWARN               ! max no. warnings
        INTEGER          NOUT                 ! tmp no. sources per time step

        LOGICAL       :: DUPFLAG  = .FALSE.  ! true: record is an actual duplicate
        LOGICAL, SAVE :: FIRSTIME = .TRUE.   ! true: first time routine called
        LOGICAL, SAVE :: HOURFLAG = .FALSE.  ! true: hour-spec
        LOGICAL, SAVE :: DFLAG    = .FALSE.  ! true: error on duplicates
        LOGICAL, SAVE :: LFLAG    = .FALSE.  ! true: iteration on special var

        CHARACTER(6  )   TYPE             ! "Hourly" or "Daily"
        CHARACTER(256)   BUFFER           ! src description buffer
        CHARACTER(300)   MESG             ! message buffer
        CHARACTER(CASLEN3) IDNAM          ! tmp Inventory Data Name
        CHARACTER(IOVLEN3) SMKNAM         ! tmp SMOKE name

        INTEGER       IDXA      ! position of first variable in source index
        INTEGER       IDXB      ! position of second in iable in source index
        INTEGER       IDXVNAM   ! position of calculated variable in source index
        INTEGER       WARNCNT   ! number of times warnings
 
        CHARACTER(16) :: PROGNAME = 'WRPDEMIS' !  program name

C***********************************************************************
C   begin body of program WRPDEMIS

C.........  If there is one or more computed output variable, get set up
C.........Allocate memory for source index 

        IF( NCOMP .GT. 0 ) THEN
            IF ( ALLOCATED( SIDX ) ) DEALLOCATE (SIDX)
            IF ( ALLOCATED( EIDX ) ) DEALLOCATE (EIDX)
            ALLOCATE( SIDX( NPDSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SIDX', PROGNAME )
            ALLOCATE( EIDX( NPDSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EIDX', PROGNAME )
            SIDX = 0 
            EIDX = 0 
        END IF

C.........  Reset FIRSTIME flag for either day-spec or hour-spec
        IF( .NOT. FIRSTIME .AND. .NOT. DAYFLAG ) FIRSTIME = .TRUE.

C.........  For the first time the routine is called...
        IF( FIRSTIME .AND. .NOT. HOURFLAG ) THEN

C.............  Get settings from the environment.
            DFLAG = ENVYN( 'RAW_DUP_CHECK',
     &                     'Error if duplicate inventory records',
     &                     .FALSE., IOS )

C.............  Get maximum number of warnings
            MXWARN = ENVINT( WARNSET , ' ', 100, I )

C.............  Allocate memory for flag for writing missing-data messages
            IF( ALLOCATED( NOMISS ) ) DEALLOCATE( NOMISS )
            ALLOCATE( NOMISS( NSRC,NVASP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'NOMISS', PROGNAME )
            NOMISS = .TRUE.  ! Array

C.............  Create reverse index for pollutants and activities
            IF (ALLOCATED (EAIDX2 )) DEALLOCATE (EAIDX2)
            MXEA = MAXVAL( EAIDX )
            ALLOCATE( EAIDX2( MXEA ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EAIDX2', PROGNAME )
            EAIDX2 = 0

            DO V = 1, NVAR
                EAIDX2( EAIDX( V ) ) = V
            END DO
 
C.............  Create reverse index for special variables
            IF (ALLOCATED (SPIDX2 )) DEALLOCATE (SPIDX2)
            ALLOCATE( SPIDX2( MXSPDAT ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SPIDX2', PROGNAME )
            SPIDX2 = 0

            DO V = 1, MXSPDAT
                K = SPIDX( V )
                IF( K .GT. 0 ) SPIDX2( K ) = V       
            END DO

C.............  Create warning array for apparent duplicates
            IF (ALLOCATED (IDNNOTE)) DEALLOCATE (IDNNOTE)
            ALLOCATE( IDNNOTE(NUNIQCAS), STAT=IOS )
            CALL CHECKMEM( IOS, 'IDNNOTE', PROGNAME )
            IDNNOTE = 0  ! Array

            FIRSTIME = .FALSE.
            IF( .NOT. DAYFLAG ) HOURFLAG = .TRUE.

        END IF

        PDIDX  = 0        ! array (index)
        PDDATA = BADVAL3  ! array (emissions/activities)
        PDTOTL = BADVAL3  ! array (total daily emissions/activities)

C.........  Sort sources & inventory data codes for current time step
C.........  Added CIDXA into sorting in SMOKE 2.3.2 to handle now using
C           inventory data codes instead of SMOKE pollutant codes only.
C           Want to prevent false duplicates caused by Inventory Table renaming 
C           multiple Inventory Data Code (i.e., CAS numbers) to the same
C           SMOKE name.  These are not really duplicates and should be
C           ignored by the warning messages below.

        CALL SORTI2( NPDPT( TIDX ), IDXSRC( 1,TIDX ), 
     &               SPDIDA( 1,TIDX ), CIDXA( 1,TIDX ) )

C.........  Store sorted records for this hour
        LS = 0  ! previous source
        LN = -9 ! previous Inventory Data Code position in UNIQCAS
        K  = 0
        DO I = 1, NPDPT( TIDX )

            J = IDXSRC( I,TIDX )
            S = SPDIDA( J,TIDX )
            V = CODEA ( J,TIDX )
            N = CIDXA ( J,TIDX )

C.............  Intialize as not a special data variable (e.g., not flow rate)
            LFLAG = .FALSE.

C.............  Check for bad index
            IF( V .LE. 0 ) THEN
                MESG = 'INTERNAL ERROR: problem indexing input '//
     &                 'pollutants to output pollutants.'
                CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )

C.............  Check for index for special variables and set V to be consistent
C               with the output structure of the file
            ELSE IF ( V .GT. CODFLAG3 ) THEN   ! CODFLAG3 = 9000 for special data types
                V2 = V - CODFLAG3      ! remove flag on index
                V = NVAR + SPIDX( V2 ) ! reset to condensed order from master
                LFLAG = .TRUE.         ! flag as a special data variable

C.............  Otherwise, set index for period-specific pollutant or activity
            ELSE
                V = EAIDX2( V )

            END IF

C.............  Initialize duplicates flag
            DUPFLAG = .FALSE.

C.............  If current source is not equal to previous source
            IF( S .NE. LS  .OR. I .EQ.  NPDPT( TIDX )) THEN
                IF( S .NE. LS ) THEN
                    K = K + 1
                    PDIDX( K ) = S
                    LS         = S
                END IF 

C...............Calculate formula if needed  
C...............Get the location of start and end index 
                IF ( NCOMP > 0 ) THEN
                    IF ( I .EQ. 1 ) THEN 
                        SIDX( K ) = I
                    ELSE IF ( I .GT. 1 .AND. I .LT. NPDPT( TIDX )) THEN 
                        SIDX( K ) = I
                        EIDX( K-1 ) = I - 1
                        CK = K-1
C...............Only specify the end index for last source 
                    ELSE IF ( I .EQ.  NPDPT( TIDX )) THEN
                        CK = K
                        EIDX( K ) = I
                    END IF 

                    IF ( I .GT. 1 ) THEN 
                      DO F = 1, NCOMP
                        IDXA = 0
                        IDXB = 0
                        IDXVNAM = 0
                        DO II = SIDX(CK), EIDX(CK)    
                            JJ = IDXSRC( II,TIDX )
                            VV = CODEA ( JJ,TIDX )
                            IF ( EANAM(VV) .EQ. VIN_A(F) ) THEN
                                IDXA = JJ
                            ELSE IF ( EANAM(VV) .EQ. VIN_B(F) ) THEN
                                IDXB = JJ
                            END IF
                            IDXVNAM = NVASP - NCOMP + F 
                        END DO    ! end of variable search loop

C...............Calculate formula variables
                        IF ( IDXA .GT. 0 .AND. IDXB .GT. 0 
     &                       .AND. IDXVNAM .GT. 0 ) THEN
                            IF ( CHKPLUS(F) )  THEN 
                                 IF ( PDDATA( CK,IDXVNAM ) .GT. 0 ) THEN
                                 PDDATA( CK,IDXVNAM ) =  PDDATA( CK,IDXVNAM ) +
     &                              EMISVA( IDXA, TIDX ) + EMISVA( IDXB, TIDX ) 
                                 ELSE
                                     PDDATA( CK,IDXVNAM ) =
     &                                  EMISVA( IDXA, TIDX ) + EMISVA( IDXB, TIDX ) 
                                 END IF
                            END IF
                            IF ( CHKMINUS(F) )  THEN 
                                IF ( PDDATA( CK,IDXVNAM ) .GT. 0 ) THEN
                                    PDDATA( CK,IDXVNAM ) = PDDATA( CK,IDXVNAM )+
     &                                  EMISVA( IDXA, TIDX ) - EMISVA( IDXB, TIDX ) 
                                ELSE
                                    PDDATA( CK,IDXVNAM ) =
     &                                  EMISVA( IDXA, TIDX ) - EMISVA( IDXB, TIDX ) 
                                END IF
                            END IF

C..........................  Check for negative values for daily value
                            IF( PDDATA( CK,IDXVNAM ) .LT. 0 )  THEN
                                WARNCNT = WARNCNT + 1
                                PDDATA( CK,IDXVNAM ) = 0.0

                                IF ( WARNCNT .LE. MXWARN ) THEN
                                    CALL FMTCSRC( CSOURC( S ), 7, BUFFER, L2 )
                                    WRITE( MESG,94020 ) 'WARNING: '//
     &                               'Resetting negative value of "'//
     &                               'average-day "'//TRIM(VNAME(F))// 
     &                               '" from',PDDATA( CK,IDXVNAM ),'to 0. for '//
     &                              'source:'//CRLF()//BLANK10//BUFFER(1:L2)
                                    MXWARN = MXWARN + 1
                                END IF 
                                CALL M3MESG( MESG )
                            END IF 

                        ELSE IF ( WARNCNT .LE. MXWARN ) THEN
                            CALL FMTCSRC( CSOURC( S ), NCHARS, BUFFER, L2 )
                            EFLAG = .TRUE.
                            MESG = 'ERROR: no sources for calculating ' //
     &                         'variable "'//TRIM(VNAME(F))//'" J '//
     &                          CRLF() // BLANK10 // BUFFER( 1:L2 )
                            CALL M3MSG2( MESG )
                            WARNCNT = WARNCNT + 1
                            CYCLE
                        END IF
                       
                      END DO    ! end of formula loop
                    END IF
                END IF   ! end of formula calculation

C............. If source is the same, look for duplicate data variables
            ELSE IF ( N .EQ. LN ) THEN
                DUPFLAG = .TRUE.      !  This iteration has a duplicate
            END IF

C.............  If emissions are not yet set for current source and variable
            IF( PDDATA( K,V ) .LT. AMISS3 ) THEN  ! PDDATA is blank

                PDDATA( K,V ) = EMISVA( J,TIDX )
                PDTOTL( K,V ) = DYTOTA( J,TIDX )

C.............  Otherwise, sum emissions, report, and set error if needed
            ELSE 

                IF( .NOT. LFLAG .AND. EMISVA( J,TIDX ) .GE. 0. )
     &              PDDATA( K,V )  = PDDATA( K,V ) + EMISVA( J,TIDX )

                IF( .NOT. LFLAG .AND. DYTOTA( J,TIDX ) .GE. 0. )
     &              PDTOTL( K,V )  = PDTOTL( K,V ) + DYTOTA( J,TIDX )

C.................  If the source is an actual duplicate (the same source
C                   and Inventory Data Name), proceed accordingly.  Do
C                   not write warnings for Inventory Data Names that
C                   are combined into the same SMOKE name.
                IF ( DUPFLAG ) THEN
                    CALL FMTCSRC( CSOURC( S ), NCHARS, BUFFER, L2 )
                    NWARN1 = NWARN1 + 1

                    IF( NWARN1 <= MXWARN ) THEN
                    IF( DFLAG ) THEN
                        EFLAG = .TRUE.
                        MESG = 'ERROR: Duplicate source in inventory:'//
     &                         CRLF() // BLANK10 // BUFFER( 1:L2 )
                        CALL M3MESG( MESG )
                        CYCLE
                    ELSE IF ( LFLAG ) THEN
                        MESG = 'WARNING: Duplicate source in ' //
     &                         'inventory. Will store only one value '//
     &                         'for '// CRLF()// BLANK10//SPDATDSC(V2)//
     &                         ':' // CRLF() // BLANK10 // BUFFER(1:L2)
                    ELSE
                        MESG = 'WARNING: Duplicate source in ' //
     &                         'inventory will have summed emissions:' 
     &                         //CRLF() // BLANK10 // BUFFER( 1:L2 )
                        CALL M3MESG( MESG )
                    END IF
                    END IF

C.................  If multiple entries for this K, V then record for later
C                   note on this Inventory Data Name.  This section is 
C                   encountered only if multiple Inventory Date Names
C                   for this source/hour stored as same SMOKE name.
                ELSE
                    IDNNOTE( N ) = V 

                END IF   !  If duplicate or not
            END IF

            LN = N  ! Set LN for next iteration

        END DO  ! End loop over data for this time step

C.........  Set tmp variable for loops
        NOUT = K

C.........  Check if there are missing values and output errors, if the
C           flag is set to treat these as errors
C.........  Also use loop to create diurnal profiles from emission values,
C           if needed.
        DO V = 1, NVASP
            DO I = 1, NOUT

                S = PDIDX( I )

C.................  Format source information
                CALL FMTCSRC( CSOURC( S ), NCHARS, BUFFER, L2 )

C.................  Check for missing values
                IF ( PDDATA( I,V ) .LT. AMISS3 ) THEN

C.....................  For fires data, reset "missing" values
C                       to zero, since Temporal will not need to
C                       know the difference between missing and zero.
C                       This is because there are no "annual" values
C                       to fill in for missing daily/hourly data for fires.
                    IF( FIREFLAG ) PDDATA( I,V ) = 0.

                    IF( NOMISS( S,V ) ) THEN
                        NWARN = NWARN + 1

                        IF( NWARN <= MXWARN ) THEN
                            IF ( V .LE. NVAR ) THEN
                                MESG = 'WARNING: Data missing for: ' //
     &                             CRLF()//BLANK10//BUFFER( 1:L2 )//
     &                             ' VAR: '// EANAM( EAIDX( V ) )
                            ELSE
                                K = V - NVAR
                                MESG = 'WARNING: Data missing for: ' //
     &                             CRLF()//BLANK10//BUFFER( 1:L2 )//
     &                             ' VAR: '//SPDATNAM( SPIDX2( K ) )
                            END IF

                            CALL M3MESG( MESG )
                        END IF
                        NOMISS( S,V ) = .FALSE.
                    END IF

                END IF

C.................  Skip next section if special data variable
                IF ( V .GT. NVAR ) CYCLE

C.................  If profiles need to be created instead of hourly emissions
C.................  Make sure totals data are available!
                IF ( PFLAG .AND. PDTOTL( I,V ) .LT. AMISS3 ) THEN
                    EFLAG = .TRUE.
                    MESG = 'ERROR: Cannot create profiles ' //
     &                     'because daily total data missing for:' //
     &                     CRLF() // BLANK10 // BUFFER( 1:L2 ) //
     &                     ' VAR: ' // EANAM( EAIDX( V ) )
                    CALL M3MESG( MESG )
                    CYCLE

C.................  Prevent divide by zero
C.................  If totals are zero, then profile values should be zero
                ELSE IF( PFLAG .AND. PDTOTL( I,V ) .GT. 0 ) THEN

                    PDDATA( I,V ) = PDDATA( I,V ) / PDTOTL( I,V )

                END IF

            END DO
        END DO
  
C.........  Write emissions for this time step

        IF ( .NOT. WRITE3( FNAME, ALLVAR3, JDATE, JTIME, PDIDX ) ) THEN

            L2   = LEN_TRIM( FNAME )
            MESG= 'Error writing output file "' // FNAME(1:L2) // '"'
            CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )

        END IF        

C.........  If this is the last time step, then write notes
        IF( LASTSTEP ) THEN
            TYPE = 'Hourly'
            IF( DAYFLAG ) TYPE = 'Daily' 
            DO I = 1, NUNIQCAS
                IF( IDNNOTE( I ) .GT. 0 ) THEN
                    IDNAM = UNIQCAS( I )
                    SMKNAM = EANAM( IDNNOTE( I ) )
                    MESG= 'NOTE: '//TRIM( TYPE )// ' inventory data "'//
     &                    TRIM( IDNAM ) // '" summed with other '//
     &                    'pollutants to compute "'//TRIM(SMKNAM)//'".'
                    CALL M3MSG2( MESG )
                END IF
            END DO
        END IF
 
        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )
94020   FORMAT( 10( A, :, E10.2, :, 1X ) )


C******************  INTERNAL SUBPROGRAMS  *****************************

        END SUBROUTINE WRPDEMIS
