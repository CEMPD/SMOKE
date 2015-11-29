
        SUBROUTINE RDCFPRO( CFDEV )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C       Reads control factor data from CFPRO file
C
C  PRECONDITIONS REQUIRED:
C       CFDEV must be opened
C
C  SUBROUTINES AND FUNCTIONS CALLED:  none
C
C  REVISION  HISTORY:
C     06/12: Created by Dongmei Yang
C     06/15: Updated by BH Baek
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

C.........  MODULES for public variables
C.........  This module contains the major data structure and control flags
        USE MODMERGE, ONLY:  NIPPA, EANAM, NMSPC, EMNAM, NSMATV,
     &                       TSVDESC

C.........  This module contains data structures and flags specific to Movesmrg
        USE MODMVSMRG, ONLY: CFPRO, EXPCFFLAG, REFCFFLAG

C.........  This module contains the lists of unique source characteristics
        USE MODLISTS, ONLY: NINVIFIP, INVCFIP, NINVSCC, INVSCC

C.........  This module is used for reference county information
        USE MODMBSET, ONLY: NREFC, MCREFIDX, NINVC, MCREFSORT

       
        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        LOGICAL       BLKORCMT
        LOGICAL       CHKINT
        LOGICAL       CHKREAL
        INTEGER       GETFLINE
        INTEGER       INDEX1
        INTEGER       FIND1
        INTEGER       FINDC
        INTEGER       STR2INT
        INTEGER       ENVINT
        REAL          STR2REAL
        CHARACTER(2)  CRLF
        
        EXTERNAL BLKORCMT, CHKINT, CHKREAL, FIND1, GETFLINE, 
     &           STR2INT, STR2REAL, CRLF, INDEX1, ENVINT

C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT (IN) :: CFDEV             ! CFPRO file unit no.

C...........   Local allocatable arrays
        CHARACTER(20)  SEGMENT( 10 )          ! parsed input line

        INTEGER,            ALLOCATABLE :: NLFIPS( : )     ! FIPS matched
        INTEGER,            ALLOCATABLE :: NLSCCS( : )     ! SCC matched
        INTEGER,            ALLOCATABLE :: NLPOLS( : )     ! POL matched
        INTEGER,            ALLOCATABLE :: NLMONS( : )     ! month matched
        INTEGER           , ALLOCATABLE :: INDXCF( : )     ! index for search matched x-ref
        
        CHARACTER(FPSLEN3+CASLEN3+2), ALLOCATABLE :: CHRT01( : )
        CHARACTER(FPSLEN3+CASLEN3+2), ALLOCATABLE :: CHRT02( : )
        CHARACTER(FPSLEN3+CASLEN3+2), ALLOCATABLE :: CHRT03( : )
        CHARACTER(FPSLEN3+CASLEN3+2), ALLOCATABLE :: CHRT04( : )
        CHARACTER(FPSLEN3+CASLEN3+2), ALLOCATABLE :: CHRT05( : )
        CHARACTER(FPSLEN3+CASLEN3+2), ALLOCATABLE :: CHRT06( : )
        CHARACTER(FPSLEN3+CASLEN3+2), ALLOCATABLE :: CHRT07( : )
        CHARACTER(FPSLEN3+CASLEN3+2), ALLOCATABLE :: CHRT08( : )
        CHARACTER(FPSLEN3+CASLEN3+2), ALLOCATABLE :: CHRT09( : )
        CHARACTER(FPSLEN3+CASLEN3+2), ALLOCATABLE :: CHRT10( : )
        CHARACTER(FPSLEN3+CASLEN3+2), ALLOCATABLE :: CHRT11( : )
        CHARACTER(FPSLEN3+CASLEN3+2), ALLOCATABLE :: CHRT12( : )
        CHARACTER(FPSLEN3+CASLEN3+2), ALLOCATABLE :: CHRT13( : )
        CHARACTER(FPSLEN3+CASLEN3+2), ALLOCATABLE :: CHRT14( : )
        CHARACTER(FPSLEN3+CASLEN3+2), ALLOCATABLE :: CHRT15( : )
        CHARACTER(FPSLEN3+CASLEN3+2), ALLOCATABLE :: CHRT16( : )
        CHARACTER(FPSLEN3+CASLEN3+2), ALLOCATABLE :: CHRT17( : )
        CHARACTER(FPSLEN3+CASLEN3+2), ALLOCATABLE :: CHRT18( : )
        CHARACTER(FPSLEN3+CASLEN3+2), ALLOCATABLE :: CHRT19( : )
        CHARACTER(FPSLEN3+CASLEN3+2), ALLOCATABLE :: CHRT20( : )
        CHARACTER(FPSLEN3+CASLEN3+2), ALLOCATABLE :: DUPCHECK( : )   ! store duplicate entries
        CHARACTER(300),     ALLOCATABLE :: CFIPSCCLIST( : )          ! store matched entries

C...........   Other local variables
        INTEGER         I, J, K, L, M, N, NS, NX, NSCC, L1, L2    ! counters and indexes
        INTEGER         IOS         ! error status
        INTEGER      :: IREC = 0    ! record counter
        INTEGER         NFIPS       ! total matched FIPS
        INTEGER         NSCCS       ! total matched SCC
        INTEGER         NPOLS       ! total matched POL 
        INTEGER         POLIDX      ! current POL index
        INTEGER         SPCIDX      ! current species index
        INTEGER         SCCIDX      ! current scc index
        INTEGER         NLINES      ! number of lines
        INTEGER         NXREF       ! max no of matached x-ref entries
        INTEGER         IFIP        ! current FIPS code
        INTEGER         MON, NMONS  ! current Month
        INTEGER         IDUM        ! tmp dummy integer
        INTEGER         MXWARN      ! maximum number of warnings
        INTEGER      :: NWARN = 0   ! current number of warnings
        INTEGER         ISTA, PSTA, NSTA
        INTEGER         NCHRT01, NCHRT02, NCHRT03, NCHRT04, NCHRT05
        INTEGER         NCHRT06, NCHRT07, NCHRT08, NCHRT09, NCHRT10
        INTEGER         NCHRT11, NCHRT12, NCHRT13, NCHRT14, NCHRT15
        INTEGER         NCHRT16, NCHRT17, NCHRT18, NCHRT19, NCHRT20
        INTEGER         FF, F1, F2, F3, F4, F5, F6, F7, F8, F9, F10
        INTEGER         F11, F12, F13, F14, F15, F16, F17, F18, F19, F20
        
        REAL            CFVAL       ! control factor value

        LOGICAL, SAVE:: FIRSTIME = .TRUE. ! true: firstime
        LOGICAL      :: EFLAG = .FALSE.   ! true: error found
        LOGICAL      :: PFLAG, SKIPREC    ! true: skip pollutant 

        CHARACTER(300)     LINE     ! line buffer
        CHARACTER(300)     MESG     ! message buffer
        CHARACTER(2)       CMON, BLKMON, FILMON
        CHARACTER(FPSLEN3+CASLEN3+2) DUPTMP      ! tmp duplicate check buffer
        CHARACTER(FIPLEN3) CFIP, BLKFIP          ! tmp (character) FIPS code : CHRT03
        CHARACTER(SCCLEN3) CSCC, BLKSCC          ! current SCC : CHRT04
        CHARACTER(IOVLEN3) POLNAM, SPCNAM, CBUF, CPOL, BLKPOL, FILPOL  ! current pollutant-species name 
        CHARACTER(FPSLEN3+CASLEN3+2) CFIPSCC, BLKFIPSCC    ! tmp FIPS code // SCC : CHRT06

        CHARACTER(16) :: PROGNAME = 'RDCFPRO'    ! program name

C***********************************************************************
C   begin body of subroutine RDCFPRO

        IF( FIRSTIME ) THEN
C.............  Get maximum number of warnings
            MXWARN = ENVINT( WARNSET, ' ', 100, IOS )

C.............  No of states
            NSTA = 0
            PSTA = 0
            DO I = 1 , NINVIFIP
                ISTA = STR2INT( INVCFIP( I )( 16:17 ) )
                IF( ISTA /= PSTA ) THEN
                    NSTA = NSTA + 1
                    PSTA = ISTA
                END IF
            END DO

C.............  Allocate locate FIPS/SCC hierarchy arrays
            NCHRT01 = NSTA               ! SCC=0,   FIP=state, Mon=0, Pol=0
            NCHRT02 = NSTA               ! SCC=0,   FIP=state, Mon=0, Pol=pol
            NCHRT03 = NSTA               ! SCC=0,   FIP=state, Mon=mon, Pol=0
            NCHRT04 = NSTA               ! SCC=0,   FIP=state, Mon=mon, Pol=pol
            NCHRT05 = NINVIFIP           ! SCC=0,   FIP=all, Mon=0, Pol=0
            NCHRT06 = NINVIFIP           ! SCC=0,   FIP=all, Mon=0, Pol=pol
            NCHRT07 = NINVIFIP           ! SCC=0,   FIP=all, Mon=mon, Pol=0
            NCHRT08 = NINVIFIP           ! SCC=0,   FIP=all, Mon=mon, Pol=pol
            NCHRT09 = NINVSCC            ! SCC=all, FIP=0, Mon=0, Pol=0
            NCHRT10 = NINVSCC            ! SCC=all, FIP=0, Mon=0, Pol=pol
            NCHRT11 = NINVSCC            ! SCC=all, FIP=0, Mon=mon, Pol=0
            NCHRT12 = NINVSCC            ! SCC=all, FIP=0, Mon=mon, Pol=pol
            NCHRT13 = NSTA * NINVSCC     ! SCC=all, FIP=state, Mon=0, Pol=0
            NCHRT14 = NSTA * NINVSCC     ! SCC=all, FIP=state, Mon=0, Pol=pol
            NCHRT15 = NSTA * NINVSCC     ! SCC=all, FIP=state, Mon=mon, Pol=0
            NCHRT16 = NSTA * NINVSCC     ! SCC=all, FIP=state, Mon=mon, Pol=pol
            NCHRT17 = NINVIFIP * NINVSCC ! SCC=all, FIP=all, Mon=0, Pol=0
            NCHRT18 = NINVIFIP * NINVSCC ! SCC=all, FIP=all, Mon=0, Pol=pol
            NCHRT19 = NINVIFIP * NINVSCC ! SCC=all, FIP=all, Mon=mon, Pol=0
            NCHRT20 = NINVIFIP * NINVSCC ! SCC=all, FIP=all, Mon=mon, Pol=pol

            ALLOCATE( CHRT01( NCHRT01 ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CHRT01', PROGNAME )
            ALLOCATE( CHRT02( NCHRT02 ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CHRT02', PROGNAME )
            ALLOCATE( CHRT03( NCHRT03 ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CHRT03', PROGNAME )
            ALLOCATE( CHRT04( NCHRT04 ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CHRT04', PROGNAME )
            ALLOCATE( CHRT05( NCHRT05 ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CHRT05', PROGNAME )
            ALLOCATE( CHRT06( NCHRT06 ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CHRT06', PROGNAME )
            ALLOCATE( CHRT07( NCHRT07 ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CHRT07', PROGNAME )
            ALLOCATE( CHRT08( NCHRT08 ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CHRT08', PROGNAME )
            ALLOCATE( CHRT09( NCHRT09 ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CHRT09', PROGNAME )
            ALLOCATE( CHRT10( NCHRT10 ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CHRT10', PROGNAME )
            ALLOCATE( CHRT11( NCHRT11 ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CHRT11', PROGNAME )
            ALLOCATE( CHRT12( NCHRT12 ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CHRT12', PROGNAME )
            ALLOCATE( CHRT13( NCHRT13 ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CHRT13', PROGNAME )
            ALLOCATE( CHRT14( NCHRT14 ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CHRT14', PROGNAME )
            ALLOCATE( CHRT15( NCHRT15 ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CHRT15', PROGNAME )
            ALLOCATE( CHRT16( NCHRT16 ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CHRT16', PROGNAME )
            ALLOCATE( CHRT17( NCHRT17 ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CHRT17', PROGNAME )
            ALLOCATE( CHRT18( NCHRT18 ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CHRT18', PROGNAME )
            ALLOCATE( CHRT19( NCHRT19 ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CHRT19', PROGNAME )
            ALLOCATE( CHRT20( NCHRT20 ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CHRT20', PROGNAME )
            CHRT01 = ''
            CHRT02 = ''
            CHRT03 = ''
            CHRT04 = ''
            CHRT05 = ''
            CHRT06 = ''
            CHRT07 = ''
            CHRT08 = ''
            CHRT09 = ''
            CHRT10 = ''
            CHRT11 = ''
            CHRT12 = ''
            CHRT13 = ''
            CHRT14 = ''
            CHRT15 = ''
            CHRT16 = ''
            CHRT17 = ''
            CHRT18 = ''
            CHRT19 = ''
            CHRT20 = ''

C.............  Build FIPS and SCC hierarchy tables
            N = 0
            BLKMON = '00'
            FILMON = 'MM'
            BLKFIP = '000000000000'
            BLKSCC = '00000000000000000000'
            BLKPOL = '0000000000000000'
            FILPOL = 'PPPPPPPPPPPPPPPP'
            BLKFIPSCC = BLKFIP // BLKSCC // BLKMON // BLKPOL

            DO I = 1, NINVSCC

                CHRT09( I ) = BLKFIP // INVSCC( I ) // BLKMON // BLKPOL
                CHRT10( I ) = BLKFIP // INVSCC( I ) // BLKMON // FILPOL
                CHRT11( I ) = BLKFIP // INVSCC( I ) // FILMON // BLKPOL
                CHRT12( I ) = BLKFIP // INVSCC( I ) // FILMON // FILPOL

                NSTA = 0
                PSTA = -1 
                DO J = 1, NINVIFIP

                    N = N + 1
                    CHRT05(J) = INVCFIP(J) // BLKSCC // BLKMON // BLKPOL
                    CHRT06(J) = INVCFIP(J) // BLKSCC // BLKMON // FILPOL
                    CHRT07(J) = INVCFIP(J) // BLKSCC // FILMON // BLKPOL
                    CHRT08(J) = INVCFIP(J) // BLKSCC // FILMON // FILPOL

                    CHRT17(N) = INVCFIP(J) // INVSCC(I) // BLKMON // BLKPOL
                    CHRT18(N) = INVCFIP(J) // INVSCC(I) // BLKMON // FILPOL
                    CHRT19(N) = INVCFIP(J) // INVSCC(I) // FILMON // BLKPOL
                    CHRT20(N) = INVCFIP(J) // INVSCC(I) // FILMON // FILPOL

                    NS = 0
                    ISTA = STR2INT( INVCFIP( J )( 16:17 ) )

                    IF( ISTA /= PSTA ) THEN

                        NSTA = NSTA + 1
                        WRITE( CHRT01( NSTA ),'(I12.12,A)' ) ISTA, BLKSCC // BLKMON // BLKPOL
                        WRITE( CHRT02( NSTA ),'(I12.12,A)' ) ISTA, BLKSCC // BLKMON // FILPOL
                        WRITE( CHRT03( NSTA ),'(I12.12,A)' ) ISTA, BLKSCC // FILMON // BLKPOL
                        WRITE( CHRT04( NSTA ),'(I12.12,A)' ) ISTA, BLKSCC // FILMON // FILPOL

                        DO K = 1, NINVSCC
                            NS = NS + 1
                            WRITE( CHRT13(NS),'(I12.12,A)' ) ISTA, INVSCC(K) // BLKMON // BLKPOL
                            WRITE( CHRT14(NS),'(I12.12,A)' ) ISTA, INVSCC(K) // BLKMON // FILPOL
                            WRITE( CHRT15(NS),'(I12.12,A)' ) ISTA, INVSCC(K) // FILMON // BLKPOL
                            WRITE( CHRT16(NS),'(I12.12,A)' ) ISTA, INVSCC(K) // FILMON // FILPOL

                        END DO

                        PSTA = ISTA

                    END IF

                END DO    ! FIPS loop

            END DO     ! SCC loop

C.............  Allocate storage based on number of FIPs and SCCs in inventory
            ALLOCATE( CFPRO( NINVIFIP, NINVSCC, NIPPA+NMSPC, 12 ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CFPRO', PROGNAME )
            CFPRO = 1.0   ! array
            
C.............  Get the number of lines in the file
            NLINES = GETFLINE( CFDEV, 'CFPRO file' )
            ALLOCATE( CFIPSCCLIST( NLINES ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CFIPSCCLIST', PROGNAME )
            ALLOCATE( INDXCF( NLINES ), STAT=IOS )
            CALL CHECKMEM( IOS, 'INDXCF', PROGNAME )
            ALLOCATE( DUPCHECK( NLINES ), STAT=IOS )
            CALL CHECKMEM( IOS, 'DUPCHECK', PROGNAME )
            ALLOCATE( NLFIPS( NINVIFIP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'NLFIPS', PROGNAME )
            ALLOCATE( NLSCCS( NINVSCC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'NLSCCS', PROGNAME )
            ALLOCATE( NLPOLS( NIPPA+NMSPC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'NLPOLS', PROGNAME )
            ALLOCATE( NLMONS( 12 ), STAT=IOS )
            CALL CHECKMEM( IOS, 'NLMONS', PROGNAME )
            CFIPSCCLIST = ''
            DUPCHECK   = ''
            INDXCF = 0
            NLFIPS = 0
            NLSCCS = 0
            NLPOLS = 0
            NLMONS = 0

            FIRSTIME = .FALSE.

        END IF

C.........  Read through file and store hourly data
        NX = 0
        DO I = 1, NLINES

            READ( CFDEV, 93000, END=999, IOSTAT=IOS ) LINE
            
            IREC = IREC + 1
            
            IF( IOS .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG, 94010 ) 'I/O error', IOS,
     &            'reading control factor file at line', IREC
                CALL M3MESG( MESG )
                CYCLE
            END IF

C.............  Skip blank or comment lines
            IF( BLKORCMT( LINE ) ) CYCLE

C.............  Parse line into fields
            CALL PARSLINE( LINE, 6, SEGMENT )
            
            CFIP   = TRIM( SEGMENT( 1 ) )
            CSCC   = TRIM( SEGMENT( 2 ) )
            POLNAM = ADJUSTL( SEGMENT( 3 ) )
            MON    = STR2INT( SEGMENT( 4 ) )

            IF( POLNAM == ' ' ) THEN
                CPOL = BLKPOL
            ELSE
                CPOL = FILPOL
            END IF

            IF( MON < 1 ) THEN
                MON = 0
                CMON = BLKMON
            ELSE
                CMON = FILMON
            END IF
            
            CBUF = ' '
            CALL FLTRXREF( CFIP, CBUF, CSCC, POLNAM, CBUF,
     &                     IDUM, IDUM, IDUM, PFLAG, SKIPREC )

            CFIPSCC = CFIP // CSCC // CMON // CPOL

            F1 = INDEX1( CFIPSCC, NCHRT01, CHRT01 )  ! SCC=0,   FIP=state, Mon=0, Pol=0
            F2 = INDEX1( CFIPSCC, NCHRT02, CHRT02 )  ! SCC=0,   FIP=state, Mon=0, Pol=pol
            F3 = INDEX1( CFIPSCC, NCHRT03, CHRT03 )  ! SCC=0,   FIP=state, Mon=mon, Pol=0
            F4 = INDEX1( CFIPSCC, NCHRT04, CHRT04 )  ! SCC=0,   FIP=state, Mon=mon, Pol=pol
            F5 = INDEX1( CFIPSCC, NCHRT05, CHRT05 )  ! SCC=0,   FIP=all, Mon=0, Pol=0
            F6 = INDEX1( CFIPSCC, NCHRT06, CHRT06 )  ! SCC=0,   FIP=all, Mon=0, Pol=pol
            F7 = INDEX1( CFIPSCC, NCHRT07, CHRT07 )  ! SCC=0,   FIP=all, Mon=mon, Pol=0
            F8 = INDEX1( CFIPSCC, NCHRT08, CHRT08 )  ! SCC=0,   FIP=all, Mon=mon, Pol=pol
            F9 = INDEX1( CFIPSCC, NCHRT09, CHRT09 )  ! SCC=all, FIP=0, Mon=0, Pol=0
            F10 = INDEX1( CFIPSCC, NCHRT10, CHRT10 ) ! SCC=all, FIP=0, Mon=0, Pol=pol
            F11 = INDEX1( CFIPSCC, NCHRT11, CHRT11 ) ! SCC=all, FIP=0, Mon=mon, Pol=0
            F12 = INDEX1( CFIPSCC, NCHRT12, CHRT12 ) ! SCC=all, FIP=0, Mon=mon, Pol=pol
            F13 = INDEX1( CFIPSCC, NCHRT13, CHRT13 ) ! SCC=all, FIP=state, Mon=0, Pol=0
            F14 = INDEX1( CFIPSCC, NCHRT14, CHRT14 ) ! SCC=all, FIP=state, Mon=0, Pol=pol
            F15 = INDEX1( CFIPSCC, NCHRT15, CHRT15 ) ! SCC=all, FIP=state, Mon=mon, Pol=0
            F16 = INDEX1( CFIPSCC, NCHRT16, CHRT16 ) ! SCC=all, FIP=state, Mon=mon, Pol=pol
            F17 = INDEX1( CFIPSCC, NCHRT17, CHRT17 ) ! SCC=all, FIP=all, Mon=0, Pol=0
            F18 = INDEX1( CFIPSCC, NCHRT18, CHRT18 ) ! SCC=all, FIP=all, Mon=0, Pol=pol
            F19 = INDEX1( CFIPSCC, NCHRT19, CHRT19 ) ! SCC=all, FIP=all, Mon=mon, Pol=0
            F20 = INDEX1( CFIPSCC, NCHRT20, CHRT20 ) ! SCC=all, FIP=all, Mon=mon, Pol=pol

            FF = -1
            IF( F20 > 0 ) FF = 20
            IF( F19 > 0 ) FF = 19
            IF( F18 > 0 ) FF = 18
            IF( F17 > 0 ) FF = 17
            IF( F16 > 0 ) FF = 16
            IF( F15 > 0 ) FF = 15
            IF( F14 > 0 ) FF = 14
            IF( F13 > 0 ) FF = 13
            IF( F12 > 0 ) FF = 12
            IF( F11 > 0 ) FF = 11
            IF( F10 > 0 ) FF = 10
            IF( F9 > 0 ) FF = 9
            IF( F8 > 0 ) FF = 8
            IF( F7 > 0 ) FF = 7
            IF( F6 > 0 ) FF = 6
            IF( F5 > 0 ) FF = 5
            IF( F4 > 0 ) FF = 4
            IF( F3 > 0 ) FF = 3
            IF( F2 > 0 ) FF = 2
            IF( F1 > 0 ) FF = 1

            IF( CFIPSCC == BLKFIPSCC ) FF = 0

            IF( FF < 0 ) CYCLE     ! skip unmatched entries
            NX = NX + 1

            WRITE( DUPTMP,'(A,I2.2)' ) CFIP//CSCC//POLNAM,MON
            L = INDEX1( DUPTMP, NLINES, DUPCHECK )
            IF( L > 0 ) THEN
                WRITE( MESG,94010 ) 'ERROR: Duplicate entry found at', IREC
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            ELSE
                DUPCHECK( NX ) = DUPTMP
            END IF    

            INDXCF( NX ) = FF
            CFIPSCCLIST( NX ) = TRIM( LINE )

C.............  No of CFPRO entries check
            IF( NX == 0 ) THEN
                MESG = 'ERROR: No applicable entry found from CFPRO input file'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 ) 
            END IF
            
        END DO

        NXREF = NX

C.........  Processing a list of valid CFPRO entries
        DO N = 0, 20     ! six hierarchy

          DO I = 1, NXREF
        
            IF( INDXCF( I ) /= N ) CYCLE   ! process low hierarchy first

            NFIPS = 0
            NSCCS = 0
            NPOLS = 0
       
C.............  Parse line into fields
            CALL PARSLINE( CFIPSCCLIST( I ), 6, SEGMENT )
            
            CFIP   = TRIM( SEGMENT( 1 ) )
            CSCC   = TRIM( SEGMENT( 2 ) )
            POLNAM = ADJUSTL( SEGMENT( 3 ) )

            CBUF = ' '
            CALL FLTRXREF( CFIP, CBUF, CSCC, POLNAM, CBUF,
     &                     IDUM, IDUM, IDUM, PFLAG, SKIPREC )

C.............  Convert FIP to integer
            K = 0
            IFIP = STR2INT( CFIP )
            IF( IFIP == 0  ) THEN 

C.....................  State-level is not applicable when REF_CFPRO_YN is set to Y
                IF( REFCFFLAG ) THEN
                    MESG = 'WARNING: FIPS code '//CFIP//' is not a reference county'
                    CALL M3MESG( MESG )
                END IF

                NFIPS = NINVIFIP
                DO J = 1, NINVIFIP
                    NLFIPS( J ) = J
                END DO

            ELSE         ! FIPS is not zero
                CFIP = ADJUSTR( SEGMENT( 1 ) )
                CALL PADZERO( CFIP )
                L1 = FINDC( CFIP, NINVIFIP, INVCFIP )

                IF( L1 > 0 ) THEN 
                    NFIPS = 1
                    NLFIPS(1) = L1

C.....................  Propagate reference-county-specific control factor to inventory counties.
                    IF( REFCFFLAG ) THEN
                        L2 = FINDC( CFIP, NREFC, MCREFIDX( :,1 ) )
                        IF( L2 < 1 ) THEN
                            MESG = 'WARNING: FIPS code '//CFIP//' is not a reference county'
                            CALL M3MESG( MESG )
                            CYCLE
                        END IF

C.........................  find ref county and apply CF to ref-inventory counties
                        DO J = 1, NINVC
                            IF( CFIP == MCREFSORT( J,2 ) ) THEN  ! found matched ref county
                                L1 = FINDC( MCREFSORT( J,1 ), NINVIFIP, INVCFIP )
                                IF( L1 > 0 ) THEN
                                    K = K + 1
                                    NLFIPS( K ) = L1
                                END IF
                            END IF
                        END DO
                        NFIPS = K
                    END IF
                ELSE     ! FIPS is country/state  

C.....................  State-level is not applicable when REF_CFPRO_YN is set to Y
                    IF( CFIP( 18:20 ) /= '000' ) THEN
                        MESG = 'WARNING: FIPS '//CFIP//' is not an inventory county'
                        CALL M3MESG( MESG )
                        CYCLE
                    END IF

                    DO J = 1, NINVIFIP 
                        IF( CFIP(16:17) == INVCFIP(J)(16:17) ) THEN
                            K = K + 1
                            NLFIPS( K ) = J
                        END IF
                    END DO
                    NFIPS = K

                END IF

            END IF

            IF( NFIPS == 0 ) CYCLE

C.............  Find SCC in inventory list
            K = 0 

            IF( CSCC == ' ' .OR. CSCC == BLKSCC ) THEN 
                NSCCS = NINVSCC
                DO J = 1, NINVSCC
                    NLSCCS(J) = J
                END DO
            ELSE
                SCCIDX = FINDC( CSCC, NINVSCC, INVSCC )
                IF( SCCIDX < 1 ) THEN
                    MESG = 'WARNING: SCC '//CSCC//' is not in the inventory.'
                    CALL M3MESG( MESG )
                    CYCLE
                END IF
                NSCCS =1
                NLSCCS(1) = SCCIDX
            END IF

C.............  Check pollutant name and mode
            K = 0
            IF( POLNAM == ' ' ) THEN
                NPOLS  = NIPPA+NMSPC
                DO J = 1, NPOLS
                    NLPOLS( J ) =  J
                END DO
            ELSE 
                NPOLS = 0
                IF ( EXPCFFLAG ) THEN
                    POLIDX = INDEX1( POLNAM, NIPPA, EANAM )
                    SPCIDX = INDEX1( POLNAM, NMSPC, EMNAM )
                    NPOLS = 1
                    IF( POLIDX > 0 ) THEN
                        NLPOLS( NPOLS ) = POLIDX
                    ELSE IF( SPCIDX > 0 ) THEN
                        NLPOLS( NPOLS ) = NIPPA + SPCIDX
                    ELSE
                        MESG = 'WARNING: '//TRIM(POLNAM)//' is not found in' //
     &                         ' both inventory pollutant and model species lists.'
                        CALL M3MESG( MESG )
                        CYCLE
                    END IF

                    IF( SPCIDX > 0 .AND. POLIDX > 0 ) THEN
                        NPOLS = 2
                        NLPOLS( 1 ) = POLIDX
                        NLPOLS( 2 ) = NIPPA + SPCIDX
                    END IF

                ELSE        ! inv poll-specific control factors application
                    POLIDX = INDEX1( POLNAM, NIPPA, EANAM )
                    IF ( POLIDX > 0) THEN
                        NPOLS = NPOLS + 1 
                        NLPOLS( NPOLS ) = POLIDX
                        DO J = 1, NSMATV
                           L1 = INDEX( TSVDESC( J ), SPJOIN )
                           L2 = LEN_TRIM( TSVDESC( J ) )
                           IF( POLNAM == TSVDESC( J )( 1:L1-1 ) ) THEN
                               SPCNAM = TSVDESC( J )( L1+1:L2 )
                               SPCIDX = INDEX1( SPCNAM, NMSPC, EMNAM )
                               NPOLS = NPOLS + 1
                               NLPOLS( NPOLS ) = NIPPA + SPCIDX
                           END IF
                        END DO
                    ELSE
                        MESG = 'WARNING: '//TRIM(POLNAM)//' is not found in' //
     &                         ' inventory pollutant list.'
                        CALL M3MESG( MESG )
                        CYCLE
                    END IF
                END IF
            END IF

C.............  Check month values
            NMONS = 0
            IF( SEGMENT( 4 ) == ' ' ) SEGMENT( 4 ) = '0'
            IF( STR2INT( SEGMENT( 4 ) ) == 0 ) THEN
                NMONS = 12
                DO J = 1, NMONS
                    NLMONS( J ) = J
                END DO
            ELSE IF ( .NOT. CHKINT( SEGMENT( 4 ) ) ) THEN
                EFLAG = .TRUE.
                MESG = 'ERROR: Bad month format '//TRIM(SEGMENT(4))
                CALL M3MESG( MESG )
                CYCLE
            ELSE         ! month is integer value
                MON = STR2INT( SEGMENT( 4 ) )
                IF( MON < 1  .OR. MON > 12 ) THEN
                    EFLAG = .TRUE.
                    MESG = 'ERROR: Bad month format '//TRIM(SEGMENT(4))
                    CALL M3MESG( MESG )
                    CYCLE
                END IF
                NMONS = 1
                NLMONS( 1 ) = MON
            END IF

C.............  check no of fips/scc/poll/mon
            IF( NFIPS < 1 .OR. NSCCS < 1 .OR. NPOLS < 1 .OR. NMONS < 1 )  CYCLE

C.............  check and get up control factor values
            IF ( .NOT. CHKREAL( SEGMENT( 5 ) ) ) THEN
                EFLAG = .TRUE.
                MESG = 'ERROR: Bad control factor value'// TRIM(SEGMENT(5))
                CALL M3MESG( MESG )
                CYCLE
            END IF
 
            CFVAL = STR2REAL( SEGMENT( 5 ) )

            DO L = 1, NFIPS
              DO J = 1, NSCCS
                DO K = 1, NPOLS
                  DO M = 1, NMONS
                    CFPRO( NLFIPS(L),NLSCCS(J),NLPOLS(K),NLMONS(M) ) = CFVAL
                  END DO
                END DO
              END DO
            END DO

          END DO  

        END DO  

        REWIND( CFDEV )
        
        IF( EFLAG ) THEN
            MESG = 'Problem found in control factor file.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF
        
        RETURN

999     MESG = 'End of file'
        MESG = 'End of file reached unexpectedly. ' //
     &         'Check format of CFPRO' // CRLF() // BLANK5 //
     &         'input file.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )  
      
C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )
94050   FORMAT(  A,I8,A,F5.2,A,F5.2,A,I6.6,5A,I2 )
        
        END SUBROUTINE RDCFPRO

