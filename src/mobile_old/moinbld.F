
C Version "@(#)$Id$ $Source$ $Date$ 

        SUBROUTINE MOINBLD( MODEV, IPARM, YEAR, TMPRTR )

C***********************************************************************
C  program body starts at line
C
C  DESCRIPTION:
C      This program accepts a parameter scheme index IPARM  
C      and an ambient temperature, and creates an input file for 
C      MOBILE5a.  The parameter x-ref arrays are provided through a common 
C      block.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Prototype 12/95 by MRH
C
C***********************************************************************

        IMPLICIT NONE

C...........   INCLUDES:
c       INCLUDE 'EFS.EXT'
        INCLUDE 'CHDIMS3.EXT'
        INCLUDE 'MBDIMS3.EXT'

C...........   EXTERNAL FUNCTIONS and their descriptions:
        INTEGER      FIND1
        LOGICAL      GETYN
        INTEGER      TRIMLEN
        EXTERNAL     GETYN, TRIMLEN

C...........   LOCAL VARIABLES and their descriptions:

C...........   Subroutine arguments
        INTEGER MODEV                    ! Unit number for MOBILE5a file
        INTEGER IPARM                    ! Scenario specific parameter index no.
        INTEGER YEAR                     ! Year for calculating EFs
        REAL    TMPRTR                   ! Ambient temperature

C...........   Normal MOBILE5a parameter x-ref variables
        INTEGER IPSI  ( NPREF )          ! Normal x-ref scheme index
        INTEGER NIDX  ( NPREF, NNFAC )   ! Normal x-ref entries
        INTEGER CIDX  ( NPREF, NCFAC )   ! Control Normal x-ref entries
        INTEGER MRNUM                    ! Actual no. of normal x-ref entries

        COMMON / MPARMS / MRNUM, IPSI, NIDX, CIDX 

C...........   Main data table indices
        INTEGER IACCUM( MXFACS )
        INTEGER IATP  ( MXFACS )
        INTEGER IAVSPD( MXFACS )
        INTEGER ICTOFF( MXFACS )
        INTEGER IDIESL( MXFACS )
        INTEGER IEFIM ( MXFACS )
        INTEGER ILAPMO( MXFACS )
        INTEGER IOMODE( MXFACS )
        INTEGER IOXYGF( MXFACS )         
        INTEGER IPREST( MXFACS )
        INTEGER IPURGT( MXFACS )
        INTEGER IRDIST( MXFACS )
        INTEGER ITAMRT( MXFACS )
        INTEGER ITPLEN( MXFACS )

        SAVE    IACCUM, IATP, IAVSPD, ICTOFF, IDIESL, IEFIM , ILAPMO, 
     &          IOMODE, IOXYGF, IPREST, IPURGT, IRDIST, ITAMRT, ITPLEN

c       INTEGER IBSEXH( )
c       INTEGER IEFLEV( )
c       INTEGER IEFRFL( )


C...........   Main data tables ! number of rows second to help with indexing

        REAL    TACCUM( 200, MXFACS )  ! Mileage accum rates: 25 yrs, 8 vtype
        REAL    TAVSPD( 8  , MXFACS )  ! Average speeds for all vehicle types
        REAL    TDIESL( 51 , MXFACS )  ! Dsl sales % for 25 years LDVs & LDTs
                                       !   Addtional space for year.
        REAL    TOMODE( 3  , MXFACS )  ! 3 operating mode parameters
        REAL    TOXYGF( 5 , MXFACS )   ! Oxygenated fuels params, 5 is integer
        REAL    TRDIST( 200, MXFACS )  ! Registration dist'n: 25 yrs, 8 vtype
        REAL    TTPLEN( 7  , MXFACS )  ! Speed, Trip lengths for 6 trip types

        SAVE    TACCUM, TAVSPD, TDIESL, TOMODE, TRDIST, TTPLEN

C. . . . . . . . LAP/Month parameters:
        INTEGER     TLAPMO( 3  , MXFACS )  ! month, RFG code, OXY code
        REAL        TREIDS( 3  , MXFACS )  ! Reid's v.p. for 3 periods
        CHARACTER*1 TFUELV( MXFACS )

        SAVE        TLAPMO, TREIDS, TFUELV

C. . . . . . . . Tampering rates: 3 time categroies, 4 vehicle types, 8 kinds
        REAL    TTRZML( 96 , MXFACS )  ! Zero-mile level
        REAL    TTRDR ( 96 , MXFACS )  ! Deterioration rate

        SAVE    TTRZML, TTRDR

C. . . . . . . . I/M program tables:
        INTEGER TEFIM ( 16, MXFACS )   !  Table of I/M integers
        INTEGER IMCX1, IMCX2, ISTAT, TTYPE
        REAL    TCTOFF( 3 , MXFACS )   !  Table of I/M transient cutoffs
        REAL    CUTOFF
        LOGICAL ONEFILE
        LOGICAL READIMC1( 0:MAXT12 ), READIMC2( 0:MAXT4P )


        INTEGER T12P1, T4PP1
        PARAMETER( T12P1 = MAXT12 + 1 )
        PARAMETER( T4PP1 = MAXT4P + 1 )

        DATA    READIMC1 / T12P1 * .TRUE. /
        DATA    READIMC2 / T4PP1 * .TRUE. /

        CHARACTER*256 FNAME
        CHARACTER*100 IMENV1, IMENV2
        CHARACTER*16  TNAME1, TNAME2
        CHARACTER*4   CCUT( 3 )

        SAVE    TEFIM, ONEFILE, TCTOFF, READIMC1, READIMC2

C. . . . . . . . Anti-Tampering program tables:
        INTEGER TATP  ( 18, MXFACS )   ! ATP years, vtypes, insp types, etc.
        INTEGER TPREST( 10, MXFACS )   ! Pressure test years, vtypes, etc.
        INTEGER TPURGT( 10, MXFACS )   ! Purge test years, vtypes, etc.

        SAVE    TATP, TPREST, TPURGT

C...........   Numbers of actual entries in tables 
        INTEGER NACCUM, NRDIST, NTAMRT, NLAPMO, NTPLEN, NAVSPD,
     &          NOMODE, NEFIM, NDIESL, NCTOFF, NOXYGF, NATP, 
     &          NPREST, NPURGT 

        SAVE    NACCUM, NRDIST, NTAMRT, NLAPMO, NTPLEN, NAVSPD,
     &          NOMODE, NEFIM, NDIESL, NCTOFF, NOXYGF, NATP,
     &          NPREST, NPURGT

        INTEGER IOXY, JOPM, JTPLEN, JLAPMO, JOXY, IRFG
 
        LOGICAL READFAC( 19 ), FIRSTIME( 19 ), LSPD
        
        DATA READFAC  /  19*.FALSE. /
        DATA FIRSTIME /  19*.TRUE.  /

        SAVE READFAC, FIRSTIME

C...........   Output variables
        
        INTEGER PROMPT, IOUNEW, TAMFLG, SPDFLG,
     &          VMFLAG, MYMRFG, NEWFLG, IMFLAG,
     &          ALHFLG, ATPFLG, RLFLAG, LOCFLG,
     &          TEMFLG, OUTFMT, PRTFLG, IDLFLG,
     &          NMHFLG, HCFLAG

        INTEGER DSFLAG, OXYFLG, RFGFLG, REGION
c       INTEGER INPSI_LAST, ICPSI_LAST

        INTEGER NSPEED, MONTH
        REAL SPEED( 8 )

        CHARACTER*60 PROJID

        DATA PROJID / 'Default name' /

C...........   Data statements

c       DATA    INPSI_LAST, ICPSI_LAST / -1, -1 /

C...........   Other local variables

        INTEGER G6, G7, G8
        INTEGER I, ICUT, IP, IPPREV, J, JDSL, JJ, JK, K
        CHARACTER*256 MESG, BUFFER

C***********************************************************************
C   begin body of program MOINBLD

C.........  Set defaults to control variables
        PROMPT = 3
        IOUNEW = 0
        TAMFLG = 1
        SPDFLG = 2
        VMFLAG = 1
        MYMRFG = 1
        NEWFLG = 1
        IMFLAG = 1
        ALHFLG = 1
        ATPFLG = 1
        RLFLAG = 1
        LOCFLG = 1
        TEMFLG = 2
        OUTFMT = 3
        PRTFLG = 4
        IDLFLG = 1
        NMHFLG = 3
        HCFLAG = 3


C.........  Find normal parameters from x-ref index
        IP = FIND1( IPARM, MRNUM, IPSI )

        IF( IP .LT. 0 ) THEN
            WRITE( MESG,94100 ) 'Parameter profile', IPARM, 
     &                          'not found in MPREF'
            CALL M3EXIT( 'MOINBLD', 0, 0, MESG, 2 )
        ENDIF

C.........  Make sure parameter set in consistent

C Notes:
C  1) G4 and G5 can't be positive if G3 is not

C.........  Determine which factors must be read in given parameters
        J = 0
        DO 22 I = 1, NNFAC
    
            J = J + 1
            IF( NIDX( IP,I ) .NE. 0 ) READFAC( J ) = .TRUE.

   22   CONTINUE

        DO 33 I = 1, NCFAC

            J = J + 1
            IF( CIDX( IP,I ) .NE. 0 ) READFAC( J ) = .TRUE.

   33   CONTINUE

C.........  Read in main factor matrices and indices, unless already have
        IF ( READFAC( 1 ) .AND. FIRSTIME( 1 ) ) THEN
            CALL GETFACR( 'MILEAGE ACCUMULATION RATES', 'MACUM',
     &                    MXFACS, NACCUM, IACCUM, 200, TACCUM )
            FIRSTIME( 1 ) = .FALSE.

        ENDIF

        IF ( READFAC( 2 ) .AND. FIRSTIME( 2 ) ) THEN
            CALL GETFACR( 'REGISTRATION DISTRIBUTIONS', 'MRGDS',
     &                    MXFACS, NRDIST, IRDIST, 200, TRDIST )
            FIRSTIME( 2 ) = .FALSE.

        ENDIF

        IF ( READFAC( 3 ) .AND. FIRSTIME( 3 ) ) THEN
            CALL GETFACRR( 'TAMPERING RATES', 'MTRATE', 
     &                    MXFACS, NTAMRT, ITAMRT, 96, TTRZML,
     &                    96, TTRDR )
            FIRSTIME( 3 ) = .FALSE.

        ENDIF

        IF ( READFAC( 4 ) .AND. FIRSTIME( 4 ) ) THEN
            CALL GETFACICR( 'LAP/MONTH', 'MLAPMO', 
     &                      MXFACS, NLAPMO, ILAPMO, 3, TLAPMO, 
     &                      1, TFUELV, 3, TREIDS )
            FIRSTIME( 4 ) = .FALSE.
        ENDIF

        IF ( READFAC( 5 ) .AND. FIRSTIME( 5 ) ) THEN
            CALL GETFACR( 'TRIP LENGTH DISTRIBUTIONS', 'MTLEN',
     &                    MXFACS, NTPLEN, ITPLEN, 7, TTPLEN )
            FIRSTIME( 5 ) = .FALSE.

        ENDIF

        IF ( READFAC( 6 ) .AND. FIRSTIME( 6 ) ) THEN
           
            CALL GETFACR( 'DIESEL FUELS SALES FRACTIONS', 'MDLSF',
     &                     MXFACS, NDIESL, IDIESL, 51, TDIESL )
            FIRSTIME( 6 ) = .FALSE.

        ENDIF

        IF ( READFAC( 7 ) .AND. FIRSTIME( 7 ) ) THEN
            CALL GETFACR( 'AVERAGE SPEEDS','MAVSPD',
     &                    MXFACS, NAVSPD, IAVSPD, 8, TAVSPD )
            FIRSTIME( 7 ) = .FALSE.
    
        ENDIF

        IF ( READFAC( 8 ) .AND. FIRSTIME( 8 ) ) THEN
            CALL GETFACR( 'OPERATING MODES','MOPMOD',
     &                    MXFACS, NOMODE, IOMODE, 3, TOMODE )
            FIRSTIME( 8 ) = .FALSE.

        ENDIF

c           IF ( READFAC( 9 ) ) THEN

c           ENDIF

c           IF ( READFAC( 11 ) .AND. FIRSTIME( 10 ) ) THEN
c               CALL GETFAC ( 'BASIC EMISSIONS RATES', 'MBSEXH',
c    &                        MXFACS, NBSEXH, IBSEXH,  ... )
c               FIRSTIME( 10 ) = .FALSE.
c           ENDIF

c           IF ( READFAC( 12 ) .AND. FIRSTIME( 11 ) ) THEN
c               CALL GETFACIR( 'NEW EVAPORATIVE EMISSIONS TESTS',
c    &                         'MNEVAP', MXFACS, NNEVAP, INEVAP,
c    &                         #, , #, )
c               FIRSTIME( 11 ) = .FALSE.
c           ENDIF

        IF ( READFAC( 13 ) .AND. FIRSTIME( 13 ) ) THEN
            CALL GETFACI ( 'I/M PROGRAMS', 'MEFIM', MXFACS, 
     &                     NEFIM, IEFIM, 16, TEFIM )
            FIRSTIME( 13 ) = .FALSE.

            MESG = 'Use default I/M Credit files throughout?'
            ONEFILE = GETYN ( MESG, .TRUE. )

        ENDIF

        IF ( READFAC( 13 ) ) THEN
C.............  Set default I/M credits file settings
            IMCX1  = 0   ! Index to default Tech I-II read indicator
            IMCX2  = 0   ! Index to default Tech IV+ read indicator
            TNAME1 = 'DEFAULT'
            TNAME2 = 'DEFAULT'

C.............  Set I/M credit file settings individually
            IF( .NOT. ONEFILE ) THEN  

                IMCX1 = CIDX( IP,4 )
                IMCX2 = CIDX( IP,5 )
                IF( CIDX( IP,4 ) .GT. 0 ) THEN 
                    TNAME1 = 'NON-DEFAULT'
                ENDIF

                IF( CIDX( IP,5 ) .GT. 0 ) THEN
                    TNAME2 = 'NON-DEFAULT'
                ENDIF

            ENDIF

C.............  Read Tech I-II file
            IF( READIMC1( IMCX1 ) ) THEN

   44           WRITE( MESG, 94100 ) 
     &                'Enter logical name of '// 
     &                TNAME1( 1:TRIMLEN( TNAME1 ) ) // ' TECH I-II '//
     &                'FILE', IMCX1

                IF( IMCX1 .GE. 0 .AND. IMCX1 .LE. 9 ) 
     &              WRITE( FNAME, '(A,I1)') 'TECH12_',IMCX1
                IF( IMCX1 .GE. 10 .AND. IMCX1 .LE. 99 ) 
     &              WRITE( FNAME, '(A,I2)') 'TECH12_',IMCX1
                IF( IMCX1 .GE. 99 .AND. IMCX1 .LE. 999 ) 
     &              WRITE( FNAME, '(A,I3)') 'TECH12_',IMCX1

                CALL GETSTR( MESG, FNAME, FNAME )
                CALL ENVSTR( FNAME, 'Non-default Tech I-II file name',
     &                       'TECH12.D', IMENV1, ISTAT )

                IF( ISTAT .LT. 0 ) THEN

                    WRITE( *,'(/5X,A,///10X,A)' )
     &                    'ERROR: input file not found: '//
     &                    FNAME( 1:TRIMLEN( FNAME ) ),
     &                    'Could not open input file "' //
     &                    FNAME( 1:TRIMLEN( FNAME ) ) // '"'

                    IF( GETYN( '     Try again', .TRUE. ) ) THEN
                        GO TO 44
                    ELSE
                        MESG = 'Ending program MOINBLD.'
                        CALL M3EXIT( 'MOINBLD', 0, 0, MESG, 2 )
                    ENDIF

                ENDIF

                READIMC1( IMCX1 ) = .FALSE.

            ENDIF

C.............  Read Tech IV+ file
            IF( READIMC2( IMCX2) ) THEN
   55           WRITE( MESG, 94100 ) 
     &                'Enter logical name of '// 
     &                TNAME2( 1:TRIMLEN( TNAME2 ) ) // ' TECH IV+ '//
     &                'FILE', IMCX2

                IF( IMCX2 .GE. 0 .AND. IMCX2 .LE. 9 ) 
     &              WRITE( FNAME, '(A,I1)') 'TECH4P_',IMCX2
                IF( IMCX2 .GE. 10 .AND. IMCX2 .LE. 99 ) 
     &              WRITE( FNAME, '(A,I2)') 'TECH4P_',IMCX2
                IF( IMCX2 .GE. 99 .AND. IMCX2 .LE. 999 ) 
     &              WRITE( FNAME, '(A,I3)') 'TECH4P_',IMCX2

                CALL GETSTR( MESG, FNAME, FNAME )
                CALL ENVSTR( FNAME, 'Non-default Tech IV+ file name',
     &                       'IMDATA.D', IMENV2, ISTAT )

                IF( ISTAT .LT. 0 ) THEN

                    WRITE( *,'(/5X,A,///10X,A)' )
     &                    'ERROR: input file not found: '//
     &                    FNAME( 1:TRIMLEN( FNAME ) ),
     &                    'Could not open input file "' //
     &                    FNAME( 1:TRIMLEN( FNAME ) ) // '"'

                    IF( GETYN( '     Try again', .TRUE. ) ) THEN
                        GO TO 55
                    ELSE
                        MESG = 'Ending program MOINBLD'
                        CALL M3EXIT( 'MOINBLD', 0, 0, MESG, 2 )
                    ENDIF

                ENDIF

                READIMC2( IMCX2 ) = .FALSE.
            ENDIF

        ENDIF

C.........  Read anti-tampering programs
        IF ( READFAC( 16 ) .AND. FIRSTIME( 16 ) ) THEN
            CALL GETFACI( 'ANTI-TAMPERING PROGRAMS','MEFATP',
     &                     MXFACS, NATP, IATP, 18, TATP )
            FIRSTIME( 16 ) = .FALSE.

        ENDIF

C.........  Read functional pressure test programs
        IF ( READFAC( 17 ) .AND. FIRSTIME( 17 ) ) THEN
            CALL GETFACI( 'FUNCTIONAL PRESSURE TEST','MPREST', MXFACS,
     &                     NPREST, IPREST, 10, TPREST )
            FIRSTIME( 17 ) = .FALSE.

        ENDIF

C.........  Read functional purge test programs
        IF ( READFAC( 18 ) .AND. FIRSTIME( 18 ) ) THEN
            CALL GETFACI( 'FUNCTIONAL PURGE TEST','MPURGT', MXFACS,
     &                     NPURGT, IPURGT, 10, TPURGT )
            FIRSTIME( 18 ) = .FALSE.

        ENDIF

C.............  Set control flags

        IF( NIDX( IP,3 ) .GT. 0 ) THEN
            TAMFLG = 2
        ENDIF

        IF(     NIDX( IP,1 ) .GT. 0 .AND. NIDX( IP,2 ) .GT. 0 ) THEN
            MYMRFG = 4
        ELSEIF( NIDX( IP,1 ) .EQ. 0 .AND. NIDX( IP,2 ) .GT. 0 ) THEN
            MYMRFG = 3
        ELSEIF( NIDX( IP,1 ) .GT. 0 .AND. NIDX( IP,2 ) .EQ. 0 ) THEN
            MYMRFG = 2
        ENDIF

        IF(     CIDX( IP,1 ) .GT. 0 .AND. CIDX( IP,2 ) .GT. 0 ) THEN
            NEWFLG = 4
        ELSEIF( CIDX( IP,1 ) .EQ. 0 .AND. CIDX( IP,2 ) .GT. 0 ) THEN
            NEWFLG = 3
        ELSEIF( CIDX( IP,1 ) .GT. 0 .AND. CIDX( IP,2 ) .EQ. 0 ) THEN
            NEWFLG = 2
        ENDIF

        IF( CIDX( IP,3 ) .GT. 0 ) THEN
            IMFLAG = 2
        ENDIF

        IF( NIDX( IP,10 ) .GT. 0 ) THEN
            ALHFLG = 3
        ENDIF

        G6 = CIDX( IP,6 )
        G7 = CIDX( IP,7 )
        G8 = CIDX( IP,8 )
        IF( G6 .GT. 0 .AND. G7 .GT. 0 .AND. G8 .GT. 0 ) THEN
            ATPFLG = 8
        ELSEIF( G6 .GT. 0 .AND. G7 .GT. 0 ) THEN
            ATPFLG = 5
        ELSEIF( G6 .GT. 0 .AND. G8 .GT. 0 ) THEN
            ATPFLG = 6
        ELSEIF( G7 .GT. 0 .AND. G8 .GT. 0 ) THEN
            ATPFLG = 7
        ELSEIF( G6 .GT. 0 ) THEN
            ATPFLG = 2
        ELSEIF( G7 .GT. 0 ) THEN
            ATPFLG = 3
        ELSEIF( G8 .GT. 0 ) THEN
            ATPFLG = 4
        ENDIF

C........  Deal with AVE SPEEDS now, because if they are all the same, need to
C               modify SPDFLG
        IF( NIDX( IP,7 ) .GT. 0 .AND. 
     &      NIDX( IP,5 ) .EQ. 0       ) THEN

            NSPEED = 8

            J = FIND1( NIDX( IP,7 ), NAVSPD, IAVSPD )
            IF( J .LT. 0 ) THEN
                WRITE( MESG,94100 ) 
     &                 'AVERAGE SPEED profile', NIDX( IP,7 ),
     &                 'from parameter scheme', IP,
     &                 'not found in average speed table.'
                CALL M3EXIT( 'MOINBLD', 0, 0, MESG, 2 )
            ENDIF

            DO 111 I = 1, 8
                SPEED( I ) = TAVSPD( I, J )
  111       CONTINUE

C.............  Make adjustment for output if speeds all the same
            LSPD = .TRUE.                          ! Speed the same
            DO 122 I = 2, 8
                IF( SPEED( I ) .NE. SPEED( 1 ) ) LSPD = .FALSE.
  122       CONTINUE

            IF( LSPD ) THEN
                NSPEED = 1                         ! Only output one speed
                SPDFLG = 1
            ENDIF

        ELSEIF( NIDX( IP,5 ) .GT. 0       ) THEN

            JTPLEN = FIND1( NIDX( IP,5 ), NTPLEN, ITPLEN )
            IF( JTPLEN .LT. 0 ) THEN
                WRITE( MESG,94100 ) 
     &                 'TRIP LENGTH profile', NIDX( IP,5 ),
     &                 'from parameter scheme', IP,
     &                 'not found in trip length table.'
                CALL M3EXIT( 'MOINBLD', 0, 0, MESG, 2 )
            ENDIF

            IF( NIDX( IP,7 ) .GT. 0 ) THEN
                WRITE( MESG,94100 ) 'AVERAGE SPEED profile', NIDX(IP,7),
     &                              'from parameter scheme', IP,
     &                              'not applied because trip length'//
     &                              'distribution', NIDX( IP,5 ), 
     &                              'applied instead.'
                CALL M3MESG( MESG )
            ENDIF

            NSPEED = 1
            SPEED( 1 ) = TTPLEN( 1, JTPLEN )
            SPDFLG = 3

        ENDIF

C.........  Write out control section  ..........................
        WRITE( MODEV, 93005 ) PROMPT, IOUNEW
        WRITE( MODEV, 93000 ) PROJID
        WRITE( MODEV, 93010 ) TAMFLG, SPDFLG,
     &        VMFLAG, MYMRFG, NEWFLG, IMFLAG,
     &        ALHFLG, ATPFLG, RLFLAG, LOCFLG,
     &        TEMFLG, OUTFMT, PRTFLG, IDLFLG,
     &        NMHFLG, HCFLAG

C.........  Write out One-time data section  ....................

C.........  Tampering rates for {Non-I/M case} OR {Non-I/M and I/M cases}
        IF( NIDX( IP,3 ) .GT. 0 ) THEN

            J = FIND1( NIDX( IP,3 ), NTAMRT, ITAMRT )
            IF( J .LT. 0 ) THEN
                WRITE( MESG,94100 )
     &                 'TAMPERING RATES profile', NIDX( IP,3 ),
     &                 'from parameter scheme', IP,
     &                 'not found in tampering rates table.'
                CALL M3EXIT( 'MOINBLD', 0, 0, MESG, 2 )
            ENDIF

            WRITE( MODEV, 93020 ) ( TTRZML( I,J ), I=1,96 ) ! Non-I/M ZML

            IF( CIDX( IP,3 ) .GT. 0 ) THEN                  ! I/M ZML 

                JJ = FIND1( CIDX( IP,3 ), NEFIM, IEFIM  )
                IF( JJ .LT. 0 ) THEN
                    WRITE( MESG,94100 ) 
     &                     'I/M PROGRAM profile', CIDX( IP,3 ),
     &                     'from parameter scheme', IP,
     &                     'not found in I/M PROGRAMS table.'
                    CALL M3EXIT( 'MOINBLD', 0, 0, MESG, 2 )
                ENDIF

                IF( TEFIM( 16,JJ ) .LE. 0 ) THEN

                     WRITE( MESG,94100 )
     &                'I/M profile',CIDX( IP,3 ),
     &                'from parameter scheme', IP,
     &                'gives illegal tampering profile',TEFIM( 16,JJ ),
     &                'to use for non-default I/M tampering rates.'
                     CALL M3EXIT( 'MOINBLD', 0, 0, MESG, 2 )
 
                ELSE
                    JK = FIND1( TEFIM( 16,JJ ), NTAMRT, ITAMRT )
                    IF( JK .LT. 0 ) THEN
                        WRITE( MESG,94100 )
     &                        'TAMPERING RATES profile', TEFIM( 16,JJ ), 
     &                        '(requested by I/M profile', CIDX( IP,3 ),
     &                        'from parameter scheme', IP,
     &                        ') not in tampering rates table.'
                        CALL M3EXIT( 'MOINBLD', 0, 0, MESG, 2 )
                    ENDIF

                ENDIF

                WRITE( MODEV, 93020 ) ( TTRZML( I,JK ), I=1,96 )

            ENDIF

            WRITE( MODEV, 93020 ) ( TTRDR ( I,J ), I=1,96 ) ! Non-I/M DR

            IF( CIDX( IP,3 ) .GT. 0 ) THEN                  ! I/M DR
                WRITE( MODEV, 93020 ) ( TTRDR ( I,JK ), I=1,96 )
            ENDIF

        ENDIF

C.........  Annual mileage accumulation rates
        IF( NIDX( IP,1 ) .GT. 0 ) THEN

            J = FIND1( NIDX( IP,1 ), NACCUM, IACCUM )
            IF( J .LT. 0 ) THEN
                WRITE( MESG,94100 ) 
     &                 'MILEAGE ACCUM profile', NIDX( IP,1 ),
     &                 'from parameter scheme', IP,
     &                 'not found in mileage accum table.'
                CALL M3EXIT( 'MOINBLD', 0, 0, MESG, 2 )
            ENDIF

            WRITE( MODEV, 93030 ) ( TACCUM( I,J ), I=1,200 )

        ENDIF

C.........  Registration distributions
        IF( NIDX( IP,2 ) .GT. 0 ) THEN

            J = FIND1( NIDX( IP,2 ), NRDIST, IRDIST )
            IF( J .LT. 0 ) THEN
                WRITE( MESG,94100 ) 
     &                 'REGISTRATION DISTN profile', NIDX( IP,2 ),
     &                 'from parameter scheme', IP,
     &                 'not found in registration distn table.'
                CALL M3EXIT( 'MOINBLD', 0, 0, MESG, 2 )
            ENDIF

            WRITE( MODEV, 93040 ) ( TRDIST( I,J ), I=1,200 )

        ENDIF

C.........  Basic emissions rates
c           IF( CIDX( IP,1 ) .GT. 0 ) THEN
C! Hard - integer and real
c           ENDIF

C.........  New Evaporative Emissions test procedure
c           IF( CIDX( IP,2 ) .GT. 0 ) THEN
C! Hard - integer and real
c           ENDIF

C.........  Inspection and Maintenence Program
        IF( CIDX( IP,3 ) .GT. 0 ) THEN

C.............  Get index of I/M program
            J = FIND1( CIDX ( IP,3 ), NEFIM , IEFIM )
            IF( J .LT. 0 ) THEN
                WRITE( MESG,94100 ) 
     &                 'I/M PROGRAM profile', CIDX( IP,3 ),
     &                 'from parameter scheme', IP,
     &                 'not found in I/M PROGRAMS table.'
                CALL M3EXIT( 'MOINBLD', 0, 0, MESG, 2 )
            ENDIF


C.............  Set defaults
            ICUT = 1    ! No non-default cutoffs

            TTYPE = TEFIM( 14,J )
C.............  Set emissions test parms and ensure internal consistency
            IF( TTYPE .EQ. 1 .OR. TTYPE .EQ. 2 .OR.
     &          TTYPE .EQ. 3                        ) THEN

C.................  Check for non-default cut-offs for types 1,2,3
                IF( TEFIM( 15,J ) .GT. 0 ) THEN
                    IF ( IMCX1 .EQ. 0 .OR. 
     &                   IMCX2 .EQ. 0      ) THEN 

                         WRITE( MESG,94100 )
     &                     'MOBILE5 parameters profile',IP,
     &                     'has I/M profile with test type',TTYPE,
     &                     'and non-default cutoffs, '//
     &                     'but no non-default credits files.'
                         CALL M3EXIT( 'MOINBLD', 0, 0, MESG, 2 )

                    ELSEIF( TTYPE .EQ. 1 ) THEN
                        MESG = 'IDLE TEST NON-DEFAULT CUTOFFS'

                    ELSEIF( TTYPE .EQ. 2 ) THEN
                        MESG = '2500/IDLE TEST NON-DEFAULT CUTOFFS'

                    ELSEIF( TTYPE .EQ. 3 ) THEN
                        MESG = 'LOADED TEST NON-DEFAULT CUTOFFS'

                    ENDIF

                ENDIF

            ELSEIF( TTYPE .EQ. 4 ) THEN

C..................  Ensure cut-offs will be read for type 4
                 IF( TEFIM( 15,J ) .LE. 0 ) THEN

                     WRITE( MESG,94100 )
     &                     'I/M Profile number',CIDX( IP,3 ),
     &                     'from parameter scheme', IP,
     &                     'has non-postive cutpoint index for'//
     &                     'transient test type.'
                     CALL M3EXIT( 'MOINBLD', 0, 0, MESG, 2 )

                 ELSE
                     MESG = 'TRANSIENT TEST CUTOFFS'
                 ENDIF

            ELSE

                 WRITE( MESG,94100 ) 
     &                 'I/M Profile number',CIDX( IP,3 ),
     &                 'from parameter scheme', IP,
     &                 'has illegal test type specifier of', TTYPE
                 CALL M3EXIT( 'MOINBLD', 0, 0, MESG, 2 )

            ENDIF

C.............  Read cutoff file for given case, if specified
            IF( TEFIM( 15,J ) .GT. 0 ) THEN

                IF( IP .NE. IPPREV ) 
     &          CALL GETFACR( MESG, 'MCTOFF', MXFACS, 
     &                        NCTOFF, ICTOFF, 3, TCTOFF )

                ICUT = 2
                JJ = FIND1( TEFIM( 15,J ), NCTOFF, ICTOFF )
                IF( JJ .LT. 0 ) THEN
                    BUFFER = MESG
                    WRITE( MESG,94100 ) MESG( 1:TRIMLEN(MESG) )
     &                                  //' profile', TEFIM( 15,J ),
     &                                  'not found in cutoffs table.'
                    CALL M3EXIT( 'MOINBLD', 0, 0, MESG, 2 )
                ENDIF


            ENDIF

C.............  Write I/M record
            IF( ICUT .EQ. 1 ) THEN
                WRITE( MODEV, 93100 ) ( TEFIM( I,J ), I=1,14 ), ICUT,
     &                                  2, 2
            ELSE

                DO 133 I = 1,3

                    CUTOFF = TCTOFF( I,JJ )
                    IF( CUTOFF .GE. 0. .AND. CUTOFF .LT. 10. ) THEN

                        WRITE( CCUT( I ),93101 ) CUTOFF

                    ELSEIF( CUTOFF .GE. 10. .AND. 
     &                      CUTOFF .LT. 100.       ) THEN

                        WRITE( CCUT( I ),93102 ) CUTOFF

                    ELSEIF( CUTOFF .GT. 100. .AND. 
     &                      CUTOFF .LT. 1000.      ) THEN

                        WRITE( CCUT( I ),93103 ) CUTOFF

                    ELSE

                        WRITE( MESG,94100 ) 'Cutoff value',I,
     &                        'from cutoff profile', JJ,
     &                        'is an illegal value.'
                        CALL M3EXIT( 'MOINBLD', 0, 0, MESG, 2 )

                    ENDIF

  133           CONTINUE

                WRITE( MODEV, 93100 ) ( TEFIM( I,J ), I=1,14 ), ICUT,
     &                                  2, 2, ( CCUT( I ), I=1,3 )

            ENDIF

C.............  Write file names for tech12 and tech4+ files
C... Note: Normal mobile5 will only accept 40 character names!
            WRITE( MODEV,93000 ) IMENV1( 1:TRIMLEN( IMENV1 ) )  
            WRITE( MODEV,93000 ) IMENV2( 1:TRIMLEN( IMENV2 ) )

        ENDIF

C.........  Anti-tampering programs
        IF( CIDX( IP,6 ) .GT. 0 ) THEN

            J = FIND1( CIDX( IP,6 ), NATP, IATP )

            IF( J .LT. 0 ) THEN
                WRITE( MESG,94100 ) 
     &                 'Anti-tampering profile ', CIDX( IP,6 ),
     &                 'from parameter scheme', IP,
     &                 'not found in ANTI-TAMPERING table.'
                CALL M3EXIT( 'MOINBLD', 0, 0, MESG, 2 )
            ENDIF

            WRITE( MODEV,93200 ) ( TATP( I,J ), I=1,18 )

        ENDIF

C.........  Functional pressure check
        IF( CIDX( IP,7 ) .GT. 0 ) THEN

            J = FIND1( CIDX( IP,7 ), NPREST, IPREST )

            IF( J .LT. 0 ) THEN
                WRITE( MESG,94100 ) 
     &                 'Pressure test profile ', CIDX( IP,7 ),
     &                 'from parameter scheme', IP,
     &                 'not found in PRESSURE TEST table.'
                CALL M3EXIT( 'MOINBLD', 0, 0, MESG, 2 )
            ENDIF

            WRITE( MODEV,93200 ) ( TPREST( I,J ), I=1,10 )

        ENDIF

C.........  Functional purge check
        IF( CIDX( IP,8 ) .GT. 0 ) THEN

            J = FIND1( CIDX( IP,8 ), NPURGT, IPURGT )

            IF( J .LT. 0 ) THEN
                WRITE( MESG,94100 ) 
     &                 'Purge test profile ', CIDX( IP,8 ),
     &                 'from parameter scheme', IP,
     &                 'not found in PURGE TEST table.'
                CALL M3EXIT( 'MOINBLD', 0, 0, MESG, 2 )
            ENDIF

            WRITE( MODEV,93200 ) ( TPURGT( I,J ), I=1,10 )

        ENDIF

C.............  Refueling controls 
c           IF( CIDX( IP,10 ) .GT. 0 ) THEN
C! Hard - integer and real
c           ENDIF

C........... Check LAP/Month code for consistencies, get index
        IF( NIDX( IP,4 ) .LE. 0 ) THEN ! Error

            WRITE( MESG, 94100 ) 
     &            'LAP/month record number', NIDX( IP,4 ), 
     &            'from parameter scheme', IP,
     &            'not allowed.  Value MUST be positive.'
            CALL M3EXIT( 'MOINBLD', 0, 0, MESG, 2 )

        ELSE

            JLAPMO = FIND1( NIDX( IP,4 ), NLAPMO, ILAPMO )
            IF( JLAPMO .LT. 0 ) THEN
                WRITE( MESG,94100 ) 
     &                 'LAP/month profile', NIDX( IP,4 ),
     &                 'from parameter scheme', IP,
     &                 'not found in LAP/month table.'
                CALL M3EXIT( 'MOINBLD', 0, 0, MESG, 2 )
            ENDIF

            MONTH  = TLAPMO( 1, JLAPMO )
            IOXY   = TLAPMO( 2, JLAPMO )
            IRFG   = TLAPMO( 3, JLAPMO )

            IF( MONTH .LT. 1 .OR. MONTH .GT. 7 .OR.
     &          ( MONTH .GT. 1 .AND. MONTH .LT. 7 ) ) THEN ! Error

                 WRITE( MESG,94100 ) 
     &            'Month in LAPMO record', NIDX( IP,4 ),
     &            'from parameter scheme', IP,
     &            'must be 1 or 7.'
                 CALL M3EXIT( 'MOINBLD', 0, 0, MESG, 2 )

            ENDIF

            IF( IRFG .GT. 0 ) THEN

                 IF( TFUELV( JLAPMO ) .LT. 'A' .OR.
     &               TFUELV( JLAPMO ) .GT. 'E'      ) THEN ! Error

                     WRITE( MESG,94100 ) 
     &                 'Fuel volatility in LAPMO record', NIDX( IP,4 ),
     &                 'from parameter scheme', IP,
     &                 'must be A, B, C, D, or E when RFG controls '//
     &                 'are present.'

                     CALL M3EXIT( 'MOINBLD', 0, 0, MESG, 2 )

                 ELSE

                     RFGFLG = 2

                 ENDIF

            ELSEIF( IRFG .EQ. 0 ) THEN

                RFGFLG = 1

c           ELSE ! Error 

            ENDIF
        ENDIF

        IF( IOXY .GT. 0 ) THEN

            OXYFLG = 2 

        ELSEIF( IOXY .EQ. 0 ) THEN

            OXYFLG = 1

        ELSE  ! Error

        ENDIF

C.........  Write out scenario section .........................

C.........  Scenario descriptive record

C........... REGION
        IF( NIDX( IP,9 ) .EQ. 1 .AND. 
     &      CIDX( IP,9 ) .EQ. 0       ) THEN

            REGION = 1

        ELSEIF( NIDX( IP,9 ) .EQ. 2 .AND.
     &          CIDX( IP,9 ) .EQ. 0       ) THEN

            REGION = 2

        ELSEIF( NIDX( IP,9 ) .EQ. 1  .AND.
     &          CIDX( IP,9 ) .EQ. -1      ) THEN

            IF( CIDX( IP,6 ) .GT. 0 .AND.
     &          CIDX( IP,3 ) .GT. 0       ) THEN

                REGION = 3
            ELSE                             ! Error

                WRITE( *, 93000 ) 
     &                 'Must use ATP and I/M programs even when ',
     &                 'using default LEV program (LEV parm=-1)'
                WRITE( MESG,94000 ) 
     &                'Illegal combination of parameters'
                CALL M3EXIT( 'MOINBLD', 0, 0, MESG, 2 )

            ENDIF

        ELSEIF( NIDX( IP,9 ) .EQ. 1 .AND.
     &          CIDX( IP,9 ) .GT. 0       ) THEN

            IF( CIDX( IP,3 ) .GT. 0 ) THEN

                REGION = 4

            ELSE ! Error
 
                WRITE( *,93000 ) 
     &                 'Must specify I/M program for LEV program'
                CALL M3EXIT( 'MOINBLD', 0, 0, MESG, 2 )
                WRITE( MESG,94000 ) 
     &                'Illegal combination of parameters'

            ENDIF

        ELSE ! Error

            WRITE( *,93000 ) 'Altitude setting',NIDX( IP,9 )
            WRITE( *,93000 ) 'LEV program selection',CIDX( IP,9 )
            WRITE( MESG,94000 ) 
     &             'Illegal combination of parameters'
            CALL M3EXIT( 'MOINBLD', 0, 0, MESG, 2 )

        ENDIF

C........... OPERATING MODE FRACTIONS
        IF( NIDX( IP,8 ) .GT. 0 ) THEN

            JOPM = FIND1( NIDX( IP,8 ), NOMODE, IOMODE )
            IF( JOPM .LT. 0 ) THEN
                WRITE( MESG,94100 )
     &                 'OPERATING MODES profile', NIDX( IP,8 ),
     &                 'from parameter scheme', IP,
     &                 'not found in operating modes table.'
                CALL M3EXIT( 'MOINBLD', 0, 0, MESG, 2 )
            ENDIF


        ELSE ! Error

            WRITE( MESG,94000 )
     &             'Unable to set operating modes '//
     &             'for parameter scheme', IP
            CALL M3EXIT( 'MOINBLD', 0, 0, MESG, 2 )

        ENDIF

C.........  Write out scenario record: Set no diurnal/diurnal to 'N'
C              
        IF( NSPEED .EQ. 1 ) THEN

            WRITE( MODEV, 93400 ) 
     &            REGION, YEAR-1900, SPEED( 1 ), TMPRTR,
     &            ( TOMODE( K,JOPM ), K=1,3 ), MONTH, 'N'

        ELSE

            WRITE( MODEV, 93410 ) 
     &              REGION, YEAR-1900, 
     &              ( SPEED ( K ), K=1,8 ), TMPRTR,
     &              ( TOMODE( K,JOPM ), K=1,3 ), MONTH, 'N'

        ENDIF


C.........  Low Emissions vehicles (LEV) record

c           IF( CIDX( IP,9 ) .GT. 0 ) THEN

c           ENDIF 

c           WRITE( MODEV, 93XXX )

C.........  Local Area Parameter Record (LAP)

C...........  DSFLAG
        IF( NIDX( IP,6 ) .EQ. 0 ) THEN

            DSFLAG = 1

        ELSEIF( NIDX( IP,6 ) .GT. 0 ) THEN

            DSFLAG = 2

        ELSE ! Error
            WRITE( MESG,94000 ) 'Error setting DSFLAG'
            CALL M3EXIT( 'MOINBLD', 0, 0, MESG, 2 )
 
        ENDIF

C...........  Find road class name by searching numbers index
c       JRCL = FIND1( RCL, NRCLAS, MROADS3 )

C...........  Blank out fuel volatility if missing
        IF( TFUELV( JLAPMO ) .EQ. 'N' ) TFUELV( JLAPMO ) = ' '

        WRITE( MODEV, 93500 ) 
     &         'Road Class Name', 
     &         TFUELV( JLAPMO ), 
     &         TMPRTR - 1., TMPRTR + 1.0,     !  Min/Max temp unimportant
     &         ( TREIDS( K,JLAPMO ), K = 1,2 ), 
     &         INT( TREIDS( 3, JLAPMO ) ),
     &         OXYFLG, DSFLAG, RFGFLG

C.........  Oxygenated Fuels Descriptive Record: index w/ NIDX(IP,4)
        IF( IOXY .GT. 0 ) THEN

            IF ( FIRSTIME( 19 ) ) THEN
                CALL GETFACR( 'OXYGENATED FUELS PROGRAM', 'MOXYF',
     &                        MXFACS, NOXYGF, IOXYGF, 5, TOXYGF )
                FIRSTIME( 19 ) = .FALSE.

            ENDIF

            JOXY = FIND1( IOXY, NOXYGF, IOXYGF )
                
            IF( JOXY.LT. 0 ) THEN
                WRITE( MESG,94100 ) 
     &                 'Oxygenated fuels profile', IOXY, 
     &                 'from parameter scheme',IP, 'is not '//
     &                 'found in the oxygenated fuels programs table.'
                CALL M3EXIT( 'MOINBLD', 0, 0, MESG, 2 )
            ENDIF

            WRITE( MODEV, 93510 ) ( TOXYGF( I,JOXY ), I = 1,4 ),
     &                              INT( TOXYGF( 5,JOXY ) )
        ENDIF

C.........  Diesel Sales Fractions
        IF ( NIDX( IP,6 ) .GT. 0 ) THEN

            JDSL = FIND1( NIDX( IP,6 ), NDIESL, IDIESL )

            IF( INT( TDIESL( 1,JDSL ) ) .NE. YEAR ) THEN
                WRITE( MESG,94100 ) 
     &                 'DIESEL SALES profile', NIDX( IP,6 ),
     &                 'from parameter scheme', IP,
     &                 'has year',INT( TDIESL( 1,JDSL ) ),
     &                 'but modeling year specified is', YEAR
                CALL M3EXIT( 'MOINBLD', 0, 0, MESG, 2 )
            ENDIF

            IF( JDSL .LT. 0 ) THEN
                WRITE( MESG,94100 ) 
     &                 'DIESEL SALES profile', NIDX( IP,6 ),
     &                 'from parameter scheme', IP,
     &                 'is not found in diesel sales fractions table.'
                CALL M3EXIT( 'MOINBLD', 0, 0, MESG, 2 )
            ENDIF

            WRITE( MODEV, 93520 ) ( TDIESL( I,JDSL ), I = 2,51 )

        ENDIF

C.........  Trip length distribution
        IF( NIDX( IP,7 ) .EQ. 0 .AND.
     &      NIDX( IP,5 ) .GT. 0       ) THEN

            WRITE( MODEV, 93530 )
     &           ( TTPLEN( I, JTPLEN ), I=2,7 )
        ENDIF

C.............  Additional Correction factors for LDGVs
c           IF( NIDX( IP,10 ) .GT. 0 ) THEN
C! Hard for many reasons
c           ENDIF

        IPPREV = IP   ! Store IP for compsrison on next call

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

93000   FORMAT( A )
93005   FORMAT( 2 (I1, 2X) )
93010   FORMAT( 25(I1, 1X) )
93020   FORMAT( 8F8.4 )
93030   FORMAT( 7( 10F7.5 / 10F7.5 / 5F7.5 / ) 10F7.5 / 10F7.5 / 5F7.5 )
93040   FORMAT( 7( 10(1X,F4.3) / 10(1X,F4.3) / 5(1X,F4.3) / ),
     &             10(1X,F4.3) / 10(1X,F4.3) / 5(1X,F4.3)      )
93100   FORMAT( 4(I2,1X),2(I2,1X),I3,1X,2(I1,1X),4I1,1X,4I1,3(1X,A4) )
93101   FORMAT( F4.2 )
93102   FORMAT( F4.1 )
93103   FORMAT( F4.0 )
93200   FORMAT( 3( I2, 1X ), 4I1, 1X, 2I1, 1X, I4, 1X, 8I1 )
93400   FORMAT( I1,1X,I2,5 (1X,F4.1),2X,I1,15X,A1 )
93410   FORMAT( I1,1X,I2,12(1X,F4.1),2X,I1,15X,A1 )
93500   FORMAT( A, 1X, A1, 2F5.0, 2F5.1, 1X, I2, 3(1X,I1) )
93510   FORMAT( 4( F4.3, 1X ), I1 )
93520   FORMAT( 20F4.3 / 20F4.3 / 10F4.3 )
93530   FORMAT( 6(1X,F4.1) )
93600   FORMAT( 8F4.3 )

94000   FORMAT( A )
94100   FORMAT( 10( A, 1X, I5, 1X ) )
        END
