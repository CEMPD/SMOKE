
        SUBROUTINE RDMTREF( XDEV, IFIP, IRCLAS, ILINK,
     &                      NMON, NWEK, NDIU, NEND,
     &                      MONREF, WEKREF, DIUREF, ENDREF,
     &                      MDEX, WDEX, DDEX )

C***********************************************************************
C  subroutine body starts at line 176
C
C  DESCRIPTION:
C       Constructs sorted, structured mobile source temporal 
C       cross reference array structure; then constructs monthly,
C	weekly, and diurnal maps MDEX, WDEX, DDEX of sources to
C	profile numbers.
C
C  PRECONDITIONS REQUIRED:  
C       temporal cross-reference file opened on XDEV and positioned
C       at beginning of file
C	[ MON, WEK, DIU ]REFs and their counts read from MTPRO.
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       Models-3 I/O; 
C       FIND1, FIND2, FIND3 GETEFILE, GETNUM, GETYN, DAYMON, 
C
C  REVISION  HISTORY:
C       Prototype  1/96 by CJC.
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1996, MCNC--North Carolina Supercomputing Center
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
C****************************************************************************

      IMPLICIT NONE


C...........   INCLUDES:
        
        INCLUDE 'TMDIMS3.EXT'   !  temporal parameters
        INCLUDE 'MBDIMS3.EXT'   !  mobile-source dimensioning parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API constants


C...........   ARGUMENTS:

        INTEGER         XDEV                   !  unit number for xref file
        INTEGER         IFIP  ( NMSRC )        ! FIP codes for sources
        INTEGER         IRCLAS( NMSRC )        ! EPA road classification code
        INTEGER         ILINK ( NMSRC )        ! link number or 0 for non-link
        INTEGER         TFLAG ( NMSRC )        ! use month, week profiles??
        INTEGER         NMON, NWEK, NDIU, NEND ! numbers of profiles read
        INTEGER         MONREF( NMCOD )     ! Monthly         XREF table
        INTEGER         WEKREF( NWCOD )     ! Weekly          XREF table
        INTEGER         DIUREF( NDCOD )     ! Weekday-diurnal XREF table
        INTEGER         ENDREF( NDCOD )     ! Weekend-diurnal XREF table

        INTEGER         MDEX ( NMSRC, NVTYPE, NETYPE ) ! monthly       profile subscript
        INTEGER         WDEX ( NMSRC, NVTYPE, NETYPE ) ! weeky         profile subscript
        INTEGER         DDEX ( NMSRC, NVTYPE, NETYPE ) ! wkday diurnal profile subscript
        INTEGER         EDEX ( NMSRC, NVTYPE, NETYPE ) ! wkend-diurnal profile subscript


C...........   EXTERNAL FUNCTIONS and their descriptions:

        INTEGER         FIND1, FIND2, FIND3
        INTEGER         INDEX1
        INTEGER         TRIMLEN
        
        EXTERNAL    FIND1, FIND2, FIND3, INDEX1, TRIMLEN

                                            
C...........   PARAMETERS and their descriptions

        CHARACTER*16    BLANK16
        PARAMETER     ( BLANK16 = ' ' )


C...........   LOCAL VARIABLES and their descriptions:
C...........   Mobile Sources input and output arrays
        
        INTEGER         MREF ( NVTYPE, NETYPE )  ! monthly       profile subscript
        INTEGER         WREF ( NVTYPE, NETYPE )  ! weeky         profile subscript
        INTEGER         DREF ( NVTYPE, NETYPE )  ! wkday diurnal profile subscript

C...........   Temporal Profile XREF tables:  give monthly, weekly, hourly 
C...........   profile-indexes for lookup into temporal profile tables.
C.......   Before sorting:

        INTEGER         NPRFA
        INTEGER         INDXA( MXREF )	!  subscript (for sorting)
        INTEGER         FIPFA( MXREF )	!  FIP code
        INTEGER         RDTFA( MXREF )	!  road type (from MROADS3)
        INTEGER         LNKFA( MXREF )	!  LINK ID, or 0
        INTEGER         VHTFA( MXREF )	!  Vehicle type (index into VTYPE3)
        INTEGER         MPRFA( MXREF )	!  month-of-year profile index
        INTEGER         WPRFA( MXREF )	!  day-of-week   profile index
        INTEGER         DPRFA( MXREF )	!  hour-of-day   profile index
        INTEGER         CODEA( MXREF )	!  emis mode code (index into EMTNAM)

C.......   After sorting:
C.......   Default no-keys profile

        INTEGER         MPR00( NVTYPE, NETYPE )  !  month-in-year
        INTEGER         WPR00( NVTYPE, NETYPE )  !  day-in-week
        INTEGER         DPR00( NVTYPE, NETYPE )  !  hour-in-day (diurnal)

C.......   No-FIP, by-RDT profiles

        INTEGER         NPRNF
        INTEGER         RDTNF( NRCLAS )
        INTEGER         MPRNF( NRCLAS, NVTYPE, NETYPE )
        INTEGER         WPRNF( NRCLAS, NVTYPE, NETYPE )
        INTEGER         DPRNF( NRCLAS, NVTYPE, NETYPE )

C.......   No-county, (by-STATE FIP code-and-RDT) profiles

        INTEGER         NPRNC
        INTEGER         FIPNC( MXSREF )
        INTEGER         RDTNC( MXSREF )
        INTEGER         MPRNC( MXSREF, NVTYPE, NETYPE )
        INTEGER         WPRNC( MXSREF, NVTYPE, NETYPE )
        INTEGER         DPRNC( MXSREF, NVTYPE, NETYPE )

C.......   FIP-RDT-specific profiles

        INTEGER         NPRFC
        INTEGER         FIPFC( MXFREF )
        INTEGER         RDTFC( MXFREF )
        INTEGER         MPRFC( MXFREF, NVTYPE, NETYPE )
        INTEGER         WPRFC( MXFREF, NVTYPE, NETYPE )
        INTEGER         DPRFC( MXFREF, NVTYPE, NETYPE )

C.......   FIP-[RDTY-]-LNK-specific profiles

        INTEGER         NPRFL
        INTEGER         FIPFL( MXLREF )
        INTEGER         RDTFL( MXLREF )
        INTEGER         LNKFL( MXLREF )
        INTEGER         MPRFL( MXLREF, NVTYPE, NETYPE )
        INTEGER         WPRFL( MXLREF, NVTYPE, NETYPE )
        INTEGER         DPRFL( MXLREF, NVTYPE, NETYPE )


C...........   Other local variables

        LOGICAL         EFLAG   !  error-flag
        LOGICAL         MFLAG   !  use monthly profiles
        LOGICAL         WFLAG   !  use weekly profiles
        INTEGER         IOS
        INTEGER         IREC
        INTEGER         I, J, K, M, F, S, V
        INTEGER         FIP, RDT, VHT, LNK, COD, EMM, IMON, IWEK, IDIU

        CHARACTER*256   MESG    !  buffer for M3EXIT()


C***********************************************************************
C   begin body of subroutine RDMTREF
        
        DO  12  J = 1, NETYPE
        DO  11  I = 1, NVTYPE
            MPR00( I,J ) = IMISS3
            WPR00( I,J ) = IMISS3
            DPR00( I,J ) = IMISS3
            DO  7  S = 1, NRCLAS
                MPRNF( S,I,J ) = IMISS3
                WPRNF( S,I,J ) = IMISS3
                DPRNF( S,I,J ) = IMISS3
7           CONTINUE
            DO  8  S = 1, MXSREF
                MPRNC( S,I,J ) = IMISS3
                WPRNC( S,I,J ) = IMISS3
                DPRNC( S,I,J ) = IMISS3
8           CONTINUE
            DO  9  S = 1, MXFREF
                MPRFC( S,I,J ) = IMISS3
                WPRFC( S,I,J ) = IMISS3
                DPRFC( S,I,J ) = IMISS3
9           CONTINUE
            DO  10  S = 1, MXLREF
                MPRFL( S,I,J ) = IMISS3
                WPRFL( S,I,J ) = IMISS3
                DPRFL( S,I,J ) = IMISS3
10          CONTINUE
11      CONTINUE
12      CONTINUE


C.......   Read temporal-profile cross-reference file.

        CALL M3MSG2( 'Reading temporal-profile XREF file...' )

        EFLAG = .FALSE.
        IREC  = 1
        K     = 0

101     CONTINUE
            
            IREC = IREC + 1
            READ( XDEV, *, END=102, IOSTAT=IOS ) 
     &              FIP, RDT, LNK, VHT, EMM, COD, IMON, IWEK, IDIU
            
            IF ( IOS .NE. 0 ) THEN

                WRITE( MESG,94010 ) 'Line number ', IREC, 
     &                              'I/O status  ', IOS
                CALL M3MSG2( MESG )
                CALL M3EXIT( 'TMPMOBIL', 0, 0, 
     &                       'Error reading CROSS-REFERENCE file', 2 )
 
            ELSE IF ( IMON .LE. 0 ) THEN
 
                WRITE( MESG,94020 ) 
     &              'Bad IMON=', IMON, 
     &              'at record', IREC,
     &              'for FIP'  , FIP, 
     &              'RCLASS'   , RDT,
     &              'ILINK'    , LNK,
     &              'VTYPE'    , VHT,
     &              'EMMODE'   , EMM,
     &              'POLCODE'  , COD
                CALL M3MESG( MESG )
                EFLAG = .TRUE.
                GO TO  101
 
            ELSE IF ( IWEK .LE. 0 ) THEN
 
                WRITE( MESG,94020 ) 
     &              'Bad IWEK=', IWEK, 
     &              'at record', IREC,
     &              'for FIP'  , FIP, 
     &              'RCLASS'   , RDT,
     &              'ILINK'    , LNK,
     &              'VTYPE'    , VHT,
     &              'EMMODE'   , EMM,
     &              'POLCODE'  , COD
                CALL M3MESG( MESG )
                EFLAG = .TRUE.
                GO TO  101
 
            ELSE IF ( IDIU .LE. 0 ) THEN
 
                WRITE( MESG,94020 ) 
     &              'Bad IDIU=', IDIU, 
     &              'at record', IREC,
     &              'for FIP'  , FIP, 
     &              'RCLASS'   , RDT,
     &              'ILINK'    , LNK,
     &              'VTYPE'    , VHT,
     &              'EMMODE'   , EMM,
     &              'POLCODE'  , COD
                CALL M3MESG( MESG )
                EFLAG = .TRUE.
                GO TO  101
 
            END IF

            IF ( K .GE. MXREF ) THEN !  FIP-independent fallback profiles

                WRITE( MESG,94020 ) 
     &              'XREF table overflow at record', IREC
                CALL M3MESG( MESG )
                EFLAG = .TRUE.
 
            ELSE
                
                K          = K + 1
                INDXA( K ) = K
                FIPFA( K ) = FIP
                RDTFA( K ) = RDT
                VHTFA( K ) = VHT
                LNKFA( K ) = LNK
                MPRFA( K ) = IMON
                WPRFA( K ) = IWEK
                DPRFA( K ) = IDIU
                IF ( EMM .EQ. 1 ) THEN		 !  "exhaust" mode:  all polls
                    I = FIND1( COD, NIPOL, PCODE )
                    IF ( I .GT. 0 ) THEN
                        CODEA( K ) = I
                    ELSE
                        K = K - 1	!  back out of this entry 
                        WRITE( MESG,94010 ) 
     &                  'Pollutant code', COD, 'not being modeled'
                        CALL M3MESG( MESG )
                    END IF
                ELSE IF ( EMM .LE. NEMODE ) THEN !  evp | rnl | rst
                    IF ( COD .EQ. 43104 ) THEN	 !  voc--ok
                        CODEA( K ) = NIPOL + EMM - 1
                    ELSE
                        K = K - 1	!  back out of this entry 
                        WRITE( MESG,94010 ) 
     &                  'Pollutant code', COD, 
     &                  'incorrect for emis mode ' // MEMODE( EMM )
                        CALL M3MESG( MESG )
                    END IF
                ELSE
                    K = K - 1	!  back out of this entry 
                    WRITE( MESG,94010 ) 
     &              'Illegal emissions mode', EMM
                    CALL M3MESG( MESG )
                END IF

            END IF      !  if i=0; else if fip=0; else...

            GO TO  101   !  to head of read-XREF loop
                            
102     CONTINUE
                    
        IF ( EFLAG ) THEN
            CALL M3EXIT( 'TMPMOBIL', 0, 0, 'Bad XREF file.', 2 )
        END IF
        
        NPRFA = K


C.......   Sort and process temporal-profile cross-reference file.

        CALL M3MSG2( 'Processing temporal-profile XREF file...' )
        CALL SORTI4( NPRFA, INDXA, FIPFA, RDTFA, LNKFA, VHTFA )

        EFLAG = .FALSE.
        IREC  = 0
        I     = 0
        NPRNF = 0
        NPRNC = 0
        NPRFC = 0
        
        DO  199  K = 1, NPRFA
            
            J    = INDXA( K )
            FIP  = FIPFA( J )
            RDT  = RDTFA( J )
            VHT  = VHTFA( J )
            LNK  = LNKFA( J )
            EMM  = CODEA( J )
            IMON = MPRFA( J )
            IWEK = WPRFA( J )
            IDIU = DPRFA( J )

            IF ( FIP .EQ. 0 ) THEN !  FIP-independent fallback profiles
                
                IF ( RDT .EQ. 0 ) THEN  !  ultimate fallback profiles
                    
                    IF ( VHT .EQ. 0 ) THEN
                        IF ( EMM .EQ. 0 ) THEN
                            DO 111  M = 1, NETYPE
                            DO 110  V = 1, NVTYPE
                                MPR00( V,M ) = IMON
                                WPR00( V,M ) = IWEK
                                DPR00( V,M ) = IDIU
110                         CONTINUE
111                         CONTINUE
                        ELSE
                            DO 113  V = 1, NVTYPE
                                MPR00( V,EMM ) = IMON
                                WPR00( V,EMM ) = IWEK
                                DPR00( V,EMM ) = IDIU
                                DO 112  M = 1, NETYPE
                                    IF ( MPR00( V,M ).EQ.IMISS3 ) 
     &                                   MPR00( V,M ) = IMON
                                    IF ( WPR00( V,M ).EQ.IMISS3 ) 
     &                                   WPR00( V,M ) = IWEK
                                    IF ( DPR00( V,M ).EQ.IMISS3 ) 
     &                                   DPR00( V,M ) = IDIU
112                             CONTINUE
113                         CONTINUE
                        END IF
                    ELSE
                        IF ( EMM .EQ. 0 ) THEN
                            DO  115  M = 1, NETYPE
                                MPR00( VHT,M ) = IMON
                                WPR00( VHT,M ) = IWEK
                                DPR00( VHT,M ) = IDIU
                                DO  114  V = 1, NVTYPE
                                    IF ( MPR00( V,M ).EQ.IMISS3 ) 
     &                                   MPR00( V,M ) = IMON
                                    IF ( WPR00( V,M ).EQ.IMISS3 ) 
     &                                   WPR00( V,M ) = IWEK
                                    IF ( DPR00( V,M ).EQ.IMISS3 ) 
     &                                   DPR00( V,M ) = IDIU
114                             CONTINUE
115                         CONTINUE
                        ELSE
                            MPR00( VHT,EMM ) = IMON
                            WPR00( VHT,EMM ) = IWEK
                            DPR00( VHT,EMM ) = IDIU
                            DO 117  M = 1, NETYPE
                            DO 116  V = 1, NVTYPE
                                IF ( MPR00( V,M ).EQ.IMISS3 ) 
     &                               MPR00( V,M ) = IMON
                                IF ( WPR00( V,M ).EQ.IMISS3 ) 
     &                               WPR00( V,M ) = IWEK
                                IF ( DPR00( V,M ).EQ.IMISS3 ) 
     &                               DPR00( V,M ) = IDIU
116                         CONTINUE
117                         CONTINUE
                        END IF
                    END IF
                    
                ELSE                    !  RDT?-dependent-only profiles
                
                    IF ( NPRNF .EQ. 0   .OR.
     &                   RDT .NE. RDTNF( NPRNF ) ) THEN
                        
                        NPRNF = NPRNF + 1
                        IF ( NPRNF .GT. NRCLAS ) THEN
                            EFLAG = .TRUE.
                            WRITE( MESG,94010 )
     &                              'Line number ', IREC,
     &                              'FIP         ', FIP,
     &                              'RDT         ', RDT,
     &                              'VHT         ', VHT
                            CALL M3MESG( MESG )
                            CALL M3MESG( 
     &                           '   Subscript NPRNF out of range' )
                            GO TO 199
                        END IF

                        RDTNF( NPRNF ) = RDT

                    END IF
                    
                    IF ( VHT .EQ. 0 ) THEN
                        IF ( EMM .EQ. 0 ) THEN
                            DO  121  M = 1, NETYPE
                            DO  120  V = 1, NVTYPE
                                MPRNF( NPRNF,V,M ) = IMON
                                WPRNF( NPRNF,V,M ) = IWEK
                                DPRNF( NPRNF,V,M ) = IDIU
120                         CONTINUE
121                         CONTINUE
                        ELSE
                            DO  123  V = 1, NVTYPE
                                MPRNF( NPRNF,V,EMM ) = IMON
                                WPRNF( NPRNF,V,EMM ) = IWEK
                                DPRNF( NPRNF,V,EMM ) = IDIU
                                DO  122  M = 1, NETYPE
                                    IF ( MPRNF( NPRNF,V,M ).EQ.IMISS3 )
     &                                   MPRNF( NPRNF,V,M ) = IMON
                                    IF ( WPRNF( NPRNF,V,M ).EQ.IMISS3 )
     &                                   WPRNF( NPRNF,V,M ) = IWEK
                                    IF ( DPRNF( NPRNF,V,M ).EQ.IMISS3 )
     &                                   DPRNF( NPRNF,V,M ) = IDIU
122                             CONTINUE
123                         CONTINUE
                        END IF
                    ELSE
                        IF ( EMM .EQ. 0 ) THEN
                            DO  125  M = 1, NETYPE
                                MPRNF( NPRNF,VHT,M ) = IMON
                                WPRNF( NPRNF,VHT,M ) = IWEK
                                DPRNF( NPRNF,VHT,M ) = IDIU
                                DO  124  V = 1, NVTYPE
                                    IF ( MPRNF( NPRNF,V,M ).EQ.IMISS3 ) 
     &                                   MPRNF( NPRNF,V,M ) = IMON
                                    IF ( WPRNF( NPRNF,V,M ).EQ.IMISS3 ) 
     &                                   WPRNF( NPRNF,V,M ) = IWEK
                                    IF ( DPRNF( NPRNF,V,M ).EQ.IMISS3 )
     &                                   DPRNF( NPRNF,V,M ) = IDIU
124                             CONTINUE
125                         CONTINUE
                        ELSE
                            MPRNF( NPRNF,VHT,EMM ) = IMON
                            WPRNF( NPRNF,VHT,EMM ) = IWEK
                            DPRNF( NPRNF,VHT,EMM ) = IDIU
                            DO  127  M = 1, NETYPE
                            DO  126  V = 1, NVTYPE
                                IF ( MPRNF( NPRNF,V,M ).EQ.IMISS3 ) 
     &                               MPRNF( NPRNF,V,M ) = IMON
                                IF ( WPRNF( NPRNF,V,M ).EQ.IMISS3 ) 
     &                               WPRNF( NPRNF,V,M ) = IWEK
                                IF ( DPRNF( NPRNF,V,M ).EQ.IMISS3 )
     &                               DPRNF( NPRNF,V,M ) = IDIU
126                         CONTINUE
127                         CONTINUE
                        END IF
                    END IF

                END IF
            
            ELSE IF ( MOD( FIP,1000 ) .EQ. 0 ) THEN     !  state-ASC dependent
            
                FIP = FIP/1000
                
                IF ( NPRNC .EQ. 0               .OR.
     &               FIP   .NE. FIPNC( NPRNC )  .OR.
     &               RDT   .NE. RDTNC( NPRNC ) ) THEN

                    NPRNC = NPRNC + 1
                    IF ( NPRNC .GT. MXSREF ) THEN
                        EFLAG = .TRUE.
                        WRITE( MESG,94010 )
     &                          'Line number ', IREC,
     &                          'FIP         ', FIP,
     &                          'RDT         ', RDT,
     &                          'VHT         ', VHT
                        CALL M3MESG( MESG )
                        CALL M3MESG(
     &                       '   Subscript NPRNC out of range' )
                        GO TO 199
                    END IF

                    FIPNC( NPRNC ) = FIP
                    RDTNC( NPRNC ) = RDT

                END IF
                    
                IF ( VHT .EQ. 0 ) THEN
                    IF ( EMM .EQ. 0 ) THEN
                        DO 131  M = 1, NETYPE
                        DO 130  I = 1, NVTYPE
                            MPRNC( NPRNC,V,M ) = IMON
                            WPRNC( NPRNC,V,M ) = IWEK
                            DPRNC( NPRNC,V,M ) = IDIU
130                     CONTINUE                    
131                     CONTINUE                    
                    ELSE
                        DO 133  I = 1, NVTYPE
                            MPRNC( NPRNC,V,EMM ) = IMON
                            WPRNC( NPRNC,V,EMM ) = IWEK
                            DPRNC( NPRNC,V,EMM ) = IDIU
                            DO  132  M = 1, NETYPE
                                IF ( MPRNC( NPRNC,V,M ).EQ.IMISS3 ) 
     &                               MPRNC( NPRNC,V,M ) = IMON
                                IF ( WPRNC( NPRNC,V,M ).EQ.IMISS3 ) 
     &                               WPRNC( NPRNC,V,M ) = IWEK
                                IF ( DPRNC( NPRNC,V,M ).EQ.IMISS3 )
     &                               DPRNC( NPRNC,V,M ) = IDIU
132                         CONTINUE                    
133                     CONTINUE                   
                                                
                    END IF
                ELSE
                    IF ( EMM .EQ. 0 ) THEN
                        DO  135  M = 1, NETYPE
                            MPRNC( NPRNC,VHT,M ) = IMON
                            WPRNC( NPRNC,VHT,M ) = IWEK
                            DPRNC( NPRNC,VHT,M ) = IDIU
                            DO 134  V = 1, NVTYPE
                                IF ( MPRNC( NPRNC,V,M ).EQ.IMISS3 ) 
     &                               MPRNC( NPRNC,V,M ) = IMON
                                IF ( WPRNC( NPRNC,V,M ).EQ.IMISS3 ) 
     &                               WPRNC( NPRNC,V,M ) = IWEK
                                IF ( DPRNC( NPRNC,V,M ).EQ.IMISS3 )
     &                               DPRNC( NPRNC,V,M ) = IDIU
134                         CONTINUE                    
135                     CONTINUE                    
                    ELSE
                        MPRNC( NPRNC,VHT,EMM ) = IMON
                        WPRNC( NPRNC,VHT,EMM ) = IWEK
                        DPRNC( NPRNC,VHT,EMM ) = IDIU
                        DO 137  M = 1, NETYPE
                        DO 136  V = 1, NVTYPE
                            IF ( MPRNC( NPRNC,V,M ).EQ.IMISS3 ) 
     &                           MPRNC( NPRNC,V,M ) = IMON
                            IF ( WPRNC( NPRNC,V,M ).EQ.IMISS3 ) 
     &                           WPRNC( NPRNC,V,M ) = IWEK
                            IF ( DPRNC( NPRNC,V,M ).EQ.IMISS3 )
     &                           DPRNC( NPRNC,V,M ) = IDIU
136                     CONTINUE                    
137                     CONTINUE                    
                    END IF
                END IF

            ELSE IF ( LNK .EQ. 0 ) THEN    !  FIP-RDT dependent profiles

                IF ( NPRFC .EQ. 0             .OR.
     &               FIP .NE. FIPFC( NPRFC )  .OR.
     &               RDT .NE. RDTFC( NPRFC ) ) THEN
                                           
                    NPRFC = NPRFC + 1
                    IF ( NPRFC .GT. MXFREF ) THEN
                        EFLAG = .TRUE.
                        WRITE( MESG,94010 )
     &                          'Line number ', IREC,
     &                          'FIP         ', FIP,
     &                          'RDT         ', RDT,
     &                          'VHT         ', VHT
                        CALL M3MESG( MESG )
                        CALL M3MESG(
     &                      '   Subscript NPRFC out of range' )
                        GO TO 199
                    END IF

                    FIPFC( NPRFC ) = FIP
                    RDTFC( NPRFC ) = RDT

                END IF

                IF ( VHT .EQ. 0 ) THEN
                    IF ( EMM .EQ. 0 ) THEN
                        DO 141  M = 1, NETYPE
                        DO 140  V = 1, NVTYPE
                            MPRFC( NPRFC,V,M ) = IMON
                            WPRFC( NPRFC,V,M ) = IWEK
                            DPRFC( NPRFC,V,M ) = IDIU
140                     CONTINUE
141                     CONTINUE
                    ELSE
                        DO 143  V = 1, NVTYPE
                            MPRFC( NPRFC,V,EMM ) = IMON
                            WPRFC( NPRFC,V,EMM ) = IWEK
                            DPRFC( NPRFC,V,EMM ) = IDIU
                            DO 142  M = 1, NETYPE
                                IF ( MPRFC( NPRFC,V,M ).EQ.IMISS3 ) 
     &                               MPRFC( NPRFC,V,M ) = IMON
                                IF ( WPRFC( NPRFC,V,M ).EQ.IMISS3 ) 
     &                               WPRFC( NPRFC,V,M ) = IWEK
                                IF ( DPRFC( NPRFC,V,M ).EQ.IMISS3 )
     &                               DPRFC( NPRFC,V,M ) = IDIU
142                         CONTINUE
143                     CONTINUE
                    END IF
                ELSE
                    IF ( EMM .EQ. 0 ) THEN
                        DO  145  M = 1, NETYPE
                            MPRFC( NPRFC,VHT,M ) = IMON
                            WPRFC( NPRFC,VHT,M ) = IWEK
                            DPRFC( NPRFC,VHT,M ) = IDIU
                            DO  144  V = 1, NVTYPE
                                IF ( MPRFC( NPRFC,V,M ).EQ.IMISS3 )
     &                               MPRFC( NPRFC,V,M ) = IMON
                                IF ( WPRFC( NPRFC,V,M ).EQ.IMISS3 ) 
     &                               WPRFC( NPRFC,V,M ) = IWEK
                                IF ( DPRFC( NPRFC,V,M ).EQ.IMISS3 )
     &                               DPRFC( NPRFC,V,M ) = IDIU
144                         CONTINUE
145                     CONTINUE
                    ELSE
                        MPRFC( NPRFC,VHT,EMM ) = IMON
                        WPRFC( NPRFC,VHT,EMM ) = IWEK
                        DPRFC( NPRFC,VHT,EMM ) = IDIU
                        DO 147  M = 1, NETYPE
                        DO 146  V = 1, NVTYPE
                            IF ( MPRFC( NPRFC,V,M ).EQ.IMISS3 )
     &                           MPRFC( NPRFC,V,M ) = IMON
                            IF ( WPRFC( NPRFC,V,M ).EQ.IMISS3 ) 
     &                           WPRFC( NPRFC,V,M ) = IWEK
                            IF ( DPRFC( NPRFC,V,M ).EQ.IMISS3 )
     &                           DPRFC( NPRFC,V,M ) = IDIU
146                     CONTINUE
147                     CONTINUE
                    END IF
                END IF

            ELSE                 !  FIP-LNK dependent profiles

                IF ( NPRFL .EQ. 0             .OR.
     &               FIP .NE. FIPFL( NPRFL )  .OR.
     &               LNK .NE. LNKFL( NPRFL ) ) THEN
                                           
                    NPRFL = NPRFL + 1
                    IF ( NPRFL .GT. MXLREF ) THEN
                        EFLAG = .TRUE.
                        WRITE( MESG,94010 )
     &                          'Line number ', IREC,
     &                          'FIP         ', FIP,
     &                          'RDT         ', RDT,
     &                          'LNK         ', LNK
                        CALL M3MESG( MESG )
                        CALL M3MESG(
     &                      '   Subscript NPRFC out of range' )
                        GO TO 199
                    END IF

                    FIPFL( NPRFL ) = FIP
                    RDTFL( NPRFL ) = RDT
                    LNKFL( NPRFL ) = LNK

                END IF

                IF ( VHT .EQ. 0 ) THEN
                    IF ( EMM .EQ. 0 ) THEN
                        DO 151  M = 1, NETYPE
                        DO 150  V = 1, NVTYPE
                            MPRFL( NPRFL,V,M ) = IMON
                            WPRFL( NPRFL,V,M ) = IWEK
                            DPRFL( NPRFL,V,M ) = IDIU
150                     CONTINUE
151                     CONTINUE
                    ELSE
                        DO 153  V = 1, NVTYPE
                            MPRFL( NPRFL,V,EMM ) = IMON
                            WPRFL( NPRFL,V,EMM ) = IWEK
                            DPRFL( NPRFL,V,EMM ) = IDIU
                            DO 152  M = 1, NETYPE
                                IF ( MPRFL( NPRFL,V,M ).EQ.IMISS3 )
     &                               MPRFL( NPRFL,V,M ) = IMON
                                IF ( WPRFL( NPRFL,V,M ).EQ.IMISS3 ) 
     &                               WPRFL( NPRFL,V,M ) = IWEK
                                IF ( DPRFL( NPRFL,V,M ).EQ.IMISS3 )
     &                               DPRFL( NPRFL,V,M ) = IDIU
152                         CONTINUE
153                     CONTINUE
                    END IF
                ELSE
                    IF ( EMM .EQ. 0 ) THEN
                        DO 155  M = 1, NETYPE
                            MPRFL( NPRFL,VHT,M ) = IMON
                            WPRFL( NPRFL,VHT,M ) = IWEK
                            DPRFL( NPRFL,VHT,M ) = IDIU
                            DO 154  V = 1, NVTYPE
                                IF ( MPRFL( NPRFL,V,M ).EQ.IMISS3 )
     &                               MPRFL( NPRFL,V,M ) = IMON
                                IF ( WPRFL( NPRFL,V,M ).EQ.IMISS3 ) 
     &                               WPRFL( NPRFL,V,M ) = IWEK
                                IF ( DPRFL( NPRFL,V,M ).EQ.IMISS3 )
     &                               DPRFL( NPRFL,V,M ) = IDIU
154                         CONTINUE
155                     CONTINUE
                    ELSE
                        MPRFL( NPRFL,VHT,EMM ) = IMON
                        WPRFL( NPRFL,VHT,EMM ) = IWEK
                        DPRFL( NPRFL,VHT,EMM ) = IDIU
                        DO 157  M = 1, NETYPE
                        DO 156  V = 1, NVTYPE
                            IF ( MPRFL( NPRFL,V,M ).EQ.IMISS3 )
     &                           MPRFL( NPRFL,V,M ) = IMON
                            IF ( WPRFL( NPRFL,V,M ).EQ.IMISS3 ) 
     &                           WPRFL( NPRFL,V,M ) = IWEK
                            IF ( DPRFL( NPRFL,V,M ).EQ.IMISS3 )
     &                           DPRFL( NPRFL,V,M ) = IDIU
156                     CONTINUE
157                     CONTINUE
                    END IF
                END IF

            END IF      !  if i=0; else if fip=0; else...

199     CONTINUE
                    
        IF ( EFLAG ) THEN
            CALL M3EXIT( 'TMPMOBIL:RDMTREF', 0, 0, 'Bad XREF file.', 2 )
        END IF


C...........   Construct mapping of source number to effective
C...........   composite temporal profile number and to vehicle mixes.
              
        CALL M3MSG2( 'Building source indexes...' )
        EFLAG = .FALSE.

        DO  299  S = 1, NMSRC
        
            MFLAG = ( MOD( TFLAG( S ), MTPRFAC ) .EQ. 0 )
            WFLAG = ( MOD( TFLAG( S ), WTPRFAC ) .EQ. 0 )
            
            IF ( .NOT. MFLAG ) THEN
                DO  202  M = 1, NETYPE
                DO  201  V = 1, NVTYPE
                    MDEX( S,V,M ) = 0
201             CONTINUE
202             CONTINUE
            END IF

            IF ( .NOT. WFLAG ) THEN
                DO  204  M = 1, NETYPE
                DO  203  V = 1, NVTYPE
                    WDEX( S,V,M ) = 0
203             CONTINUE
204             CONTINUE
            END IF

            FIP = IFIP  ( S )
            RDT = IRCLAS( S )
            LNK = ILINK ( S )

C...........   Search for FIP-LNK match; else for FIP-RDT match;
C...........   else for STATE-RDT match; else for     RDT match; 
C...........   else use ultimate fall back default.
            
            F = FIND3( FIP, RDT, LNK, NPRFL, FIPFL, RDTFL, LNKFL ) ! FIP-RDT-LNK

            IF ( F .GT. 0 ) THEN

                DO  212  M = 1, NETYPE
                DO  211  V = 1, NVTYPE
                    MREF( V,M ) = MPRFL( F,V,M )
                    WREF( V,M ) = WPRFL( F,V,M )
                    DREF( V,M ) = DPRFL( F,V,M )
211             CONTINUE
212             CONTINUE

                GO TO  266       !  to end of TREF lookups

            END IF      !  if find3() FIP-RDT-LNK-VHT succeeded

            F = FIND2( FIP, RDT, NPRFC, FIPFC, RDTFC ) ! FIP-LNK

            IF ( F .GT. 0 ) THEN

                DO  222  M = 1, NETYPE
                DO  221  V = 1, NVTYPE
                    MREF( V,M ) = MPRFC( F,V,M )
                    WREF( V,M ) = WPRFC( F,V,M )
                    DREF( V,M ) = DPRFC( F,V,M )
221             CONTINUE
222             CONTINUE

                GO TO  266       !  to end of TREF lookups

            END IF      !  if find3() FIP-and-RDT succeeded

            F = FIND2( FIP/1000, RDT, NPRNC, FIPNC, RDTNC )  ! STATE-RDT
  
            IF ( F .GT. 0 ) THEN

                DO  233  M = 1, NETYPE
                DO  232  V = 1, NVTYPE
                    MREF( V,M ) = MPRNC( F,V,M )
                    WREF( V,M ) = WPRNC( F,V,M )
                    DREF( V,M ) = DPRNC( F,V,M )
232             CONTINUE
233             CONTINUE

                GO TO  266       !  to end of TREF lookups

            END IF      !  if find3() STATE-and-RDT succeeded

            F = FIND1( RDT, NPRNF, RDTNF )    ! RDT search
  
            IF ( F .GT. 0 ) THEN

                DO  244  M = 1, NETYPE
                DO  243  V = 1, NVTYPE
                    MREF( V,M ) = MPRNF( F,V,M )
                    WREF( V,M ) = WPRNF( F,V,M )
                    DREF( V,M ) = DPRNF( F,V,M )
243             CONTINUE
244             CONTINUE

            ELSE                        !  fallback to RDT-independent default

                DO  255  M = 1, NETYPE
                DO  254  V = 1, NVTYPE
                    MREF( V,M ) = MPR00( V,M )
                    WREF( V,M ) = WPR00( V,M )
                    DREF( V,M ) = DPR00( V,M )
254             CONTINUE
255             CONTINUE

            END IF      !  if find2() RDT (only) succeeded, or not

266         CONTINUE    !  check index results:

            DO  277  M = 1, NETYPE
            DO  276  V = 1, NVTYPE

                IF ( MFLAG ) THEN
                    MDEX( S,V,M ) = FIND1( MREF( V,M ), NMON, MONREF )
                    IF( MDEX( S,V,M ) .LT. 0 ) THEN
                        WRITE( MESG,94010 )
     &                  'Cross-reference index not found for MREF=',
     &                   MREF( V,M )
                        CALL M3MESG( MESG )
                        EFLAG = .TRUE.
                    END IF
                END IF

                IF ( WFLAG ) THEN
                    WDEX( S,V,M ) = FIND1( WREF( V,M ), NWEK, WEKREF )
                    IF( WDEX( S,V,M ) .LT. 0 ) THEN
                        WRITE( MESG,94010 )
     &                  'Cross-reference index not found for WREF=',
     &                   WREF( V,M )
                        CALL M3MESG( MESG )
                        EFLAG = .TRUE.
                    END IF
                END IF

                DDEX( S,V,M ) = FIND1( DREF( V,M ), NDIU, DIUREF )
                IF( DDEX( S,V,M ) .LT. 0 .AND. NDIU .GT. 0 ) THEN
                    WRITE( MESG,94010 )
     &              'Cross-reference index not found for DREF=',
     &               DREF( V,M )
                    CALL M3MESG( MESG )
                    EFLAG = .TRUE.
                END IF

                EDEX( S,V,M ) = FIND1( DREF( V,M ), NEND, ENDREF )
                IF( EDEX( S,V,M ) .LT. 0 .AND. NEND .GT. 0 ) THEN
                    WRITE( MESG,94010 )
     &              'Cross-reference index not found for EREF=',
     &               DREF( V,M )
                    CALL M3MESG( MESG )
                    EFLAG = .TRUE.
                END IF

276         CONTINUE!  end index lookup loop on species V
277         CONTINUE!  end index lookup loop on emissions type M

299     CONTINUE        !  end loop on sources S, constructing TNDX

        
C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )

92010   FORMAT( 5X, A, I10 )


C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10 ( A, :, I5, :, 2X ) )

94020   FORMAT( A, I5, 1X, A, I7, 1X, 10( A, I7, :, 1X ) )

        END

