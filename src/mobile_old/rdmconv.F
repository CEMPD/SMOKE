
        SUBROUTINE RDMCONV( KDEV, TDEV,
     &                      FAC0, FIP1, RDT1, FAC1, N1,
     &                            FIP2, RDT2, FAC2, N2,
     &                                  RDT3, FAC3, N3  )

C***********************************************************************
C  function body starts at line  124
C
C  FUNCTION:
C       Reads the indicated volatile organic conversion file,
C       sorts it, and returns the number of profile entries found.
C       In case of overflow, returns negative of number of profile
C       entries found.
C
C  PRECONDITIONS REQUIRED:
C       volatile conversion file already opened on unit KDEV
C       speciation list file already opened on unit TDEV
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C
C  REVISION  HISTORY:
C       prototype 7/96 by MRH
C
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

      INCLUDE 'CHDGEN3.EXT'     ! emis chem info (both inventory and model)
      INCLUDE 'MBDIMS3.EXT'     ! mobile-specific information
      INCLUDE 'PARMS3.EXT'      ! I/O API constants
      INCLUDE 'FDESC3.EXT'      ! I/O API file description data structure
      INCLUDE 'IODECL3.EXT'     ! I/O API function declarations

C...........   ARGUMENTS and their descriptions:

        INTEGER     KDEV           !  unit number for speciation profile file
        INTEGER     TDEV           !  unit number for speciation list file
        INTEGER     N1             !  number of entries in group 1
        INTEGER     N2             !  number of entries in group 2
        INTEGER     N3             !  number of entries in group 3
        INTEGER     NRDT           !  indicator for valid road class
        INTEGER     FIP1( MXSCNV ) ! sorted FIP codes in group 1
        INTEGER     FIP2( MXSCNV ) ! sorted FIP codes in group 1
        INTEGER     RDT1( MXSCNV ) ! sorted Road Class codes in group 1
        INTEGER     RDT2( MXSCNV ) ! sorted Road Class codes in group 2
        INTEGER     RDT3( NRCLAS ) ! sorted Road Class codes in group 2

        REAL  FAC0(         NVTYPE, MXEPRC ) ! sorted conv factors in group 0
        REAL  FAC1( MXSCNV, NVTYPE, MXEPRC ) ! sorted conv factors in group 1
        REAL  FAC2( MXSCNV, NVTYPE, MXEPRC ) ! sorted conv factors in group 2
        REAL  FAC3( NRCLAS, NVTYPE, MXEPRC ) ! sorted conv factors in group 3

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER*2 CRLF
        INTEGER     CSVNAME
        INTEGER     FIND1
        INTEGER     INDEX1
        INTEGER     STR2INT
        REAL        STR2REAL

        EXTERNAL    CRLF, CSVNAME, FIND1, INDEX1, STR2INT, STR2REAL

C...........   PARAMETERS

        CHARACTER*5 BLANK5
        PARAMETER ( BLANK5 = ' ' )

C...........   SCRATCH LOCAL VARIABLES and their descriptions:

        INTEGER         I, I1, I2, I3, J, S ! Counters

        INTEGER         EM               ! Temporary emission mode index
        INTEGER         FIP              ! Temporary FIP code
        INTEGER         IOS              ! I/O Status code
        INTEGER         IREC             ! Record number of input file
        INTEGER         NEMODE           ! Number of emissions modes
        INTEGER         NTOT             ! Total items read in
        INTEGER         RDT              ! Temporary road class
        INTEGER         VHT              ! Temporary vehicle type

        INTEGER         EMA ( MXSCNV )   ! unsorted emissions modes index
        INTEGER         FIPA( MXSCNV )   ! unsorted FIPS
        INTEGER         INDX( MXSCNV )   ! index for sorting
        INTEGER         RDTA( MXSCNV )   ! unsorted road classes
        INTEGER         VHTA( MXSCNV )   ! unsorted vehicle types
  
        REAL            FAC              ! Temporary conversion factor
        REAL            FACA( MXSCNV )   ! unsorted conversion factors

        LOGICAL         EFLAG            ! Error flag

        CHARACTER*8     MEMODE  ( MXEPRC )        !  Emission modes
        CHARACTER*128   LINE
        CHARACTER*256   MESG

C***********************************************************************
C   begin body of function  RDMCONV

        NEMODE = CSVNAME( TDEV, .FALSE., .FALSE., .TRUE. , MEMODE )

C.........  Read volatile pollutants conversion factors file
        IREC  =  0
        I     =  0
        EFLAG = .FALSE.
77      CONTINUE        !  head of the KDEV-read loop

            READ( KDEV, 93010, END=99, IOSTAT=IOS ) LINE

            IREC = IREC + 1

            IF ( IOS .NE. 0 ) THEN

                EFLAG = .TRUE.
                WRITE( MESG, 94010 )
     &              'I/O error', IOS, 
     &              'reading VOLATILE CONVERSION line', IREC
                CALL M3MESG( MESG )

            ELSE

                FIP = STR2INT ( LINE(  1:5  ) )
                RDT = STR2INT ( LINE( 11:13 ) )
                VHT = STR2INT ( LINE(  9:10 ) )
                EM  = INDEX1  ( LINE( 14:16 ), NEMODE, MEMODE )
                FAC = STR2REAL( LINE( 19:25 ) )

                IF( RDT .NE. 0 ) THEN  ! Make sure we keep the defaults
                    NRDT = FIND1( RDT, NRCLAS, MROADS3 )
                ELSE
                    NRDT = 1
                ENDIF

C.................  NOTE: roundabout gotos prevent HP compiler bug
                IF( NRDT .LE. 0 ) THEN
                    GO TO 88   !  to head of the KDEV-read loop

                ELSEIF( EM .LE. 0 ) THEN
                    GO TO 88   !  to head of the KDEV-read loop

                ELSEIF( FAC .LT. 0. ) THEN
                    EFLAG = .TRUE.

                    WRITE( MESG,94020 )
     &                'Bad factor at line', IREC,
     &                'of VOLATILE CONVERSION FACTOR file =', FAC
                    CALL M3MSG2( MESG )

                    GO TO 88  ! to head of KDEV-read loop
                ENDIF

                I = I + 1
                IF ( I .LE. MXSCNV ) THEN

                    INDX( I ) = I
                    FIPA( I ) = FIP
                    RDTA( I ) = RDT
                    VHTA( I ) = VHT
                    EMA ( I ) = EM
                    FACA( I ) = FAC

                END IF

            END IF

88          GO TO  77   !  to head of loop

99      CONTINUE        !  end of the KDEV-read loop

        NTOT = I

        IF( NTOT .GT. MXSCNV ) THEN                
            EFLAG = .TRUE.
            WRITE( MESG, 94010 ) 
     &        'Max VOLATILE CONVERSION FACTOR table size (MXSCNV)=', 
     &        MXSCNV, CRLF() // BLANK5 // 'Actual number = ', NTOT
            CALL M3MESG( MESG )
        ENDIF

        IF ( EFLAG ) THEN
            CALL M3EXIT( 'RDMCONV', 0, 0, 
     &           'Error reading VOLATILE CONVERSION file', 2 )

        ENDIF

        CALL SORTI2( NTOT, INDX, FIPA, RDTA )

C.......   Initialize category-specific conversion factors

        DO 107 J = 1, MXEPRC
            DO 105 I = 1, NVTYPE

                FAC0( I,J ) = 1.0

                DO 101 S = 1, NRCLAS
                    FAC3( S,I,J ) = 1.0
  101           CONTINUE

                DO 103 S = 1, MXSCNV
                    FAC2( S,I,J ) = 1.0
                    FAC1( S,I,J ) = 1.0
  103           CONTINUE

  105       CONTINUE
  107   CONTINUE

        I1 =  0
        I2 =  0
        I3 =  0
        DO 201 I = 1, NTOT

            J   = INDX( I )
            FIP = FIPA( J )
            RDT = RDTA( J )
            VHT = VHTA( J )
            EM  = EMA ( J )
            FAC = FACA( J )

            IF ( FIP .EQ. 0 ) THEN !  FIP-independent fallback profiles
 
                IF ( RDT .EQ. 0 ) THEN  !  ultimate fallback profiles

                    CALL GRPFAC( VHT, EM, FAC, 
     &                           FAC0, 1, NVTYPE, MXEPRC,
     &                           1, NVTYPE, NEMODE )
 
                ELSE                    !  RDT?-dependent-only profiles

                    IF( I3  .EQ. 0          .OR.        ! First instance
     &                  RDT .NE. RDT3( I3 )      ) THEN ! Not repeat

                        I3 = I3 + 1
                        IF ( I3 .LE. NRCLAS ) THEN

                            RDT3( I3 ) = RDT
 
                            CALL GRPFAC( VHT, EM, FAC,
     &                                   FAC3, NRCLAS, NVTYPE, MXEPRC,
     &                                   I3  , NVTYPE, NEMODE )

                        END IF
                    END IF
                END IF

            ELSE IF ( MOD( FIP,1000 ) .EQ. 0 ) THEN     ! State dependent
 
                FIP = ( FIP/1000 ) * 1000

                IF( I2  .EQ. 0          .OR.         ! First instance
     &              FIP .NE. FIP2( I2 ) .OR.         ! Not previous
     &              RDT .NE. RDT2( I2 )      ) THEN  ! Not previous
 
                    I2 = I2 + 1
                    IF ( I2 .LE. MXSCNV ) THEN
 
                        FIP2( I2 ) = FIP
                        RDT2( I2 ) = RDT

                        CALL GRPFAC( VHT, EM, FAC, 
     &                               FAC2, MXSREF, NVTYPE, MXEPRC,
     &                               I2  , NVTYPE, NEMODE )

                    END IF
                END IF

            ELSE                                  !  FIP-RDT dependent profiles

                IF ( I1  .EQ. 0          .OR.        ! First instance
     &               FIP .NE. FIP1( I1 ) .OR.        ! Not previous
     &               RDT .NE. RDT1( I1 )      ) THEN ! Not previous
 
                    I1 = I1 + 1
                    IF ( I1 .LE. MXSCNV ) THEN
 
                        FIP1( I1 ) = FIP
                        RDT1( I1 ) = RDT

                        CALL GRPFAC( VHT, EM, FAC,
     &                               FAC1, MXFREF, NVTYPE, MXEPRC,
     &                               I1  , NVTYPE, NEMODE )

                    END IF
                END IF

            ENDIF

201     CONTINUE

        N1 = I1
        N2 = I2
        N3 = I3
        
        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xx

93010   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10 ( A, :, I5, :, 2X ) )

94020   FORMAT( A, I8, 2X, 10 ( A, :, E10.3, :, 1X ) )


        END

