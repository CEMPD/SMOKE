 
        SUBROUTINE RDACONV( KDEV, FIPST, ASC7S, ASC3S, FACST, NST,
     &                      ASC7N, ASC3N, FACNN, NNN )

C***********************************************************************
C  subroutine body starts at line  105
C
C  DESCRIPTION:
C       Reads the indicated volatile organic conversion file,
C	filters out point-source and mobile-source records,
C       sorts it, and returns the number of profile entries found.
C       In case of overflow, returns negative of number of 
C	area source profile entries found.
C
C  PRECONDITIONS REQUIRED:
C       volatile conversion file already opened on unit KDEV
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C	SORTI3()
C
C  REVISION  HISTORY:
C       Adapted 8/96 by CJC from point-source rdconv.F
C
C***********************************************************************
C PROJECT TITLE: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File Version @(#)$Id$
C Pathname:    $Source$
C Last updated: $Date$ 
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
C***************************************************************************/

      IMPLICIT NONE

C...........   INCLUDES:

      INCLUDE 'CHDIMS3.EXT'     ! emis chem info (both inventory and model)
      INCLUDE 'PARMS3.EXT'      ! I/O API constants
      INCLUDE 'FDESC3.EXT'      ! I/O API file description data structure
      INCLUDE 'IODECL3.EXT'     ! I/O API function declarations

C...........   ARGUMENTS and their descriptions:

        INTEGER     KDEV           !  unit number for speciation profile file
        INTEGER     NST             !  number of entries in group 1
        INTEGER     NNN             !  number of entries in group 2
        INTEGER     FIPST( MXSCNV ) ! sorted state codes in group 1
        INTEGER     ASC7S( MXSCNV ) ! sorted state-specific ASC7 codes 
        INTEGER     ASC3S( MXSCNV ) ! sorted state-specific ASC3 codes 
        REAL        FACST( MXSCNV ) ! sorted state-specific conversion factors
        INTEGER     ASC7N( MXSCNV ) ! sorted fallback ASC7 codes 
        INTEGER     ASC3N( MXSCNV ) ! sorted fallback ASC3 codes 
        REAL        FACNN( MXSCNV ) ! sorted fallback conversion factors

C...........   EXTERNAL FUNCTIONS and their descriptions:

        INTEGER     LBLANK
        INTEGER     STR2INT
        REAL        STR2REAL
        INTEGER     TRIMLEN

        EXTERNAL    LBLANK, STR2INT, STR2REAL, TRIMLEN

C...........   SCRATCH LOCAL VARIABLES and their descriptions:

        INTEGER         FIP              ! Temporary FIP code
        INTEGER         I, I1, I2, J     ! Counters
        INTEGER         IOS              ! I/O Status code
        INTEGER         IREC             ! Record number of input file
        INTEGER         NTOT             ! Total items read in
        INTEGER         ASC7, ASC3       ! Temporary ASC code
        INTEGER         SID              ! Temporary state ID
        INTEGER         COL, SIZ

        INTEGER         FIPA( MXSCNV )   ! unsorted FIPS
        INTEGER         INDX( MXSCNV )   ! index for sorting
        INTEGER         AS7A( MXSCNV )   ! unsorted SCCS
        INTEGER         AS3A( MXSCNV )   ! unsorted SCCS

        REAL            FAC              ! Temporary conversion factor
        REAL            FACA( MXSCNV )   ! unsorted conversion factors

        LOGICAL         EFLAG            ! Error flag

        CHARACTER*1     LETR
        CHARACTER*128   LINE
        CHARACTER*256   MESG

C***********************************************************************
C   begin body of function  RDACONV

C.........  Read volatile pollutants conversion factors file

        IREC  =  0
        I     =  0
        EFLAG = .FALSE.

11      CONTINUE        !  head of the KDEV-read loop

            IREC = IREC + 1
            READ( KDEV, 93010, END=99, IOSTAT=IOS ) LINE

            IF ( IOS .NE. 0 ) THEN

                EFLAG = .TRUE.
                WRITE( MESG, 94010 )
     &              'I/O error', IOS, 
     &              'reading VOLATILE CONVERSION line', IREC
                CALL M3MESG( MESG )
                GO TO  11		!  to head of read-loop

            END IF		!  if error during read

            IF ( INDEX( LINE( 1:2 ), '.' ) .EQ. 0 ) THEN
                SID = MAX( STR2INT( LINE( 1:2 ) ), 0 )
            ELSE
                SID = 0		!  "fallback" state placeholder
            END IF

            IF ( INDEX( LINE( 3 : 17 ), '.' ) .EQ. 0 ) THEN

                COL  = 3 + LBLANK( LINE( 3 : 17 ) )
                LETR = LINE( COL : COL )

                IF ( LETR .LT. 'A' .OR. LETR .GT. 'Z' ) THEN ! no MV sources

                    SIZ = TRIMLEN( LINE( COL : 17 ) )

                    IF ( SIZ .GE. 10 )  THEN ! exclude point sources

                        ASC7 = STR2INT( LINE( COL   : COL+6 ) )
                        ASC3 = STR2INT( LINE( COL+7 : COL+9 ) )

                    ELSE
                        GO TO 11
                    ENDIF
                ELSE
                    GO TO 11
                ENDIF

            ELSE		!  else "fallback" ASC

                ASC7 = 0
                ASC3 = 0

            END IF		!  if "sas-default" ASC value or not

            FAC = STR2REAL( LINE( 19:23 ) )
            IF( FAC .LT. 0.0 ) THEN

                EFLAG = .TRUE.
                WRITE( MESG,94020 )
     &            'Bad factor at line', IREC,
     &            'of VOLATILE CONVERSION FACTOR file =', FAC
                CALL M3MESG( MESG )
                GO TO 11  ! to head of read loop

            END IF		!  if fac invalid

            I = I + 1
            IF ( I .LE. MXSCNV ) THEN

                INDX( I ) = I
                FIPA( I ) = SID * 1000
                AS7A( I ) = ASC7
                AS3A( I ) = ASC3
                FACA( I ) = FAC

            END IF		!  if I in bounds

            GO TO  11   !  to head of loop

99      CONTINUE        !  end of the KDEV-read loop

        IF ( I .GT. MXSCNV ) THEN
            EFLAG = .TRUE.
            WRITE( MESG, 94010 )
     &          'Actual CONVERSION table size:', I,
     &          'Dimensioned (MXSCNV):', MXSCNV
            CALL M3MSG2( MESG )
        END IF		!  if I not in bounds

        IF ( EFLAG ) THEN
            CALL M3EXIT( 'RDACONV', 0, 0, 
     &           'Error reading VOLATILE CONVERSION file', 2 )

        ENDIF

C...........   Indirect sort on <FIPA, AS7A, AS3A>, then select 
C...........   into state-specific and state-independent tables.
C...........   (Overflow impossible, since both table-sizes are
C...........   the same as the unsorted input tables)

        NTOT = I
                
        CALL SORTI3( NTOT, INDX, FIPA, AS7A, AS3A )

        I1 =  0
        I2 =  0
        DO 101 I = 1, NTOT

            J    = INDX( I )
            FIP  = FIPA( J )
            ASC7 = AS7A( J )
            ASC3 = AS3A( J )
            FAC  = FACA( J )

            IF( FIP .GT. 0 ) THEN 

                I1 = I1 + 1
                IF ( I1 .LE. MXSCNV ) THEN
                    FIPST( I1 ) = FIP
                    ASC7S( I1 ) = ASC7
                    ASC3S( I1 ) = ASC3
                    FACST( I1 ) = FAC
                END IF

            ELSE

                I2 = I2 + 1
                IF ( I2 .LE. MXSCNV ) THEN
                    ASC7N( I2 ) = ASC7
                    ASC3N( I2 ) = ASC3
                    FACNN( I2 ) = FAC
                END IF

            END IF

101     CONTINUE

        IF ( I1 .LE. MXSCNV ) THEN
            NST =  I1
        ELSE
            NST = -I1
            WRITE( MESG,94010 )
     &          'FIP-dependent table overflow in RDACONV(); ' //
     &          'Dimensioned', MXSCNV,
     &          'actual',      I1
        END IF

        IF ( I2 .LE. MXSCNV ) THEN
            NNN =  I2
        ELSE
            NNN = -I2
            WRITE( MESG,94010 )
     &          'FIP-independent table overflow in RDACONV(); ' //
     &          'Dimensioned', MXSCNV,
     &          'actual',      I2
        END IF
        
        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xx

93010   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10 ( A, :, I5, :, 2X ) )

94020   FORMAT( A, I8, 2X, 10 ( A, :, E10.3, :, 1X ) )


        END

