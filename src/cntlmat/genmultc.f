
        SUBROUTINE GENMULTC( ADEV, CDEV, GDEV, LDEV, NCPE, PYEAR,
     &                       ENAME, MNAME, CFLAG, GFLAG, LFLAG, SFLAG )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      This subroutine computes the multiplicative control factors
C      and writes out the multiplicative control matrix.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     
C
C***************************************************************************
C
C Project Title: Sparse Matrix Operator Kernel gsions (SMOKE) Modeling
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
C***************************************************************************

C.........  MODULES for public variables
C.........  This module contains the inventory arrays
        USE MODSOURC

C.........  This module contains the control packet data and control matrices
        USE MODCNTRL

C.........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  i/o api parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.
        INCLUDE 'SETDECL.EXT'   !  FileSetAPI variables and functions
        INCLUDE 'FLTERR.EXT'    !  functions for comparing two numbers

C...........   EXTERNAL FUNCTIONS and their descriptions:
        LOGICAL   ENVYN
        INTEGER   GETEFILE
        INTEGER   INDEX1
        INTEGER   PROMPTFFILE
        REAL      YR2DAY

        EXTERNAL  ENVYN, GETEFILE, INDEX1, PROMPTFFILE, YR2DAY

C...........   SUBROUTINE ARGUMENTS

        INTEGER     , INTENT (IN) :: ADEV   ! file unit no. for tmp ADD file
        INTEGER     , INTENT (IN) :: CDEV   ! file unit no. for tmp CTL file 
        INTEGER     , INTENT (IN) :: GDEV   ! file unit no. for tmp CTG file
        INTEGER     , INTENT (IN) :: LDEV   ! file unit no. for tmp ALW file
        INTEGER     , INTENT (IN) :: NCPE   ! no. of control packet entries
        INTEGER     , INTENT (IN) :: PYEAR  ! projected year, or missing
        CHARACTER*16, INTENT (IN) :: ENAME  ! logical name for i/o api 
                                            ! inventory input file
        CHARACTER*16, INTENT (IN OUT) :: MNAME  ! logical name for mult. cntl. mat.
        LOGICAL     , INTENT (IN) :: CFLAG  ! true = apply CTL controls
        LOGICAL     , INTENT (IN) :: GFLAG  ! true = apply CTG controls
        LOGICAL     , INTENT (IN) :: LFLAG  ! true = apply ALW controls
        LOGICAL     , INTENT (IN) :: SFLAG  ! true = apply EMS_CTL controls

C...........   Local allocatable arrays
        INTEGER, ALLOCATABLE :: ALWINDX ( :,: ) ! indices to ALW controls table
        INTEGER, ALLOCATABLE :: CTGINDX ( :,: ) ! indices to CTG controls table
        INTEGER, ALLOCATABLE :: CTLINDX ( :,: ) ! indices to CTL controls table
        INTEGER, ALLOCATABLE :: GRPINDX ( : )   ! index from sources to groups
        INTEGER, ALLOCATABLE :: GRPSTIDX( : )   ! sorting index

        REAL   , ALLOCATABLE :: BACKOUT ( : )   ! factor used to account for pol
                                                ! specific control info that is
                                                ! already in the inventory
        REAL   , ALLOCATABLE :: CTLEFF  ( : )   ! control efficiency
        REAL   , ALLOCATABLE :: EMIS    ( : )   ! base inventory emissions
        REAL   , ALLOCATABLE :: FACTOR  ( : )   ! multiplicative controls
        REAL   , ALLOCATABLE :: RULEFF  ( : )   ! rule effectiveness
        REAL   , ALLOCATABLE :: RULPEN  ( : )   ! rule penetration
        REAL   , ALLOCATABLE :: GRPINEM ( :,: ) ! initial emissions
        REAL   , ALLOCATABLE :: GRPOUTEM( :,: ) ! controlled emissions

        LOGICAL, ALLOCATABLE :: GRPFLAG ( : )   ! true: group controlled

        CHARACTER(LEN=STALEN3+SCCLEN3), ALLOCATABLE :: GRPCHAR( : ) ! group chars

C.........   Local arrays
        INTEGER                 OUTTYPES( NVCMULT,6 ) ! var type:int/real

        CHARACTER(LEN=IOVLEN3)  OUTNAMES( NVCMULT,6 ) ! var names
        CHARACTER(LEN=IOULEN3)  OUTUNITS( NVCMULT,6 ) ! var units
        CHARACTER(LEN=IODLEN3)  OUTDESCS( NVCMULT,6 ) ! var descriptions

C...........   Other local variables
        INTEGER          I, J, K, S  ! counters and indices

        INTEGER          IDX      ! group index
        INTEGER          IOS      ! input/output status
        INTEGER          NGRP     ! number of reporting groups
        INTEGER          PIDX     ! previous IDX
        INTEGER          RDEV     ! Report unit number
        INTEGER          SCCBEG   ! begining of SCC in CSOURC string
        INTEGER          SCCEND   ! end of SCC in CSOURC string

        REAL             ALWFAC   ! allowable control factor
        REAL             ALWEMIS  ! allowable emissions
        REAL             CAP      ! emissions cap
        REAL             CTGFAC   ! control technology control factor
        REAL             CTGFAC2  ! MAXACF or RSNACF
        REAL             CUTOFF   ! CTG cutoff for application of control
        REAL             DENOM    ! denominator of control back-out factor
        REAL             E1, E2   ! tmp emissions values
        REAL             FAC      ! tmp control factor
        REAL             MACT     ! max. achievable cntrl tech. cntrl factor
        REAL             RACT     ! reasonably achiev. cntrl tech. cntrl factor
        REAL             REPLACE  ! replacement emissions

        LOGICAL          LO3SEAS  ! true: use ozone-season emissions
        LOGICAL, SAVE :: APPLFLAG = .FALSE. ! true: something has been applied
        LOGICAL, SAVE :: OPENFLAG = .FALSE. ! true: output file has been opened

        CHARACTER*100          OUTFMT     ! header format buffer
        CHARACTER*300          MESG       ! message buffer
        CHARACTER(LEN=FPLLEN3) CPLT       ! tmp point src info through plant
        CHARACTER(LEN=FPLLEN3) PPLT       ! previous CPLT
        CHARACTER(LEN=STALEN3) CSTA       ! tmp char state
        CHARACTER(LEN=STALEN3) PSTA       ! previous char state
        CHARACTER(LEN=SCCLEN3) TSCC       ! tmp SCC
        CHARACTER(LEN=SCCLEN3) PSCC       ! previous SCC
        CHARACTER(LEN=SRCLEN3) CSRC       ! tmp source chars
        CHARACTER(LEN=IOVLEN3) PNAM       ! tmp pollutant name

        CHARACTER*16  :: PROGNAME = 'GENMULTC' ! program name

C***********************************************************************
C   begin body of subroutine GENMULTC

C.........  Get environment variables that control program behavior
        MESG = 'Use annual or ozone season emissions'
        LO3SEAS = ENVYN( 'SMK_O3SEASON_YN', MESG, .FALSE., IOS )

C.........  Open reports file
        RPTDEV( 1 ) = PROMPTFFILE( 
     &                'Enter logical name for MULTIPLICATIVE ' //
     &                'CONTROLS REPORT',
     &                .FALSE., .TRUE., CRL // 'CREP', PROGNAME )
        RDEV = RPTDEV( 1 )

C.........  Allocate index to reporting groups
        ALLOCATE( GRPINDX( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'GRPINDX', PROGNAME )
        ALLOCATE( GRPSTIDX( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'GRPSTIDX', PROGNAME )
        ALLOCATE( GRPCHAR( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'GRPCHAR', PROGNAME )
        GRPINDX  = 0  ! array

C.........  Get set up for group reporting...
        IF( CATEGORY .EQ. 'POINT' ) THEN

C.............  Count the number of groups in the inventory
            PPLT = ' '
            NGRP = 0
            DO S = 1, NSRC
                CPLT = CSOURC( S )( 1:FPLLEN3 )
                IF( CPLT .NE. PPLT ) THEN
                    NGRP = NGRP + 1
                    PPLT = CPLT
                END IF
                GRPINDX ( S ) = NGRP
                GRPSTIDX( S ) = S     ! Needed for loops, but not used to sort
            END DO

        ELSE 

            IF( CATEGORY .EQ. 'AREA' ) THEN
                SCCBEG = ARBEGL3( 2 )
                SCCEND = ARENDL3( 2 )
            ELSE            ! MOBILE
                SCCBEG = MBBEGL3( 4 )
                SCCEND = MBENDL3( 4 )
            END IF

C.............  Build and sort source array for SCC-state grouping
            DO S = 1, NSRC
                CSTA = CSOURC( S )( 1     :STALEN3 )
                TSCC = CSOURC( S )( SCCBEG:SCCEND  )

                GRPSTIDX( S ) = S  
                GRPCHAR ( S ) = CSTA // TSCC
            END DO

            CALL SORTIC( NSRC, GRPSTIDX, GRPCHAR )

C.............  Count the number of state/SCCs in the domain
            PSTA = ' '
            PSCC = ' '
            SCCBEG = STALEN3 + 1
            SCCEND = STALEN3 + SCCLEN3
            DO S = 1, NSRC
                J = GRPSTIDX( S )
                CSTA = GRPCHAR( J )( 1     :STALEN3 )
                TSCC = GRPCHAR( J )( SCCBEG:SCCEND  )
                IF( CSTA .NE. PSTA .OR. TSCC .NE. PSCC ) THEN
                    NGRP = NGRP + 1
                    PSTA = CSTA
                    PSCC = TSCC
                END IF
                GRPINDX( J ) = NGRP
            END DO

        END IF

C...........  Allocate memory for the number of groups for storing emissions
          ALLOCATE( GRPFLAG( NGRP ), STAT=IOS )
          CALL CHECKMEM( IOS, 'GRPFLAG', PROGNAME )
          ALLOCATE( GRPINEM( NGRP, NVCMULT ), STAT=IOS )
          CALL CHECKMEM( IOS, 'GRPINEM', PROGNAME )
          ALLOCATE( GRPOUTEM( NGRP, NVCMULT ), STAT=IOS )
          CALL CHECKMEM( IOS, 'GRPOUTEM', PROGNAME )

C...........  Initialize
          GRPINEM  = 0. ! array
          GRPOUTEM = 0. ! array
          GRPFLAG  = .FALSE.  ! array

C...........  Allocate local memory
          ALLOCATE( ALWINDX( NSRC, NVCMULT ), STAT=IOS )
          CALL CHECKMEM( IOS, 'ALWINDX', PROGNAME )
          ALLOCATE( CTGINDX( NSRC, NVCMULT ), STAT=IOS )
          CALL CHECKMEM( IOS, 'CTGINDX', PROGNAME )
          ALLOCATE( CTLINDX( NSRC, NVCMULT ), STAT=IOS )
          CALL CHECKMEM( IOS, 'CTLINDX', PROGNAME )
          ALLOCATE( BACKOUT( NSRC ), STAT=IOS )
          CALL CHECKMEM( IOS, 'BACKOUT', PROGNAME )
          ALLOCATE( CTLEFF( NSRC ), STAT=IOS )
          CALL CHECKMEM( IOS, 'CTLEFF', PROGNAME )
          ALLOCATE( EMIS( NSRC ), STAT=IOS )
          CALL CHECKMEM( IOS, 'EMIS', PROGNAME )
          ALLOCATE( FACTOR( NSRC ), STAT=IOS )
          CALL CHECKMEM( IOS, 'FACTOR', PROGNAME )
          ALLOCATE( RULEFF( NSRC ), STAT=IOS )
          CALL CHECKMEM( IOS, 'RULEFF', PROGNAME )
          ALLOCATE( RULPEN( NSRC ), STAT=IOS )
          CALL CHECKMEM( IOS, 'RULPEN', PROGNAME )

C...........  For each pollutant that receives controls, obtain variable
C             names for control efficiency, rule effectiveness, and, in the
C             case of AREA sources, rule penetration. These variable names
C             will be used in reading the inventory file.

        CALL BLDENAMS( CATEGORY, NVCMULT, 6, PNAMMULT, OUTNAMES,
     &                 OUTUNITS, OUTTYPES, OUTDESCS )

C...........  Read in indices from temporary files. No error checking is
C             performed because it is assumed that the program has already
C             successfully written the temporary files.

        DO I = 1, NVCMULT

           IF( PCTLFLAG( I, 1 ) ) THEN
              DO S = 1, NSRC
                 READ(CDEV,*) CTLINDX( S, I )
              END DO
           END IF

           IF( PCTLFLAG( I, 2 ) ) THEN
              DO S = 1, NSRC
                READ(GDEV,*) CTGINDX( S, I )
              END DO
           END IF

           IF( PCTLFLAG( I, 3 ) ) THEN
              DO S = 1, NSRC
                 READ(LDEV,*) ALWINDX( S, I )
              END DO
           END IF

        END DO ! end pollutant loop

        REWIND( CDEV )
        REWIND( GDEV )
        REWIND( LDEV )

C...........  Fractionalize control-packet information

        IF ( CFLAG ) THEN
            FACCEFF = FACCEFF*0.01  ! array
            FACREFF = FACREFF*0.01  ! array
            FACRLPN = FACRLPN*0.01  ! array
        END IF

        IF( SFLAG ) THEN
            BASCEFF = BASCEFF*0.01  ! array
            BASREFF = BASREFF*0.01  ! array
            BASRLPN = BASRLPN*0.01  ! array
            EMSCEFF = EMSCEFF*0.01  ! array
            EMSREFF = EMSREFF*0.01  ! array
            EMSRLPN = EMSRLPN*0.01  ! array
        END IF

C...........  Loop through pollutants that receive controls

        DO I = 1, NVCMULT

C............  Set tmp pollutant name 
           PNAM = PNAMMULT( I )

C............  Initialize control factor array
           FACTOR = 1.0  ! array

C...........  Read in emissions data from inventory file
C...........  Ozone-season emissions
           IF( LO3SEAS ) THEN

               IF ( .NOT. 
     &              READSET( ENAME, OUTNAMES(I,2), 1, -1, 
     &                       0, 0, EMIS ) ) THEN
                   CALL WRITE_MESG_EXIT( OUTNAMES(I,1), PROGNAME )
               END IF

C...........  Annual emissions
           ELSE 

               IF ( .NOT. 
     &              READSET( ENAME, OUTNAMES(I,1), 1, -1,
     &                       0, 0, EMIS )  ) THEN
                   CALL WRITE_MESG_EXIT( OUTNAMES(I,1), PROGNAME )
               END IF

C...............  Divide annual emissions to get average day
              FAC = YR2DAY( BYEAR )
              EMIS = EMIS * FAC      ! array

            END IF

            DO S = 1, NSRC

                IF ( EMIS( S ) .LT. AMISS3 ) THEN
                    EMIS( S ) = 0.0
                END IF

C.................  Compute group emissions before controls
                J = GRPINDX( S )
                GRPINEM ( J,I ) = GRPINEM ( J,I ) + EMIS( S )

            END DO ! end source loop

C...........  If CONTROL packet is present: For the current pollutant, read
C             in control efficiency, rule effectiveness, and, in the case of 
C             AREA sources, rule penetration.
           IF ( CFLAG ) THEN

              SELECT CASE( CATEGORY )

C.............  Area sources...
              CASE( 'AREA' )

              IF ( .NOT. READSET( ENAME, OUTNAMES(I,4), 1, -1,
     &                          0, 0, CTLEFF ) ) THEN
                 CALL WRITE_MESG_EXIT( OUTNAMES(I,4), PROGNAME )
              END IF

              IF ( .NOT. READSET( ENAME, OUTNAMES(I,5), 1, -1,
     &                          0, 0, RULEFF ) ) THEN
                 CALL WRITE_MESG_EXIT( OUTNAMES(I,5), PROGNAME )
              END IF

              IF ( .NOT. READSET( ENAME, OUTNAMES(I,6), 1, -1,
     &                          0, 0, RULPEN ) ) THEN
                 CALL WRITE_MESG_EXIT( OUTNAMES(I,6), PROGNAME )
              END IF

C.............  Mobile sources...
              CASE( 'MOBILE' ) 

              CTLEFF = 0.   ! array
              RULEFF = 0.   ! array
              RULPEN = 0.   ! array

C.............  Point sources...
              CASE( 'POINT' )

              IF ( .NOT. READSET( ENAME, OUTNAMES(I,3), 1, -1,
     &                          0, 0, CTLEFF ) ) THEN
                 CALL WRITE_MESG_EXIT( OUTNAMES(I,3), PROGNAME )
              END IF

              IF ( .NOT. READSET( ENAME, OUTNAMES(I,4), 1, -1,
     &                          0, 0, RULEFF ) ) THEN
                 CALL WRITE_MESG_EXIT( OUTNAMES(I,4), PROGNAME )
              END IF

              RULPEN = 100.0  ! array

              CASE DEFAULT
                  MESG = 'Case ' // CATEGORY // ' not supported.'
                  CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

              END SELECT      ! end select on category

C...........  Fractionalize all inventory control information and set
C             missing values to one, then calculate the factor which will
C             be used to account for control information already in the
C             inventory

              DO S = 1, NSRC

                 IF ( CTLEFF( S ) .LT. AMISS3 ) THEN
                    CTLEFF( S ) = 0.0
                 ELSE
                    CTLEFF( S ) = CTLEFF( S )*0.01
                 END IF

                 IF ( RULEFF( S ) .LT. AMISS3 ) THEN
                    RULEFF( S ) = 1.0
                 ELSE
                    RULEFF( S ) = RULEFF( S )*0.01
                 END IF

                 IF ( RULPEN( S ) .LT. AMISS3 ) THEN
                    RULPEN( S ) = 1.0
                 ELSE
                    RULPEN( S ) = RULPEN( S )*0.01
                 END IF

C..................  Perform division by zero check. 
                 DENOM = ( 1.0 - CTLEFF(S)*RULEFF(S)*RULPEN(S) )

                 IF ( FLTERR( DENOM, 0.0 ) ) THEN
                    BACKOUT( S ) = 1.0/DENOM
                 ELSE
                    BACKOUT( S ) = 0.0
                 END IF

              END DO ! end source loop

C...........  If EMS-95 CONTROL packet is present:
           ELSE IF ( SFLAG ) THEN

C..................  Perform division by zero check. 
              DO J = 1, NCPE

                 DENOM = ( 1.0 - BASCEFF(J)*BASREFF(J)*BASRLPN(J) )

                 IF ( FLTERR( DENOM, 0.0 ) ) THEN
                    BACKOUT( J ) = 1.0/DENOM
                 ELSE
                    BACKOUT( J ) = 0.0
                 END IF

              END DO ! end source loop

           END IF

C.............................................................................
C...........  Apply /CONTROL/ packet controls if present for the current 
C             pollutant
C.............................................................................
           IF ( CFLAG .AND. PCTLFLAG( I, 1 ) ) THEN

C...............  Loop through sources
              DO S = 1, NSRC

                 E1 = EMIS( S )

C..................  Get index to /CONTROL/ packet from tmp file and compute
C                    control factor
                 K = CTLINDX( S, I )
                 IF ( K .GT. 0 ) THEN
                    CTLEFF( S ) = FACCEFF( K )
                    RULEFF( S ) = FACREFF( K )
                    RULPEN( S ) = FACRLPN( K )
                    FACTOR( S ) = BACKOUT( S )*
     &                          ( 1.0 - CTLEFF(S)*RULEFF(S)*RULPEN(S) )

C.....................  Overwrite temporary file line with new info
                    E2 = E1 * FACTOR( S )
                    WRITE( CDEV,93300 ) 1, PNAM, E1, E2, FACTOR( S )
                    APPLFLAG = .TRUE.

C..................  Overwrite temporary file line with new info
                 ELSE
                    E2 = E1
                    WRITE( CDEV,93300 ) 0, 'D', 0., 0., 0.

                 END IF

              END DO ! end source loop

           END IF

C.............................................................................
C...........  Compute CTL factor using EMS-95 inputs
C.............................................................................
C...........  NOTE - SFLAG and CFLAG cannot both be true
           IF ( SFLAG .AND. PCTLFLAG( I, 1 ) ) THEN

              DO S = 1, NSRC

                 E1 = EMIS( S )

                 K = CTLINDX( S, I )
                 IF ( K .GT. 0 ) THEN
                    IF( EMSTOTL( K ) .NE. 0. ) THEN
                       FACTOR( S ) = EMSTOTL( K )

                    ELSE
                       CTLEFF( S ) = EMSCEFF( K )
                       RULEFF( S ) = EMSREFF( K )
                       RULPEN( S ) = EMSRLPN( K )
                       FACTOR( S ) = BACKOUT( K ) * EMSPTCF( K ) * 
     &                               (1.- CTLEFF(S)*RULEFF(S)*RULPEN(S))

                    END IF

C.....................  Overwrite temporary file line with new info
                    E2 = E1 * FACTOR( S )
                    WRITE( CDEV,93300 ) 1, PNAM, E1, E2, FACTOR( S )
                    APPLFLAG = .TRUE.

C..................  Overwrite temporary file line with new info
                 ELSE
                    E2 = E1
                    WRITE( CDEV,93300 ) 0, 'D', 0., 0., 0.

                 END IF

              END DO ! end source loop

           END IF

C.............................................................................
C............  Apply /CTG/ packet
C.............................................................................
           IF ( GFLAG .AND. PCTLFLAG( I, 2 ) ) THEN

C...............  Compute CTG factor
              DO S = 1, NSRC

                 E1  = EMIS( S ) * FACTOR( S )
                 FAC = 1.

                 K = CTGINDX( S, I ) 
                 IF ( K .GT. 0 ) THEN
                    CUTOFF = CUTCTG ( K )
                    CTGFAC = FACCTG ( K )
                    MACT   = FACMACT( K )
                    RACT   = FACRACT( K )

C.....................  Check to see if emissions exceed cutoff and if 
C                       necessary, apply controls
C.....................  The comparison emission value already has controls
C                       from the /CONTROL/ packet
                    IF ( E1 .GT. CUTOFF ) THEN

C........................  Initialize CTG factor with base factor from packet
                       FAC = CTGFAC

C........................  Compute output emissions with /CONTROL/ and base CTG 
                       E2  = E1*CTGFAC

C........................  If emissions still exceed cutoff, apply second
C                          CTG factor
                       IF ( E2 .GT. CUTOFF ) THEN

C...........................  Use MACT factor if it is defined 
                          IF ( MACT .GT. 0 ) THEN
                             CTGFAC2 = MACT

C...........................  Otherwise, use RACT factor 
                          ELSE IF ( RACT .GT. 0 ) THEN
                             CTGFAC2 = RACT

C...........................  Otherwise, set to cutoff value 
                          ELSE
                             CTGFAC2 = CUTOFF / E2

                          END IF

                          FAC = FAC * CTGFAC2
                          E2  = E2  * CTGFAC2

                       END IF

C........................  Compute aggregate factor for current source
                       FACTOR( S ) = FACTOR( S ) * FAC 

C........................  Overwrite temporary file line with new info
                       WRITE( GDEV,93300 ) 1, PNAM, E1, E2, FAC
                       APPLFLAG = .TRUE.

C.....................  If no controls, then overwrite temporary line only
                    ELSE
                       WRITE( GDEV,93300 ) 0, 'D', 0., 0., 0.

                    END IF

                 ELSE ! If K = 0 (for sources not having applied "Controls"

                     WRITE( GDEV,93300 ) 0, 'D', 0., 0., 0.

                 END IF

              END DO ! end source loop

           END IF

           IF ( LFLAG .AND. PCTLFLAG( I, 3 ) ) THEN

C...........  Process ALW packet

              DO S = 1, NSRC

                 E1 = EMIS( S )

                 K = ALWINDX( S, I ) 
                 IF ( K .GT. 0 ) THEN
                    ALWFAC  = FACALW  ( K )
                    CAP     = EMCAPALW( K )
                    REPLACE = EMREPALW( K )

C.....................  Both Cap value and Replace value are defined, then
C                       compare emissions to Cap and set factor with Replace.
                    IF ( CAP .GE. 0 .AND. REPLACE .GE. 0 ) THEN

                       IF ( EMIS( S ) .GT. CAP ) THEN
                          FACTOR( S ) = REPLACE/EMIS( S )
                       END IF

C.....................  Only Cap value is defined, then compare emissions to Cap
C                       set factor with Cap
                    ELSE IF ( CAP .GE. 0 .AND. REPLACE .LT. 0 ) THEN

                       IF ( EMIS( S ) .GT. CAP ) THEN
                          FACTOR( S ) = CAP/EMIS( S )
                       END IF

C.....................  Only Replace value is defined, then set factor with
C                       Replace.
                    ELSE IF ( CAP .LT. 0 .AND. REPLACE .GE. 0 ) THEN

                       IF ( EMIS( S ) .GT. REPLACE ) THEN
                          FACTOR( S ) = REPLACE/EMIS( S )
                       END IF
                       
                    END IF

C.....................  Overwrite temporary file line with new info
                    E2 = E1 * FACTOR( S )
                    WRITE( LDEV,93300 ) 1, PNAM, E1, E2, FACTOR( S )
                    APPLFLAG = .TRUE.

C..................  If no controls, then overwrite temporary line only
                 ELSE
                     WRITE( LDEV,93300 ) 0, 'D', 0., 0., 0.

                 END IF

                END DO ! end source loop

            END IF

C.............  Store output emissions for groups
C.............  This must be in a separate loop to account for all possible
C               combinations of packets
            DO S = 1, NSRC

                E1 = EMIS( S )
                E2 = E1 * FACTOR( S )
                J  = GRPINDX( S )
                GRPOUTEM( J,I ) = GRPOUTEM( J,I ) + E2

C.................  Flag group if emissions are different
                IF( E1 .NE. E2 ) GRPFLAG( J ) = .TRUE.

            END DO

C.............  Open control matrix, if needed and if not opened before
            IF( APPLFLAG .AND. .NOT. OPENFLAG ) THEN 
                CALL OPENCMAT( ENAME, 'MULTIPLICATIVE', MNAME )
                OPENFLAG = .TRUE.
            END IF

C.............  Write multiplicative controls for current pollutant
            IF( OPENFLAG ) THEN
                IF( .NOT. WRITE3( MNAME, PNAM, 0, 0, FACTOR ) ) THEN
                    MESG = 'Failed to write multiplicative control ' // 
     &                     'factors for pollutant ' // PNAM
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF
            END IF

        END DO ! end pollutant loop

C.........  Write out controlled facilities report for point sources
        IF( APPLFLAG ) THEN
            WRITE( RDEV, 93000 ) 'Processed as '// CATDESC// ' sources'
            IF ( PYEAR .GT. 0 ) THEN
                WRITE( RDEV, 93390 ) 'Projected inventory year ', PYEAR
            ELSE
                WRITE( RDEV, 93390 ) 'Base inventory year ', BYEAR
            END IF
            IF( CFLAG ) WRITE( RDEV, 93000 ) 
     &                  'Controls applied with /CONTROL/ packet'
            IF( SFLAG ) WRITE( RDEV, 93000 ) 
     &                  'Controls applied with /EMS_CONTROL/ packet'
            IF( GFLAG ) WRITE( RDEV, 93000 ) 
     &                  'Controls applied with /CTG/ packet'
            IF( LFLAG ) WRITE( RDEV, 93000 ) 
     &                  'Controls applied with /ALLOWABLE/ packet'
            IF( LO3SEAS ) THEN
                WRITE( RDEV,93000 ) 'Ozone-season data basis in report'
            ELSE
                WRITE( RDEV,93000 ) 'Annual total data basis in report'
            END IF

            IF ( CATEGORY .EQ. 'POINT' ) THEN
                WRITE( RDEV, 93000 ) 
     &             'Emissions by controlled facility before and ' //
     &             'after controls'
            ELSE
                WRITE( RDEV, 93000 )
     &             'Emissions by controlled state/SCC before and ' //
     &             'after controls'
            END IF

            IF ( CATEGORY .EQ. 'POINT' ) THEN
                WRITE( RDEV, 93400 ) 
     &               ( PNAMMULT( I ), PNAMMULT( I ), I=1,NVCMULT )
            ELSE
                WRITE( RDEV, 93402 ) 
     &               ( PNAMMULT( I ), PNAMMULT( I ), I=1,NVCMULT )
            END IF

            PNAM = '[tons/day]'
            WRITE( RDEV, 93405 ) ( PNAM, PNAM, I=1,NVCMULT )
            
            I = 26 + NVCMULT * 43
            WRITE( RDEV,93000 ) REPEAT( '-', I )
 
            PIDX = 0
            DO S = 1, NSRC
                K   = GRPSTIDX( S )
                IDX = GRPINDX ( K )

                IF( IDX .NE. PIDX .AND.
     &              GRPFLAG( IDX )      ) THEN

                    J = FIPLEN3 + 1
                    SELECT CASE( CATEGORY )
                    CASE ( 'AREA' , 'MOBILE' )
                        CSRC = GRPCHAR( S )
                        WRITE( RDEV, 93410 ) 
     &                      CSRC( 1:STALEN3 ), CSRC( SCCBEG:SCCEND ),
     &                      ( GRPINEM ( IDX,I ), GRPOUTEM( IDX,I ),
     &                        I = 1, NVCMULT )
                    CASE ( 'POINT' )                 
                        CSRC = CSOURC( S )
                        WRITE( RDEV, 93412 ) 
     &                      CSRC( 1:FIPLEN3 ), CSRC( J:FPLLEN3 ),
     &                      ( GRPINEM ( IDX,I ), GRPOUTEM( IDX,I ),
     &                        I = 1, NVCMULT )
                    END SELECT

                    PIDX = IDX
                END IF
            END DO

        END IF

        IF( .NOT. APPLFLAG ) THEN

            MESG = 'WARNING: No CONTROL, EMS_CONTROL, CTG, or ' //
     &             'ALLOWABLE packet entries match inventory.'
            CALL M3MSG2( MESG )

            MESG = 'WARNING: Multiplicative control will not be output!'
            CALL M3MSG2( MESG )

            WRITE( RDEV, 93000 ) 'No CONTROL, EMS_CONTROL, CTG, or ' //
     &             'ALLOWABLE packet entries matched inventory.'

        END IF

C.........  Deallocate local memory
        DEALLOCATE( ALWINDX, CTGINDX, CTLINDX, BACKOUT, CTLEFF,
     &              EMIS, FACTOR, RULEFF, RULPEN )

        DEALLOCATE( GRPINDX, GRPFLAG, GRPINEM, GRPOUTEM )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

93300   FORMAT( I2, 1X, '"', A, '"', 3( 1X, E12.5 ) )

93390   FORMAT( A, I4.4 )

93400   FORMAT( ' Region;', 5X, 'Facility ID;', 1X, 
     &          100( '  In ', A16, ';', 1X, 'Out ', A16, :, ';' ) )

93402   FORMAT( '  State;', 5X, '        SCC;', 1X, 
     &          100( '  In ', A16, ';', 1X, 'Out ', A16, :, ';' ) )

93405   FORMAT( 7X, ';', 16X, ';', 1X
     &          100( 2X, A16, 3X, ';', 1X, A16, :, 4X, ';' ))

93410   FORMAT( 4X, A3, ';', 6X, A10, ';', 1X, 
     &          100( 10X, E11.4, :, ';' ))

93412   FORMAT( 1X, A6, ';', 1X, A15, ';', 1X, 
     &          100( 10X, E11.4, :, ';' ))

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

C******************  INTERNAL SUBPROGRAMS  *****************************

        CONTAINS

C.............  This internal subprogram writes and error message
C               and then terminates program execution if an error
C               is encountered reading control information from the
C               inventory file
            SUBROUTINE WRITE_MESG_EXIT( OUTNAME, PROGNAME )

C.............  Subprogram arguments
            CHARACTER*(*), INTENT (IN) :: OUTNAME   ! name of inventory
                                                    ! variable that generated
                                                    ! the error
            CHARACTER*16,  INTENT (IN) :: PROGNAME  ! name of calling subroutine

C.............  Local variables
            CHARACTER* 300   MESG                   ! message buffer

C----------------------------------------------------------------------

            MESG = 'Error reading ' // OUTNAME // ' from inventory file'

            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            END SUBROUTINE WRITE_MESG_EXIT

C----------------------------------------------------------------------

        END SUBROUTINE GENMULTC
