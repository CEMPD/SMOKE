
        SUBROUTINE GENMULTC( ADEV, CDEV, GDEV, LDEV, RDEV, NCPE, ENAME, 
     &                       MNAME, CFLAG, GFLAG, LFLAG, SFLAG )

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
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2000, MCNC--North Carolina Supercomputing Center
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
        INCLUDE 'FLTERR.EXT'    !  functions for comparing two numbers

C...........   EXTERNAL FUNCTIONS and their descriptions:
        LOGICAL   ENVYN
        INTEGER   GETEFILE

        EXTERNAL  ENVYN, GETEFILE

C...........   SUBROUTINE ARGUMENTS

        INTEGER     , INTENT (IN) :: ADEV   ! file unit no. for tmp ADD file
        INTEGER     , INTENT (IN) :: CDEV   ! file unit no. for tmp CTL file 
        INTEGER     , INTENT (IN) :: GDEV   ! file unit no. for tmp CTG file
        INTEGER     , INTENT (IN) :: LDEV   ! file unit no. for tmp ALW file
        INTEGER     , INTENT (IN) :: RDEV   ! file unit no. for output report
        INTEGER     , INTENT (IN) :: NCPE   ! no. of control packet entries
        CHARACTER*16, INTENT (IN) :: ENAME  ! logical name for i/o api 
                                            ! inventory input file
        CHARACTER*16, INTENT (IN) :: MNAME  ! logical name for mult. cntl. mat.
        LOGICAL     , INTENT (IN) :: CFLAG  ! true = apply CTL controls
        LOGICAL     , INTENT (IN) :: GFLAG  ! true = apply CTG controls
        LOGICAL     , INTENT (IN) :: LFLAG  ! true = apply ALW controls
        LOGICAL     , INTENT (IN) :: SFLAG  ! true = apply EMS_CTL controls

C...........   Local allocatable arrays
        INTEGER, ALLOCATABLE :: ALWINDX ( :,: ) ! indices to ALW controls table
        INTEGER, ALLOCATABLE :: CTGINDX ( :,: ) ! indices to CTG controls table
        INTEGER, ALLOCATABLE :: CTLINDX ( :,: ) ! indices to CTL controls table
        INTEGER, ALLOCATABLE :: PLTINDX ( : )   ! index from sources to plants

        REAL   , ALLOCATABLE :: BACKOUT ( : )   ! factor used to account for pol
                                                ! specific control info that is
                                                ! already in the inventory
        REAL   , ALLOCATABLE :: CTLEFF  ( : )   ! control efficiency
        REAL   , ALLOCATABLE :: EMIS    ( : )   ! base inventory emissions
        REAL   , ALLOCATABLE :: FACTOR  ( : )   ! multiplicative controls
        REAL   , ALLOCATABLE :: RULEFF  ( : )   ! rule effectiveness
        REAL   , ALLOCATABLE :: RULPEN  ( : )   ! rule penetration
        REAL   , ALLOCATABLE :: PLTINEM ( :,: ) ! initial emissions
        REAL   , ALLOCATABLE :: PLTOUTEM( :,: ) ! controlled emissions

        LOGICAL, ALLOCATABLE :: PLTFLAG ( : )   ! true: plant controlled

C.........  Local arrays
        INTEGER                 OUTTYPES( NVCMULT,6 ) ! var type:int/real

        CHARACTER(LEN=IOVLEN3)  OUTNAMES( NVCMULT,6 ) ! var names
        CHARACTER(LEN=IOULEN3)  OUTUNITS( NVCMULT,6 ) ! var units
        CHARACTER(LEN=IODLEN3)  OUTDESCS( NVCMULT,6 ) ! var descriptions

C...........   Other local variables
        INTEGER          I, J, K, S  ! counters and indices

        INTEGER          IDX      ! plant index
        INTEGER          IOS      ! input/output status
        INTEGER          NPLT     ! number of plants
        INTEGER          PIDX     ! previous IDX

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

        CHARACTER*100          OUTFMT     ! header format buffer
        CHARACTER*300          MESG       ! message buffer
        CHARACTER(LEN=FPLLEN3) CPLT       ! tmp point src info through plant
        CHARACTER(LEN=FPLLEN3) PPLT       ! previous CPLT
        CHARACTER(LEN=SRCLEN3) CSRC       ! tmp source chars
        CHARACTER(LEN=IOVLEN3) PNAM       ! tmp pollutant name

        CHARACTER*16  :: PROGNAME = 'GENMULTC' ! program name

C***********************************************************************
C   begin body of subroutine GENMULTC

C...........  Get environment variables that control program behavior
          MESG = 'Use annual or ozone season emissions'
          LO3SEAS = ENVYN( 'SMK_O3SEASON_YN', MESG, .FALSE., IOS )

C...........  Get set up for reporting...
          IF( CATEGORY .EQ. 'POINT' ) THEN

C...............  Count the number of plants in the domain
              PPLT = ' '
              NPLT = 0
              DO S = 1, NSRC
                  CPLT = CSOURC( S )( 1:FPLLEN3 )
                  IF( CPLT .NE. PPLT ) NPLT = NPLT + 1
                  PPLT = CPLT
              END DO

C...............  Allocate memory for the number of plants for storing emissions
              ALLOCATE( PLTINDX( NSRC ), STAT=IOS )
              CALL CHECKMEM( IOS, 'PLTINDX', PROGNAME )
              ALLOCATE( PLTFLAG( NPLT ), STAT=IOS )
              CALL CHECKMEM( IOS, 'PLTFLAG', PROGNAME )
              ALLOCATE( PLTINEM( NPLT, NVCMULT ), STAT=IOS )
              CALL CHECKMEM( IOS, 'PLTINEM', PROGNAME )
              ALLOCATE( PLTOUTEM( NPLT, NVCMULT ), STAT=IOS )
              CALL CHECKMEM( IOS, 'PLTOUTEM', PROGNAME )
 
C...............  Initialize
              PLTINDX  = 0  ! array
              PLTINEM  = 0. ! array
              PLTOUTEM = 0. ! array
              PLTFLAG  = .FALSE.  ! array

C...............  Create index from sources to plants
              PPLT = ' '
              NPLT = 0
              DO S = 1, NSRC
                  CPLT = CSOURC( S )( 1:FPLLEN3 )
                  IF( CPLT .NE. PPLT ) THEN
                      NPLT = NPLT + 1
                      PPLT = CPLT
                  END IF
                  PLTINDX( S ) = NPLT
              END DO

          END IF

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

           IF( LO3SEAS ) THEN

               IF ( .NOT. 
     &              READ3( ENAME, OUTNAMES(I,2), 1, 0, 0, EMIS )  ) THEN
                   CALL WRITE_MESG_EXIT( OUTNAMES(I,1), PROGNAME )
               END IF

           ELSE 

               IF ( .NOT. 
     &              READ3( ENAME, OUTNAMES(I,1), 1, 0, 0, EMIS )  ) THEN
                   CALL WRITE_MESG_EXIT( OUTNAMES(I,1), PROGNAME )
               END IF

           END IF

           DO S = 1, NSRC

              IF ( EMIS( S ) .LT. AMISS3 ) THEN
                 EMIS( S ) = 0.0
              END IF

           END DO ! end source loop

C...........  If CONTROL packet is present: For the current pollutant, read
C             in control efficiency, rule effectiveness, and, in the case of 
C             AREA sources, rule penetration.
           IF ( CFLAG ) THEN

              SELECT CASE( CATEGORY )

              CASE( 'AREA' )

              IF ( .NOT. READ3( ENAME, OUTNAMES(I,4), 1, 0, 0, 
     &        CTLEFF ) ) THEN
                 CALL WRITE_MESG_EXIT( OUTNAMES(I,4), PROGNAME )
              END IF

              IF ( .NOT. READ3( ENAME, OUTNAMES(I,5), 1, 0, 0, 
     &        RULEFF ) ) THEN
                 CALL WRITE_MESG_EXIT( OUTNAMES(I,5), PROGNAME )
              END IF

              IF ( .NOT. READ3( ENAME, OUTNAMES(I,6), 1, 0, 0, 
     &        RULPEN ) ) THEN
                 CALL WRITE_MESG_EXIT( OUTNAMES(I,6), PROGNAME )
              END IF

              CASE( 'POINT' )

              IF ( .NOT. READ3( ENAME, OUTNAMES(I,3), 1, 0, 0, 
     &        CTLEFF ) ) THEN
                 CALL WRITE_MESG_EXIT( OUTNAMES(I,3), PROGNAME )
              END IF

              IF ( .NOT. READ3( ENAME, OUTNAMES(I,4), 1, 0, 0, 
     &        RULEFF ) ) THEN
                 CALL WRITE_MESG_EXIT( OUTNAMES(I,4), PROGNAME )
              END IF

              RULPEN = 100.0  ! array

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

C..................  Overwrite temporary file line with new info
                 ELSE
                    E2 = E1
                    WRITE( CDEV,93300 ) 0, 'D', 0., 0., 0.

                 END IF

C..................  Store input and output emissions for plant
                 IF( CATEGORY .EQ. 'POINT' ) THEN
                     J = PLTINDX( S )
                     PLTINEM ( J,I ) = PLTINEM ( J,I ) + E1
                     PLTOUTEM( J,I ) = PLTOUTEM( J,I ) + E2

C......................  Flag plant if emissions are different
                     IF( E1 .NE. E2 ) PLTFLAG( J ) = .TRUE.

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

C..................  Overwrite temporary file line with new info
                 ELSE
                    E2 = E1
                    WRITE( CDEV,93300 ) 0, 'D', 0., 0., 0.

                 END IF

C..................  Store input and output emissions for plant
                 IF( CATEGORY .EQ. 'POINT' ) THEN
                     J = PLTINDX( S )
                     PLTINEM ( J,I ) = PLTINEM ( J,I ) + E1
                     PLTOUTEM( J,I ) = PLTOUTEM( J,I ) + E2

C......................  Flag plant if emissions are different
                     IF( E1 .NE. E2 ) PLTFLAG( J ) = .TRUE.

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

C..................  If no controls, then overwrite temporary line only
                 ELSE
                     WRITE( LDEV,93300 ) 0, 'D', 0., 0., 0.

                 END IF

              END DO ! end source loop

           END IF

C...........  Write multiplicative controls for current pollutant

            IF( .NOT. WRITE3( MNAME, PNAM, 0, 0, FACTOR ) ) THEN
                MESG = 'Failed to write multiplicative control ' // 
     &                 'factors for pollutant ' // PNAM
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

        END DO ! end pollutant loop

C.........  Write out controlled facilities report for point sources
        IF( CATEGORY .EQ. 'POINT' ) THEN
            WRITE( RDEV, 93000 ) 'Controlled facilities report'
            WRITE( RDEV, 93000 ) 'Adjustments from /CONTROL/ or ' //
     &                           '/EMS CONTROL/ packets only'
            WRITE( RDEV, 93400 ) ( PNAMMULT( I ), I=1,NVCMULT )
            WRITE( OUTFMT, '( A, I3.3, A )' ) 
     &             '(30X, ', NVCMULT, '( "Before", 7X, "After", 6X ) )'
            WRITE( RDEV, OUTFMT )

            PIDX = 0
            DO S = 1, NSRC
        	 IDX = PLTINDX( S )
 
        	 IF( IDX .NE. PIDX .AND.
     &               PLTFLAG( IDX )      ) THEN

                     J = FIPLEN3 + 1
                     CSRC = CSOURC( S )                    
                     WRITE( RDEV, 93410 ) 
     &                      CSRC( 1:FIPLEN3 ), CSRC( J:FPLLEN3 ),
     &                      ( PLTINEM ( IDX,I ), PLTOUTEM( IDX,I ),
     &                        I = 1, NVCMULT )

                     PIDX = IDX
        	 END IF
            END DO
        END IF

C.........  Deallocate local memory
        DEALLOCATE( ALWINDX, CTGINDX, CTLINDX, BACKOUT, CTLEFF,
     &              EMIS, FACTOR, RULEFF, RULPEN )

        IF( CATEGORY .EQ. 'POINT' ) 
     &      DEALLOCATE( PLTINDX, PLTFLAG, PLTINEM, PLTOUTEM )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

93300   FORMAT( I2, 1X, '"', A, '"', 3( 1X, E12.5 ) )

93400   FORMAT( 'Co/St/Cy', 4X, 'Facility ID', 13X, 100( A16, :, 8X ) )

93410   FORMAT( 1X, A6, 2X, A15, 1X, 100( F11.2, :, 1X ) )

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
