
        SUBROUTINE GENMULTC( ADEV, CDEV, GDEV, LDEV, NCPE, ENAME, 
     &                       MNAME, CFLAG, GFLAG, LFLAG )

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
C****************************************************************************/
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1999, MCNC--North Carolina Supercomputing Center
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
        INTEGER         GETEFILE
        INTEGER         PROMPTFFILE

        EXTERNAL   GETEFILE, PROMPTFFILE

C...........   SUBROUTINE ARGUMENTS

        INTEGER     , INTENT (IN) ::  ADEV   ! file unit no. for tmp ADD file
        INTEGER     , INTENT (IN) ::  CDEV   ! file unit no. for tmp CTL file 
        INTEGER     , INTENT (IN) ::  GDEV   ! file unit no. for tmp CTG file
        INTEGER     , INTENT (IN) ::  LDEV   ! file unit no. for tmp ALW file
        INTEGER     , INTENT (IN) ::  NCPE   ! no. of control packet entries
        CHARACTER*16, INTENT (IN) ::  ENAME  ! logical name for i/o api 
                                             ! inventory input file
        CHARACTER*16, INTENT (IN) ::  MNAME  ! logical name for mult. cntl. mat.
        LOGICAL     , INTENT (IN) ::  CFLAG  ! true = apply CTL controls
        LOGICAL     , INTENT (IN) ::  GFLAG  ! true = apply CTG controls
        LOGICAL     , INTENT (IN) ::  LFLAG  ! true = apply ALW controls

C...........   Local temporary arrays
        INTEGER      ALWINDX  ( NSRC,NVCMULT ) ! indices to ALW controls table
        INTEGER      CTGINDX  ( NSRC,NVCMULT ) ! indices to CTG controls table
        INTEGER      CTLINDX  ( NSRC,NVCMULT ) ! indices to CTL controls table
        INTEGER      OUTTYPES ( NVCMULT,6 )    ! var type:int/real


        REAL         BACKOUT( NSRC )  ! factor used to account for pollutant
                                      ! specific control information that is
                                      ! already in the inventory
        REAL         CTLEFF ( NSRC )  ! control efficiency
        REAL         EMIS   ( NSRC )  ! base inventory emissions
        REAL         FACTOR ( NSRC )  ! multiplicative controls
        REAL         RULEFF ( NSRC )  ! rule effectiveness
        REAL         RULPEN ( NSRC )  ! rule penetration

        CHARACTER(LEN=IOVLEN3)    OUTNAMES( NVCMULT,6 ) ! var names
        CHARACTER(LEN=IOULEN3)    OUTUNITS( NVCMULT,6 ) ! var units
        CHARACTER(LEN=IODLEN3)    OUTDESCS( NVCMULT,6 ) ! var descriptions

C...........   Other local variables
        INTEGER          I,J,K,S  ! counters and indices
        INTEGER          IOS      ! input/output status

        REAL             ALWFAC   ! allowable control factor
        REAL             ALWEMIS  ! allowable emissions
        REAL             CAP      ! emissions cap
        REAL             CTGFAC   ! control technology control factor
        REAL             CTGFAC2  ! MAXACF or RSNACF
        REAL             CUTOFF   ! CTG cutoff for application of control
        REAL             DENOM    ! denominator of control back-out factor
        REAL             MACT     ! max. achievable cntrl tech. cntrl factor
        REAL             RACT     ! reasonably achievable cntrl tech. cntrl
                                  ! factor
        REAL             REPLACE  ! replacement emissions

        CHARACTER        FILENM                ! file name
        CHARACTER, SAVE  ::  PATHNM            ! path name for tmp file
        CHARACTER*300    MESG                  ! message buffer
        CHARACTER*16  :: PROGNAME = 'GENMULTC' ! program name

C***********************************************************************
C   begin body of subroutine GENMULTC

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

C...........  Fractionalize control-packet information

        IF ( CFLAG ) THEN

           DO S = 1, NCPE
              FACCEFF( S ) = FACCEFF( S )/100.0
              FACREFF( S ) = FACREFF( S )/100.0
              FACRLPN( S ) = FACRLPN( S )/100.0
           END DO ! end source loop

        END IF

C...........  Loop through pollutants that receive controls

        DO I = 1, NVCMULT

C...........  Initialize control factor array

           FACTOR = 1.0  ! array

C...........  Read in emissions data from inventory file

           IF ( .NOT. READ3( ENAME, OUTNAMES(I,1), 1, 0, 0, 
     &     EMIS ) ) THEN
              CALL WRITE_MESG_EXIT( OUTNAMES(I,1), PROGNAME )
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

              RULPEN = 100.0

              END SELECT ! end select on category

C...........  Fractionalize all inventory control information and set
C             missing values to one, then calculate the factor which will
C             be used to account for control information already in the
C             inventory

              DO S = 1, NSRC

                 IF ( CTLEFF( S ) .LT. AMISS3 ) THEN
                    CTLEFF( S ) = 0.0
                 ELSE
                    CTLEFF( S ) = CTLEFF( S )/100.0
                 END IF

                 IF ( RULEFF( S ) .LT. AMISS3 ) THEN
                    RULEFF( S ) = 1.0
                 ELSE
                    RULEFF( S ) = RULEFF( S )/100.0
                 END IF

                 IF ( RULPEN( S ) .LT. AMISS3 ) THEN
                    RULPEN( S ) = 1.0
                 ELSE
                    RULPEN( S ) = RULPEN( S )/100.0
                 END IF

C...........  Perform division by zero check.

                 DENOM = ( 1.0 - CTLEFF(S)*RULEFF(S)*RULPEN(S) )
                 IF ( FLTERR( DENOM, 0.0 ) ) THEN
                    BACKOUT( S ) = 1.0/DENOM
                 ELSE
                    BACKOUT( S ) = 0.0
                 END IF

              END DO ! end source loop

C...........  If CTG packet is present:

           ELSE IF ( GFLAG ) THEN

           END IF

           IF ( CFLAG .AND. PCTLFLAG( I, 1 ) ) THEN

C...........  Compute CTL factor

              DO S = 1, NSRC

                 K = CTLINDX( S, I )
                 IF ( K .GT. 0 ) THEN
                    CTLEFF( S ) = FACCEFF( K )
                    RULEFF( S ) = FACREFF( K )
                    RULPEN( S ) = FACRLPN( K )
                    FACTOR( S ) = BACKOUT( S )*
     &                          ( 1.0 - CTLEFF(S)*RULEFF(S)*RULPEN(S) )
                 ELSE
                    FACTOR( S ) = 1.0
                 END IF

                 EMIS( S ) = EMIS( S )*FACTOR( S ) ! apply controls

              END DO ! end source loop

           END IF

           IF ( GFLAG .AND. PCTLFLAG( I, 2 ) ) THEN

C...........  Compute CTG factor

              DO S = 1, NSRC

                 K = CTGINDX( S, I ) 
                 IF ( K .GT. 0 ) THEN
                    CUTOFF = CUTCTG ( K )
                    CTGFAC = FACCTG ( K )
                    MACT   = FACMACT( K )
                    RACT   = FACRACT( K )

C...........  Determine if MACT or RACT is to be used as the second
C             CTG factor

                    IF ( MACT .GE. 0 .AND. RACT .GE. 0 ) THEN
                       CTGFAC2 = MACT
                    ELSE IF ( MACT .GE. 0 .AND. RACT .LT. 0 ) THEN
                       CTGFAC2 = MACT
                    ELSE IF ( MACT .LT. 0 .AND. RACT .GE. 0 ) THEN
                       CTGFAC2 = RACT
                    ELSE
                       CTGFAC2 = 1.0
                    END IF

C...........  Check to see if emissions exceed cutoff and if necessary,
C             apply controls

                    IF ( EMIS( S ) .GT. CUTOFF ) THEN

                       FACTOR( S ) = FACTOR( S )*CTGFAC
                       EMIS( S )   = EMIS( S )*CTGFAC

C...........  If emissions still exceed cutoff, apply second CTG factor

                       IF ( EMIS( S ) .GT. CUTOFF ) THEN
                          FACTOR( S ) = FACTOR( S )*CTGFAC2
                       END IF

                    END IF

                 END IF

              END DO ! end source loop

           END IF

           IF ( LFLAG .AND. PCTLFLAG( I, 3 ) ) THEN

C...........  Process ALW packet

              DO S = 1, NSRC

                 K = ALWINDX( S, I ) 
                 IF ( K .GT. 0 ) THEN
                    ALWFAC  = FACALW  ( K )
                    CAP     = EMCAPALW( K )
                    REPLACE = EMREPALW( K )

C...........  Determine whether CAP or REPLACE is to be used as the value for
C             for allowable emissions. Check to see if emission exceed this
C             value.  Compute ALW control factor, if necessary.

                    IF ( CAP .GE. 0 .AND. REPLACE .GE. 0 ) THEN

                       IF ( EMIS( S ) .GT. CAP ) THEN
                          FACTOR( S ) = REPLACE/EMIS( S )
                       END IF

                    ELSE IF ( CAP .GE. 0 .AND. REPLACE .LT. 0 ) THEN

                       IF ( EMIS( S ) .GT. CAP ) THEN
                          FACTOR( S ) = CAP/EMIS( S )
                       END IF

                    ELSE IF ( CAP .LT. 0 .AND. REPLACE .GE. 0 ) THEN

                       IF ( EMIS( S ) .GT. REPLACE ) THEN
                          FACTOR( S ) = REPLACE/EMIS( S )
                       END IF

                    ELSE
                       
                    END IF

                 END IF

              END DO ! end source loop

           END IF

C...........  Write multiplicative controls for current pollutant

        IF( .NOT. WRITE3( MNAME, PNAMMULT( I ), 0, 0, FACTOR ) ) THEN
            MESG = 'Failed to write multiplicative control factors 
     &              for pollutant ' // PNAMMULT( I )
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        END DO ! end pollutant loop

        RETURN


C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

93000   FORMAT( A )

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
