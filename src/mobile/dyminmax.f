
        SUBROUTINE DYMINMAX( NSRC, JTIME, LUPDOUT, DAYBEGT, DAYENDT, 
     &                       VALBYSRC, MINBYSRC, MAXBYSRC, 
     &                       MINOUT, MAXOUT )

C***********************************************************************
C  subroutine DYMINMAX body starts at line < >
C
C  DESCRIPTION:
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION HISTORY:
C
C***************************************************************************
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
C****************************************************************************

        IMPLICIT NONE

C...........   INCLUDES:
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters

C...........   SUBROUTINE ARGUMENTS
        INTEGER     , INTENT    (IN) :: NSRC             ! no. sources
        INTEGER     , INTENT    (IN) :: JTIME            ! HHMMSS
        LOGICAL     , INTENT    (IN) :: LUPDOUT          ! true: force output update
        INTEGER     , INTENT    (IN) :: DAYBEGT ( NSRC ) ! begin. time for day
        INTEGER     , INTENT    (IN) :: DAYENDT ( NSRC ) ! ending time for day
        REAL        , INTENT    (IN) :: VALBYSRC( NSRC ) ! per-source values
        REAL        , INTENT(IN OUT) :: MINBYSRC( NSRC ) ! min (over day) 
        REAL        , INTENT(IN OUT) :: MAXBYSRC( NSRC ) ! max (over day)
        REAL        , INTENT(IN OUT) :: MINOUT  ( NSRC ) ! min for output
        REAL        , INTENT(IN OUT) :: MAXOUT  ( NSRC ) ! max for output

C...........   Other local variables
        INTEGER     S           ! counters and indices

        REAL        VAL         ! tmp value

        LOGICAL, SAVE :: IFLAG   = .FALSE. ! true: set initial to false
        LOGICAL, SAVE :: INITIAL = .TRUE.  ! true: 'til min/max by src initlzd

        CHARACTER*16 :: PROGNAME = 'DYMINMAX' ! program name

C***********************************************************************
C   begin body of subroutine DYMINMAX

C.........  For the first time, initialize all entries to missing
        IF( INITIAL ) THEN
            MINBYSRC = BADVAL3  ! array
            MAXBYSRC = BADVAL3  ! array
        END IF

C.........  Loop through sources
C.........  Process according to the source's time of day
        DO S = 1, NSRC

            VAL = VALBYSRC( S )

C.............  Initialize the min/max by source when the time is equal to the
C               start of that source's day or the routine is being called for
C               the first time

            IF( VAL .GT. AMISS3 .AND. 
     &        ( INITIAL .OR. JTIME .EQ. DAYBEGT( S ) ) )THEN

                MINBYSRC( S ) = VAL
                MAXBYSRC( S ) = VAL

C.................  Turn off setting after min/max have been initialized
                IFLAG = .TRUE.

C.............  Update the minimum and maximum by source for every time step
            ELSE IF ( VAL .GT. AMISS3 ) THEN

                IF( VAL .LT. MINBYSRC( S ) ) MINBYSRC( S ) = VAL
                IF( VAL .GT. MAXBYSRC( S ) ) MAXBYSRC( S ) = VAL

            END IF

C.............  Update the output min/max by source at the end of each
C               source's day, or if an update is being forced by the 
C               subroutine argument.
            IF( LUPDOUT .OR. JTIME .EQ. DAYENDT( S ) ) THEN
                MINOUT( S ) = MINBYSRC( S )
                MAXOUT( S ) = MAXBYSRC( S )
            END IF

        END DO

        IF ( IFLAG ) INITIAL = .FALSE.

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I9, :, 1X ) )
 
        END SUBROUTINE DYMINMAX
