
        SUBROUTINE EXPNDEFT
   
C***********************************************************************
C  subroutine EXPNDEFT body starts at line < >
C
C  DESCRIPTION:
C     Expand the records of the emission-factor index and min/max 
C     temperature combinations so that pure PSIs contain the temperature
C     indices from the combination factors that use them.  Also expand the
C     list so that any scenario group root PSIs that have been added will
C     have the temperature indices of the PSIs that are grouped with them.
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
C****************************************************************************
 
C...........   MODULES for public variables
C.........  This module contains emission factor tables and related
        USE MODEMFAC

        IMPLICIT NONE

C...........   INCLUDES:
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS
        CHARACTER*2   CRLF
        INTEGER       FIND1

        EXTERNAL      CRLF, FIND1

C...........   Local parameters
        INTEGER, PARAMETER :: MXCOL = 5   ! maximum number of input columns

C...........   Local allocatable arrays
        INTEGER, ALLOCATABLE :: NEWPSI( : )  ! tmp array of PSIs
        INTEGER, ALLOCATABLE :: NEWTMI( : )  ! tmp array of min/max tmpr indices
        INTEGER, ALLOCATABLE :: NEWIDX( : )  ! tmp sort array

C...........   Local variables
        INTEGER         I, J, JJ, K, KK, M, N, NN

        INTEGER         CINDX     ! tmp index to contributing PSIs table
        INTEGER         CINDXB    ! tmp index to contributing PSIs table
        INTEGER         CNTCOMBO  ! tmp count of contributing PSIs in combo PSI
        INTEGER         CPSI      ! tmp contributing PSI
        INTEGER         ICNT      ! tmp counter for new array entries
        INTEGER         INDXTMP   ! tmp pos of root or contrib PSI in TMMI list
        INTEGER         IOS       ! i/o status
        INTEGER         IREC      ! record counter
        INTEGER         KB, KE    ! begining and end of K loop
        INTEGER         MXINCRSE  ! maximum increase for NMMTEF
        INTEGER         MXMMTEF   ! maximum new NMMTEF before remove duplicates
        INTEGER         MXNEW     ! maximum total new allocation
        INTEGER         NTMMI     ! no. TMMIs in list for root or contrib PSI
        INTEGER         PPSI      ! PSI of previous iteration
        INTEGER         PSI       ! tmp parameter scheme index
        INTEGER         PSIROOT   ! tmp 1st PSI of a PSI group
        INTEGER         PTMMI     ! TMMI of previous iteration
        INTEGER         TMMI      ! tmp min/max temperature index

        LOGICAL      :: EFLAG = .FALSE.  ! true: error found
        LOGICAL      :: WFLAG = .TRUE.   ! true: write out warning

        CHARACTER*300          MESG      ! message buffer
        CHARACTER(LEN=IOVLEN3) ACT

        CHARACTER*16 :: PROGNAME = 'EXPNDEFT' ! program name

C***********************************************************************
C   begin body of subroutine EXPNDEFT

C.........  Determine the new size by looking at the number of temperature
C           indexes in each combination EF.  If there are no combination
C           EFs, then the total will be 0 and the routine can abort.

C.........  The PDATPNTR points to either the table of contributing PSIs for
C           combination PSIs, or it is the position of the PSI in a multi-
C           scenario group.

C.........  Loop through PSIs actually being used
        MXINCRSE = 0
        DO I = 1, NPSIALL

            PSI = PSIALL( I )

C............. Find PSI in emission factors data tables
            J = FIND1( PSI, NPSIDAT, PSIDAT )
            K = PDATINDX( J )

            CNTCOMBO = PDATTYPE( K ) ! pure or combo
            CINDX    = PDATPNTR( K ) ! pointer (see above)

C.............  If not a combination emission factor, and not a root PSI, 
C               increase max size by the number of tmpr combos for this EF
            IF( CNTCOMBO .LE. 0 .AND. CINDX .GT. 1 ) THEN
                MXINCRSE = MXINCRSE + MAX( 0,MMTEFNUM( I ) )

C.............  If a combination emission factor, increase max size by the
C               number of min/max tmprs that the combo PSI uses times the
C               number of contributing PSIs times 2 (in case the contributing
C               PSIs need their root PSIs from a group updated as well.
            ELSE IF( CNTCOMBO .GT. 0 ) THEN
                MXINCRSE = MXINCRSE + 
     &                     MAX( 0,MMTEFNUM( I ) ) * CNTCOMBO * 2

            END IF

        END DO

C.........  If there will be no increase, leave routine now
        IF( MXINCRSE .EQ. 0 ) RETURN

C.........  Allocate local arrays for storing new tables
        MXNEW = NMMTEF + MXINCRSE
        ALLOCATE( NEWPSI( MXNEW ),STAT=IOS )
        CALL CHECKMEM( IOS, 'NEWPSI', PROGNAME )
        ALLOCATE( NEWTMI( MXNEW ),STAT=IOS )
        CALL CHECKMEM( IOS, 'NEWTMI', PROGNAME )
        ALLOCATE( NEWIDX( MXNEW ),STAT=IOS )
        CALL CHECKMEM( IOS, 'NEWIDX', PROGNAME )

C.........  Initialize local arrays with existing arrays
        NEWPSI( 1:NMMTEF ) = MMTEFPSI( 1:NMMTEF )
        NEWTMI( 1:NMMTEF ) = MMTEFIDX( 1:NMMTEF )

        DO I = 1, NMMTEF
            NEWIDX( I ) = I
        END DO

C.........  Store local tables with unsorted new temperature indexes for 
C           PSIs that contribute to combination PSIs.
        ICNT = NMMTEF
        DO I = 1, NPSIALL

            PSI = PSIALL( I )

C............. Find PSI in emission factors data tables
            J = FIND1( PSI, NPSIDAT, PSIDAT )
            K = PDATINDX( J )  ! sorted position

            CNTCOMBO = PDATTYPE( K ) ! pure or combo
            PSIROOT  = PDATROOT( K ) ! root PSI for multi-scenario group
            CINDX    = PDATPNTR( K ) ! pointer (see above)

C.............  If not a combination emission factor, and is in a EF group, 
C               but is not a root PSI, add the temperature indices to the
C               list for the root PSI
            IF( CNTCOMBO .LE. 0 .AND. CINDX .GT. 1 ) THEN

                N = FIND1( PSIROOT, NPSIALL, PSIALL )
                IF( N .LE. 0 ) THEN
                    EFLAG = .TRUE.
                    CALL PSI_ERROR( PSIROOT )
                    CYCLE
                END IF
                INDXTMP = MMTEFPTR( N )
                NTMMI  = MMTEFNUM( N )

C.................  Check to make sure the root PSI *has* any factors, and
C                   if not, skip it because it will be updated with the 
C                   combination factors
                IF( INDXTMP .LE. 0 ) CYCLE

C.................  For each temperature index of current PSI, search index
C                   of root PSI and if not there, add to unsorted list
                KB = MMTEFPTR( I )
                KE = KB + MMTEFNUM( I ) - 1

                DO KK = KB, KE

                    TMMI = MMTEFIDX( KK )
                    M    = FIND1( TMMI, NTMMI, MMTEFIDX( INDXTMP ) )

                    IF( M .LE. 0 ) THEN
                        ICNT = ICNT + 1

                        IF( ICNT .LE. MXNEW ) THEN
                            NEWIDX( ICNT ) = ICNT
                            NEWPSI( ICNT ) = PSIROOT
                            NEWTMI( ICNT ) = TMMI
                        END IF  ! end dimension check
                    END IF      ! end for TMMI not found

                END DO          ! end loop on temperature indices

C.............  If a combination emission factor...
            ELSE IF( CNTCOMBO .GT. 0 ) THEN

C.................  Loop through contributing PSIs
        	DO J = 1, CNTCOMBO

                    CPSI = CMBOPSI( J,CINDX )  ! Obtain contributing PSI

C.....................  Find contributing PSI in list of all PSIs
                    N = FIND1( CPSI, NPSIALL, PSIALL )
                    IF( N .LE. 0 ) THEN
                        EFLAG = .TRUE.
                        CALL PSI_ERROR( CPSI ) 
                        CYCLE
                    END IF

                    INDXTMP = MMTEFPTR( N )
                    NTMMI   = MMTEFNUM( N )

C.....................  For each temperature index of combo PSI, search indices
C                       for contributing PSI and if not there, add to unsorted 
C                       list
C.....................  Must account for the case where there are no temperature
C                       indices for the contributing PSI by checking the
C                       value of INDXTMP for values < 0
                    KB = MMTEFPTR( I )
                    KE = KB + MMTEFNUM( I ) - 1

                    DO KK = KB, KE  ! Loop through indices of combo PSI

                        TMMI = MMTEFIDX( KK )
                        M    = -1
                        IF( INDXTMP .GT. 0 ) THEN
                            M = FIND1( TMMI, NTMMI, MMTEFIDX(INDXTMP) )
                        END IF

                	IF( M .LE. 0 ) THEN
                            ICNT = ICNT + 1

                            IF( ICNT .LE. MXNEW ) THEN
                        	NEWIDX( ICNT ) = ICNT
                        	NEWPSI( ICNT ) = CPSI
                        	NEWTMI( ICNT ) = TMMI
                            END IF  ! end dimension check
                	END IF  ! end for TMMI not found

C.........................  Now, if the contributing PSI is a member of an EF
C                           group, update the table for the root PSI of the
C                           group also.
                        JJ = FIND1( CPSI, NPSIDAT, PSIDAT )
        		NN = PDATINDX( JJ )  ! sorted position

        		PSIROOT  = PDATROOT( NN ) ! root PSI 
                        CINDXB   = PDATPNTR( NN ) ! pointer (see above)

                        IF( CINDXB .GT. 1 ) THEN
                            ICNT = ICNT + 1
                            IF( ICNT .LE. MXNEW ) THEN
                        	NEWIDX( ICNT ) = ICNT
                        	NEWPSI( ICNT ) = PSIROOT
                        	NEWTMI( ICNT ) = TMMI
                            END IF  ! end dimension check
                        END IF

                    END DO      ! end loop on temperature indices

        	END DO          ! end loop on contributing PSIs for combo PSI
            END IF              ! end section for combination PSIs

        END DO                  ! end loop on all PSIs
 
        MXMMTEF = ICNT

        IF( MXMMTEF .GT. MXNEW ) THEN
            WRITE( MESG, 94010 )
     &             'INTERNAL ERROR: dimension for storing ' //
     &             'new EF-ref/temperature table' // CRLF()// BLANK10// 
     &             'was', MXNEW, 'but actually needed', MXMMTEF
            CALL M3MSG2( MESG )
            CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )
        END IF

C.........  Only error is PSI error from PSI_ERROR call (see below)
        IF( EFLAG ) THEN

            MESG = 'Update your inventory to include the missing ' //
     &             'sources, or update your MPREF file to ' //
     &             CRLF() // BLANK10 // 'to not use the listed PSIs.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF

C.........  Sort the unsorted new arrays
        CALL SORTI2( MXMMTEF, NEWIDX, NEWPSI, NEWTMI )

C.........  Deallocate original tables and reallocate for new size
        DEALLOCATE( MMTEFPSI, MMTEFIDX )

        ALLOCATE( MMTEFPSI( MXMMTEF ),STAT=IOS )
        CALL CHECKMEM( IOS, 'MMTEFPSI', PROGNAME )
        ALLOCATE( MMTEFIDX( MXMMTEF ),STAT=IOS )
        CALL CHECKMEM( IOS, 'MMTEFIDX', PROGNAME )

C.........  Update the EF-ref/temperature arrays
C.........  NOTE - there could be duplicates in the new array because we
C           were not updating as we went along. So, make sure these are 
C           eliminated by checking previous values.
        PPSI  = -9
        PTMMI = -9
        K     = 0
        DO I = 1, MXMMTEF

            J    = NEWIDX( I )
            PSI  = NEWPSI( J )
            TMMI = NEWTMI( J )

            IF( PSI .NE. PPSI .OR. TMMI .NE. PTMMI ) THEN
                K = K + 1
                MMTEFPSI( K ) = PSI
                MMTEFIDX( K ) = TMMI
            END IF

            PPSI  = PSI
            PTMMI = TMMI

        END DO

        NMMTEF = K

C.........  Initilialize before upcoming loop
        MMTEFPTR = 0   ! array
        MMTEFNUM = 0   ! array
        MMTEFEND = 0   ! array

C.........  Set the pointers to the position in the temperature index
c           array of the first and last occurance of each PSI
        PPSI = -9
        DO I = 1, NMMTEF

            PSI = MMTEFPSI( I )

            IF( PSI .NE. PPSI ) THEN

                K = FIND1( PSI, NPSIALL, PSIALL )

C.................  Report internal error if not found
                IF( K .LE. 0 ) THEN
                    WRITE( MESG,94010 ) 'INTERNAL ERROR: PSI', PSI,
     &                     'from EF-ref/temperature array not found '//
     &                     'master PSI list'
                    CALL M3MSG2( MESG )
                    CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )
                ELSE
                    MMTEFPTR( K ) = I

                END IF

                PPSI = PSI

            END IF

            MMTEFNUM( K ) = I - MMTEFPTR( K ) + 1
            MMTEFEND( K ) = MMTEFPTR( K ) + MMTEFNUM( K ) - 1

        END DO

C.........  Deallocate local memory from the routine
        DEALLOCATE( NEWIDX, NEWPSI, NEWTMI )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I9, :, 1X ) )

C******************  INTERNAL SUBPROGRAMS  *****************************

        CONTAINS

C..............  This internal subprogram provides an error when the PSI
C                is not found in the PSIs for sources.
            SUBROUTINE PSI_ERROR( PSI )

C.............  Subroutine arguments
            INTEGER     , INTENT (IN):: PSI

C.............  Local variables
            CHARACTER*300 MESG

C......................................................................

            WRITE( MESG,94010 ) 'ERROR: PSI ', PSI, 
     &             'is not found in MEFTEMP file. This can occur'//
     &             CRLF() // BLANK10 // 'when the sources using this '//
     &             'PSI are not in the inventory.'
            CALL M3MSG2( MESG )

C--------------------  FORMAT  STATEMENTS   ----------------------------

94010       FORMAT( 10( A, :, I8, :, 1X ) )

            END SUBROUTINE PSI_ERROR

        END SUBROUTINE EXPNDEFT
