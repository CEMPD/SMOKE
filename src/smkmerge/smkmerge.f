
        PROGRAM SMKMERGE

C***********************************************************************
C  program SMKMERGE body starts at line
C
C  DESCRIPTION:
C      The purpose of this program is to merge the inventory or hourly
C      emissions files from the tmp* programs with gridding matrices and with
C      optionally any combination of speciation matrices and 3 control
C      matrices (different types).  The program can operate on from 1 to 4 
C      source categories (area, biogenic, mobile, or point sources), or any 
C      combination of these.  If a layer fractions file is input, then the 
C      output file is 3-d.  This program is not used for the MPS/MEPSE files 
C      for CMAQ.
C
C  PRECONDITIONS REQUIRED:  
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C       Copied from csgldaymrg.F version 1.7 by M Houyoux 2/99
C
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
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
C Pathname: $Source$
C Last updated: $Date$ 
C
C****************************************************************************

C.........  MODULES for public variables
C.........  This module contains the major data structure and control flags
        USE MODMERGE

C.........  This module contains the control packet data and control matrices
        USE MODCNTRL

        IMPLICIT NONE

C...........   INCLUDES:
        
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file desc. data structures

C...........   EXTERNAL FUNCTIONS and their descriptions:
        
        CHARACTER*2     CRLF
        CHARACTER*10    HHMMSS
        CHARACTER*14    MMDDYY
        INTEGER         WKDAY

        EXTERNAL    CRLF, HHMMSS, MMDDYY, WKDAY

C.........  LOCAL PARAMETERS and their descriptions:

        CHARACTER*50, PARAMETER :: SCCSW = '@(#)$Id$'

C...........   LOCAL VARIABLES and their descriptions:

C...........   Logical names and unit numbers (not in MODMERGE)
        INTEGER         LDEV
     
C...........   Other local variables
    
        INTEGER          J, K, L1, L2, L3, N, V, T ! counters and indices

        INTEGER          IOS           ! tmp I/O status
        INTEGER          JDATE         ! Julian date (YYYYDDD)
        INTEGER          JTIME         ! time (HHMMSS)
        INTEGER          K1, K2, K3    ! tmp indices for valid speciation matrix
        INTEGER          K4, K5, K6    ! tmp indices for valid reactivity matrix
        INTEGER          LDATE         ! Julian date from previous iteration
        INTEGER          MXGRP         ! max no. of pollutant groups
        INTEGER          MXPOLPGP      ! max no. of pols per group
        INTEGER          MXSPPOL       ! max no. of spcs per pol/group
        INTEGER          NGRP          ! actual no. of pollutant groups
        INTEGER          PGID          ! previous iteration group ID no.
        INTEGER          NPPGP         ! tmp actual no. pols per group

        REAL          :: RDUM = 0      ! dummy real value

        CHARACTER(LEN=PLSLEN3) SVBUF   ! pol to species description buffer
        CHARACTER(LEN=IOVLEN3) SBUF    ! tmp species name
        CHARACTER(LEN=IOVLEN3) PBUF    ! tmp pollutant name
        CHARACTER*300    MESG

        CHARACTER*16  :: PROGNAME = 'SMKMERGE' ! program name

C***********************************************************************
C   begin body of program SMKMERGE
        
        LDEV = INIT3()

C.........  Write out copywrite, version, web address, header info, and prompt
C           to continue running the program.
        CALL INITEM( LDEV, SCCSW, PROGNAME )

C.........  Retrieve control environment variables and set logical control
C           flags. Use a local module to pass the control flags.
        CALL GETMRGEV
        
C.........  Open input files and retrieve episode information
        CALL OPENMRGIN

C.........  Determine number of states and counties for each source category
!        CALL 
! NOTE: add later, when adding state/county reports

C.........  Create arrays of sorted unique pol-to-species
C.........  Create arrays of sorted unique pollutants
C.........  Create arrays of sorted unique species
        CALL MRGVNAMS

C.........  Allocate memory for fixed-size arrays by source category...
C.........  Area source fixed-size arrays
        CALL ALLOCMRG( MXGRP, MXPOLPGP, MXSPPOL )

C.........  Read reactivity matrices
        IF( ARFLAG ) CALL RDRMAT( ARNAME, ANSREAC, ARNMSPC, ACRREPEM, 
     &                            ACRIDX, ACRPRJFC, ACRMKTPN, ACRFAC )

        IF( MRFLAG ) CALL RDRMAT( MRNAME, MNSREAC, MRNMSPC, MCRREPEM, 
     &                            MCRIDX, MCRPRJFC, MCRMKTPN, MCRFAC )

        IF( PRFLAG ) CALL RDRMAT( PRNAME, PNSREAC, PRNMSPC, PCRREPEM, 
     &                            PCRIDX, PCRPRJFC, PCRMKTPN, PCRFAC )

C.........  Read gridding matrices (note: must do through subroutine because of
C           needing contiguous allocation for integer and reals)
        IF( AFLAG ) CALL RDGMAT( AGNAME, NGRID, ANGMAT, ANGMAT,
     &                           AGMATX( 1 ), AGMATX( NGRID + 1 ),
     &                           AGMATX( NGRID + ANGMAT + 1 ) )

        IF( MFLAG ) CALL RDGMAT( MGNAME, NGRID, MNGMAT, MNGMAT,
     &                           MGMATX( 1 ), MGMATX( NGRID + 1 ),
     &                           MGMATX( NGRID + MNGMAT + 1 ) )

        IF( PFLAG ) THEN

            PGMATX = 1.  ! initialize array b/c latter part not in file
            CALL RDGMAT( PGNAME, NGRID, NPSRC, 1,
     &                   PGMATX( 1 ), PGMATX( NGRID + 1 ), RDUM )
        END IF

C.........  Build indicies for pollutant/species groups
        CALL BLDMRGIDX( MXGRP, MXPOLPGP, MXSPPOL, NGRP )

C.........  Open NetCDF output files, open ASCII report files, and write headers
        CALL OPENMRGOUT
C NOTE: This is partially written.  Enough to write out the gridded output file.

C.........  Loop through processing groups (if speciation, this will be specia-
C           tion groups, but if no speciation, this will be pollutant groups,  
C           for purposes of memory usage if many pollutants and/or species)
        PGID = IMISS3
        DO N = 1, NGRP

            NPPGP = PLCNT( N )

C.............  If new pollutant group, read pollutant-specific control matrices
C.............  Note that only the pollutants in this group that are actually
C               in the control matrices are stored, and the index that says
C               which are valid is *U_EXIST and *A_EXIST
            IF( IDPGP( N ) .NE. PGID ) THEN

                IF( AUFLAG )
     &              CALL RD3MASK( AUNAME, 0, 0, NASRC, NPPGP, 
     &                      PLNAMES( 1,N ), AU_EXIST( 1,N ), ACUMATX )

                IF( MUFLAG )
     &              CALL RD3MASK( MUNAME, 0, 0, NMSRC, NPPGP, 
     &                      PLNAMES( 1,N ), MU_EXIST( 1,N ), MCUMATX )

                IF( PUFLAG )
     &              CALL RD3MASK( PUNAME, 0, 0, NPSRC, NPPGP, 
     &                      PLNAMES( 1,N ), PU_EXIST( 1,N ), PCUMATX )

                IF( AAFLAG )
     &              CALL RD3MASK( AANAME, 0, 0, NASRC, NPPGP, 
     &                      PLNAMES( 1,N ), AA_EXIST( 1,N ), ACAMATX )

                IF( MAFLAG )
     &              CALL RD3MASK( MANAME, 0, 0, NMSRC, NPPGP, 
     &                      PLNAMES( 1,N ), MU_EXIST( 1,N ), MCAMATX )

                IF( PAFLAG )
     &              CALL RD3MASK( PANAME, 0, 0, NPSRC, NPPGP, 
     &                      PLNAMES( 1,N ), PA_EXIST( 1,N ), PCAMATX )

            END IF

C.............  For speciation, read speciation matrix entries that are in
C               this group.
C.............  Also, Print message about which species are being processed
            IF( SFLAG ) THEN

                MESG = ' '
                L3 = 0
                DO V = 1, NPPGP

                    DO J = 1, NSMPPG( V,N )

                        K = SMINDEX( J,V,N )
                        IF( K .LE. 0 ) THEN
                            EXIT   ! End this loop
                        ELSE
                            SVBUF = TSVDESC( SMINDEX( J,V,N ) )
                            SBUF  = EMNAM  ( SPINDEX( J,V,N ) )
                            L1 = LEN_TRIM( MESG )
                            L2 = LEN_TRIM( SBUF )
                            L3 = L3 + L2 + 3
                            IF( L3 .GT. 60 ) THEN
                                MESG = MESG( 1:L1 )// CRLF()// BLANK10// 
     &                                 ' "' // SBUF( 1:L2 ) // '"'
                                L3 = 0
                            ELSE
                                MESG = MESG( 1:L1 ) // ' "' // 
     &                                 SBUF( 1:L2 ) // '"'
                            END IF
                        END IF

C.........................  Set position for input of speciation matrix
                        IF( AFLAG ) K1 = AS_EXIST( J,V,N ) 
                        IF( MFLAG ) K2 = MS_EXIST( J,V,N ) 
                        IF( PFLAG ) K3 = PS_EXIST( J,V,N ) 

C.........................  Read speciation matrix for current variable and
C                           position
                        IF ( K1 .GT. 0 )
     &                    CALL RDSMAT( ASNAME, SVBUF, ASMATX( 1,K1 ) )
                        IF ( K2 .GT. 0 )
     &                    CALL RDSMAT( MSNAME, SVBUF, MSMATX( 1,K2 ) )
                        IF ( K3 .GT. 0 )
     &                    CALL RDSMAT( PSNAME, SVBUF, PSMATX( 1,K3 ) )

C.........................  Set position for application of reactivity matrix
                        IF( ARFLAG ) K4 = AR_EXIST( J,V,N ) 
                        IF( MRFLAG ) K5 = MR_EXIST( J,V,N ) 
                        IF( PRFLAG ) K6 = PR_EXIST( J,V,N ) 

C.........................  Adjust speciation matrix, if needed, with reactivty
C                           matrix, and store old speciation (for use in 
C                           creating state/county totals)
C.........................  Also, update speciation "EXIST" arrays
                        IF ( K4 .GT. 0 ) THEN
                          CALL APPLRSPC( ANSREAC, LMKTPON, ACRIDX,
     &                                   ACRMKTPN, ACRFAC(1,K4), 
     &                                   ASMATSAV(1,K1), ASMATX(1,K1) )
                          AS_EXIST( J,V,N ) = K4
                        END IF

                        IF ( K5 .GT. 0 ) THEN
                          CALL APPLRSPC( MNSREAC, LMKTPON, MCRIDX,
     &                                   MCRMKTPN, MCRFAC(1,K5),
     &                                   MSMATSAV(1,K1), MSMATX(1,K1) )
                          MS_EXIST( J,V,N ) = K5
                        END IF

                        IF ( K6 .GT. 0 ) THEN
                          CALL APPLRSPC( PNSREAC, LMKTPON, PCRIDX,
     &                                   PCRMKTPN, PCRFAC(1,K6),
     &                                   PSMATSAV(1,K1), PSMATX(1,K1) )
                          PS_EXIST( J,V,N ) = K6
                        END IF

                    END DO  ! End pol-to-species loop

                END DO      ! End pollutant loop

                MESG = 'Processing species:'// CRLF()// BLANK10// MESG
                CALL M3MSG2( MESG )

C.............  Otherwise, print message about which pollutants are being 
C               processed
            ELSE

                MESG = ' '
                L3 = 0
                DO V = 1, NPPGP
                    PBUF = PLNAMES( V,N )
                    L1 = LEN_TRIM( MESG )
                    L2 = LEN_TRIM( PBUF )
                    L3 = L3 + L2 + 3
                    IF( L3 .GT. 60 ) THEN
                        MESG = MESG( 1:L1 ) // CRLF() // BLANK10 // 
     &                         ' "'// PBUF( 1:L2 ) // '"'
                        L3 = 0
                    ELSE
                        MESG = MESG( 1:L1 )// ' "'// PBUF( 1:L2 )// '"'
                    END IF
                END DO

                MESG = 'Processing pollutants:'//CRLF()//BLANK10//MESG
                CALL M3MSG2( MESG )

            ENDIF

C.............  Initialize various arrays that need to start as zero before 
C               time loop
! NOTE: Does anything need to go here?
 
! NOTE: Need to add the capability to output future year stuff based on
C current-year days of the week.  I think this needs to be handled in the new
C program that will merge the temporal emissions files in various ways and in
C the temporal programs themselves.  If the emissions file is an hourly file,
C then the output date will be from the dates in the file. If the emissions file
C is an inventory file, then the output date will be future year date from the
C FDESC3D (if it is there), otherwise, it will be the base-year date from the
C FDESC3D packet.  I should print a warning if INVYR has multiple values

C.............  Loop through output time steps
            JDATE = SDATE
            JTIME = STIME
            LDATE = 0
            DO T = 1, NSTEPS   ! at least once for time-independent

C................. For time-dependent processing, write out a few messages...
                IF( TFLAG ) THEN
                    
C.....................  Write out message for new day.  NOTE: For time-
C                       independent, LDATE and JDATE will both be zero.
                    IF( JDATE .NE. LDATE ) THEN

                        J = WKDAY( JDATE )
                        MESG = 'Processing ' // DAYS( J ) // 
     &                         MMDDYY( JDATE )
                        CALL M3MSG2( MESG )

                    ENDIF

C.....................  For new hour...
C.....................  Write to screen because WRITE3 only writes to LDEV
                    WRITE( *, 93020 ) HHMMSS( JTIME )

                END IF

C.................  If area sources, read inventory emissions for this time 
C                   step for all area-source pollutants in current pol group
                IF( AFLAG )
     &              CALL RD3MASK( ATNAME, JDATE, JTIME, NASRC, NPPGP,
     &                      PLNAMES( 1,N ), A_EXIST( 1,N ), AEMSRC   )

C.................  If mobile sources, read inventory emissions for this time 
C                   step for all area-source pollutants in current pol group
                IF( MFLAG ) 
     &              CALL RD3MASK( MTNAME, JDATE, JTIME, NMSRC, NPPGP, 
     &                      PLNAMES( 1,N ), M_EXIST( 1,N ), MEMSRC   )

C.................  If point sources, read inventory emissions for this time 
C                   step for all area-source pollutants in current pol group
                IF( PFLAG )
     &              CALL RD3MASK( PTNAME, JDATE, JTIME, NPSRC, NPPGP, 
     &                      PLNAMES( 1,N ), P_EXIST( 1,N ), PEMSRC   )

C.................  If layer fractions, read them for this time step
                IF( LFLAG ) THEN

                    IF( .NOT. READ3( PLNAME, 'LFRAC', ALLAYS3, 
     &                               JDATE, JTIME, LFRAC      ) ) THEN

                        MESG = 'Could not read LFRAC from ' // PLNAME
                        CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )

                    END IF   ! if read3() failed

                END IF

C.................  Initialize arrays for this time step that need to be zero


C.................  Loop through pollutants in pollutant group 
                DO V = 1, NPPGP

C.....................  If speciation, processing accordingly...
                    IF( SFLAG ) THEN

C.........................  Loop through pol-to-species for this pollutant/group
                        DO J = 1, NSMPPG( V,N )

C............................. Initialized gridded, merged emissions
                            TEMGRD = 0.  ! array

C.............................  For area sources, loop through sources, apply
C                               whatever matrices are valid, and store
                            IF( AFLAG ) THEN

                                K1 = A_EXIST ( V,N )
                                K2 = AU_EXIST( V,N )
                                K3 = AA_EXIST( V,N )
                                K4 = AS_EXIST( J,V,N )
                                K5 = NGRID + ANGMAT + 1
                                CALL MRGMULT( NASRC, NGRID, 1, ANGMAT, 
     &                               ANGMAT, K1, K2, K3, K4, 
     &                               AEMSRC, ACUMATX, ACAMATX, 
     &                               ASMATX, AGMATX(1), 
     &                               AGMATX(NGRID+1), AGMATX(K5), 
     &                               AEMGRD, TEMGRD )
                            END IF
                            
C.............................  For biogenic sources, read gridded emissions,
C                               add to totals and store
                            IF( BFLAG ) THEN
! NOTE: Fill in later
!                                K4 = BS_EXIST( J,V,N )
!                                CALL MRGBIO( NGRID, K4, BEMGRD, TEMGRD )

                            END IF
                            
C.............................  For mobile sources, loop through sources, apply
C                               whatever matrices are valid, and store
                            IF( MFLAG ) THEN

                                K1 = M_EXIST ( V,N )
                                K2 = MU_EXIST( V,N )
                                K3 = MA_EXIST( V,N )
                                K4 = MS_EXIST( J,V,N )
                                K5 = NGRID + MNGMAT + 1
                                CALL MRGMULT( NMSRC, NGRID, 1, MNGMAT,
     &                               MNGMAT, K1, K2, K3, K4, 
     &                               MEMSRC, MCUMATX, MCAMATX, 
     &                               MSMATX, MGMATX(1), 
     &                               MGMATX(NGRID+1), MGMATX(K5),  
     &                               MEMGRD, TEMGRD )
                            END IF
                            
C.............................  For point sources, loop through sources, apply
C                               whatever matrices are valid, and store
                            IF( PFLAG ) THEN

                                K1 = P_EXIST ( V,N )
                                K2 = PU_EXIST( V,N )
                                K3 = PA_EXIST( V,N )
                                K4 = PS_EXIST( J,V,N )
                                K5 = NGRID + NPSRC + 1
                                CALL MRGMULT( NPSRC, NGRID, EMLAYS, 
     &                               NPSRC, NPSRC, K1, K2, K3, K4,
     &                               PEMSRC, PCUMATX, PCAMATX, 
     &                               PSMATX, PGMATX(1), 
     &                               PGMATX(NGRID+1), PGMATX(K5), 
     &                               PEMGRD, TEMGRD )
                            END IF
                            
C.............................  Call subroutine to write speciated emissions
                            SBUF = EMNAM( SPINDEX( J,V,N ) )
                            CALL WRMRGGRD( SBUF, JDATE, JTIME )
            
                        END DO   ! End loop on pol-to-species

C.....................  If NOT speciation, processing accordingly...
C NOTE: Need to alloc speciation matrices even if not using because of EMMULT

                    ELSE

C......................... Initialized gridded, merged emissions
                        TEMGRD = 0.  ! array

C.........................  For area sources, loop through sources, apply
C                           whatever matrices are valid, and store
                        IF( AFLAG ) THEN

                            K1 = A_EXIST ( V,N )
                            K2 = AU_EXIST( V,N )
                            K3 = AA_EXIST( V,N )
                            K4 = 0
                            K5 = NGRID + ANGMAT + 1
                            CALL MRGMULT( NASRC, NGRID, 1, ANGMAT, 
     &                               ANGMAT, K1, K2, K3, K4,
     &                               AEMSRC, ACUMATX, ACAMATX,
     &                               ASMATX, AGMATX(1), 
     &                               AGMATX(NGRID+1), AGMATX(K5), 
     &                               AEMGRD, TEMGRD )
                        END IF
                            
C.........................  For mobile sources, loop through sources, apply
C                           whatever matrices are valid, and store
                        IF( MFLAG ) THEN

                            K1 = M_EXIST ( V,N )
                            K2 = MU_EXIST( V,N )
                            K3 = MA_EXIST( V,N )
                            K4 = 0
                            K5 = NGRID + MNGMAT + 1
                            CALL MRGMULT( NMSRC, NGRID, 1, MNGMAT,
     &                               MNGMAT, K1, K2, K3, K4,
     &                               MEMSRC, MCUMATX, MCAMATX, 
     &                               MSMATX, MGMATX(1), 
     &                               MGMATX(NGRID+1), MGMATX(K5),  
     &                               MEMGRD, TEMGRD )
                        END IF
                            
C.........................  For point sources, loop through sources, apply
C                           whatever matrices are valid, and store
                        IF( PFLAG ) THEN

                            K1 = P_EXIST ( V,N )
                            K2 = PU_EXIST( V,N )
                            K3 = PA_EXIST( V,N )
                            K4 = 0
                            K5 = NGRID + NPSRC + 1
                            CALL MRGMULT( NPSRC, NGRID, EMLAYS, NPSRC,
     &                               NPSRC, K1, K2, K3, K4,
     &                               PEMSRC, PCUMATX, PCAMATX, 
     &                               PSMATX, PGMATX(1), 
     &                               PGMATX(NGRID+1), PGMATX(K5), 
     &                               PEMGRD, TEMGRD )
                        END IF

C.........................  Call subroutine to write inventory data
                        CALL WRMRGGRD( PLNAMES( V,N ), JDATE, JTIME )

                    END IF  ! End speciation or not

                END DO      ! End loop on pollutants in group

                CALL NEXTIME( JDATE, JTIME, TSTEP )     !  update model clock

C.................  Write report data for previous day...
                IF( JDATE .NE. LDATE .OR. T .EQ. NSTEPS ) THEN
!                    CALL WRMRGREP( )
! Add report writing later
                END IF

                LDATE = JDATE

            END DO          ! End loop on time steps

        END DO   ! End of loop on pollutant/pol-to-spcs groups

999     CONTINUE          !  exit program

C.........  Successful completion of program
        CALL M3EXIT( PROGNAME, 0, 0, ' ', 0 )


C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93020   FORMAT( 8X, 'at time ', A8 )

C...........   Internal buffering formats............ 94xxx

94000   FORMAT( A )

94010   FORMAT( 10 ( A, :, I10, :, 2X ) )

94020   FORMAT( A, I4, 2X, 10 ( A, :, 1PG14.6, :, 2X ) )

94030   FORMAT( 8X, 'at time ', A8 )

        END PROGRAM SMKMERGE

