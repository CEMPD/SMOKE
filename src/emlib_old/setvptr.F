
        SUBROUTINE SETVPTR( MXFNAM, NDINDESC, CASEDESC, MXDREF, NDREF,
     &                      NDIN, VNAMREF, VNAMIN, IVPTR )

C***********************************************************************
C  subroutine body starts at line 79
C
C  DESCRIPTION:
C     Subroutine SETVPTR compares a list of names to another list of
C     names, and outputs a pointer from the input list to the reference
C     list.  It also checks the actual dimension to the maximum dimension
C
C  PRECONDITIONS REQUIRED:  
C     Correct calling parameters
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       CRLF, INDEX1, TRIMLEN
C
C  REVISION  HISTORY:
C       Initial code  8/98 by M Houyoux
CC
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
C*************************************************************************

C.........  Subroutine paramaters and their descriptions
        CHARACTER*8   MXFNAM           ! name of fortran variable from caller
        CHARACTER*256 NDINDESC         ! description of NDIN
        CHARACTER*256 CASEDESC         ! desc of the file/var being checked
        INTEGER       MXDREF           ! maximum number of reference names
        INTEGER       NDREF            ! actual number of reference names
        INTEGER       NDIN             ! actual number of test names
        CHARACTER*(*) VNAMREF( NDREF ) ! reference names
        CHARACTER*(*) VNAMIN ( NDIN )  ! test names
        INTEGER       IVPTR  ( NDIN )  ! output pointer from test to reference

C.........  External functions
        CHARACTER*2   CRLF
        INTEGER       INDEX1
        INTEGER       TRIMLEN

        EXTERNAL      CRLF, INDEX1, TRIMLEN

C.........  Local parameters
        CHARACTER*5   BLANK5
        PARAMETER   ( BLANK5 = ' ' )

C.........  Local variables
        INTEGER       L, L1, L2, L3
        INTEGER       J, V

        LOGICAL       EFLAG

        CHARACTER*16  VBUF       
        CHARACTER*256 MESG       

C***********************************************************************
C   begin body of subroutine SETVPTR

         L1 = TRIMLEN( MXFNAM )
         L2 = TRIMLEN( NDINDESC )
         L3 = TRIMLEN( CASEDESC )

        IF( NDIN .GT. MXDREF ) THEN

            WRITE( MESG,94010 ) 'Dimension mismatch for number of ' // 
     &             NDINDESC( 1:L2 ) // 's' //
     &             CRLF() // BLANK5 // CASEDESC( 1:L3 ) // ':',
     &             NDIN, 'program (' // MXFNAM( 1:L1 ) // '):',
     &             MXDREF
            CALL M3EXIT( 'SETVPTR', 0, 0, MESG, 2 )

        ENDIF

        EFLAG = .FALSE.
        DO 11 V = 1, NDIN

            VBUF = VNAMIN( V )
            L = TRIMLEN( VBUF )
            J = INDEX1 ( VBUF( 1:L ), NDREF, VNAMREF )

            IF( J .GT. 0 ) THEN
               IVPTR( V ) = J 

            ELSE
                EFLAG = .TRUE.
                MESG = 'Error: ' // NDINDESC // ' "' // VBUF( 1:L ) // 
     &                 '" not found in available list'
                CALL M3MSG2( MESG )

            ENDIF

11      CONTINUE

        IF( EFLAG ) THEN
            MESG = 'Error checking ' // NDINDESC // 's.'
            CALL M3EXIT( 'SETVPTR', 0, 0, MESG, 2 )
        ENDIF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats.............94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END
