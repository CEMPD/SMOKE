OPTIONS LINESIZE=110;
/* This program calculates projection factors based on a set of
   growth factors and control info by device and scc.  Since the
   projection factors are based on stid cyid pltid and scc.  The
   device control factors are used by implementing a fuel-use
   weighted factor for each device at a pltid.  

Prototype : Jeff Vukovich 10/97
*/
/*
C***************************************************************************
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
*/


/* read in filenames */

DATA _NULL_;
CALL SYMPUT ('filnm',((SCAN("&sysparm",1,";"))));
CALL SYMPUT ('filnm2',((SCAN("&sysparm",2,";"))));
CALL SYMPUT ('filnm3',((SCAN("&sysparm",3,";"))));
CALL SYMPUT ('filnm4',((SCAN("&sysparm",4,";"))));
stop;
RUN;

/* read control data by scc */

DATA scc (keep=stid cyid scc cescc );

   INFILE "&filnm" DLM=',';
   INPUT  stid   2.       /* STID  code */
          cyid   3.     /* CYID  code */
          scc    $8.     /* Round3 zone*/
          cescc   7.3;
RUN;

PROC SORT DATA=scc;
BY stid cyid scc;
RUN;

/* read control data by device */

DATA device (keep=stid cyid pltid scc dvid cedev);

   INFILE "&filnm2" DLM=',';
   INPUT  stid   2.       /* STID  code */
          cyid   3.     /* CYID  code */
          pltid    $15.     /* Round3 zone*/
          dvid     $15.     /* Round3 zone*/
          stkid    $12.     /* Round3 zone*/
          oris      $6.     /* Round3 zone*/
          blrid     $5.     /* Round3 zone*/
          seg       $2.     /* Round3 zone*/
          plt      $40.     /* Round3 zone*/
          plt2     $40.     /* Round3 zone*/
          scc      $8.     /* Round3 zone*/
          cedev    9.3;    /* Round3 zone*/
RUN;

PROC SORT DATA=device;
BY stid cyid pltid scc dvid ;
RUN;

/* read in growth factors */

DATA grw (keep=stid cyid pltid scc maygr jungr julgr auggr sepgr)
;

   INFILE "&filnm3" DLM=',';
   INPUT  stid   2.       /* STID  code */
          cyid   3.     /* CYID  code */
          pltid    $15.     /* Round3 zone*/
          scc      $8.     /* Round3 zone*/
          maygr   10.5
          jungr   10.5
          julgr   10.5
          auggr   10.5
          sepgr   10.5 ;
RUN;

PROC SORT DATA= grw;
BY stid cyid scc;
RUN;

/* read in the projection year's fuel use for fuel use weighting */

DATA fu95 (keep= stid cyid pltid blrid dvid scc may5 jun5 jul5 aug5 sep5);

   INFILE "&filnm4" DLM=',';
   INPUT   stid 2.
           cyid 3.
           pltid $15.
           blrid $5.
           dvid  $12.
           scc   $8.
           may5   9.2
           jun5   9.2
           jul5   9.2
           aug5   9.2
           sep5   9.2;
RETURN;


PROC SORT DATA=fu95;
BY stid cyid pltid scc dvid ;
RUN;

/* calculate fuel use weighting factors */

DATA fu95  ;
 MERGE fu95 (IN=a) device (IN=b);
 BY stid cyid pltid scc dvid ;
 if may5 = '.' then delete;
 tdv = may5 + jun5 + jul5 + aug5 + sep5;
RUN;

DATA fu95;
 set fu95;
 if cedev = '.' then delete;
RUN;


PROC SUMMARY DATA=fu95;
BY stid cyid pltid scc ;
VAR tdv;
OUTPUT out=totplt sum=tfuse;
RUN;

/* apply weighting factors */

DATA cewts  ;
 MERGE totplt (IN=a) fu95 (IN=b);
 BY stid cyid pltid scc;
 wfac = tdv/tfuse;
 cewtd = cedev * wfac;
RUN;


PROC SORT DATA=cewts;
BY stid cyid pltid scc;
RUN;


PROC SUMMARY DATA=cewts;
BY stid cyid pltid scc;
VAR cewtd;
OUTPUT out = cefuwtd sum= ce_fuwt;
RUN;
 
PROC SORT DATA=cefuwtd;
BY stid cyid pltid scc;
RUN;

/* merge growth and scc controls */ 

DATA fu95y  ;
 MERGE grw (IN=a) scc (IN=b);
 by stid cyid scc;
 if maygr eq '.' then delete;
RUN;

PROC SORT DATA=fu95y;
by stid cyid pltid scc;
RUN;

/* merge growth scc and fuel use weighted info */
 
DATA fu95z(keep = stid cyid pltid scc maygc jungc julgc auggc sepgc ce cescc ce_fuwt)   ;
 MERGE cefuwtd (IN=a) fu95y (IN=b);
 by stid cyid pltid scc;
 ce = ce_fuwt;
 if ce_fuwt eq '.' then ce = cescc;
 if ce = '.' then ce = 0.0;
 maygc = maygr * (1 - (ce/100));
 jungc = jungr * (1 - (ce/100));
 julgc = julgr * (1 - (ce/100));
 auggc = auggr * (1 - (ce/100));
 sepgc = sepgr * (1 - (ce/100));
 if maygr eq '.' then delete;
RUN;


PROC SORT DATA=fu95z NODUPREC FORCE;
BY stid cyid pltid scc maygc jungc julgc auggc sepgc;
RUN;

DATA _null_;
 SET fu95z END=EOF;
 FILE print   NOTITLES;

  PUT   @1 stid 2.
        @3 cyid 3.
        @6 pltid $15.
        @21 scc   8.
        @29 maygc   8.4
        @37 jungc   8.4
        @45 julgc   8.4
        @53 auggc   8.4
        @61 sepgc   8.4
        @70 ce      5.2
        @76 cescc   5.2
        @82 ce_fuwt 5.2; 
RETURN;


 
