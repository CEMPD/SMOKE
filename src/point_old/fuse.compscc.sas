OPTIONS LINESIZE=110;
/* This program takes two different years of fuel use data with OTAG id
   info and computes growth factors based on the fuel use data.  The
   growth factors are by stid cyid pltid and scc.
 
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
stop;
RUN;

/* read in fuel use data for a particular year with OTAG id info */

DATA fu90 (keep= stid cyid pltid dvid scc may jun jul aug sep);

   INFILE "&filnm" DLM=',';
   INPUT   stid 2.
           cyid 3. 
           pltid $15.
           dvid $12.
           scc   8.
           may  10.2
           jun  10.2
           jul  10.2
           aug  10.2
           sep  10.2;
RETURN;

PROC SORT DATA=fu90;
BY stid cyid pltid scc dvid;
RUN;


/* read in fuel use data for another particular year with OTAG id info */

DATA fu95 (keep= stid cyid pltid dvid scc may5 jun5 jul5 aug5 sep5);

   INFILE "&filnm2" DLM=',';
   INPUT   stid 2.
           cyid 3.  
           pltid $15.
           dvid $12.
           scc   8. 
           may5  10.2
           jun5  10.2
           jul5  10.2
           aug5  10.2
           sep5  10.2;
RETURN;


PROC SORT DATA=fu95;
BY stid cyid pltid scc dvid ;
RUN;

/* merge the two years and find remaining matches.  Must have
   data for both years to be a match */

DATA fu95x ;
MERGE fu95 (IN=a) fu90 (IN=b);
   by stid cyid pltid scc dvid;
   if a and not b then delete;
   if b and not a then delete;
RUN;


/* sum fuel use by month and scc */

PROC SUMMARY DATA=fu95x;
BY stid cyid pltid scc ;
VAR may jun jul aug sep may5 jun5 jul5 aug5 sep5;
OUTPUT out=total sum=may jun jul aug sep may5 jun5 jul5 aug5 sep5;
RUN;

PROC SORT DATA=total;
BY stid cyid pltid scc ;
RUN;

/* calculate growth factors */

DATA total ;
set total;

   sum95 = jun5 + jul5 + aug5;
   sum90 = jun + jul + aug;
   avg95 = sum95/3;
   maygr = (sum95/sum90) * (may5/avg95);
   jungr = (sum95/sum90) * (jun5/avg95); 
   julgr = (sum95/sum90) * (jul5/avg95); 
   auggr = (sum95/sum90) * (aug5/avg95); 
   sepgr = (sum95/sum90) * (sep5/avg95); 

RUN;



PROC SORT DATA=total (keep = stid cyid pltid scc maygr jungr julgr
  auggr sepgr);
BY stid cyid pltid scc ;
RUN;

/* print out growth factors */

DATA _null_;
 SET total END=EOF;
 FILE print   NOTITLES;

  PUT   @1 stid 2.
        @3 cyid 3.
        @6 pltid $15.
        @21 scc   8.
        @29 maygr  10.5
        @39 jungr  10.5
        @49 julgr  10.5
        @59 auggr  10.5
        @69 sepgr  10.5;
RETURN;

 
