/*
  Copyright (c) 2003 by Stefan Kurtz and The Institute for
  Genomic Research.  This is OSI Certified Open Source Software.
  Please see the file LICENSE for licensing information and
  the file ACKNOWLEDGEMENTS for names of contributors to the
  code base.
*/

//\IgnoreLatex{

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "types.h"
#include "megabytes.h"
#include "maxmatdef.h"


//}

/*EE
  This module contains the main function of maxmatch3. It calls
  the following three functions in an appropriate order and with
  proper arguments.
*/

/*EE
  The following function is imported form \texttt{maxmatopt.c}.
*/

Sint parsemaxmatoptions (MMcallinfo *maxmatcallinfo,
                         Argctype argc,
                         char **argv);

/*EE
  The following function is imported form \texttt{maxmatinp.c}.
*/

Sint getmaxmatinput (Multiseq *subjectmultiseq,
                     bool matchnucleotidesonly,
                     char *subjectfile);

/*EE
  The following function is imported form \texttt{procmaxmat.c}.
*/

Sint procmaxmatches(MMcallinfo *mmcallinfo,
                    Multiseq *subjectmultiseq);

//\IgnoreLatex{

/*
  This is the main function.
*/

MAINFUNCTION {
    Sint retcode;
    MMcallinfo mmcallinfo;
    Multiseq subjectmultiseq;

    initclock();
    retcode = parsemaxmatoptions (&mmcallinfo, argc, argv);
    if (retcode < 0) {
        STANDARDMESSAGE;  // return error code and show message
    }
    if (retcode == 1) { // program was called with option -help
        return EXIT_SUCCESS;
    }
		
    if (getmaxmatinput (&subjectmultiseq,
    mmcallinfo.matchnucleotidesonly,
    &mmcallinfo.subjectfile[0]) != 0) {
        STANDARDMESSAGE;
    }
    if(procmaxmatches(&mmcallinfo,&subjectmultiseq) != 0) {
        STANDARDMESSAGE;
    }
    freemultiseq (&subjectmultiseq);
    fprintf(stderr,"# COMPLETETIME %s %s %.2f\n",
    argv[0],&mmcallinfo.subjectfile[0],
    getruntime());
    fprintf(stderr,"# SPACE %s %s %.2f\n",argv[0],
    &mmcallinfo.subjectfile[0],
    MEGABYTES(getspacepeak()+mmgetspacepeak()));
    return EXIT_SUCCESS;
}

//}
