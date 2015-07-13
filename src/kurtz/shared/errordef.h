/*
  Copyright (c) 2003 by Stefan Kurtz and The Institute for
  Genomic Research.  This is OSI Certified Open Source Software.
  Please see the file LICENSE for licensing information and
  the file ACKNOWLEDGEMENTS for names of contributors to the
  code base.
*/

//\Ignore{

#ifndef ERRORDEF_H
#define ERRORDEF_H
#include <stdio.h>
#include <stdlib.h>
#include "types.h"


char *messagespace(void);
Sint maxerrormsg(void);


//}

/*
  This file contains some macros to write error messages into a
  buffer returned by the function \texttt{messagespace}.
*/


/*
  The following is a macro to show the usage line for all programs
  which have options and an indexname as the last argument.
*/

#define USAGEOUT\
    printf("Usage: %s options indexname\n"\
         "%s -help shows possible options",\
        argv[0],argv[0]);

/*
  The following is the standard message in the main function. It shows the
  program name and the error message as returned by the function
  \texttt{messagespace}.
*/

#define STANDARDMESSAGE\
    fprintf(stderr,"%s: %s\n",argv[0],messagespace());\
    return EXIT_FAILURE

#define SIMPLESTANDARDMESSAGE\
    fprintf(stderr,"%s\n",messagespace());\
    return EXIT_FAILURE

/*
  The following is a generell error message which leads to a termination
  of the program.
*/

#define NOTSUPPOSED\
    fprintf(stderr,"%s: line %lu: This case is not supposed to .cppur\n",\
             __FILE__,(Showuint) __LINE__);\
    exit(EXIT_FAILURE)


#endif

//}
