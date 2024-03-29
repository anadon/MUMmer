/*
  Copyright (c) 2003 by Stefan Kurtz and The Institute for
  Genomic Research.  This is OSI Certified Open Source Software.
  Please see the file LICENSE for licensing information and
  the file ACKNOWLEDGEMENTS for names of contributors to the
  code base.
*/

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/mman.h>
#include "args.h"
#include "types.h"
#include "types.h"
#include "streemac.h"
#include "streeacc.h"

MAINFUNCTION {
    Uchar *text;
    Uint textlen;
    Suffixtree stree;


    CHECKARGNUM(2,"filename");
    text = (Uchar *) CREATEMEMORYMAP(argv[1],False,&textlen);
    if(text == NULL) {
        STANDARDMESSAGE;
    }
    CONSTRUCTSTREE(&stree,text,textlen,return EXIT_FAILURE);
		
    freestree(&stree);
    if(munmap(text, 0) != 0) {
        STANDARDMESSAGE;
    }
    return EXIT_SUCCESS;
}
