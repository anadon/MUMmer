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
#include "streemac.h"
#include "streetyp.h"
#include "construct.h"

/*@unused@*/ static void progresswithdot(/*@unused@*/ Uint nextstep,
        /*@unused@*/ void *info) {
    (void) putc('.',stderr);
    (void) fflush(stderr);
}

/*@unused@*/ static void finalprogress(/*@unused@*/ void *info) {
    (void) putc('\n',stderr);
}

/*EE
  The following function constructs the suffix tree.
*/

int main(int argc, char** argv){
    Uchar *text;
    size_t textlen;
		FILE *fd;
    Suffixtree stree;
    char *filename;


    CHECKARGNUM(2,"filename");
		
    fd = fopen(argv[1], "r");
		fseek(fd, 0, SEEK_END);
		textlen = ftell(fd);
		rewind(fd);
    text = (Uchar *) mmap(NULL, textlen, PROT_READ, MAP_FIXED, fileno(fd), 0);
    if(text == NULL) {
        fprintf(stderr,"%s: cannot open file \"%s\" ",argv[0],filename);
        fprintf(stderr,"or file \"%s\" is empty\n",filename);
        return EXIT_FAILURE;
    }
    if(textlen == 0) {
        fprintf(stderr,"%s: file \"%s\" is empty\n",argv[0],filename);
        return EXIT_FAILURE;
    }
    fprintf(stderr,"# construct suffix tree for sequence of length %lu\n",
    (Showuint) textlen);
		
    //if(constructprogressstree(&stree,text,textlen,NULL,NULL,NULL) != 0) {
    //    fprintf(stderr,"%s %s: %s\n",argv[0],filename,messagespace());
    //    return EXIT_FAILURE;
    //}
    /*
      addleafcountsstree(&stree);
    */
    if(munmap(text, textlen) != 0) {
        //STANDARDMESSAGE;
    }
		fclose(fd);
    freestree(&stree);
		
    return EXIT_SUCCESS;
}
