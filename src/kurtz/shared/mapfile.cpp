/*
  Copyright (c) 2003 by Stefan Kurtz and The Institute for
  Genomic Research.  This is OSI Certified Open Source Software.
  Please see the file LICENSE for licensing information and
  the file ACKNOWLEDGEMENTS for names of contributors to the
  code base.
*/

//\Ignore{

#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/uio.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <unistd.h>
#include "types.h"
#include "megabytes.h"

#define MAXMAPPEDFILES 32

//}

/*EE
  This file contains functions to map files to memory, to store pointers to
  the corresponding secondary memory, and to maintain the size of the mapped
  memory. */

/*
  Pointers to the memory areas, initialized to NULL, i.e. not occupied.
*/

static void *memoryptr[MAXMAPPEDFILES] = {NULL};

static Uint currentspace = 0,              // currently mapped num of bytes
            spacepeak = 0,                 // maximally mapped num of bytes
            mappedbytes[MAXMAPPEDFILES] = {0};  // size of the memory map

/*
  The following two tables store important information to
  generate meaningfull error messages.
*/


/*
  file where the mmap was created
*/

static char *filemapped[MAXMAPPEDFILES] = {NULL};

/*
  line where the mmap was created
*/

static Uint linemapped[MAXMAPPEDFILES] = {0};

/*
  The following two functions \texttt{mmaddspace} and \texttt{mmsubtractspace}
  maintain the variables \texttt{currentspace} and \texttt{spacepeak}.
*/



/*EE
  The following function opens a file, and stores the size of the
  file in \texttt{numofbytes}. It returns the file descriptor of the
  file, or a negative error code if something went wrong.
*/

int simplefileOpen(char *filename) {
    return open(filename,O_RDONLY);
}

/*EE
  The following function creates a memory map for a given file
  descriptor \texttt{fd}. \texttt{writemap} is true iff the map
  should be writable, and \texttt{numofbytes} is the
  size of the file to be mapped.
*/

/*EE
  The following function returns a memory map for a given filename, or
  \texttt{NULL} if something went wrong.
*/

void *creatememorymap(char *filename, bool writemap, size_t numofbytes) {
    return mmap(0, numofbytes, writemap ? (PROT_READ | PROT_WRITE) : PROT_READ,
                        MAP_PRIVATE, open(filename,O_RDONLY), (off_t) 0);
}

/*EE
  The following function unmaps the memorymap referenced by
  \texttt{mappedfile}. It returns a negative value if the
  \texttt{mappedfile} is \texttt{NULL}, or the corresponding
  filedescriptor cannot be found, or the \texttt{munmap} operation
  fails.
*/

Sint deletememorymap(char *file,Uint line,void *mappedfile) {
    Filedesctype fd;

    if(mappedfile == NULL) {
        return -1;
    }
    for(fd=0; fd<MAXMAPPEDFILES; fd++) {
        if(memoryptr[fd] == mappedfile) {
            break;
        }
    }
    if(fd == MAXMAPPEDFILES) {
        return -2;
    }
    if(munmap(memoryptr[fd],(size_t) mappedbytes[fd]) != 0) {
        return -3;
    }
    memoryptr[fd] = NULL;
    mappedbytes[fd] = 0;
    if(close(fd) != 0) {
        return -4;
    }
    return 0;
}


/*EE
  The following function frees the space for all memory maps
  which have not already been unmapped.
*/

Sint mmwrapspace(void) {
    Filedesctype fd;

    for(fd=0; fd<MAXMAPPEDFILES; fd++) {
        if(memoryptr[fd] != NULL) {
            if(munmap(memoryptr[fd],(size_t) mappedbytes[fd]) != 0) {
                return -1;
            }
            memoryptr[fd] = NULL;
            if(close(fd) != 0) {
                return -4;
            }
            memoryptr[fd] = NULL;
            mappedbytes[fd] = 0;
            filemapped[fd] = NULL;
            linemapped[fd] = 0;
        }
    }
    return 0;
}

/*EE
  The following function shows the space peak in megabytes on \texttt{stderr}.
*/

void mmshowspace(void) {
    fprintf(stderr,"# mmap space peak in megabytes: %.2f\n",MEGABYTES(spacepeak));
}

/*EE
  The following function returns the space peak in bytes.
*/

Uint mmgetspacepeak(void) {
    return spacepeak;
}
