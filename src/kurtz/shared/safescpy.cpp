/*
  Copyright (c) 2003 by Stefan Kurtz and The Institute for
  Genomic Research.  This is OSI Certified Open Source Software.
  Please see the file LICENSE for licensing information and
  the file ACKNOWLEDGEMENTS for names of contributors to the
  code base.
*/

//\Ignore{

#include <string.h>
#include "errordef.h"

//}

/*EE
  The following function copies the 0-terminated string pointed to by
  \texttt{source} to the memory area pointed to by \texttt{dest}, provided
  \texttt{source} is shorter than \texttt{maxlen}. If this is not
  true, then the function returns a negative error code. Otherwise
  the return value is 0.
*/
