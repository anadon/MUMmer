/*
  Copyright (c) 2003 by Stefan Kurtz and The Institute for
  Genomic Research.  This is OSI Certified Open Source Software.
  Please see the file LICENSE for licensing information and
  the file ACKNOWLEDGEMENTS for names of contributors to the
  code base.
*/

#include "types.h"
#include "streemac.h"
#include "streeacc.h"
#include "streetyp.h"

bool exactlytwoleavesstree(Suffixtree *stree,PairUint *twoleaves,Bref start) {
    bool firstleaffound = False;
    Uint tmpval, node;

    node = GETCHILD(start);
    while(True) {
        if(ISLEAF(node)) {
            if(firstleaffound) {
                twoleaves->uint1 = GETLEAFINDEX(node);
                if(twoleaves->uint0 > twoleaves->uint1) {
                    tmpval = twoleaves->uint1;
                    twoleaves->uint1 = twoleaves->uint0;
                    twoleaves->uint0 = tmpval;
                }
								
                if(NILPTR(LEAFBROTHERVAL(stree->leaftab[GETLEAFINDEX(node)]))) {
                    return True;
                }
                return False;
            }
            twoleaves->uint0 = GETLEAFINDEX(node);
            firstleaffound = True;
            node = LEAFBROTHERVAL(stree->leaftab[GETLEAFINDEX(node)]);
        } else {
            return False;
        }
    }
}
