/*
  Copyright (c) 2003 by Stefan Kurtz and The Institute for
  Genomic Research.  This is OSI Certified Open Source Software.
  Please see the file LICENSE for licensing information and
  the file ACKNOWLEDGEMENTS for names of contributors to the
  code base.
*/

//}

//\FILEINFO{construct.c}{Stefan Kurtz}{November 1999}

//\Ignore{

#include "intbits.h"
#include "types.h"
#include "streemac.h"
#include "streeacc.h"
#include "streetyp.h"


#define VALIDINIT 0

//}

/*
 This file contains code for the improved linked list implementation,
 as described in \cite{KUR:1998,KUR:BAL:1999}. It can be
 compiled with two options:
 \begin{itemize}
 \item
 for short strings of length \(\leq 2^{21}-1=2\) megabytes,
 we recommend the option
 \texttt{STREESMALL}. This results in a representation which requires
 \(2\) integers for each small node, and three integers for each large node.
 See also the header file \texttt{streesmall.h}.
 \item
 for long strings of length \(\leq 2^{27}-1=12\) megabytes, we
 recommend the option
 \texttt{STREELARGE}. This results in a representation which requires
 \(2\) integers for each small node, and four integers for each large node.
 See also the header file \texttt{streelarge.h}.
 \end{itemize}
*/

//\subsection{Space Management}

/*
 For a string of length \(n\) we initially allocate space for
 \(\texttt{STARTFACTOR}\cdot\texttt{SMALLINTS}\cdot n\) integers to store
 the branching nodes. This usually suffices for most cases. In case we need
 more integers, we allocate space for \(\texttt{ADDFACTOR}\cdot n\)
 (at least 16) extra branching nodes.
*/

#ifndef STARTFACTOR
#define STARTFACTOR 0.5
#endif

#define ADDFACTOR   0.05
#define MINEXTRA    16

/*
 Before a new node is stored, we check if there is enough space available.
 If not, the space is enlarged by a small amount. Since some global pointers
 directly refer into the table, these have to be adjusted after reallocation.
*/

static void spaceforbranchtab(Suffixtree *stree) {

    if(stree->nextfreebranch >= stree->firstnotallocated) {
        Uint tmpheadnode,
             tmpchainstart = 0,
             extra = (Uint) (ADDFACTOR * (MULTBYSMALLINTS(stree->textlen+1)));
        if(extra < MINEXTRA) {
            extra = MULTBYSMALLINTS(MINEXTRA);
        }
        stree->currentbranchtabsize += extra;
        tmpheadnode = BRADDR2NUM(stree,stree->headnode);
        if(stree->chainstart != NULL) {
            tmpchainstart = BRADDR2NUM(stree,stree->chainstart);
        }

        stree->branchtab = (Uint*) malloc(sizeof(Uint) * stree->currentbranchtabsize);
        stree->nextfreebranch = stree->branchtab + stree->nextfreebranchnum;
        stree->headnode = stree->branchtab + tmpheadnode;
        if(stree->chainstart != NULL) {
            stree->chainstart = stree->branchtab + tmpchainstart;
        }
        stree->firstnotallocated
            = stree->branchtab + stree->currentbranchtabsize - LARGEINTS;
    }
}

//\subsection{Initializing and Retrieving Headpositions, Depth, and Suffixlinks}

/*
 We have three functions to initialize and retrieve head positions, depth, and
 suffix links. The implementation depends on the bit layout.
 \begin{enumerate}
 \item
 The function \emph{setdepthnum} stores the \emph{depth} and the
 \emph{head position} of a new large node.
 \item
 The function \emph{setsuffixlink} stores the \emph{suffixlink}
 of a new large node.
 \item
 The function \emph{getlargelinkconstruction} retrieves the \emph{suffixlink}
 of a large node, which is referenced by \emph{headposition}.
 \end{enumerate}
*/


static Uint getlargelinkconstruction(Suffixtree *stree) {
    SYMBOL secondchar;

    if(stree->headnodedepth == 1) {
        return 0;        // link refers to root
    }
    if(stree->headnodedepth == 2) { // determine second char of egde
        if(stree->headend == NULL) {
            secondchar = *(stree->tailptr-1);
        } else {
            secondchar = *(stree->tailptr - (stree->headend - stree->headstart + 2));
        }
        return stree->rootchildren[(Uint) secondchar];
    }
    return *(stree->headnode+4);
}

//\subsection{Insertion of Nodes}

/*
  The function \emph{insertleaf} inserts a leaf and a corresponding leaf
  edge outgoing from the current \emph{headnode}.
  \emph{insertprev} refers to the node to the left of the leaf to be inserted.
  If the leaf is the first child, then \emph{insertprev} is
  \texttt{NULL}.
*/

static void insertleaf(Suffixtree *stree) {
    Uint *ptr, newleaf;

    newleaf = MAKELEAF(stree->nextfreeleafnum);
    if(stree->headnodedepth == 0) {              // head is the root
        if(stree->tailptr != stree->sentinel) {    // no \$-edge initially
            stree->rootchildren[(Uint) *(stree->tailptr)] = newleaf;
            *(stree->nextfreeleafptr) = VALIDINIT;
        }
    } else {
        if (stree->insertprev == 0) { // newleaf = first child
            *(stree->nextfreeleafptr) = GETCHILD(stree->headnode);
            SETCHILD(stree->headnode,newleaf);
        } else {
            if(ISLEAF(stree->insertprev)) { // previous node is leaf
                ptr = stree->leaftab + GETLEAFINDEX(stree->insertprev);
                *(stree->nextfreeleafptr) = LEAFBROTHERVAL(*ptr);
                SETLEAFBROTHER(ptr,newleaf);
            } else { // previous node is branching node
                ptr = stree->branchtab + GETBRANCHINDEX(stree->insertprev);
                *(stree->nextfreeleafptr) = GETBROTHER(ptr);
                SETBROTHER(ptr,newleaf);
            }
        }
    }
    RECALLSUCC(newleaf);     // recall node on s.cppessor path of \emph{headnode}
    stree->nextfreeleafnum++;
    stree->nextfreeleafptr++;
}

/*
  The function \emph{insertbranch} inserts a branching node and splits
  the appropriate edges, .cppording to the canonical location of the current
  head. \emph{insertprev} refers to the node to the left of the branching
  node to be inserted. If the branching node is the first child, then
  \emph{insertprev} is \texttt{NULL}. The edge to be split ends
  in the node referred to by \emph{insertnode}.
*/

static void insertbranchnode(Suffixtree *stree) {
    Uint *ptr, *insertnodeptr, *insertleafptr, insertnodeptrbrother;

    spaceforbranchtab(stree);
    if(stree->headnodedepth == 0) {    // head is the root
        stree->rootchildren[(Uint) *(stree->headstart)]
            = MAKEBRANCHADDR(stree->nextfreebranchnum);
        *(stree->nextfreebranch+1) = VALIDINIT;
    } else {
        if(stree->insertprev == 0) { // new branch = first child
            SETCHILD(stree->headnode,MAKEBRANCHADDR(stree->nextfreebranchnum));
        } else {
            if(ISLEAF(stree->insertprev)) { // new branch = right brother of leaf
                ptr = stree->leaftab + GETLEAFINDEX(stree->insertprev);
                SETLEAFBROTHER(ptr,MAKEBRANCHADDR(stree->nextfreebranchnum));
            } else {                   // new branch = brother of branching node
                SETBROTHER(stree->branchtab + GETBRANCHINDEX(stree->insertprev),
                           MAKEBRANCHADDR(stree->nextfreebranchnum));
            }
        }
    }
    if(ISLEAF(stree->insertnode)) { // split edge is leaf edge
        insertleafptr = stree->leaftab + GETLEAFINDEX(stree->insertnode);
        if (stree->tailptr == stree->sentinel ||
                *(stree->headend+1) < *(stree->tailptr)) {
            SETNEWCHILDBROTHER(MAKELARGE(stree->insertnode),  // first child=oldleaf
                               LEAFBROTHERVAL(*insertleafptr));  // inherit brother
            RECALLNEWLEAFADDRESS(stree->nextfreeleafptr);
            SETLEAFBROTHER(insertleafptr,                     // new leaf =
                           MAKELEAF(stree->nextfreeleafnum)); // right brother of old leaf
        } else {
            SETNEWCHILDBROTHER(MAKELARGELEAF(stree->nextfreeleafnum),  // first child=new leaf
                               LEAFBROTHERVAL(*insertleafptr));  // inherit brother
            *(stree->nextfreeleafptr) = stree->insertnode;  // old leaf = right brother of of new leaf
            RECALLLEAFADDRESS(insertleafptr);
        }
    } else { // split edge leads to branching node
        insertnodeptr = stree->branchtab + GETBRANCHINDEX(stree->insertnode);
        insertnodeptrbrother = GETBROTHER(insertnodeptr);
        if (stree->tailptr == stree->sentinel ||
                *(stree->headend+1) < *(stree->tailptr)) {
            SETNEWCHILDBROTHER(MAKELARGE(stree->insertnode), // first child new branch
                               insertnodeptrbrother);        // inherit right brother
            RECALLNEWLEAFADDRESS(stree->nextfreeleafptr);
            SETBROTHER(insertnodeptr,MAKELEAF(stree->nextfreeleafnum)); // new leaf = brother of old branch
        } else {
            SETNEWCHILDBROTHER(MAKELARGELEAF(stree->nextfreeleafnum), // first child is new leaf
                               insertnodeptrbrother);        // inherit brother
            *(stree->nextfreeleafptr) = stree->insertnode;   // new branch is brother of new leaf
            RECALLBRANCHADDRESS(insertnodeptr);
        }
    }
    SETNILBIT;
    RECALLSUCC(MAKEBRANCHADDR(stree->nextfreebranchnum)); // node on s.cpp. path
    stree->currentdepth = stree->headnodedepth + (Uint) (stree->headend-stree->headstart+1);
    SETDEPTHEADPOS(stree->currentdepth,stree->nextfreeleafnum);
    SETMAXBRANCHDEPTH(stree->currentdepth);
    stree->nextfreeleafnum++;
    stree->nextfreeleafptr++;
}

//\subsection{Finding the Head-Locations}

/*
  The function \emph{rescan} finds the location of the current head.
  In order to scan down the tree, it suffices to look at the first
  character of each edge.
*/

static void rescan (Suffixtree *stree) {
    Uint *nodeptr, *largeptr = NULL, distance = 0, node, prevnode,
                    nodedepth, edgelen, wlen, leafindex, headposition;
    SYMBOL headchar, edgechar;

    if(stree->headnodedepth == 0) { // head is the root
        headchar = *(stree->headstart);  // headstart is assumed to be not empty
        node = stree->rootchildren[(Uint) headchar];
        if(ISLEAF(node)) { // stop if s.cppessor is leaf
            stree->insertnode = node;
            return;
        }
        nodeptr = stree->branchtab + GETBRANCHINDEX(node);
        GETONLYDEPTH(nodedepth,nodeptr);
        wlen = (Uint) (stree->headend - stree->headstart + 1);
        if(nodedepth > wlen) {  // cannot reach the s.cppessor node
            stree->insertnode = node;
            return;
        }
        stree->headnode = nodeptr;        // go to s.cppessor node
        stree->headnodedepth = nodedepth;
        if(nodedepth == wlen) {           // location has been scanned
            stree->headend = NULL;
            return;
        }
        (stree->headstart) += nodedepth;
    }
    while(True) { // \emph{headnode} is not the root
        headchar = *(stree->headstart);  // \emph{headstart} is assumed to be nonempty
        prevnode = 0;
        node = GETCHILD(stree->headnode);
        while(True) {           // traverse the list of s.cppessors
            if(ISLEAF(node)) { // s.cppessor is leaf
                leafindex = GETLEAFINDEX(node);
                edgechar = stree->text[stree->headnodedepth + leafindex];
                if(edgechar == headchar) {  // correct edge found
                    stree->insertnode = node;
                    stree->insertprev = prevnode;
                    return;
                }
                prevnode = node;
                node = LEAFBROTHERVAL(stree->leaftab[leafindex]);
            } else { // s.cppessor is branch node
                nodeptr = stree->branchtab + GETBRANCHINDEX(node);
                GETONLYHEADPOS(headposition,nodeptr);
                edgechar = stree->text[stree->headnodedepth + headposition];
                if(edgechar == headchar) { // correct edge found
                    break;
                }
                prevnode = node;
                node = GETBROTHER(nodeptr);
            }
        }

        GETDEPTHAFTERHEADPOS(nodedepth,nodeptr);     // get info about s.cpp node
        edgelen = nodedepth - stree->headnodedepth;
        wlen = (Uint) (stree->headend - stree->headstart + 1);
        if(edgelen > wlen) {   // cannot reach the s.cpp node
            stree->insertnode = node;
            stree->insertprev = prevnode;
            return;
        }
        stree->headnode = nodeptr;    // go to the s.cppessor node
        stree->headnodedepth = nodedepth;
        if(edgelen == wlen) {                  // location is found
            stree->headend = NULL;
            return;
        }
        (stree->headstart) += edgelen;
    }
}

/*
  The function \emph{taillcp} computes the length of the longest common prefix
  of two strings. The first string is between pointers \emph{start1} and
  \emph{end1}. The second string is the current tail, which is between  the
  pointers \emph{tailptr} and \emph{sentinel}.
*/

static Uint taillcp(Suffixtree *stree,SYMBOL *start1, SYMBOL *end1) {
    SYMBOL *ptr1 = start1, *ptr2 = stree->tailptr + 1;
    while(ptr1 <= end1 && ptr2 < stree->sentinel && *ptr1 == *ptr2) {
        ptr1++;
        ptr2++;
    }
    return (Uint) (ptr1-start1);
}

/*
 The function \emph{scanprefix} scans a prefix of the current tail
 down from a given node.
*/

static void scanprefix(Suffixtree *stree) {
    Uint *nodeptr = NULL, *largeptr = NULL, leafindex, nodedepth, edgelen, node,
          distance = 0, prevnode, prefixlen, headposition;
    SYMBOL *leftborder = (SYMBOL *) NULL, tailchar, edgechar = 0;

    if(stree->headnodedepth == 0) { // headnode is root
        if(stree->tailptr == stree->sentinel) { // there is no \$-edge
            stree->headend = NULL;
            return;
        }
        tailchar = *(stree->tailptr);
        if((node = stree->rootchildren[(Uint) tailchar]) == 0) {
            stree->headend = NULL;
            return;
        }
        if(ISLEAF(node)) { // s.cppessor edge is leaf, compare tail and leaf edge label
            leftborder = stree->text + GETLEAFINDEX(node);
            prefixlen = 1 + taillcp(stree,leftborder+1,stree->sentinel-1);
            (stree->tailptr) += prefixlen;
            stree->headstart = leftborder;
            stree->headend = leftborder + (prefixlen-1);
            stree->insertnode = node;
            return;
        }
        nodeptr = stree->branchtab + GETBRANCHINDEX(node);
        GETBOTH(nodedepth,headposition,nodeptr);  // get info for branch node
        leftborder = stree->text + headposition;
        prefixlen = 1 + taillcp(stree,leftborder+1,leftborder + nodedepth - 1);
        (stree->tailptr)+= prefixlen;
        if(nodedepth > prefixlen) { // cannot reach the s.cppessor, fall out of tree
            stree->headstart = leftborder;
            stree->headend = leftborder + (prefixlen-1);
            stree->insertnode = node;
            return;
        }
        stree->headnode = nodeptr;
        stree->headnodedepth = nodedepth;
    }
    while(True) { // \emph{headnode} is not the root
        prevnode = 0;
        node = GETCHILD(stree->headnode);
        if(stree->tailptr == stree->sentinel) { //  \$-edge
            do { // there is no \$-edge, so find last s.cppessor, of which it becomes right brother
                prevnode = node;
                if(ISLEAF(node)) {
                    node = LEAFBROTHERVAL(stree->leaftab[GETLEAFINDEX(node)]);
                } else {
                    node = GETBROTHER(stree->branchtab + GETBRANCHINDEX(node));
                }
            } while(!NILPTR(node));
            stree->insertnode = NILBIT;
            stree->insertprev = prevnode;
            stree->headend = NULL;
            return;
        }
        tailchar = *(stree->tailptr);

        do { // find s.cppessor edge with firstchar = tailchar
            if(ISLEAF(node)) { // s.cppessor is leaf
                leafindex = GETLEAFINDEX(node);
                leftborder = stree->text + (stree->headnodedepth + leafindex);
                if((edgechar = *leftborder) >= tailchar) { // edge will not come later
                    break;
                }
                prevnode = node;
                node = LEAFBROTHERVAL(stree->leaftab[leafindex]);
            } else { // s.cppessor is branch node
                nodeptr = stree->branchtab + GETBRANCHINDEX(node);
                GETONLYHEADPOS(headposition,nodeptr);
                leftborder = stree->text + (stree->headnodedepth + headposition);
                if((edgechar = *leftborder) >= tailchar) { // edge will not come later
                    break;
                }
                prevnode = node;
                node = GETBROTHER(nodeptr);
            }
        } while(!NILPTR(node));
        if(NILPTR(node) || edgechar > tailchar) { // edge not found
            stree->insertprev = prevnode;   // new edge will become brother of this
            stree->headend = NULL;
            return;
        }
        if(ISLEAF(node)) { // correct edge is leaf edge, compare its label with tail
            prefixlen = 1 + taillcp(stree,leftborder+1,stree->sentinel-1);
            (stree->tailptr) += prefixlen;
            stree->headstart = leftborder;
            stree->headend = leftborder + (prefixlen-1);
            stree->insertnode = node;
            stree->insertprev = prevnode;
            return;
        }
        GETDEPTHAFTERHEADPOS(nodedepth,nodeptr); // we already know headposition
        edgelen = nodedepth - stree->headnodedepth;
        prefixlen = 1 + taillcp(stree,leftborder+1,leftborder + edgelen - 1);
        (stree->tailptr) += prefixlen;
        if(edgelen > prefixlen) { // cannot reach next node
            stree->headstart = leftborder;
            stree->headend = leftborder + (prefixlen-1);
            stree->insertnode = node;
            stree->insertprev = prevnode;
            return;
        }
        stree->headnode = nodeptr;
        stree->headnodedepth = nodedepth;
    }
}

//\subsection{Completion and Initialization}

/*
  The function \emph{completelarge} is called whenever a large node
  is inserted. It basically sets the appropriate distance values of the small
  nodes of the current chain.
*/

static void completelarge(Suffixtree *stree) {
    Uint distance, *backwards;

    if(stree->smallnotcompleted > 0) {
        backwards = stree->nextfreebranch;
        for(distance = 1; distance <= stree->smallnotcompleted; distance++) {
            backwards -= SMALLINTS;
            SETDISTANCE(backwards,distance);
        }
        stree->smallnotcompleted = 0;
        stree->chainstart = NULL;
    }
    stree->nextfreebranch += LARGEINTS;
    stree->nextfreebranchnum += LARGEINTS;
    stree->largenode++;
}

/*
  The function \emph{linkrootchildren} constructs the s.cppessor chain
  for the children of the root. This is done at the end of the algorithm
  in one sweep over table \emph{rootchildren}.
*/

static void linkrootchildren(Suffixtree *stree) {
    Uint *rcptr, *prevnodeptr, prev = 0;

    stree->alphasize = 0;
    for(rcptr = stree->rootchildren;
            rcptr <= stree->rootchildren + LARGESTCHARINDEX; rcptr++) {
        if(*rcptr != 0) {
            stree->alphasize++;
            if(prev == 0) {
                SETCHILD(stree->branchtab,MAKELARGE(*rcptr));
            } else {
                if(ISLEAF(prev)) {
                    stree->leaftab[GETLEAFINDEX(prev)] = *rcptr;
                } else {
                    prevnodeptr = stree->branchtab + GETBRANCHINDEX(prev);
                    SETBROTHER(prevnodeptr,*rcptr);
                }
            }
            prev = *rcptr;
        }
    }
    if(ISLEAF(prev)) {
        stree->leaftab[GETLEAFINDEX(prev)] = MAKELEAF(stree->textlen);
    } else {
        prevnodeptr = stree->branchtab + GETBRANCHINDEX(prev);
        SETBROTHER(prevnodeptr,MAKELEAF(stree->textlen));
    }
    stree->leaftab[stree->textlen] = NILBIT;
}

/*
  \newpage
  \emph{initSuffixtree} allocates and initializes the data structures for
  McCreight's Algorithm.
*/

static void initSuffixtree(Suffixtree *stree,SYMBOL *text,Uint textlen) {
    Uint i, *ptr;

    stree->currentbranchtabsize
        = (Uint) (STARTFACTOR * MULTBYSMALLINTS(textlen+1));
    if(stree->currentbranchtabsize < MINEXTRA) {
        stree->currentbranchtabsize = MULTBYSMALLINTS(MINEXTRA);
    }
    stree->leaftab = (Uint*) malloc(sizeof(Uint) * (textlen+2) );
    stree->rootchildren = (Uint*) malloc(sizeof(Uint) * (LARGESTCHARINDEX + 1) );
    stree->branchtab = (Uint*) malloc(sizeof(Uint) * stree->currentbranchtabsize);

    stree->text = stree->tailptr = text;
    stree->textlen = textlen;
    stree->sentinel = text + textlen;
    stree->firstnotallocated = stree->branchtab + stree->currentbranchtabsize - LARGEINTS;
    stree->headnode = stree->nextfreebranch = stree->branchtab;
    stree->headend = NULL;
    stree->headnodedepth = stree->maxbranchdepth = 0;
    for(ptr=stree->rootchildren; ptr<=stree->rootchildren+LARGESTCHARINDEX; ptr++) {
        *ptr = 0;
    }
    for(i=0; i<LARGEINTS; i++) {
        stree->branchtab[i] = 0;
    }
    stree->nextfreebranch = stree->branchtab;
    stree->nextfreebranchnum = 0;
    SETDEPTHEADPOS(0,0);
    SETNEWCHILDBROTHER(MAKELARGELEAF(0),0);
    SETBRANCHNODEOFFSET;
    stree->rootchildren[(Uint) *text] = MAKELEAF(0);
    stree->leaftab[0] = VALIDINIT;
    stree->leafcounts = NULL;
    stree->nextfreeleafnum = 1;
    stree->nextfreeleafptr = stree->leaftab + 1;
    stree->nextfreebranch = stree->branchtab + LARGEINTS;
    stree->nextfreebranchnum = LARGEINTS;
    stree->insertnode = stree->insertprev = 0;
    stree->smallnotcompleted = 0;
    stree->chainstart = NULL;
    stree->largenode = stree->smallnode = 0;

}

void freestree(Suffixtree *stree) {
    free(stree->leaftab);
    free(stree->rootchildren);
    free(stree->branchtab);
    if(stree->nonmaximal != NULL)
        free(stree->nonmaximal);
				
    if(stree->leafcounts != NULL)
        free(stree->leafcounts);
				
}

//\subsection{Computing the Suffix Tree}

/*
  \emph{constructstree} implements McCreight Algorithm to compute the suffix
  tree for a \texttt{text} of length \texttt{textlen}. For explanations, see
  \cite{KUR:1998}. The number \((i)\) refers to the cases of Section 6 in
  \cite{KUR:1998}.
*/

