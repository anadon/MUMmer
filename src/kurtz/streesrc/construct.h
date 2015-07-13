#ifndef STREEFILEDOC_H
#define STREEFILEDOC_H

Sint depthfirststree(Suffixtree *stree,Reference *startnode,
                     Sint (*processleaf)(Uint,Bref,void *),
                     bool (*processbranch1)(Bref,void *),
                     Sint (*processbranch2)(Bref,void *),
                     bool (*stoptraversal)(void *),
                     void *stopinfo,void *info);
										 
Sint constructprogressstree(Suffixtree *stree,SYMBOL *text,
                            Uint textlen,
                            void (*progress)(Uint,void *),
                            void (*finalprogress)(void *),
                            void *info);

void freestree(Suffixtree *stree);

#endif