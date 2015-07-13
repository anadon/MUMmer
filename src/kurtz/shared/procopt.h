#ifndef PROCMAXMAT_H
#define PROCMAXMAT_H

void showoptionswithoutexclude(FILE *outfp,char *program,
                               OptionDescription *opt,
                               Sint *excludetab,Uint numofopt);

Sint addoption(OptionDescription *options,Uint numofoptions,
               Uint optnum,char *optname,char *optdesc);
							 
Sint procoption(OptionDescription *opt,Uint numofopt,char *optstring);

#endif