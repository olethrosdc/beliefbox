#ifndef DEBUG_H
#define DEBUG_H

#define Serror   printf("# ERROR (%s - %s:  %d) ", __FUNCTION__,__FILE__, __LINE__); printf
#define Swarning  printf("# WARNING (%s - %s:  %d) ", __FUNCTION__,__FILE__, __LINE__); printf
//#define logmsg (void)
//#define dbgmsg (void)
#define logmsg printf ("# "); printf
#define dbgmsg printf ("# "); printf


/* Some cool error codes to use */
#define DEC_MEM_ERROR    0x0001
#define DEC_RES_INVALID  0x0002
#define DEC_SANITY_ERROR 0x0004
#define DEC_ARG_INVALID  0x1000
#define DEC_ARG_CORRUPT  0x2000


#define DEC_MASK 0xffff


int ShowError(int ec);


long FileLength(char* fname);

#endif
