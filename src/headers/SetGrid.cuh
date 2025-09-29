#ifndef _SETGRIDCUH
#define _SETGRIDCUH

const uint32_t NX = 128;
const uint32_t NY = 128;
const uint32_t NZ = 50000;
#ifdef DISPERSION
const uint32_t NT = 128;
#else
const uint32_t NT = 1;
#endif
const uint32_t BLKT = 16;
const uint32_t BLKX = 16;
const uint32_t BLKY = 16;
const uint32_t SIZE = NX*NY*NT;


#endif // _SETGRIDCUH
