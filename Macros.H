#ifndef MACROS_H
#define MACROS_H

#define USE_OMP
#define REAL double
#define GRIDDIM 2
#define SPACEDIM 2

#if GRIDDIM == 1
    #define GRIDDIM_TERM(a, b, c) a
    #define GRIDDIM_DECL(a, b, c) a
#elif GRIDDIM == 2
    #define GRIDDIM_TERM(a, b, c) a b
    #define GRIDDIM_DECL(a, b, c) a, b
#else  // GRIDDIM == 3
    #define GRIDDIM_TERM(a, b, c) a b c
    #define GRIDDIM_DECL(a, b, c) a, b, c
#endif

#if SPACEDIM == 1
    #define SPACEDIM_TERM(a, b, c) a
    #define SPACEDIM_DECL(a, b, c) a
#elif SPACEDIM == 2
    #define SPACEDIM_TERM(a, b, c) a b
    #define SPACEDIM_DECL(a, b, c) a, b
#else  // SPACEDIM == 3
    #define SPACEDIM_TERM(a, b, c) a b c
    #define SPACEDIM_DECL(a, b, c) a, b, c
#endif

#endif
