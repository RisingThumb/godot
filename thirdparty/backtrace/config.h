/* config.h.cmake */
/* ELF size: 32 or 64 */
/* #undef BACKTRACE_ELF_SIZE */
/* XCOFF size: 32 or 64 */
/* #undef BACKTRACE_XCOFF_SIZE */
/* Define to 1 if you have the __atomic functions */
/* #undef HAVE_ATOMIC_FUNCTIONS */
/* Define to 1 if you have the `clock_gettime' function. */
/* #undef HAVE_CLOCK_GETTIME */
/* Define to 1 if you have the declaration of `strnlen', and to 0 if you
   don't. */
/* #undef HAVE_DECL_STRNLEN */
/* Define to 1 if you have the <dlfcn.h> header file. */
/* #undef HAVE_DLFCN_H */
/* Define if dl_iterate_phdr is available. */
/* #undef HAVE_DL_ITERATE_PHDR */
/* Define to 1 if you have the fcntl function */
/* #undef HAVE_FCNTL */
/* Define if getexecname is available. */
/* #undef HAVE_GETEXECNAME */
/* Define if _Unwind_GetIPInfo is available. */
/* #undef HAVE_GETIPINFO */
/* Define to 1 if you have the <inttypes.h> header file. */
/* #undef HAVE_INTTYPES_H */
/* Define to 1 if you have the `z' library (-lz). */
/* #undef HAVE_LIBZ */
/* Define to 1 if you have the <link.h> header file. */
/* #undef HAVE_LINK_H */
/* Define if AIX loadquery is available. */
/* #undef HAVE_LOADQUERY */
/* Define to 1 if you have the `lstat' function. */
/* #undef HAVE_LSTAT */
/* Define to 1 if you have the <memory.h> header file. */
/* #undef HAVE_MEMORY_H */
/* Define to 1 if you have the `readlink' function. */
/* #undef HAVE_READLINK */
/* Define to 1 if you have the <stdint.h> header file. */
/* #undef HAVE_STDINT_H */
/* Define to 1 if you have the <stdlib.h> header file. */
/* #undef HAVE_STDLIB_H */
/* Define to 1 if you have the <strings.h> header file. */
/* #undef HAVE_STRINGS_H */
/* Define to 1 if you have the <string.h> header file. */
/* #undef HAVE_STRING_H */
/* Define to 1 if you have the __sync functions */
/* #undef HAVE_SYNC_FUNCTIONS */
/* Define to 1 if you have the <sys/ldr.h> header file. */
/* #undef HAVE_SYS_LDR_H */
/* Define to 1 if you have the <sys/mman.h> header file. */
/* #undef HAVE_SYS_MMAN_H */
/* Define to 1 if you have the <sys/stat.h> header file. */
/* #undef HAVE_SYS_STAT_H */
/* Define to 1 if you have the <sys/types.h> header file. */
/* #undef HAVE_SYS_TYPES_H */
/* Define to 1 if you have the <unistd.h> header file. */
/* #undef HAVE_UNISTD_H */
/* Define if -lz is available. */
/* #undef HAVE_ZLIB */

/* Enable extensions on AIX 3, Interix.  */
#ifndef _ALL_SOURCE
/* #undef _ALL_SOURCE */
#endif
/* Enable GNU extensions on systems that have them.  */
#ifndef _GNU_SOURCE
/* #undef _GNU_SOURCE */
#endif
/* Enable threading extensions on Solaris.  */
#ifndef _POSIX_PTHREAD_SEMANTICS
/* #undef _POSIX_PTHREAD_SEMANTICS */
#endif
/* Enable extensions on HP NonStop.  */
#ifndef _TANDEM_SOURCE
/* #undef _TANDEM_SOURCE */
#endif
/* Enable general extensions on Solaris.  */
#ifndef __EXTENSIONS__
/* #undef __EXTENSIONS__ */
#endif
