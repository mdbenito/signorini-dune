/* config.h.  Generated from config.h.in by configure.  */
/* config.h.in.  Generated from configure.ac by autoheader.  */

/* Alberta version found by configure, either 0x200 for 2.0 or 0x300 for 3.0
   */
#define DUNE_ALBERTA_VERSION 0x300
#define ENABLE_ALBERTA 1
#define ALBERTA_DIM 2

/* If this is set, the member 'size' of FieldVector is a method rather than an
   enum */
#define DUNE_COMMON_FIELDVECTOR_SIZE_IS_METHOD 1

/* Define to the version of dune-common */
#define DUNE_COMMON_VERSION "2.2.0"

/* Define to the major version of dune-common */
#define DUNE_COMMON_VERSION_MAJOR 2

/* Define to the minor version of dune-common */
#define DUNE_COMMON_VERSION_MINOR 2

/* Define to the revision of dune-common */
#define DUNE_COMMON_VERSION_REVISION 0

/* Define to the version of dune-geometry */
#define DUNE_GEOMETRY_VERSION "2.2.0"

/* Define to the major version of dune-geometry */
#define DUNE_GEOMETRY_VERSION_MAJOR 2

/* Define to the minor version of dune-geometry */
#define DUNE_GEOMETRY_VERSION_MINOR 2

/* Define to the revision of dune-geometry */
#define DUNE_GEOMETRY_VERSION_REVISION 0

/* If this is set, public access to the implementation of facades like Entity,
   Geometry, etc. is granted. */
/* #undef DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS */

/* Define to the version of dune-grid */
#define DUNE_GRID_VERSION "2.2.0"

/* Define to the major version of dune-grid */
#define DUNE_GRID_VERSION_MAJOR 2

/* Define to the minor version of dune-grid */
#define DUNE_GRID_VERSION_MINOR 2

/* Define to the revision of dune-grid */
#define DUNE_GRID_VERSION_REVISION 0

/* Define to the version of dune-istl */
#define DUNE_ISTL_VERSION "2.2.0"

/* Define to the major version of dune-istl */
#define DUNE_ISTL_VERSION_MAJOR 2

/* Define to the minor version of dune-istl */
#define DUNE_ISTL_VERSION_MINOR 2

/* Define to the revision of dune-istl */
#define DUNE_ISTL_VERSION_REVISION 0

/* Define to the version of dune-localfunctions */
#define DUNE_LOCALFUNCTIONS_VERSION "2.2.0"

/* Define to the major version of dune-localfunctions */
#define DUNE_LOCALFUNCTIONS_VERSION_MAJOR 2

/* Define to the minor version of dune-localfunctions */
#define DUNE_LOCALFUNCTIONS_VERSION_MINOR 2

/* Define to the revision of dune-localfunctions */
#define DUNE_LOCALFUNCTIONS_VERSION_REVISION 0

/* Standard debug streams with a level below will collapse to doing nothing */
#define DUNE_MINIMAL_DEBUG_LEVEL 4

/* Define to dummy `main' function (if any) required to link to the Fortran
   libraries. */
/* #undef F77_DUMMY_MAIN */

/* Define to dummy `main' function (if any) required to link to the Fortran
   libraries. */
/* #undef FC_DUMMY_MAIN */

/* Define if F77 and FC dummy `main' functions are identical. */
/* #undef FC_DUMMY_MAIN_EQ_F77 */

/* Define to a macro mangling the given C identifier (in lower and upper
   case), which must not contain underscores, for linking with Fortran. */
#define FC_FUNC(name,NAME) name ## _

/* As FC_FUNC, but for C identifiers containing underscores. */
#define FC_FUNC_(name,NAME) name ## _

/* Define to the version of Signorini */
#define GRIDHOWTO_VERSION "0.1"

/* Define to the major version of Signorini */
#define GRIDHOWTO_VERSION_MAJOR 0

/* Define to the minor version of Signorini */
#define GRIDHOWTO_VERSION_MINOR 1

/* Define to the revision of Signorini */
#define GRIDHOWTO_VERSION_REVISION 0

/* does the compiler support __attribute__((deprecated))? */
#define HAS_ATTRIBUTE_DEPRECATED 1

/* does the compiler support __attribute__((deprecated("message"))? */
#define HAS_ATTRIBUTE_DEPRECATED_MSG 1

/* does the compiler support __attribute__((unused))? */
#define HAS_ATTRIBUTE_UNUSED 1

/* This is only true if alberta-library was found by configure _and_ if the
   application uses the ALBERTA_CPPFLAGS */
  //#define HAVE_ALBERTA ENABLE_ALBERTA

/* Was AlgLib for DUNE found and ALGLIB_CPPFLAGS used? */
/* #undef HAVE_ALGLIB */

/* This is only true if alugrid-library was found by configure _and_ if the
   application uses the ALUGRID_CPPFLAGS */
#define HAVE_ALUGRID ENABLE_ALUGRID

/* Define to 1 if you have the <alugrid_parallel.h> header file. */
#define HAVE_ALUGRID_PARALLEL_H 1

/* Define to 1 if you have the <alugrid_serial.h> header file. */
#define HAVE_ALUGRID_SERIAL_H 1

/* Define to 1 if amiramesh-library is found */
/* #undef HAVE_AMIRAMESH */

/* Use the Apple OpenGL framework. */
/* #undef HAVE_APPLE_OPENGL_FRAMEWORK */

/* Define to 1 if the <array> C++0x is available and support array::fill */
/* #undef HAVE_ARRAY */

/* Define if you have a BLAS library. */
#define HAVE_BLAS 1

/* Define to ENABLE_BOOST if the Boost library is available */
/* #undef HAVE_BOOST */

/* define if the Boost::Fusion headers are available */
/* #undef HAVE_BOOST_FUSION */

/* Define to 1 if you have <boost/make_shared.hpp>. */
/* #undef HAVE_BOOST_MAKE_SHARED_HPP */

/* Define to 1 if you have the <boost/shared_ptr.hpp> header file. */
/* #undef HAVE_BOOST_SHARED_PTR_HPP */

/* Define to 1 if you have the <dlfcn.h> header file. */
#define HAVE_DLFCN_H 1

/* Define to 1 if dune-common was found */
#define HAVE_DUNE_COMMON 1

/* Define to 1 if dune-geometry was found */
#define HAVE_DUNE_GEOMETRY 1

/* Define to 1 if dune-grid was found */
#define HAVE_DUNE_GRID 1

/* Define to 1 if dune-istl was found */
#define HAVE_DUNE_ISTL 1

/* Define to 1 if dune-localfunctions was found */
#define HAVE_DUNE_LOCALFUNCTIONS 1

/* Was GMP found and GMP_CPPFLAGS used? */
/* #undef HAVE_GMP */

/* This is only true if grape-library was found by configure _and_ if the
   application uses the GRAPE_CPPFLAGS */
/* #undef HAVE_GRAPE */

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define if you have LAPACK library. */
#define HAVE_LAPACK 1

/* Define to 1 if you have the `gmp' library (-L$with_gmp/lib -lgmp). */
/* #undef HAVE_LIBGMP */

/* Define to 1 if you have the `m' library (-lm). */
#define HAVE_LIBM 1

/* Define to 1 if SHARED_PTR_NAMESPACE::make_shared is usable. */
/* #undef HAVE_MAKE_SHARED */

/* Define to 1 if your system has a GNU libc compatible `malloc' function, and
   to 0 otherwise. */
#define HAVE_MALLOC 1

/* Define to 1 if you have the <malloc.h> header file. */
/* #undef HAVE_MALLOC_H */

/* Define to 1 if you have the <memory> header file. */
#define HAVE_MEMORY 1

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* Define to 1 if `expansions' is a member of `mem_usage_t'. */
/* #undef HAVE_MEM_USAGE_T_EXPANSIONS */

/* Define if you have METIS library */
/* #undef HAVE_METIS */

/* Define if you have the MPI library. This is only true if MPI was found by
   configure _and_ if the application uses the DUNEMPICPPFLAGS (or the
   deprecated MPI_CPPFLAGS) */
#define HAVE_MPI ENABLE_MPI

/* Define to 1 if nullptr is supported */
/* #undef HAVE_NULLPTR */

/* Define if you have a PARDISO library. */
/* #undef HAVE_PARDISO */

/* Define if you have the Parmetis library. This is only true if MPI was found
   by configure _and_ if the application uses the PARMETIS_CPPFLAGS */
/* #undef HAVE_PARMETIS */

/* Define to 1 if psurface-library is found */
/* #undef HAVE_PSURFACE */

/* Define if you have POSIX threads libraries and header files. */
#define HAVE_PTHREAD 1

/* Define to 1 if you have the <rpc/rpc.h> header file. */
#define HAVE_RPC_RPC_H 1

/* Define to 1 if you have the `sqrt' function. */
#define HAVE_SQRT 1

/* Define to 1 if static_assert is supported */
/* #undef HAVE_STATIC_ASSERT */

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the `strchr' function. */
#define HAVE_STRCHR 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to ENABLE_SUPERLU if SUPERLU is found */
/* #undef HAVE_SUPERLU */

/* Define to 1 if SUPERLU_DIST is found */
/* #undef HAVE_SUPERLU_DIST */

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <tr1/memory> header file. */
  //#define HAVE_TR1_MEMORY 1

/* Define to 1 if you have the <tr1/tuple> header file. */
  //#define HAVE_TR1_TUPLE 1

/* Define to 1 if you have the <tr1/type_traits> header file. */
  //#define HAVE_TR1_TYPE_TRAITS 1

/* Define to 1 if you have the <tuple> header file. */
/* #undef HAVE_TUPLE */

/* Define to 1 if you have the <type_traits> header file. */
/* #undef HAVE_TYPE_TRAITS */

/* This is only true if UG was found by configure _and_ if the application
   uses the UG_CPPFLAGS */
/* #undef HAVE_UG */

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Define to 1 if you have the <windows.h> header file. */
/* #undef HAVE_WINDOWS_H */

/* Define to the sub-directory in which libtool stores uninstalled libraries.
   */
#define LT_OBJDIR ".libs/"

/* Define to 1 MPI supports MPI-2 */
#define MPI_2 1

/* Name of package */
#define PACKAGE "Signorini"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "nonick@8027.org"

/* Define to the full name of this package. */
#define PACKAGE_NAME "Signorini"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "Signorini 0.1"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "Signorini"

/* Define to the home page for this package. */
#define PACKAGE_URL ""

/* Define to the version of this package. */
#define PACKAGE_VERSION "0.1"

/* The namespace prefix of the psurface library */
/* #undef PSURFACE_NAMESPACE */

/* Define to necessary symbol if this constant uses a non-standard name on
   your system. */
/* #undef PTHREAD_CREATE_JOINABLE */

/* The header in which SHARED_PTR can be found */
#define SHARED_PTR_HEADER <tr1/memory>

/* The namespace in which SHARED_PTR can be found */
#define SHARED_PTR_NAMESPACE std::tr1

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* define to 1 if there SLU_DOUBLE imported by header slu_ddefs.h from SuperLU
   */
/* #undef SUPERLU_MIN_VERSION_4_3 */

/* define to 1 if there is a header slu_ddefs.h in SuperLU */
/* #undef SUPERLU_POST_2005_VERSION */

/* Define to 1 if your <sys/time.h> declares `struct tm'. */
/* #undef TM_IN_SYS_TIME */

/* use UG LGM domain */
/* #undef UG_LGMDOMAIN */

/* Version number of package */
#define VERSION "0.1"

/* Define to 1 if the X Window System is missing or not being used. */
/* #undef X_DISPLAY_MISSING */

/* Define to empty if `const' does not conform to ANSI C. */
/* #undef const */

/* Define to `__inline__' or `__inline' if that's what the C compiler
   calls it, or to nothing if 'inline' is not supported under any name.  */
#ifndef __cplusplus
/* #undef inline */
#endif

/* Define to rpl_malloc if the replacement function should be used. */
/* #undef malloc */

/* Define to `unsigned int' if <sys/types.h> does not define. */
/* #undef size_t */

#include <dune/common/deprecated.hh>

/* add GRIDTYPE typedef for grid implementation Dune::YaspGrid< dimgrid >:
    defining YASPGRID during compilation typedefs this grid implementation as GridType
    in namespace Dune::GridSelector;
    also integer constants dimgrid and dimworld are set in this namespace.
    The required headers for this grid implementation are also included.
  */
 #if defined YASPGRID && ! defined USED_YASPGRID_GRIDTYPE
  #if HAVE_GRIDTYPE
   #error "Ambiguous definition of GRIDTYPE."
  #endif 

  #ifndef WORLDDIM
    #define WORLDDIM GRIDDIM
  #endif
  #if not (WORLDDIM >= GRIDDIM)
    #error "WORLDDIM < GRIDDIM does not make sense."
  #endif

  #if ! (GRIDDIM == WORLDDIM)
    #error "Preprocessor assertion GRIDDIM == WORLDDIM failed."
  #endif

  #include <dune/grid/yaspgrid.hh>
  #include <dune/grid/io/file/dgfparser/dgfyasp.hh>

  namespace Dune
  {
    namespace GridSelector
    {
      const int dimgrid = GRIDDIM;
      const int dimworld = WORLDDIM;
      typedef Dune::YaspGrid< dimgrid > GridType;
    }
  }
  #define HAVE_GRIDTYPE 1
  #define USED_YASPGRID_GRIDTYPE 1
#endif // #if defined YASPGRID

#include <dune/common/unused.hh>

/* add GRIDTYPE typedef for grid implementation Dune::AlbertaGrid< dimgrid >:
    defining ALBERTAGRID during compilation typedefs this grid implementation as GridType
    in namespace Dune::GridSelector;
    also integer constants dimgrid and dimworld are set in this namespace.
    The required headers for this grid implementation are also included.
  */
 #if defined ALBERTAGRID && ! defined USED_ALBERTAGRID_GRIDTYPE
  #if HAVE_GRIDTYPE
   #error "Ambiguous definition of GRIDTYPE."
  #endif 

  #ifndef WORLDDIM
    #define WORLDDIM GRIDDIM
  #endif
  #if not (WORLDDIM >= GRIDDIM)
    #error "WORLDDIM < GRIDDIM does not make sense."
  #endif

  #if ! (WORLDDIM == ALBERTA_DIM)
    #error "Preprocessor assertion WORLDDIM == ALBERTA_DIM failed."
  #endif

  #include <dune/grid/albertagrid.hh>
  #include <dune/grid/albertagrid/dgfparser.hh>

  namespace Dune
  {
    namespace GridSelector
    {
      const int dimgrid = GRIDDIM;
      const int dimworld = WORLDDIM;
      typedef Dune::AlbertaGrid< dimgrid > GridType;
    }
  }
  #define HAVE_GRIDTYPE 1
  #define USED_ALBERTAGRID_GRIDTYPE 1
#endif // #if defined ALBERTAGRID

/* add GRIDTYPE typedef for grid implementation Dune::UGGrid< dimgrid >:
    defining UGGRID during compilation typedefs this grid implementation as GridType
    in namespace Dune::GridSelector;
    also integer constants dimgrid and dimworld are set in this namespace.
    The required headers for this grid implementation are also included.
  */
 #if defined UGGRID && ! defined USED_UGGRID_GRIDTYPE
  #if HAVE_GRIDTYPE
   #error "Ambiguous definition of GRIDTYPE."
  #endif 

  #ifndef WORLDDIM
    #define WORLDDIM GRIDDIM
  #endif
  #if not (WORLDDIM >= GRIDDIM)
    #error "WORLDDIM < GRIDDIM does not make sense."
  #endif

  #if ! (GRIDDIM == WORLDDIM)
    #error "Preprocessor assertion GRIDDIM == WORLDDIM failed."
  #endif

  #include <dune/grid/uggrid.hh>
  #include <dune/grid/io/file/dgfparser/dgfug.hh>

  namespace Dune
  {
    namespace GridSelector
    {
      const int dimgrid = GRIDDIM;
      const int dimworld = WORLDDIM;
      typedef Dune::UGGrid< dimgrid > GridType;
    }
  }
  #define HAVE_GRIDTYPE 1
  #define USED_UGGRID_GRIDTYPE 1
#endif // #if defined UGGRID

/* add GRIDTYPE typedef for grid implementation Dune::ALUConformGrid< dimgrid, dimworld >:
    defining ALUGRID_CONFORM during compilation typedefs this grid implementation as GridType
    in namespace Dune::GridSelector;
    also integer constants dimgrid and dimworld are set in this namespace.
    The required headers for this grid implementation are also included.
  */
 #if defined ALUGRID_CONFORM && ! defined USED_ALUGRID_CONFORM_GRIDTYPE
  #if HAVE_GRIDTYPE
   #error "Ambiguous definition of GRIDTYPE."
  #endif 

  #ifndef WORLDDIM
    #define WORLDDIM GRIDDIM
  #endif
  #if not (WORLDDIM >= GRIDDIM)
    #error "WORLDDIM < GRIDDIM does not make sense."
  #endif

  #include <dune/grid/alugrid.hh>
  #include <dune/grid/io/file/dgfparser/dgfalu.hh>

  namespace Dune
  {
    namespace GridSelector
    {
      const int dimgrid = GRIDDIM;
      const int dimworld = WORLDDIM;
      typedef Dune::ALUConformGrid< dimgrid, dimworld > GridType;
    }
  }
  #define HAVE_GRIDTYPE 1
  #define USED_ALUGRID_CONFORM_GRIDTYPE 1
#endif // #if defined ALUGRID_CONFORM

/* add GRIDTYPE typedef for grid implementation Dune::ALUCubeGrid< dimgrid, dimworld >:
    defining ALUGRID_CUBE during compilation typedefs this grid implementation as GridType
    in namespace Dune::GridSelector;
    also integer constants dimgrid and dimworld are set in this namespace.
    The required headers for this grid implementation are also included.
  */
 #if defined ALUGRID_CUBE && ! defined USED_ALUGRID_CUBE_GRIDTYPE
  #if HAVE_GRIDTYPE
   #error "Ambiguous definition of GRIDTYPE."
  #endif 

  #ifndef WORLDDIM
    #define WORLDDIM GRIDDIM
  #endif
  #if not (WORLDDIM >= GRIDDIM)
    #error "WORLDDIM < GRIDDIM does not make sense."
  #endif

  #include <dune/grid/alugrid.hh>
  #include <dune/grid/io/file/dgfparser/dgfalu.hh>

  namespace Dune
  {
    namespace GridSelector
    {
      const int dimgrid = GRIDDIM;
      const int dimworld = WORLDDIM;
      typedef Dune::ALUCubeGrid< dimgrid, dimworld > GridType;
    }
  }
  #define HAVE_GRIDTYPE 1
  #define USED_ALUGRID_CUBE_GRIDTYPE 1
#endif // #if defined ALUGRID_CUBE

/* add GRIDTYPE typedef for grid implementation Dune::ALUSimplexGrid< dimgrid, dimworld >:
    defining ALUGRID_SIMPLEX during compilation typedefs this grid implementation as GridType
    in namespace Dune::GridSelector;
    also integer constants dimgrid and dimworld are set in this namespace.
    The required headers for this grid implementation are also included.
  */
 #if defined ALUGRID_SIMPLEX && ! defined USED_ALUGRID_SIMPLEX_GRIDTYPE
  #if HAVE_GRIDTYPE
   #error "Ambiguous definition of GRIDTYPE."
  #endif 

  #ifndef WORLDDIM
    #define WORLDDIM GRIDDIM
  #endif
  #if not (WORLDDIM >= GRIDDIM)
    #error "WORLDDIM < GRIDDIM does not make sense."
  #endif

  #include <dune/grid/alugrid.hh>
  #include <dune/grid/io/file/dgfparser/dgfalu.hh>

  namespace Dune
  {
    namespace GridSelector
    {
      const int dimgrid = GRIDDIM;
      const int dimworld = WORLDDIM;
      typedef Dune::ALUSimplexGrid< dimgrid, dimworld > GridType;
    }
  }
  #define HAVE_GRIDTYPE 1
  #define USED_ALUGRID_SIMPLEX_GRIDTYPE 1
#endif // #if defined ALUGRID_SIMPLEX

/* add GRIDTYPE typedef for grid implementation Dune::OneDGrid:
    defining ONEDGRID during compilation typedefs this grid implementation as GridType
    in namespace Dune::GridSelector;
    also integer constants dimgrid and dimworld are set in this namespace.
    The required headers for this grid implementation are also included.
  */
 #if defined ONEDGRID && ! defined USED_ONEDGRID_GRIDTYPE
  #if HAVE_GRIDTYPE
   #error "Ambiguous definition of GRIDTYPE."
  #endif 

  #ifndef WORLDDIM
    #define WORLDDIM GRIDDIM
  #endif
  #if not (WORLDDIM >= GRIDDIM)
    #error "WORLDDIM < GRIDDIM does not make sense."
  #endif

  #if ! ((GRIDDIM == 1) && (WORLDDIM == 1))
    #error "Preprocessor assertion (GRIDDIM == 1) && (WORLDDIM == 1) failed."
  #endif

  #include <dune/grid/onedgrid.hh>
  #include <dune/grid/io/file/dgfparser/dgfoned.hh>

  namespace Dune
  {
    namespace GridSelector
    {
      const int dimgrid = GRIDDIM;
      const int dimworld = WORLDDIM;
      typedef Dune::OneDGrid GridType;
    }
  }
  #define HAVE_GRIDTYPE 1
  #define USED_ONEDGRID_GRIDTYPE 1
#endif // #if defined ONEDGRID

/* add GRIDTYPE typedef for grid implementation Dune::SGrid< dimgrid, dimworld >:
    defining SGRID during compilation typedefs this grid implementation as GridType
    in namespace Dune::GridSelector;
    also integer constants dimgrid and dimworld are set in this namespace.
    The required headers for this grid implementation are also included.
  */
 #if defined SGRID && ! defined USED_SGRID_GRIDTYPE
  #if HAVE_GRIDTYPE
   #error "Ambiguous definition of GRIDTYPE."
  #endif 

  #ifndef WORLDDIM
    #define WORLDDIM GRIDDIM
  #endif
  #if not (WORLDDIM >= GRIDDIM)
    #error "WORLDDIM < GRIDDIM does not make sense."
  #endif

  #include <dune/grid/sgrid.hh>
  #include <dune/grid/io/file/dgfparser/dgfs.hh>

  namespace Dune
  {
    namespace GridSelector
    {
      const int dimgrid = GRIDDIM;
      const int dimworld = WORLDDIM;
      typedef Dune::SGrid< dimgrid, dimworld > GridType;
    }
  }
  #define HAVE_GRIDTYPE 1
  #define USED_SGRID_GRIDTYPE 1
#endif // #if defined SGRID
