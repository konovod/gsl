/* unistd.h - minimal compatibility shim for the MSVC build.
 *
 * Some test programs include <unistd.h> unconditionally for access().
 * MSVC provides it through <io.h>.  This header is only added to the
 * include path on Windows.
 */

#ifndef GSL_MSVC_UNISTD_H
#define GSL_MSVC_UNISTD_H

#include <io.h>

#endif /* GSL_MSVC_UNISTD_H */
