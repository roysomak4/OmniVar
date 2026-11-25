// Linked into every executable of a OMNIVAR_SANITIZE build (CMakeLists.txt).
//
// With libc++ 19 as packaged by Ubuntu/Debian, std::runtime_error's message string is allocated
// by libc++ (operator new) and freed by the system libc++abi (free). AddressSanitizer reports
// that as alloc-dealloc-mismatch for every caught exception, including the RecoverableError the
// caller throws for recoverable errors (reproduced with a four-line program that throws and catches
// a std::runtime_error). The check is switched off; every other ASan check stays on.

extern "C" const char* __asan_default_options();
extern "C" const char* __asan_default_options() { return "alloc_dealloc_mismatch=0"; }
