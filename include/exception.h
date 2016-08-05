#ifndef __EXCEPTION_H__
#define __EXCEPTION_H__

#include "config/config.h"

#include <string>
#include <exception>
#include <iostream>
#include <sstream>

class Exception : public std::exception {
protected:
    std::string msg;
public:
    // Constructor
    Exception(const std::string& _msg = "Unspecified error occurred") throw() : msg(_msg) { }

    virtual ~Exception() throw() { }

    // Gets message text
    const char * what() const throw () { return msg.c_str(); }
};

#define THROW throw(Exception)
#define THROW_EXCEPTION(what) \
{ \
    std::ostringstream os; \
    os << "EXCEPTION (" << __FILE__ << ":" << __LINE__ << ", " << __func__ << "): " << what << std::endl; \
    throw Exception(os.str()); \
}

#define unlikely(expr) __builtin_expect(!!(expr), 0)
#define likely(expr)   __buildin_expect(!!(expr), 1)

#ifndef NO_ASSERT
#define ASSERT(expression, what) \
{ \
    if (unlikely(!(expression))) \
    THROW_EXCEPTION(what); \
}
#else
#define ASSERT(expression,what)
#endif

#define ASSERT_SIZE(got, expected) ASSERT(got == expected, "Wrong size: " << got << " (expected " << expected << ")")
#define ASSERT_SIZE_LESS(got,expected) ASSERT(uint(got) < uint(expected), "Index is out of boundaries: " << got << " (max = " << expected << ")")


#endif // __EXCEPTION_H__
