#ifndef __EXCEPTION_H__
#define __EXCEPTION_H__

#include "config/config.h"

#include "include/tools.h"

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

#ifndef NO_ASSERT
#define ASSERT(expression, what) \
{ \
    if (!(expression)) \
	THROW_EXCEPTION(what); \
}
#else
#define ASSERT(expression,what)
#endif


#endif // __EXCEPTION_H__
