#ifndef __LOGGER_H__
#define __LOGGER_H__

#include "config/config.h"

#include <vector>
#include <map>
#include <set>
#include <stack>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cstdarg>
#include <cstdio>

#ifndef NO_LOGGER

#include <log4cxx/basicconfigurator.h>
#include <log4cxx/logger.h>


namespace log4cxx {
    class SysLogger;

    /* ostringstream wrapper for SysLogger */
    class LoggerStream {
    private:
	SysLogger* logger;
	LevelPtr priority;
	std::ostringstream* buf;
    public:
	LoggerStream(const LoggerStream& ls);
	LoggerStream(const LevelPtr& p, SysLogger* l);
	~LoggerStream();

	/* To output class must support smth like ostringstream& operator<<(ostringstream&) */
	template<class T>
	LoggerStream& operator<<(const T& t) {
	    if (buf)
		(*buf) << t;
	    return (*this);
	}

	/* Manipulators std::hex, etc */
	template<typename X, typename Y>
	LoggerStream& operator<<(X& (*__pf)(Y&)) {
	    if (buf)
		(*buf) << __pf;
	    return (*this);
	}
    };


    class SysLogger : public LoggerPtr {
    private:
	void valog(const LevelPtr& p, const char *fmt, va_list vl);
    public:
	void logf(const LevelPtr& p, const char* fmt, ...);

	/* macro for Logger.info(char* fmt, ...), etc */
#define __LOG_FUNC(fn_name, fn_prio) \
	void fn_name(char* fmt, ...) \
	{ \
	    va_list vl; va_start(vl, fmt); \
	    valog(Level::get##fn_prio(), fmt, vl); \
	    va_end(vl); \
	}
	__LOG_FUNC(fatal, Fatal);
	__LOG_FUNC(error, Error);
	__LOG_FUNC(warn,  Warn);
	__LOG_FUNC(info,  Info);
	__LOG_FUNC(debug, Debug);
#undef __LOG_FUNC

	SysLogger(const LoggerPtr& b) : LoggerPtr(b) {
	}

	LevelPtr getLogLevel() const {
	    return (*this)->getLevel();
	}

	void setLogLevel(const LevelPtr& level) {
	    (*this)->setLevel(level);
	}

	void addAppender(AppenderPtr newAppender) {
	    (*this)->addAppender(newAppender);
	}

	void removeAllAppenders() {
	    (*this)->removeAllAppenders();
	}

	LoggerStream operator<<(const LevelPtr& p) {
	    return LoggerStream(p, this);
	}
    };


    inline void SysLogger::logf(const LevelPtr& p, const char* fmt, ...) {
	va_list vl; va_start(vl, fmt);
	valog(p, fmt, vl);
	va_end(vl);
    }

    /// First it dumps fmt+vl to some string. But it is only for standart objects, all implemented classes should suppor <<
    inline void SysLogger::valog(const LevelPtr& p_level, const char* fmt, va_list vl) {
	int size = 1024;
	char* buf = new char[size];

	while (vsnprintf(buf,size,fmt,vl) >= size-16) {
	    delete buf;
	    size += 1024;
	    buf = new char[size];
	}

	(*this)->log(p_level, std::string(buf));
	delete buf;
    }

    inline LoggerStream::LoggerStream(const LoggerStream& ls) {
	logger = ls.logger;
	priority = ls.priority;
	if ((*logger)->isEnabledFor(priority))
	    buf = new std::ostringstream;
	else buf = NULL;
    }

    inline LoggerStream::LoggerStream(const LevelPtr& p, SysLogger* l) : logger(l), priority(p), buf(0) {
	if ((*logger)->isEnabledFor(priority))
	    buf = new std::ostringstream;
	else buf = NULL;
    }

    inline LoggerStream::~LoggerStream() {
	if (buf) {
	    (*logger)->log(priority, (const std::string)(buf->str().c_str()));
	    delete buf;
	    buf = NULL;
	}
    }
} /* namespace log4cxx */

/* Some useful macro for the case we have one instance of logger for the module */
#define DEFINE_LOGGER(l)	static log4cxx::SysLogger logger = log4cxx::Logger::getLogger(l)
#define _LOG_DEBUG(lvl,v)	logger << log4cxx::Level::getDebug() << __func__ << " " << v
#define _LOG_INFO(lvl,v)	logger << log4cxx::Level::getInfo () << __func__ << " " << v
#define _LOG_WARN(lvl,v)	logger << log4cxx::Level::getWarn () << __func__ << " " << v
#define _LOG_ERROR(lvl,v)	logger << log4cxx::Level::getError() << __func__ << " " << v
#define _LOG_FATAL(lvl,v)	logger << log4cxx::Level::getFatal() << __func__ << " " << v

#define _LOG_DEBUG_P(fmt,...)	logger.logf(log4cxx::Level::getDebug(), fmt, __VA_ARGS__)

#else // #ifndef NO_LOGGER

#define DEFINE_LOGGER(l)	static int ____logger____
#define _LOG(lvl,v)		std::cout << #lvl " : " << __func__ << " : " << v << std::endl
#define _LOG_DEBUG(lvl,v)	std::cout << "DEBUG : " << __func__ << " " << v
#define _LOG_INFO(lvl,v)	std::cout << "INFO : " << __func__ << " " << v
#define _LOG_WARN(lvl,v)	std::cout << "WARN : " << __func__ << " " << v
#define _LOG_ERROR(lvl,v)	std::cout << "ERROR : " << __func__ << " " << v
#define _LOG_FATAL(lvl,v)	std::cout << "FATAL : " << __func__ << " " << v

#define _LOG_DEBUG_P(fmt,...)	printf("DEBUG : " fmt, __VA_ARGS__)

#endif // #ifndef NO_LOGGER

#define LOG_DEBUG(v)		_LOG_DEBUG(DEBUG,v)
#define LOG_INFO(v)		_LOG_INFO(INFO,v)
#define LOG_WARN(v)		_LOG_WARN(WARN,v)
#define LOG_ERROR(v)		_LOG_ERROR(ERROR,v)
#define LOG_FATAL(v)		_LOG_FATAL(FATAL,v)

#define LOG_DEBUG_P(fmt,...)	_LOG_DEBUG_P(fmt, __VA_ARGS__)

#define LOG_VAR(v) LOG_DEBUG(#v " = " << std::scientific << (v))

#ifndef NO_LOGGER
#define LLL_INFO(v)  LOG_INFO(v),  std::cout << "INFO : " << v << std::endl
#define LLL_DEBUG(v) LOG_DEBUG(v), std::cout << "DEBUG : " << v << std::endl
#define LLL_WARN(v)  LOG_WARN(v),  std::cout << "WARN : " << v << std::endl
#else
#define LLL_INFO(v)  LOG_INFO(v)
#define LLL_DEBUG(v) LOG_DEBUG(v)
#define LLL_WARN(v)  LOG_WARN(v)
#endif

/* Logging of some common classes */
template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v) {
    os << " size = " << v.size() << std::endl;;
    for (typename std::vector<T>::const_iterator it = v.begin(); it != v.end(); it++)
	os << " " << it - v.begin() << ": " << *it << std::endl;
    return os;
}

template<typename T1,typename T2>
std::ostream& operator<<(std::ostream& os, const std::pair<T1,T2>& p) {
    return os << "(" << p.first << "," << p.second << ")";
}

template<typename T1, typename T2>
std::ostream& operator<<(std::ostream& os, const std::map<T1,T2>& p) {
    for (typename std::map<T1,T2>::const_iterator it = p.begin(); it != p.end(); it++)
	os << *it << std::endl;
    return os;
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const std::set<T>& s) {
    for (typename std::set<T>::const_iterator it = s.begin(); it != s.end(); it++)
	os << " " << *it;
    return os;
}

template<typename T>
std::ostream& operator<<(std::ostream& os, std::stack<T> s) {
    std::vector<T> v;
    while (!s.empty()) {
	v.push_back(s.top());
	s.pop();
    }
    os << " size = " << v.size() << std::endl;;
    for (int i = v.size()-1; i >= 0; i--)
	os << " " << v.size()-1-i << ": " << v[i] << std::endl;
    return os;
}

#endif // __LOGGER_H__
