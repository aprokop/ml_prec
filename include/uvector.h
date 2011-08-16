#ifndef __UVECTOR_H__
#define __UVECTOR_H__

#include <cassert>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <vector>

#include "include/exception.h"

/* uvector works only for scalar types (int, unsigned, double) */

#if 0
#define uvector std::vector

#else
template<typename T>
class uvector {
private:
    T*	start;
    T*	finish;
    T*	end_of_storage;

    T* allocate(size_t n) {
	return (n != 0 ? (T*)malloc(n*sizeof(T)) : 0);
    }
    T* reallocate(T* p, size_t n) {
	return (p != 0 || n != 0) ? (T*)realloc(p, n*sizeof(T)) : 0;
    }
    void deallocate(T* p) {
	if (p)
	    free(p);
	p = 0;
    }

    typedef uvector<T>					vector_type;

public:
    typedef T						value_type;
    typedef T*						pointer;
    typedef const T*					const_pointer;
    typedef T&						reference;
    typedef const T&					const_reference;

    typedef __gnu_cxx::__normal_iterator<pointer, vector_type> iterator;
    typedef __gnu_cxx::__normal_iterator<const_pointer, vector_type>
	    const_iterator;
    typedef std::reverse_iterator<const_iterator>	const_reverse_iterator;
    typedef std::reverse_iterator<iterator>		reverse_iterator;

    typedef size_t					size_type;
    typedef ptrdiff_t					difference_type;

public:
    uvector(size_type n = 0) {
	start = allocate(n);
	finish = start + n;
	end_of_storage = start + n;
    }
    uvector(size_type n, const value_type& val) {
	start = allocate(n);
	finish = start + n;
	end_of_storage = start + n;
	for (size_type i = 0; i < n; i++)
	    *(start + i) = val;
    }

    uvector(const uvector& v) {
	size_type n = v.size();
	start = allocate(n);
	finish = start + n;
	end_of_storage = start + n;
	memcpy(start, v.start, n*sizeof(value_type));
    }
    ~uvector() {
	deallocate(start);
    }

    /* Iterators */
    iterator			begin()			{ return iterator(start); }
    const_iterator		begin() const		{ return const_iterator(start); }
    iterator			end()			{ return iterator(finish); }
    const_iterator		end() const		{ return const_iterator(finish); }
    reverse_iterator		rbegin()		{ return reverse_iterator(end()); }
    const_reverse_iterator	rbegin() const		{ return const_reverse_iterator(end()); }
    reverse_iterator		rend()			{ return reverse_iterator(begin()); }
    const_reverse_iterator	rend() const		{ return const_reverse_iterator(begin()); }

    size_type size() const {
	return finish - start;
    }

    /*
     * This function does not resize to increasing size (only if start = 0)
     * When I tried doing that it seems that compiler stops inlining it
     * and the call becomes significantly slower
     */
    void resize(size_type new_size) {
	if (start + new_size <= end_of_storage) {
	    finish = start + new_size;
	} else {
#if 0
	    /* Useful for checking if we resize arrays from
	     * already nonzero somewhere */
	    assert(!start);
#endif
	    start = reallocate(start, new_size);
	    finish = start + new_size;
	    end_of_storage = finish;
	}
    }
    void resize(size_type new_size, value_type val) {
	unsigned old_size = size();
	resize(new_size);
	for (unsigned i = old_size; i < new_size; i++)
	    *(start + i) = val;
    }

    bool empty() const {
	return start == finish;
    }

    void reserve(size_type new_size) {
	if (start + new_size > end_of_storage) {
	    size_type old_size = size();

	    start = reallocate(start, new_size);
	    finish = start + old_size;
	    end_of_storage = start + new_size;
	}
    }

    const uvector& operator=(const uvector& x) {
	deallocate(start);

	size_type n = x.size();
	start = allocate(n);
	finish = start + n;
	end_of_storage = start + n;

	if (n)
	    memcpy(start, x.start, n*sizeof(T));
	return (*this);
    }

    reference operator[](size_type i) {
	ASSERT(i < size(), "Index is out of bundaries: i = " << i << ", n = " << size());
	return *(start + i);
    }

    const_reference operator[](size_type i) const {
	ASSERT(i < size(), "Index is out of bundaries: i = " << i << ", n = " << size());
	return *(start + i);
    }

    reference back() {
	return *(end()-1);
    }
    const_reference back() const {
	return *(end()-1);
    }

    pointer data() {
	return pointer(start);
    }
    const_pointer data() const {
	return const_pointer(start);
    }

    void push_back(value_type x) {
	if (finish != end_of_storage) {
	    *(finish++) = x;
	} else {
	    const size_type old_size = size();
	    const size_type len = old_size != 0 ? 2 * old_size : 1;
	    start = reallocate(start, len);
	    finish = start + old_size;
	    end_of_storage = start + len;
	    *(finish++) = x;
	}
    }
    void swap(uvector& x) {
	std::swap(start, x.start);
	std::swap(finish, x.finish);
	std::swap(end_of_storage, x.end_of_storage);
    }

    void clear() {
	finish = start;
    }
};
#endif

template<typename T1, typename T2>
bool operator==(const uvector<T1>& x1, const uvector<T2>& x2) {
    if (x1.size() != x2.size())
	return false;
    uint n = x1.size();
    for (uint i = 0; i < n; i++)
	if (x1[i] != x2[i])
	    return false;
    return true;
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const uvector<T>& v) {
    os << " size = " << v.size() << std::endl;
    for (typename uvector<T>::const_iterator it = v.begin(); it != v.end(); it++)
	os << " " << it - v.begin() << ": " << *it << std::endl;
    return os;
}

#endif // #define __UVECTOR_H__
