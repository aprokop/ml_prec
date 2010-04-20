#ifndef __SVECTOR_H__
#define __SVECTOR_H__

#include <cassert>
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <vector>

#include "include/exception.h"

/*
 * Replacement for std::set<T> with T being a scalar type (int, unsigned, double)
 * Typical use case: few insertions, many searches
 */

template<typename T>
class svector {
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

    typedef svector<T>					vector_type;

public:
    typedef T*						pointer;
private:
    typedef T&						reference;
    typedef __gnu_cxx::__normal_iterator<pointer, vector_type> iterator;
    typedef std::reverse_iterator<iterator>		reverse_iterator;

public:
    typedef T						value_type;
    typedef const T*					const_pointer;
    typedef const T&					const_reference;

    typedef __gnu_cxx::__normal_iterator<const_pointer, vector_type> const_iterator;
    typedef std::reverse_iterator<const_iterator>	const_reverse_iterator;

    typedef size_t					size_type;
    typedef ptrdiff_t					difference_type;

private:
    /* Iterators */
    // iterator		begin()	    { return iterator(start); }
    // iterator		end()	    { return iterator(finish); }
    // reverse_iterator	rbegin()    { return reverse_iterator(end()); }
    // reverse_iterator	rend()	    { return reverse_iterator(begin()); }

    /* Insert before the element at pos */
    void insert_(iterator& pos, value_type x) {
	unsigned new_size = size() + 1;

	if (start + new_size > end_of_storage) {
	    unsigned alloc_size = 2*new_size;
	    int shift = pos - begin();

	    start = reallocate(start, alloc_size);
	    end_of_storage = start + alloc_size;

	    pos = iterator(start) + shift;
	}
	finish = start + new_size;

	reverse_iterator rit = reverse_iterator(iterator(finish)), rlim = reverse_iterator(pos+1);
	for (; rit != rlim; rit++)
	    *rit = *(rit+1);
	*rit = x;
    }

public:
    svector() {
#if 1
	unsigned n = 5; /* we want to allocate a bit */
	start = allocate(n);
	finish = start;
	end_of_storage = start + n;
#else
	start = finish = end_of_storage = NULL;
#endif
    }

    svector(const svector& v) {
	size_type n = v.size();
	start = allocate(n);
	finish = start + n;
	end_of_storage = start + n;
	memcpy(start, v.start, n*sizeof(value_type));
    }
    ~svector() {
	deallocate(start);
    }

    /* Iterators */
    const_iterator		begin() const		{ return const_iterator(start); }
    const_iterator		end() const		{ return const_iterator(finish); }
    const_reverse_iterator	rbegin() const		{ return const_reverse_iterator(end()); }
    const_reverse_iterator	rend() const		{ return const_reverse_iterator(begin()); }

    size_type size() const {
	return finish - start;
    }

    /* This function does resizes only to <= size */
    // void resize(size_type new_size) {
	// unsigned old_size = size();
	// ASSERT(new_size <= old_size, "Attempt to resize to a bigger size");
	// finish = start + new_size;
    // }

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

    const svector& operator=(const svector& x) {
	deallocate(start);

	size_type n = x.size();
	start = allocate(n);
	finish = start + n;
	end_of_storage = start + n;

	if (n)
	    memcpy(start, x.start, n*sizeof(T));
	return (*this);
    }

    const_reference operator[](size_type i) const {
	ASSERT(i < size(), "Index is out of bundaries: i = " << i << ", n = " << size());
	return *(start + i);
    }

    const_reference back() const {
	return *(end()-1);
    }

    void swap(svector& x) {
	std::swap(start, x.start);
	std::swap(finish, x.finish);
	std::swap(end_of_storage, x.end_of_storage);
    }

    void insert(value_type x) {
	iterator it = std::lower_bound(iterator(start), iterator(finish), x);
	insert_(it, x);
    }

    const_iterator find(value_type x) const {
	if (!size())
	    return end();

	const_iterator it = lower_bound(x);
	return (*it == x) ? it : end();
    }

    const_iterator lower_bound(value_type x) const {
	return std::lower_bound(begin(), end(), x);
    }

    const_iterator upper_bound(value_type x) const {
	return std::upper_bound(begin(), end(), x);
    }

#if 0
    template<typename Iterator>
    void insert(const Iterator& begin, const Iterator& end) {
	THROW_EXCEPTION("Not implemented");
    }
#endif

    void clear() {
	finish = start;
    }
};

#if 0
template<typename T>
std::ostream& operator<<(std::ostream& os, const svector<T>& v) {
    os << " size = " << v.size() << std::endl;
    for (typename svector<T>::const_iterator it = v.begin(); it != v.end(); it++)
	os << " " << it - v.begin() << ": " << *it << std::endl;
    return os;
}
#else
template<typename T>
std::ostream& operator<<(std::ostream& os, const svector<T>& v) {
    for (typename svector<T>::const_iterator it = v.begin(); it != v.end(); it++)
	os << " " << *it;
    return os;
}
#endif

#endif // #define __UVECTOR_H__
