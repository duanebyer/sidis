#ifndef SIDIS_PIMPL_HPP
#define SIDIS_PIMPL_HPP

namespace sidis {
namespace pimpl {

/**
 * Convenience class for using the PIMPL idiom easily. This is basically a
 * unique_ptr, but with a few minor modifications for the specific use case.
 */
template<typename T>
class Pimpl final {
	T* _impl;
	// TODO: There's a warning emitted by clang because of the noexcept
	// specifier. The specifier could just be commented out, but we would like
	// to document that the deleter can't throw exceptions.
	void (*_deleter)(T*) noexcept;
	using Deleter = decltype(_deleter);
	Pimpl(T* impl, Deleter deleter) noexcept :
		_impl(impl),
		_deleter(deleter) { }

public:
	Pimpl(T const& other) = delete;
	Pimpl(T&& other) noexcept :
			_impl(other._impl),
			_deleter(other._deleter) {
		other._impl = nullptr;
	}
	Pimpl<T>& operator=(Pimpl<T> other) noexcept {
		swap(*this, other);
		return *this;
	}
	~Pimpl() noexcept {
		_deleter(_impl);
		_impl = nullptr;
	}

	T const& operator*() const { return *_impl; }
	T& operator*() { return *_impl; }
	T const* operator->() const { return _impl; }
	T* operator->() { return _impl; }

	friend void swap(Pimpl<T>& a, Pimpl<T>& b) noexcept {
		T* c_impl = a._impl;
		a._impl = b._impl;
		b._impl = c_impl;
		Deleter c_deleter = a._deleter;
		a._deleter = b._deleter;
		b._deleter = c_deleter;
	}

	template<typename T1, typename... Args>
	friend Pimpl<T1> make_pimpl(Args... args);
};

/**
 * Similar to Pimpl, but for types that can be copied.
 */
template<typename T>
class PimplCopy final {
	T* _impl;
	T* (*_copyer)(T*);
	void (*_deleter)(T*) noexcept;
	using Copyer = decltype(_copyer);
	using Deleter = decltype(_deleter);
	PimplCopy(T* impl, Copyer copyer, Deleter deleter) noexcept :
		_impl(impl),
		_copyer(copyer),
		_deleter(deleter) { }

public:
	PimplCopy(T const& other) :
		_impl(other._copyer(other._impl)),
		_copyer(other._copyer),
		_deleter(other._deleter) { }
	PimplCopy(T&& other) noexcept :
			_impl(other._impl),
			_copyer(other._copyer),
			_deleter(other._deleter) {
		other._impl = nullptr;
	}
	PimplCopy<T>& operator=(PimplCopy<T> other) noexcept {
		swap(*this, other);
		return *this;
	}
	~PimplCopy() noexcept {
		_deleter(_impl);
		_impl = nullptr;
	}

	T const& operator*() const { return *_impl; }
	T& operator*() { return *_impl; }
	T const* operator->() const { return _impl; }
	T* operator->() { return _impl; }

	friend void swap(PimplCopy<T>& a, PimplCopy<T>& b) noexcept {
		T* c_impl = a._impl;
		a._impl = b._impl;
		b._impl = c_impl;
		Copyer c_copyer = a._copyer;
		a._copyer = b._copyer;
		b._copyer = c_copyer;
		Deleter c_deleter = a._deleter;
		a._deleter = b._deleter;
		b._deleter = c_deleter;
	}

	template<typename T1, typename... Args>
	friend PimplCopy<T1> make_pimpl_copy(Args... args);
};

template<typename T>
T* default_copyer(T* val) {
	if (val == nullptr) {
		return nullptr;
	} else {
		return new T(*val);
	}
}

template<typename T>
void default_deleter(T* val) noexcept {
	if (val != nullptr) {
		delete val;
	}
}

template<typename T, typename... Args>
Pimpl<T> make_pimpl(Args... args) {
	return Pimpl<T>(new T(args...), default_deleter<T>);
}

template<typename T, typename... Args>
PimplCopy<T> make_pimpl_copy(Args... args) {
	return PimplCopy<T>(new T(args...), default_copyer<T>, default_deleter<T>);
}

}
}

#endif

