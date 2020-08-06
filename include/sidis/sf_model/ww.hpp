#ifndef __SIDIS_SF_MODEL_WW_HPP__
#define __SIDIS_SF_MODEL_WW_HPP__

#include <stdexcept>

#include "sidis/numeric.hpp"
#include "sidis/structure_function.hpp"

namespace sidis {
namespace sf {
namespace model {

/**
 * Wandzura-Wilczek approximation to the structure functions [2].
 */
class WW final {
	struct Impl;
	Impl* _impl;

public:
	struct EnvironmentInitException : public std::runtime_error {
		EnvironmentInitException();
	};
	struct LinkInitException : public std::runtime_error {
		int code;
		LinkInitException(int code);
	};
	struct TimeoutInitException : public std::runtime_error {
		int code;
		TimeoutInitException(int code);
	};
	struct LibraryLoadException : public std::runtime_error {
		char const* filepath;
		LibraryLoadException(char const* filepath);
	};
	struct SendPacketException : public std::runtime_error {
		int code;
		SendPacketException(int code);
	};
	struct ReceivePacketException : public std::runtime_error {
		int code;
		ReceivePacketException(int code);
	};
	struct UnexpectedPacketException : public std::runtime_error {
		int code;
		UnexpectedPacketException(int code);
	};
	struct UnexpectedPacketContentsException : public std::runtime_error {
		int code;
		UnexpectedPacketContentsException(int code);
	};
	struct NextPacketFailureException : public std::runtime_error {
		int code;
		NextPacketFailureException(int code);
	};

	WW(unsigned timeout_packet = 1000, unsigned timeout_init = 60000);
	~WW();

	WW(WW const& other) = delete;
	WW(WW&& other) noexcept;
	WW& operator=(WW const& other) = delete;
	WW& operator=(WW&& other) noexcept;

	SfUU sf_uu(Real x, Real z, Real Q_sq, Real ph_t) const;
	SfUL sf_ul(Real x, Real z, Real Q_sq, Real ph_t) const;
	SfUT sf_ut(Real x, Real z, Real Q_sq, Real ph_t) const;
	SfLU sf_lu(Real x, Real z, Real Q_sq, Real ph_t) const;
	SfLL sf_ll(Real x, Real z, Real Q_sq, Real ph_t) const;
	SfLT sf_lt(Real x, Real z, Real Q_sq, Real ph_t) const;

	SfXU sf_xu(Real x, Real z, Real Q_sq, Real ph_t) const {
		return { sf_uu(x, z, Q_sq, ph_t), sf_lu(x, z, Q_sq, ph_t) };
	}
	SfXL sf_xl(Real x, Real z, Real Q_sq, Real ph_t) const {
		return { sf_ul(x, z, Q_sq, ph_t), sf_ll(x, z, Q_sq, ph_t) };
	}
	SfXT sf_xt(Real x, Real z, Real Q_sq, Real ph_t) const {
		return { sf_ut(x, z, Q_sq, ph_t), sf_lt(x, z, Q_sq, ph_t) };
	}
	SfUX sf_ux(Real x, Real z, Real Q_sq, Real ph_t) const {
		return {
			sf_uu(x, z, Q_sq, ph_t),
			sf_ul(x, z, Q_sq, ph_t),
			sf_ut(x, z, Q_sq, ph_t),
		};
	}
	SfLX sf_lx(Real x, Real z, Real Q_sq, Real ph_t) const {
		return {
			sf_lu(x, z, Q_sq, ph_t),
			sf_ll(x, z, Q_sq, ph_t),
			sf_lt(x, z, Q_sq, ph_t),
		};
	}

	Sf sf(Real x, Real z, Real Q_sq, Real ph_t) const {
		return { sf_ux(x, z, Q_sq, ph_t), sf_lx(x, z, Q_sq, ph_t) };
	}
};

}
}
}

#endif

