#ifndef SIDIS_SF_SET_STF_HPP
#define SIDIS_SF_SET_STF_HPP

#include "sidis/structure_function.hpp"
#include "sidis/tmd.hpp"

namespace sidis {

namespace constant {
	enum class Hadron;
}

namespace sf {
namespace set {

/**
 * Use data files from STFLIB to calculate TMDs and FFs.
 */
class StfTmdSet final : public GaussianTmdSet {
private:
	struct Impl;
	Impl* _impl;

public:
	StfTmdSet();
	StfTmdSet(StfTmdSet const&) = delete;
	StfTmdSet(StfTmdSet&&) noexcept;
	StfTmdSet& operator=(StfTmdSet const&) = delete;
	StfTmdSet& operator=(StfTmdSet&& other) noexcept;
	virtual ~StfTmdSet();

	Real charge(unsigned fl) const;

	Real xf1(unsigned fl, Real x, Real Q_sq) const override;
	Real xf1Tperp(unsigned fl, Real x, Real Q_sq) const override;
	Real xh1(unsigned fl, Real x, Real Q_sq) const override;

	Real D1(part::Hadron h, unsigned fl, Real z, Real Q_sq) const override;
	Real H1perp(part::Hadron h, unsigned fl, Real z, Real Q_sq) const override;
};

}
}
}

#endif

