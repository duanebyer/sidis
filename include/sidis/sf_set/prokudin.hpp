#ifndef SIDIS_PROKUDIN_HPP
#define SIDIS_PROKUDIN_HPP

#include "sidis/tmd.hpp"

namespace sidis {

namespace constant {
	enum class Hadron;
}

namespace sf {
namespace model {

// Use data files from [2] to calculate TMDs and FFs.
class ProkudinTmdSet final : public GaussianWwTmdSet {
private:
	struct Impl;
	Impl* _impl;

public:
	ProkudinTmdSet();
	ProkudinTmdSet(ProkudinTmdSet&&) noexcept;
	ProkudinTmdSet& operator=(ProkudinTmdSet&& other) noexcept;
	~ProkudinTmdSet();

	Real charge(unsigned fl) const;

	Real xf1(unsigned fl, Real x, Real Q_sq) const override;
	Real xf1Tperp(unsigned fl, Real x, Real Q_sq) const override;
	Real xg1(unsigned fl, Real x, Real Q_sq) const override;
	Real xg1Tperp(unsigned fl, Real x, Real Q_sq) const override;
	Real xh1(unsigned fl, Real x, Real Q_sq) const override;
	Real xh1perp(unsigned fl, Real x, Real Q_sq) const override;
	Real xh1Lperp(unsigned fl, Real x, Real Q_sq) const override;
	Real xh1Tperp(unsigned fl, Real x, Real Q_sq) const override;

	Real D1(constant::Hadron h, unsigned fl, Real z, Real Q_sq) const override;
	Real H1perp(constant::Hadron h, unsigned fl, Real z, Real Q_sq) const override;
};

}
}
}

#endif

