#ifndef SIDIS_SF_SET_TEST_HPP
#define SIDIS_SF_SET_TEST_HPP

#include "sidis/structure_function.hpp"

namespace sidis {
namespace sf {
namespace set {

/**
 * Simple (and inaccurate) parameterization of structure functions, useful for
 * testing purposes.
 *
 * The form of the parameterization is
 *
 * \f$ N Q^{-q} \left[ L(x;\alpha,\beta) + f L(x;\gamma,\delta) \right] L(z;\eta,\zeta) p_{ht}^{\rho} \exp{-\frac{p_{ht}^2}{2 \sigma^2}} \f$
 * \f$ L(x;\alpha,\beta) = \frac{\left(\alpha + \beta\right)^{\alpha + \beta}}{\alpha^\alpha \beta^\beta} x^\alpha (1 - x)^\beta \f$
 *
 * The parameterization is chosen because it is
 * * valid over arbitrary kinematics (no grids)
 * * fast to compute
 * * simple to debug
 * * easy to implement (for comparisons to other SIDIS software)
 *
 * This class is initialized by default to a simple fit, but the parameters can
 * be changed manually.
 *
 * \sa TestSfSet
 */
struct TestSfSetParams {
	struct Params {
		Real q;
		Real N;
		Real alpha;
		Real beta;
		Real f;
		Real gamma;
		Real delta;
		Real eta;
		Real zeta;
		Real rho;
		Real sigma;

		/// Evaluate the structure function described by the parameters.
		Real operator()(Real x, Real z, Real Q_sq, Real ph_t_sq) const;
	};

	//                                 q   |   N    | alph | beta |   f   | gamma | delt |  eta  | zeta | rho | sigma
	//                              -------|--------|------|------|-------|-------|------|-------|------|-----|-------
	Params F_UUL                 = { 1.23,   0.004,   1.1,   4.2,   0.,     0.,     0.,    0.94,   2.4,   0.,   0.33 }; // Completely made up.
	Params F_UUT                 = { 0.20,   1.200,   0.85,  3.0,   0.90,  -0.11,   9.5,  -0.80,   2.7,   0.,   0.33 };
	Params F_UU_cos_phih         = { 1.21,  -0.710,   1.0,   3.2,   0.95,  -0.15,   6.4,   0.14,   2.5,   1.,   0.33 };
	Params F_UU_cos_2phih        = { 0.26,  -0.064,   1.1,   5.3,   0.,     0.,     0.,    1.06,   2.1,   2.,   0.31 };

	Params F_UL_sin_phih         = { 1.11,   0.035,   2.3,   7.8,   0.,     0.,     0.,    0.15,   2.6,   1.,   0.31 };
	Params F_UL_sin_2phih        = { 0.21,  -0.013,   2.3,   7.8,   0.,     0.,     0.,    1.00,   2.9,   2.,   0.33 };

	Params F_UTL_sin_phih_m_phis = { 1.21,   0.005,   1.5,   5.3,   0.,     0.,     0.,    1.05,   2.1,   1.,   0.32 }; // Completely made up.
	Params F_UTT_sin_phih_m_phis = { 0.23,   0.110,   0.7,   4.5,   0.,     0.,     0.,    0.18,   2.0,   1.,   0.32 };
	Params F_UT_sin_2phih_m_phis = { 1.30,  -0.022,   1.2,   3.1,   0.,     0.,     0.,    0.85,   1.7,   2.,   0.32 };
	Params F_UT_sin_3phih_m_phis = { 0.53,   0.008,   2.9,   4.9,   0.,     0.,     0.,    1.62,   2.0,   3.,   0.31 };
	Params F_UT_sin_phis         = { 1.18,   0.007,   1.8,   6.1,   0.,     0.,     0.,    0.84,   2.0,   0.,   0.25 };
	Params F_UT_sin_phih_p_phis  = { 0.24,   0.111,   1.7,   6.2,   0.,     0.,     0.,    0.22,   2.9,   1.,   0.32 };

	Params F_LU_sin_phih         = { 1.19,   0.006,   1.4,   4.1,   0.,     0.,     0.,    0.56,   2.2,   1.,   0.33 }; // Completely made up.

	Params F_LL                  = { 0.26,   0.485,   1.3,   3.4,   0.,     0.,     0.,   -0.80,   2.3,   0.,   0.33 };
	Params F_LL_cos_phih         = { 1.13,  -0.285,   1.3,   3.3,   0.,     0.,     0.,    0.15,   2.3,   1.,   0.33 };

	Params F_LT_cos_phih_m_phis  = { 0.24,   0.218,   1.6,   4.7,   0.,     0.,     0.,    0.11,   2.3,   1.,   0.33 };
	Params F_LT_cos_2phih_m_phis = { 1.42,  -0.034,   1.7,   4.6,   0.,     0.,     0.,    1.06,   2.4,   2.,   0.34 };
	Params F_LT_cos_phis         = { 1.14,  -0.104,   1.6,   4.6,   0.,     0.,     0.,   -0.78,   2.1,   0.,   0.33 };
};

/**
 * Structure function set that uses the TestSfSetParams parameterization. Note
 * that this parameterization is extremely inaccurate, and shouldn't be used for
 * any actual physics purposes.
 *
 * \sa TestSfSetParams
 */
class TestSfSet : public SfSet {
	TestSfSetParams _params;
public:
	TestSfSet(TestSfSetParams params=TestSfSetParams()) : SfSet(part::Nucleus::P), _params(params) { }
	TestSfSet(TestSfSet const& other) : TestSfSet(other._params) { }
	TestSfSet& operator=(TestSfSet const& other) = delete;
	virtual ~TestSfSet() = default;

	Real F_UUL(part::Hadron, Real x, Real z, Real Q_sq, Real ph_t_sq) const override {
		return _params.F_UUL(x, z, Q_sq, ph_t_sq);
	}
	Real F_UUT(part::Hadron, Real x, Real z, Real Q_sq, Real ph_t_sq) const override {
		return _params.F_UUT(x, z, Q_sq, ph_t_sq);
	}
	Real F_UU_cos_phih(part::Hadron, Real x, Real z, Real Q_sq, Real ph_t_sq) const override {
		return _params.F_UU_cos_phih(x, z, Q_sq, ph_t_sq);
	}
	Real F_UU_cos_2phih(part::Hadron, Real x, Real z, Real Q_sq, Real ph_t_sq) const override {
		return _params.F_UU_cos_2phih(x, z, Q_sq, ph_t_sq);
	}

	Real F_UL_sin_phih(part::Hadron, Real x, Real z, Real Q_sq, Real ph_t_sq) const override {
		return _params.F_UL_sin_phih(x, z, Q_sq, ph_t_sq);
	}
	Real F_UL_sin_2phih(part::Hadron, Real x, Real z, Real Q_sq, Real ph_t_sq) const override {
		return _params.F_UL_sin_2phih(x, z, Q_sq, ph_t_sq);
	}

	Real F_UTL_sin_phih_m_phis(part::Hadron, Real x, Real z, Real Q_sq, Real ph_t_sq) const override {
		return _params.F_UTL_sin_phih_m_phis(x, z, Q_sq, ph_t_sq);
	}
	Real F_UTT_sin_phih_m_phis(part::Hadron, Real x, Real z, Real Q_sq, Real ph_t_sq) const override {
		return _params.F_UTT_sin_phih_m_phis(x, z, Q_sq, ph_t_sq);
	}
	Real F_UT_sin_2phih_m_phis(part::Hadron, Real x, Real z, Real Q_sq, Real ph_t_sq) const override {
		return _params.F_UT_sin_2phih_m_phis(x, z, Q_sq, ph_t_sq);
	}
	Real F_UT_sin_3phih_m_phis(part::Hadron, Real x, Real z, Real Q_sq, Real ph_t_sq) const override {
		return _params.F_UT_sin_3phih_m_phis(x, z, Q_sq, ph_t_sq);
	}
	Real F_UT_sin_phis(part::Hadron, Real x, Real z, Real Q_sq, Real ph_t_sq) const override {
		return _params.F_UT_sin_phis(x, z, Q_sq, ph_t_sq);
	}
	Real F_UT_sin_phih_p_phis(part::Hadron, Real x, Real z, Real Q_sq, Real ph_t_sq) const override {
		return _params.F_UT_sin_phih_p_phis(x, z, Q_sq, ph_t_sq);
	}

	Real F_LU_sin_phih(part::Hadron, Real x, Real z, Real Q_sq, Real ph_t_sq) const override {
		return _params.F_LU_sin_phih(x, z, Q_sq, ph_t_sq);
	}

	Real F_LL(part::Hadron, Real x, Real z, Real Q_sq, Real ph_t_sq) const override {
		return _params.F_LL(x, z, Q_sq, ph_t_sq);
	}
	Real F_LL_cos_phih(part::Hadron, Real x, Real z, Real Q_sq, Real ph_t_sq) const override {
		return _params.F_LL_cos_phih(x, z, Q_sq, ph_t_sq);
	}

	Real F_LT_cos_phih_m_phis(part::Hadron, Real x, Real z, Real Q_sq, Real ph_t_sq) const override {
		return _params.F_LT_cos_phih_m_phis(x, z, Q_sq, ph_t_sq);
	}
	Real F_LT_cos_2phih_m_phis(part::Hadron, Real x, Real z, Real Q_sq, Real ph_t_sq) const override {
		return _params.F_LT_cos_2phih_m_phis(x, z, Q_sq, ph_t_sq);
	}
	Real F_LT_cos_phis(part::Hadron, Real x, Real z, Real Q_sq, Real ph_t_sq) const override {
		return _params.F_LT_cos_phis(x, z, Q_sq, ph_t_sq);
	}
};

}
}
}

#endif

