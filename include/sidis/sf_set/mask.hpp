#ifndef SIDIS_SF_SET_MASK_HPP
#define SIDIS_SF_SET_MASK_HPP

#include <cstddef>
#include <memory>
#include <utility>

#include "sidis/structure_function.hpp"

namespace sidis {
namespace sf {
namespace set {

std::size_t const NUM_SF = 18;

/// Mask that selects only leading structure functions.
bool const MASK_LEADING[NUM_SF] {
	0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0,
};
/// Mask that selects only subleading structure functions.
bool const MASK_SUBLEADING[NUM_SF] {
	1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1,
};
/// Masks for different structure function polarizations.
/// \{
bool const MASK_UU[NUM_SF] {
	1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
};
bool const MASK_UL[NUM_SF] {
	0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
};
bool const MASK_UT[NUM_SF] {
	0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
};
bool const MASK_UP[NUM_SF] {
	0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
};
bool const MASK_UX[NUM_SF] {
	1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
};
bool const MASK_LU[NUM_SF] {
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
};
bool const MASK_LL[NUM_SF] {
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0,
};
bool const MASK_LT[NUM_SF] {
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1,
};
bool const MASK_LP[NUM_SF] {
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1,
};
bool const MASK_LX[NUM_SF] {
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1,
};
bool const MASK_XU[NUM_SF] {
	1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
};
bool const MASK_XL[NUM_SF] {
	0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0,
};
bool const MASK_XT[NUM_SF] {
	0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1,
};
bool const MASK_XP[NUM_SF] {
	0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1,
};
bool const MASK_XX[NUM_SF] {
	1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
};
/// \}


/**
 * Wrapper around another structure function set that only calculates certain
 * structure functions (for example, leading order only). Any of the 18 leading
 * and sub-leading structure functions can be enabled/disabled using this
 * wrapper. The indices for each structure function are as follows:
 *
 * | index | structure function                  | leading order? |
 * |------:|:------------------------------------|:---------------|
 * |     0 | \f$F_{UUL}\f$                       | no             |
 * |     1 | \f$F_{UUT}\f$                       | **yes**        |
 * |     2 | \f$F_{UU}^{\cos{\phi_h}}\f$         | no             |
 * |     3 | \f$F_{UU}^{\cos{2\phi_h}}\f$        | **yes**        |
 * |     4 | \f$F_{UL}^{\sin{\phi_h}}\f$         | no             |
 * |     5 | \f$F_{UL}^{\sin{2\phi_h}}\f$        | **yes**        |
 * |     6 | \f$F_{UTL}^{\sin{\phi_h-\phi_S}}\f$ | no             |
 * |     7 | \f$F_{UTT}^{\sin{\phi_h-\phi_S}}\f$ | **yes**        |
 * |     8 | \f$F_{UT}^{\sin{2\phi_h-\phi_S}}\f$ | no             |
 * |     9 | \f$F_{UT}^{\sin{3\phi_h-\phi_S}}\f$ | **yes**        |
 * |    10 | \f$F_{UT}^{\sin{\phi_S}}\f$         | no             |
 * |    11 | \f$F_{UT}^{\sin{\phi_h+\phi_S}}\f$  | **yes**        |
 * |    12 | \f$F_{LU}^{\sin{\phi_h}}\f$         | no             |
 * |    13 | \f$F_{LL}\f$                        | **yes**        |
 * |    14 | \f$F_{LL}^{\cos{\phi_h}}\f$         | no             |
 * |    15 | \f$F_{LT}^{\cos{\phi_h-\phi_S}}\f$  | **yes**        |
 * |    16 | \f$F_{LT}^{\cos{2\phi_h-\phi_S}}\f$ | no             |
 * |    17 | \f$F_{LT}^{\cos{\phi_S}}\f$         | no             |
 */
class MaskSfSet final : public SfSet {
	bool _mask[NUM_SF];
	std::unique_ptr<SfSet> _sf;

public:
	/// Constructs an SfSet for \p target with a \p mask determining which
	/// structure functions will be non-zero. See the class documentation for a
	/// table describing which structure functions correspond to which indices.
	MaskSfSet(const bool (&mask)[NUM_SF], std::unique_ptr<SfSet>&& sf) noexcept :
			SfSet(sf->target),
			_mask(),
			_sf(std::move(sf)) {
		// This is a dumb way of avoiding importing `size_t` header.
		for (decltype(sizeof(bool)) idx = 0; idx < sizeof(mask) / sizeof(bool); ++idx) {
			_mask[idx] = mask[idx];
		}
	}
	/// Constructs a MaskSfSet, taking ownership of the pointer \p sf.
	MaskSfSet(const bool (&mask)[NUM_SF], SfSet* sf) :
			SfSet(sf->target),
			_mask(),
			_sf(sf) {
		for (decltype(sizeof(bool)) idx = 0; idx < sizeof(mask) / sizeof(bool); ++idx) {
			_mask[idx] = mask[idx];
		}
	}
	MaskSfSet(MaskSfSet const& other) = delete;
	MaskSfSet(MaskSfSet&& other) noexcept : MaskSfSet(other._mask, std::move(other._sf)) { }
	MaskSfSet& operator=(MaskSfSet const& other) = delete;
	MaskSfSet& operator=(MaskSfSet&& other) = delete;
	virtual ~MaskSfSet() = default;

	Real F_UUL(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override {
		if (_mask[0]) {
			return _sf->F_UUL(h, x, z, Q_sq, ph_t_sq);
		} else {
			return 0.;
		}
	}
	Real F_UUT(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override {
		if (_mask[1]) {
			return _sf->F_UUT(h, x, z, Q_sq, ph_t_sq);
		} else {
			return 0.;
		}
	}
	Real F_UU_cos_phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override {
		if (_mask[2]) {
			return _sf->F_UU_cos_phih(h, x, z, Q_sq, ph_t_sq);
		} else {
			return 0.;
		}
	}
	Real F_UU_cos_2phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override {
		if (_mask[3]) {
			return _sf->F_UU_cos_2phih(h, x, z, Q_sq, ph_t_sq);
		} else {
			return 0.;
		}
	}

	Real F_UL_sin_phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override {
		if (_mask[4]) {
			return _sf->F_UL_sin_phih(h, x, z, Q_sq, ph_t_sq);
		} else {
			return 0.;
		}
	}
	Real F_UL_sin_2phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override {
		if (_mask[5]) {
			return _sf->F_UL_sin_2phih(h, x, z, Q_sq, ph_t_sq);
		} else {
			return 0.;
		}
	}

	Real F_UTL_sin_phih_m_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override {
		if (_mask[6]) {
			return _sf->F_UTL_sin_phih_m_phis(h, x, z, Q_sq, ph_t_sq);
		} else {
			return 0.;
		}
	}
	Real F_UTT_sin_phih_m_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override {
		if (_mask[7]) {
			return _sf->F_UTT_sin_phih_m_phis(h, x, z, Q_sq, ph_t_sq);
		} else {
			return 0.;
		}
	}
	Real F_UT_sin_2phih_m_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override {
		if (_mask[8]) {
			return _sf->F_UT_sin_2phih_m_phis(h, x, z, Q_sq, ph_t_sq);
		} else {
			return 0.;
		}
	}
	Real F_UT_sin_3phih_m_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override {
		if (_mask[9]) {
			return _sf->F_UT_sin_3phih_m_phis(h, x, z, Q_sq, ph_t_sq);
		} else {
			return 0.;
		}
	}
	Real F_UT_sin_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override {
		if (_mask[10]) {
			return _sf->F_UT_sin_phis(h, x, z, Q_sq, ph_t_sq);
		} else {
			return 0.;
		}
	}
	Real F_UT_sin_phih_p_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override {
		if (_mask[11]) {
			return _sf->F_UT_sin_phih_p_phis(h, x, z, Q_sq, ph_t_sq);
		} else {
			return 0.;
		}
	}

	Real F_LU_sin_phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override {
		if (_mask[12]) {
			return _sf->F_LU_sin_phih(h, x, z, Q_sq, ph_t_sq);
		} else {
			return 0.;
		}
	}

	Real F_LL(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override {
		if (_mask[13]) {
			return _sf->F_LL(h, x, z, Q_sq, ph_t_sq);
		} else {
			return 0.;
		}
	}
	Real F_LL_cos_phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override {
		if (_mask[14]) {
			return _sf->F_LL_cos_phih(h, x, z, Q_sq, ph_t_sq);
		} else {
			return 0.;
		}
	}

	Real F_LT_cos_phih_m_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override {
		if (_mask[15]) {
			return _sf->F_LT_cos_phih_m_phis(h, x, z, Q_sq, ph_t_sq);
		} else {
			return 0.;
		}
	}
	Real F_LT_cos_2phih_m_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override {
		if (_mask[16]) {
			return _sf->F_LT_cos_2phih_m_phis(h, x, z, Q_sq, ph_t_sq);
		} else {
			return 0.;
		}
	}
	Real F_LT_cos_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override {
		if (_mask[17]) {
			return _sf->F_LT_cos_phis(h, x, z, Q_sq, ph_t_sq);
		} else {
			return 0.;
		}
	}
};

}
}
}

#endif

