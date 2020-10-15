#include "sidis/structure_function.hpp"

#include "sidis/constant.hpp"

using namespace sidis;
using namespace constant;
using namespace sidis::sf;

Real Model::F_UUL(Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	static_cast<void>(h);
	static_cast<void>(x);
	static_cast<void>(z);
	static_cast<void>(Q_sq);
	static_cast<void>(ph_t_sq);
	return 0.;
}
Real Model::F_UUT(Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	static_cast<void>(h);
	static_cast<void>(x);
	static_cast<void>(z);
	static_cast<void>(Q_sq);
	static_cast<void>(ph_t_sq);
	return 0.;
}
Real Model::F_UU_cos_phih(Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	static_cast<void>(h);
	static_cast<void>(x);
	static_cast<void>(z);
	static_cast<void>(Q_sq);
	static_cast<void>(ph_t_sq);
	return 0.;
}
Real Model::F_UU_cos_2phih(Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	static_cast<void>(h);
	static_cast<void>(x);
	static_cast<void>(z);
	static_cast<void>(Q_sq);
	static_cast<void>(ph_t_sq);
	return 0.;
}

Real Model::F_UL_sin_phih(Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	static_cast<void>(h);
	static_cast<void>(x);
	static_cast<void>(z);
	static_cast<void>(Q_sq);
	static_cast<void>(ph_t_sq);
	return 0.;
}
Real Model::F_UL_sin_2phih(Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	static_cast<void>(h);
	static_cast<void>(x);
	static_cast<void>(z);
	static_cast<void>(Q_sq);
	static_cast<void>(ph_t_sq);
	return 0.;
}

Real Model::F_UTL_sin_phih_m_phis(Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	static_cast<void>(h);
	static_cast<void>(x);
	static_cast<void>(z);
	static_cast<void>(Q_sq);
	static_cast<void>(ph_t_sq);
	return 0.;
}
Real Model::F_UTT_sin_phih_m_phis(Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	static_cast<void>(h);
	static_cast<void>(x);
	static_cast<void>(z);
	static_cast<void>(Q_sq);
	static_cast<void>(ph_t_sq);
	return 0.;
}
Real Model::F_UT_sin_2phih_m_phis(Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	static_cast<void>(h);
	static_cast<void>(x);
	static_cast<void>(z);
	static_cast<void>(Q_sq);
	static_cast<void>(ph_t_sq);
	return 0.;
}
Real Model::F_UT_sin_3phih_m_phis(Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	static_cast<void>(h);
	static_cast<void>(x);
	static_cast<void>(z);
	static_cast<void>(Q_sq);
	static_cast<void>(ph_t_sq);
	return 0.;
}
Real Model::F_UT_sin_phis(Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	static_cast<void>(h);
	static_cast<void>(x);
	static_cast<void>(z);
	static_cast<void>(Q_sq);
	static_cast<void>(ph_t_sq);
	return 0.;
}
Real Model::F_UT_sin_phih_p_phis(Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	static_cast<void>(h);
	static_cast<void>(x);
	static_cast<void>(z);
	static_cast<void>(Q_sq);
	static_cast<void>(ph_t_sq);
	return 0.;
}

Real Model::F_LU_sin_phih(Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	static_cast<void>(h);
	static_cast<void>(x);
	static_cast<void>(z);
	static_cast<void>(Q_sq);
	static_cast<void>(ph_t_sq);
	return 0.;
}

Real Model::F_LL(Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	static_cast<void>(h);
	static_cast<void>(x);
	static_cast<void>(z);
	static_cast<void>(Q_sq);
	static_cast<void>(ph_t_sq);
	return 0.;
}
Real Model::F_LL_cos_phih(Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	static_cast<void>(h);
	static_cast<void>(x);
	static_cast<void>(z);
	static_cast<void>(Q_sq);
	static_cast<void>(ph_t_sq);
	return 0.;
}

Real Model::F_LT_cos_phih_m_phis(Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	static_cast<void>(h);
	static_cast<void>(x);
	static_cast<void>(z);
	static_cast<void>(Q_sq);
	static_cast<void>(ph_t_sq);
	return 0.;
}
Real Model::F_LT_cos_2phih_m_phis(Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	static_cast<void>(h);
	static_cast<void>(x);
	static_cast<void>(z);
	static_cast<void>(Q_sq);
	static_cast<void>(ph_t_sq);
	return 0.;
}
Real Model::F_LT_cos_phis(Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	static_cast<void>(h);
	static_cast<void>(x);
	static_cast<void>(z);
	static_cast<void>(Q_sq);
	static_cast<void>(ph_t_sq);
	return 0.;
}

SfUU Model::sf_uu(Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return {
		F_UUL(h, x, z, Q_sq, ph_t_sq),
		F_UUT(h, x, z, Q_sq, ph_t_sq),
		F_UU_cos_phih(h, x, z, Q_sq, ph_t_sq),
		F_UU_cos_2phih(h, x, z, Q_sq, ph_t_sq),
	};
}

SfUL Model::sf_ul(Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return {
		F_UL_sin_phih(h, x, z, Q_sq, ph_t_sq),
		F_UL_sin_2phih(h, x, z, Q_sq, ph_t_sq),
	};
}

SfUT Model::sf_ut(Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return {
		F_UTL_sin_phih_m_phis(h, x, z, Q_sq, ph_t_sq),
		F_UTT_sin_phih_m_phis(h, x, z, Q_sq, ph_t_sq),
		F_UT_sin_2phih_m_phis(h, x, z, Q_sq, ph_t_sq),
		F_UT_sin_3phih_m_phis(h, x, z, Q_sq, ph_t_sq),
		F_UT_sin_phis(h, x, z, Q_sq, ph_t_sq),
		F_UT_sin_phih_p_phis(h, x, z, Q_sq, ph_t_sq),
	};
}

SfLU Model::sf_lu(Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return {
		F_LU_sin_phih(h, x, z, Q_sq, ph_t_sq),
	};
}

SfLL Model::sf_ll(Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return {
		F_LL(h, x, z, Q_sq, ph_t_sq),
		F_LL_cos_phih(h, x, z, Q_sq, ph_t_sq),
	};
}

SfLT Model::sf_lt(Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return {
		F_LT_cos_phih_m_phis(h, x, z, Q_sq, ph_t_sq),
		F_LT_cos_2phih_m_phis(h, x, z, Q_sq, ph_t_sq),
		F_LT_cos_phis(h, x, z, Q_sq, ph_t_sq),
	};
}

SfXU Model::sf_xu(Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return {
		sf_uu(h, x, z, Q_sq, ph_t_sq),
		sf_lu(h, x, z, Q_sq, ph_t_sq),
	};
}

SfXL Model::sf_xl(Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return {
		sf_ul(h, x, z, Q_sq, ph_t_sq),
		sf_ll(h, x, z, Q_sq, ph_t_sq),
	};
}

SfXT Model::sf_xt(Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return {
		sf_ut(h, x, z, Q_sq, ph_t_sq),
		sf_lt(h, x, z, Q_sq, ph_t_sq),
	};
}

SfUP Model::sf_up(Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return {
		sf_ul(h, x, z, Q_sq, ph_t_sq),
		sf_ut(h, x, z, Q_sq, ph_t_sq),
	};
}

SfLP Model::sf_lp(Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return {
		sf_ll(h, x, z, Q_sq, ph_t_sq),
		sf_lt(h, x, z, Q_sq, ph_t_sq),
	};
}

SfUX Model::sf_ux(Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return {
		sf_uu(h, x, z, Q_sq, ph_t_sq),
		sf_ul(h, x, z, Q_sq, ph_t_sq),
		sf_ut(h, x, z, Q_sq, ph_t_sq),
	};
}

SfLX Model::sf_lx(Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return {
		sf_lu(h, x, z, Q_sq, ph_t_sq),
		sf_ll(h, x, z, Q_sq, ph_t_sq),
		sf_lt(h, x, z, Q_sq, ph_t_sq),
	};
}

SfXP Model::sf_xp(Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return {
		sf_up(h, x, z, Q_sq, ph_t_sq),
		sf_lp(h, x, z, Q_sq, ph_t_sq),
	};
}

SfXX Model::sf_xx(Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return {
		sf_ux(h, x, z, Q_sq, ph_t_sq),
		sf_lx(h, x, z, Q_sq, ph_t_sq),
	};
}

SfXX Model::sf(Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return sf_xx(h, x, z, Q_sq, ph_t_sq);
}

