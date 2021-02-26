#ifndef SIDIS_CROSS_SECTION_HPP
#define SIDIS_CROSS_SECTION_HPP

#include "sidis/constant.hpp"
#include "sidis/numeric.hpp"

namespace sidis {

namespace math {
	struct Vec3;
	struct Transform3;
}
namespace kin {
	struct Kinematics;
	struct KinematicsRad;
}
namespace lep {
	struct LepBornUU;
	struct LepBornUP;
	struct LepBornUX;
	struct LepBornLU;
	struct LepBornLP;
	struct LepBornLX;
	struct LepBornXX;

	struct LepAmmUU;
	struct LepAmmUP;
	struct LepAmmUX;
	struct LepAmmLU;
	struct LepAmmLP;
	struct LepAmmLX;
	struct LepAmmXX;

	struct LepNradUU;
	struct LepNradUP;
	struct LepNradUX;
	struct LepNradLU;
	struct LepNradLP;
	struct LepNradLX;
	struct LepNradXX;

	struct LepRadUU;
	struct LepRadUP;
	struct LepRadUX;
	struct LepRadLU;
	struct LepRadLP;
	struct LepRadLX;
	struct LepRadXX;
}
namespace had {
	struct HadUU;
	struct HadUL;
	struct HadUT;
	struct HadUP;
	struct HadLU;
	struct HadLL;
	struct HadLT;
	struct HadLP;
	struct HadXX;

	struct HadRadUU;
	struct HadRadUP;
	struct HadRadLU;
	struct HadRadLP;
	struct HadRadXX;

	struct HadRadFUU;
	struct HadRadFUP;
	struct HadRadFLU;
	struct HadRadFLP;
	struct HadRadFXX;
}
namespace sf {
	struct SfSet;
}

namespace xs {

Real born(kin::Kinematics const& kin, sf::SfSet const& sf, Real lambda_e, math::Vec3 eta);
Real amm(kin::Kinematics const& kin, sf::SfSet const& sf, Real lambda_e, math::Vec3 eta);
Real nrad_ir(kin::Kinematics const& kin, sf::SfSet const& sf, Real lambda_e, math::Vec3 eta, Real k_0_bar=INF);
Real nrad(kin::Kinematics const& kin, sf::SfSet const& sf, Real lambda_e, math::Vec3 eta, Real k_0_bar=INF, unsigned max_evals=1000000, Real prec=1e-6);
Real rad_f(kin::KinematicsRad const& kin, sf::SfSet const& sf, Real lambda_e, math::Vec3 eta);
Real rad(kin::KinematicsRad const& kin, sf::SfSet const& sf, Real lambda_e, math::Vec3 eta);

Real rad_f_integ(kin::Kinematics const& kin, sf::SfSet const& sf, Real lambda_e, math::Vec3 eta, Real k_0_bar=INF, unsigned max_evals=1000000, Real prec=1e-6);
Real rad_integ(kin::Kinematics const& kin, sf::SfSet const& sf, Real lambda_e, math::Vec3 eta, Real k_0_bar=INF, unsigned max_evals=1000000, Real prec=1e-6);

// Various corrections to Born cross-section.
Real delta_vert_rad_ir(kin::Kinematics const& kin, Real k_0_bar=INF);
Real delta_rad_ir_hard(kin::Kinematics const& kin, Real k_0_bar=INF);
Real delta_vac_lep(kin::Kinematics const& kin);
Real delta_vac_had(kin::Kinematics const& kin);

// Base functions, for calculating cross-sections directly from leptonic and
// hadronic coefficients.
struct Born {
	Real coeff;
	explicit Born(kin::Kinematics const& kin);
};

Real born_uu_base(Born const& b, lep::LepBornUU const& lep, had::HadUU const& had);
Real born_ul_base(Born const& b, lep::LepBornUP const& lep, had::HadUL const& had);
Real born_ut1_base(Born const& b, lep::LepBornUP const& lep, had::HadUT const& had);
Real born_ut2_base(Born const& b, lep::LepBornUU const& lep, had::HadUT const& had);
math::Vec3 born_up_base(Born const& b, lep::LepBornUX const& lep, had::HadUP const& had);
Real born_lu_base(Born const& b, lep::LepBornLU const& lep, had::HadLU const& had);
Real born_ll_base(Born const& b, lep::LepBornLP const& lep, had::HadLL const& had);
Real born_lt1_base(Born const& b, lep::LepBornLP const& lep, had::HadLT const& had);
Real born_lt2_base(Born const& b, lep::LepBornLU const& lep, had::HadLT const& had);
math::Vec3 born_lp_base(Born const& b, lep::LepBornLX const& lep, had::HadLP const& had);

struct Amm {
	Real coeff;
	explicit Amm(kin::Kinematics const& kin);
};

Real amm_uu_base(Amm const& b, lep::LepAmmUU const& lep, had::HadUU const& had);
Real amm_ul_base(Amm const& b, lep::LepAmmUP const& lep, had::HadUL const& had);
Real amm_ut1_base(Amm const& b, lep::LepAmmUP const& lep, had::HadUT const& had);
Real amm_ut2_base(Amm const& b, lep::LepAmmUU const& lep, had::HadUT const& had);
math::Vec3 amm_up_base(Amm const& b, lep::LepAmmUX const& lep, had::HadUP const& had);
Real amm_lu_base(Amm const& b, lep::LepAmmLU const& lep, had::HadLU const& had);
Real amm_ll_base(Amm const& b, lep::LepAmmLP const& lep, had::HadLL const& had);
Real amm_lt1_base(Amm const& b, lep::LepAmmLP const& lep, had::HadLT const& had);
Real amm_lt2_base(Amm const& b, lep::LepAmmLU const& lep, had::HadLT const& had);
math::Vec3 amm_lp_base(Amm const& b, lep::LepAmmLX const& lep, had::HadLP const& had);

struct Nrad {
	Real coeff_born;
	Real coeff_amm;
	explicit Nrad(kin::Kinematics const& kin, Real k_0_bar);
};

Real nrad_ir_uu_base(Nrad const& b, lep::LepNradUU const& lep, had::HadUU const& had);
Real nrad_ir_ul_base(Nrad const& b, lep::LepNradUP const& lep, had::HadUL const& had);
Real nrad_ir_ut1_base(Nrad const& b, lep::LepNradUP const& lep, had::HadUT const& had);
Real nrad_ir_ut2_base(Nrad const& b, lep::LepNradUU const& lep, had::HadUT const& had);
math::Vec3 nrad_ir_up_base(Nrad const& b, lep::LepNradUX const& lep, had::HadUP const& had);
Real nrad_ir_lu_base(Nrad const& b, lep::LepNradLU const& lep, had::HadLU const& had);
Real nrad_ir_ll_base(Nrad const& b, lep::LepNradLP const& lep, had::HadLL const& had);
Real nrad_ir_lt1_base(Nrad const& b, lep::LepNradLP const& lep, had::HadLT const& had);
Real nrad_ir_lt2_base(Nrad const& b, lep::LepNradLU const& lep, had::HadLT const& had);
math::Vec3 nrad_ir_lp_base(Nrad const& b, lep::LepNradLX const& lep, had::HadLP const& had);

struct Rad {
	Real coeff;
	Real R;
	explicit Rad(kin::KinematicsRad const& kin);
};

Real rad_uu_base(Rad const& b, lep::LepRadUU const& lep, had::HadRadUU const& had);
math::Vec3 rad_up_base(Rad const& b, lep::LepRadUX const& lep, had::HadRadUP const& had);
Real rad_lu_base(Rad const& b, lep::LepRadLU const& lep, had::HadRadLU const& had);
math::Vec3 rad_lp_base(Rad const& b, lep::LepRadLX const& lep, had::HadRadLP const& had);

Real rad_f_uu_base(Rad const& b, lep::LepRadUU const& lep, had::HadRadFUU const& had);
math::Vec3 rad_f_up_base(Rad const& b, lep::LepRadUX const& lep, had::HadRadFUP const& had);
Real rad_f_lu_base(Rad const& b, lep::LepRadLU const& lep, had::HadRadFLU const& had);
math::Vec3 rad_f_lp_base(Rad const& b, lep::LepRadLX const& lep, had::HadRadFLP const& had);

}
}

#endif

