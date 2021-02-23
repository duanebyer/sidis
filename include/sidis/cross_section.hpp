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
	struct LepBornLU;
	struct LepBornLP;
	struct LepBornXX;

	struct LepAmmUU;
	struct LepAmmUP;
	struct LepAmmLU;
	struct LepAmmLP;
	struct LepAmmXX;

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
	struct HadLU;
	struct HadLL;
	struct HadLT;
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

Real born(Real lambda_e, math::Vec3 eta, kin::Kinematics const& kin, sf::SfSet const& sf); Real amm(Real lambda_e, math::Vec3 eta, kin::Kinematics const& kin, sf::SfSet const& sf);
Real nrad_ir(Real lambda_e, math::Vec3 eta, kin::Kinematics const& kin, sf::SfSet const& sf, Real k_0_bar=INF);
Real nrad(Real lambda_e, math::Vec3 eta, kin::Kinematics const& kin, sf::SfSet const& sf, Real k_0_bar=INF);
Real rad_f(Real lambda_e, math::Vec3 eta, kin::KinematicsRad const& kin, sf::SfSet const& sf);
Real rad(Real lambda_e, math::Vec3 eta, kin::KinematicsRad const& kin, sf::SfSet const& sf);

// TODO: Include a precision argument, which determines how accurately the
// integrated cross-section is to be calculated.
Real rad_f_integ(Real lambda_e, math::Vec3 eta, kin::Kinematics const& kin, sf::SfSet const& sf, Real k_0_bar=INF);
Real rad_integ(Real lambda_e, math::Vec3 eta, kin::Kinematics const& kin, sf::SfSet const& sf, Real k_0_bar=INF);

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

Real born_xx_base(Real lambda_e, math::Vec3 eta, Born const& b, lep::LepBornXX const& lep, had::HadXX const& had);
Real born_uu_base(Born const& b, lep::LepBornUU const& lep, had::HadUU const& had);
Real born_ul_base(Born const& b, lep::LepBornUP const& lep, had::HadUL const& had);
Real born_ut1_base(Born const& b, lep::LepBornUP const& lep, had::HadUT const& had);
Real born_ut2_base(Born const& b, lep::LepBornUU const& lep, had::HadUT const& had);
Real born_lu_base(Born const& b, lep::LepBornLU const& lep, had::HadLU const& had);
Real born_ll_base(Born const& b, lep::LepBornLP const& lep, had::HadLL const& had);
Real born_lt1_base(Born const& b, lep::LepBornLP const& lep, had::HadLT const& had);
Real born_lt2_base(Born const& b, lep::LepBornLU const& lep, had::HadLT const& had);

struct Amm {
	Real coeff;
	explicit Amm(kin::Kinematics const& kin);
};

Real amm_xx_base(Real lambda_e, math::Vec3 eta, Amm const& b, lep::LepAmmXX const& lep, had::HadXX const& had);
Real amm_uu_base(Amm const& b, lep::LepAmmUU const& lep, had::HadUU const& had);
Real amm_ul_base(Amm const& b, lep::LepAmmUP const& lep, had::HadUL const& had);
Real amm_ut1_base(Amm const& b, lep::LepAmmUP const& lep, had::HadUT const& had);
Real amm_ut2_base(Amm const& b, lep::LepAmmUU const& lep, had::HadUT const& had);
Real amm_lu_base(Amm const& b, lep::LepAmmLU const& lep, had::HadLU const& had);
Real amm_ll_base(Amm const& b, lep::LepAmmLP const& lep, had::HadLL const& had);
Real amm_lt1_base(Amm const& b, lep::LepAmmLP const& lep, had::HadLT const& had);
Real amm_lt2_base(Amm const& b, lep::LepAmmLU const& lep, had::HadLT const& had);

struct NRadIR {
	Real coeff_born;
	Real coeff_amm;
	explicit NRadIR(kin::Kinematics const& kin, Real k_0_bar);
};

Real nrad_ir_xx_base(Real lambda_e, math::Vec3 eta, NRadIR const& b, lep::LepBornXX const& lep_born, lep::LepAmmXX const& lep_amm, had::HadXX const& had);
Real nrad_ir_uu_base(NRadIR const& b, lep::LepBornUU const& lep_born, lep::LepAmmUU const& lep_amm, had::HadUU const& had);
Real nrad_ir_ul_base(NRadIR const& b, lep::LepBornUP const& lep_born, lep::LepAmmUP const& lep_amm, had::HadUL const& had);
Real nrad_ir_ut1_base(NRadIR const& b, lep::LepBornUP const& lep_born, lep::LepAmmUP const& lep_amm, had::HadUT const& had);
Real nrad_ir_ut2_base(NRadIR const& b, lep::LepBornUU const& lep_born, lep::LepAmmUU const& lep_amm, had::HadUT const& had);
Real nrad_ir_lu_base(NRadIR const& b, lep::LepBornLU const& lep_born, lep::LepAmmLU const& lep_amm, had::HadLU const& had);
Real nrad_ir_ll_base(NRadIR const& b, lep::LepBornLP const& lep_born, lep::LepAmmLP const& lep_amm, had::HadLL const& had);
Real nrad_ir_lt1_base(NRadIR const& b, lep::LepBornLP const& lep_born, lep::LepAmmLP const& lep_amm, had::HadLT const& had);
Real nrad_ir_lt2_base(NRadIR const& b, lep::LepBornLU const& lep_born, lep::LepAmmLU const& lep_amm, had::HadLT const& had);

struct Rad {
	Real coeff;
	Real R;
	explicit Rad(kin::KinematicsRad const& kin);
};

Real rad_xx_base(Real lambda_e, math::Vec3 eta, Rad const& b, lep::LepRadXX const& lep, had::HadRadXX const& had);
Real rad_uu_base(Rad const& b, lep::LepRadUU const& lep, had::HadRadUU const& had);
math::Vec3 rad_up_base(Rad const& b, lep::LepRadUX const& lep, had::HadRadUP const& had);
Real rad_lu_base(Rad const& b, lep::LepRadLU const& lep, had::HadRadLU const& had);
math::Vec3 rad_lp_base(Rad const& b, lep::LepRadLX const& lep, had::HadRadLP const& had);

Real rad_f_xx_base(Real lambda_e, math::Vec3 eta, Rad const& b, lep::LepRadXX const& lep, had::HadRadFXX const& had);
Real rad_f_uu_base(Rad const& b, lep::LepRadUU const& lep, had::HadRadFUU const& had);
math::Vec3 rad_f_up_base(Rad const& b, lep::LepRadUX const& lep, had::HadRadFUP const& had);
Real rad_f_lu_base(Rad const& b, lep::LepRadLU const& lep, had::HadRadFLU const& had);
math::Vec3 rad_f_lp_base(Rad const& b, lep::LepRadLX const& lep, had::HadRadFLP const& had);

}
}

#endif

