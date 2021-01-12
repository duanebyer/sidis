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

Real born(Real lambda_e, math::Vec3 eta, kin::Kinematics kin, sf::SfSet const& model); Real amm(Real lambda_e, math::Vec3 eta, kin::Kinematics kin, sf::SfSet const& model);
Real nrad_ir(Real lambda_e, math::Vec3 eta, kin::Kinematics kin, sf::SfSet const& model, Real k0_cut=constant::INF);
Real nrad(Real lambda_e, math::Vec3 eta, kin::Kinematics kin, sf::SfSet const& model, Real k0_cut=constant::INF);
Real rad_f(Real lambda_e, math::Vec3 eta, kin::KinematicsRad kin, sf::SfSet const& model);
Real rad(Real lambda_e, math::Vec3 eta, kin::KinematicsRad kin, sf::SfSet const& model);

// TODO: Include a precision argument, which determines how accurately the
// integrated cross-section is to be calculated.
Real rad_f_integ(Real lambda_e, math::Vec3 eta, kin::Kinematics kin, sf::SfSet const& model, Real k0_cut=constant::INF);
Real rad_integ(Real lambda_e, math::Vec3 eta, kin::Kinematics kin, sf::SfSet const& model, Real k0_cut=constant::INF);

// Various corrections to Born cross-section.
Real delta_vert_rad_ir(kin::Kinematics kin, Real k0_cut=constant::INF);
Real delta_rad_ir_hard(kin::Kinematics kin, Real k0_cut=constant::INF);
Real delta_vac_lep(kin::Kinematics kin);
Real delta_vac_had(kin::Kinematics kin);

// Base functions, for calculating cross-sections directly from leptonic and
// hadronic coefficients.
struct Born {
	Real coeff;
	explicit Born(kin::Kinematics kin);
};

Real born_xx_base(Real lambda_e, math::Vec3 eta, Born b, lep::LepBornXX lep, had::HadXX had);
Real born_uu_base(Born b, lep::LepBornUU lep, had::HadUU had);
Real born_ul_base(Born b, lep::LepBornUP lep, had::HadUL had);
Real born_ut1_base(Born b, lep::LepBornUP lep, had::HadUT had);
Real born_ut2_base(Born b, lep::LepBornUU lep, had::HadUT had);
Real born_lu_base(Born b, lep::LepBornLU lep, had::HadLU had);
Real born_ll_base(Born b, lep::LepBornLP lep, had::HadLL had);
Real born_lt1_base(Born b, lep::LepBornLP lep, had::HadLT had);
Real born_lt2_base(Born b, lep::LepBornLU lep, had::HadLT had);

struct Amm {
	Real coeff;
	explicit Amm(kin::Kinematics kin);
};

Real amm_xx_base(Real lambda_e, math::Vec3 eta, Amm b, lep::LepAmmXX lep, had::HadXX had);
Real amm_uu_base(Amm b, lep::LepAmmUU lep, had::HadUU had);
Real amm_ul_base(Amm b, lep::LepAmmUP lep, had::HadUL had);
Real amm_ut1_base(Amm b, lep::LepAmmUP lep, had::HadUT had);
Real amm_ut2_base(Amm b, lep::LepAmmUU lep, had::HadUT had);
Real amm_lu_base(Amm b, lep::LepAmmLU lep, had::HadLU had);
Real amm_ll_base(Amm b, lep::LepAmmLP lep, had::HadLL had);
Real amm_lt1_base(Amm b, lep::LepAmmLP lep, had::HadLT had);
Real amm_lt2_base(Amm b, lep::LepAmmLU lep, had::HadLT had);

struct NRadIR {
	Real coeff_born;
	Real coeff_amm;
	explicit NRadIR(kin::Kinematics kin, Real k0_cut);
};

Real nrad_ir_xx_base(Real lambda_e, math::Vec3 eta, NRadIR b, lep::LepBornXX lep_born, lep::LepAmmXX lep_amm, had::HadXX had);
Real nrad_ir_uu_base(NRadIR b, lep::LepBornUU lep_born, lep::LepAmmUU lep_amm, had::HadUU had);
Real nrad_ir_ul_base(NRadIR b, lep::LepBornUP lep_born, lep::LepAmmUP lep_amm, had::HadUL had);
Real nrad_ir_ut1_base(NRadIR b, lep::LepBornUP lep_born, lep::LepAmmUP lep_amm, had::HadUT had);
Real nrad_ir_ut2_base(NRadIR b, lep::LepBornUU lep_born, lep::LepAmmUU lep_amm, had::HadUT had);
Real nrad_ir_lu_base(NRadIR b, lep::LepBornLU lep_born, lep::LepAmmLU lep_amm, had::HadLU had);
Real nrad_ir_ll_base(NRadIR b, lep::LepBornLP lep_born, lep::LepAmmLP lep_amm, had::HadLL had);
Real nrad_ir_lt1_base(NRadIR b, lep::LepBornLP lep_born, lep::LepAmmLP lep_amm, had::HadLT had);
Real nrad_ir_lt2_base(NRadIR b, lep::LepBornLU lep_born, lep::LepAmmLU lep_amm, had::HadLT had);

struct Rad {
	Real coeff;
	Real R;
	explicit Rad(kin::KinematicsRad kin);
};

Real rad_xx_base(Real lambda_e, math::Vec3 eta, Rad b, lep::LepRadXX lep, had::HadRadXX had);
Real rad_uu_base(Rad b, lep::LepRadUU lep, had::HadRadUU had);
math::Vec3 rad_up_base(Rad b, lep::LepRadUX lep, had::HadRadUP had);
Real rad_lu_base(Rad b, lep::LepRadLU lep, had::HadRadLU had);
math::Vec3 rad_lp_base(Rad b, lep::LepRadLX lep, had::HadRadLP had);

Real rad_f_xx_base(Real lambda_e, math::Vec3 eta, Rad b, lep::LepRadXX lep, had::HadRadFXX had);
Real rad_f_uu_base(Rad b, lep::LepRadUU lep, had::HadRadFUU had);
math::Vec3 rad_f_up_base(Rad b, lep::LepRadUX lep, had::HadRadFUP had);
Real rad_f_lu_base(Rad b, lep::LepRadLU lep, had::HadRadFLU had);
math::Vec3 rad_f_lp_base(Rad b, lep::LepRadLX lep, had::HadRadFLP had);

}
}

#endif

