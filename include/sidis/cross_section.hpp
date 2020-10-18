#ifndef SIDIS_CROSS_SECTION_HPP
#define SIDIS_CROSS_SECTION_HPP

#include "sidis/numeric.hpp"

namespace sidis {

namespace math {
	struct Vec3;
	struct Transform3;
}
namespace kin {
	struct Kinematics;
	struct KinematicsRad;
	struct KinematicsEx;
}
namespace lep {
	struct LepBornUU;
	struct LepBornUP;
	struct LepBornLU;
	struct LepBornLP;
	struct LepAmmUU;
	struct LepAmmUP;
	struct LepAmmLU;
	struct LepAmmLP;
	struct LepRadUU;
	struct LepRadUP;
	struct LepRadUX;
	struct LepRadLU;
	struct LepRadLP;
	struct LepRadLX;
}
namespace had {
	struct HadUU;
	struct HadUL;
	struct HadUT1;
	struct HadUT2;
	struct HadUX;
	struct HadUP;
	struct HadLU;
	struct HadLL;
	struct HadLT1;
	struct HadLT2;
	struct HadLX;
	struct HadLP;
}
namespace sf {
	struct Model;
}

namespace xs {

struct Born {
	Real coeff;
	explicit Born(kin::Kinematics kin);
};

struct Amm {
	Real coeff;
	explicit Amm(kin::Kinematics kin);
};

struct Rad {
	Real coeff;
	Real R;
	explicit Rad(kin::KinematicsRad kin);
};

Real born(Real lambda_e, math::Vec3 eta, kin::Kinematics kin, sf::Model const& model);
Real born_rad_factor(kin::Kinematics kin);
Real delta_vr(kin::Kinematics kin);
Real delta_vac_lep(kin::Kinematics kin);
Real delta_vac_had(kin::Kinematics kin);
Real amm(Real lambda_e, math::Vec3 eta, kin::Kinematics kin, sf::Model const& model);
Real nrad(Real lambda_e, math::Vec3 eta, kin::Kinematics kin, sf::Model const& model);

Real rad(Real lambda_e, math::Vec3 eta, kin::KinematicsRad kin, sf::Model const& model);
// TODO: Include a precision argument, which determines how accurately the
// integrated cross-section is to be calculated.
Real rad_integ(Real lambda_e, math::Vec3 eta, kin::Kinematics kin, sf::Model const& model);
Real rad_hard(Real lambda_e, math::Vec3 eta, kin::KinematicsRad kin, sf::Model const& model);
Real rad_soft(Real lambda_e, math::Vec3 eta, kin::KinematicsRad kin, sf::Model const& model);
Real rad_soft_0(Real lambda_e, math::Vec3 eta, kin::KinematicsRad kin, sf::Model const& model);

Real ex(Real lambda_e, math::Vec3 eta, kin::KinematicsEx kin, sf::Model const& model);

Real born_xu(Real lambda_e, kin::Kinematics kin, sf::Model const& model);
Real born_xl(Real lambda_e, kin::Kinematics kin, sf::Model const& model);
Real born_xt(Real lambda_e, Real phi_s, kin::Kinematics kin, sf::Model const& model);
Real born_xt1(Real lambda_e, kin::Kinematics kin, sf::Model const& model);
Real born_xt2(Real lambda_e, kin::Kinematics kin, sf::Model const& model);
Real born_lx(math::Vec3 eta, kin::Kinematics kin, sf::Model const& model);
Real born_ux(math::Vec3 eta, kin::Kinematics kin, sf::Model const& model);

Real born_uu(kin::Kinematics kin, sf::Model const& model);
Real born_ul(kin::Kinematics kin, sf::Model const& model);
Real born_ut(Real phi_s, kin::Kinematics kin, sf::Model const& model);
Real born_ut1(kin::Kinematics kin, sf::Model const& model);
Real born_ut2(kin::Kinematics kin, sf::Model const& model);
Real born_lu(kin::Kinematics kin, sf::Model const& model);
Real born_ll(kin::Kinematics kin, sf::Model const& model);
Real born_lt(Real phi_s, kin::Kinematics kin, sf::Model const& model);
Real born_lt1(kin::Kinematics kin, sf::Model const& model);
Real born_lt2(kin::Kinematics kin, sf::Model const& model);

Real amm_uu(kin::Kinematics kin, sf::Model const& model);
Real amm_ul(kin::Kinematics kin, sf::Model const& model);
Real amm_ut(Real phi_s, kin::Kinematics kin, sf::Model const& model);
Real amm_ut1(kin::Kinematics kin, sf::Model const& model);
Real amm_ut2(kin::Kinematics kin, sf::Model const& model);
Real amm_lu(kin::Kinematics kin, sf::Model const& model);
Real amm_ll(kin::Kinematics kin, sf::Model const& model);
Real amm_lt(Real phi_s, kin::Kinematics kin, sf::Model const& model);
Real amm_lt1(kin::Kinematics kin, sf::Model const& model);
Real amm_lt2(kin::Kinematics kin, sf::Model const& model);

// Base sets for cross-section calculations.
Real born_uu(Born b, lep::LepBornUU lep, had::HadUU had);
Real born_ul(Born b, lep::LepBornUP lep, had::HadUL had);
Real born_ut1(Born b, lep::LepBornUP lep, had::HadUT1 had);
Real born_ut2(Born b, lep::LepBornUU lep, had::HadUT2 had);
Real born_lu(Born b, lep::LepBornLU lep, had::HadLU had);
Real born_ll(Born b, lep::LepBornLP lep, had::HadLL had);
Real born_lt1(Born b, lep::LepBornLP lep, had::HadLT1 had);
Real born_lt2(Born b, lep::LepBornLU lep, had::HadLT2 had);

Real amm_uu(Amm b, lep::LepAmmUU lep, had::HadUU had);
Real amm_ul(Amm b, lep::LepAmmUP lep, had::HadUL had);
Real amm_ut1(Amm b, lep::LepAmmUP lep, had::HadUT1 had);
Real amm_ut2(Amm b, lep::LepAmmUU lep, had::HadUT2 had);
Real amm_lu(Amm b, lep::LepAmmLU lep, had::HadLU had);
Real amm_ll(Amm b, lep::LepAmmLP lep, had::HadLL had);
Real amm_lt1(Amm b, lep::LepAmmLP lep, had::HadLT1 had);
Real amm_lt2(Amm b, lep::LepAmmLU lep, had::HadLT2 had);

Real rad_hard_uu(Rad b, lep::LepRadUU lep, had::HadUU shift_had);
math::Vec3 rad_hard_up(Rad b, lep::LepRadUX lep, had::HadUP shift_had, math::Transform3 shift_rot);
Real rad_hard_lu(Rad b, lep::LepRadLU lep, had::HadLU shift_had);
math::Vec3 rad_hard_lp(Rad b, lep::LepRadLX lep, had::HadLP shift_had, math::Transform3 shift_rot);

Real rad_soft_uu_R(Rad b, lep::LepRadUU lep, had::HadUU shift_had);
math::Vec3 rad_soft_up_R(Rad b, lep::LepRadUX lep, had::HadUP shift_had, math::Transform3 shift_rot);
Real rad_soft_lu_R(Rad b, lep::LepRadLU lep, had::HadLU shift_had);
math::Vec3 rad_soft_lp_R(Rad b, lep::LepRadLX lep, had::HadLP shift_had, math::Transform3 shift_rot);

Real rad_soft_uu_R0(Rad b, lep::LepRadUU lep, had::HadUU had);
Real rad_soft_ul_R0(Rad b, lep::LepRadUP lep, had::HadUL had);
Real rad_soft_ut1_R0(Rad b, lep::LepRadUP lep, had::HadUT1 had);
Real rad_soft_ut2_R0(Rad b, lep::LepRadUU lep, had::HadUT2 had);
Real rad_soft_lu_R0(Rad b, lep::LepRadLU lep, had::HadLU had);
Real rad_soft_ll_R0(Rad b, lep::LepRadLP lep, had::HadLL had);
Real rad_soft_lt1_R0(Rad b, lep::LepRadLP lep, had::HadLT1 had);
Real rad_soft_lt2_R0(Rad b, lep::LepRadLU lep, had::HadLT2 had);

}
}

#endif

