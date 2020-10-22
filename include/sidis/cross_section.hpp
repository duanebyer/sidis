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

// Total SIDIS cross-section, including radiative corrections from all sources.
Real rc(Real lambda_e, math::Vec3 eta, kin::Kinematics kin, sf::Model const& model);

Real born(Real lambda_e, math::Vec3 eta, kin::Kinematics kin, sf::Model const& model);
Real amm(Real lambda_e, math::Vec3 eta, kin::Kinematics kin, sf::Model const& model);
Real nrad(Real lambda_e, math::Vec3 eta, kin::Kinematics kin, sf::Model const& model);
// TODO: Include a precision argument, which determines how accurately the
// integrated cross-section is to be calculated.
Real rad_integ(Real lambda_e, math::Vec3 eta, kin::Kinematics kin, sf::Model const& model);

Real rad(Real lambda_e, math::Vec3 eta, kin::KinematicsRad kin, sf::Model const& model);
Real rad_hard(Real lambda_e, math::Vec3 eta, kin::KinematicsRad kin, sf::Model const& model);
Real rad_soft(Real lambda_e, math::Vec3 eta, kin::KinematicsRad kin, sf::Model const& model);
Real rad_soft_0(Real lambda_e, math::Vec3 eta, kin::KinematicsRad kin, sf::Model const& model);

// Radiative corrections to Born cross-section.
Real born_rad_factor(kin::Kinematics kin);
Real delta_vr(kin::Kinematics kin);
Real delta_vac_lep(kin::Kinematics kin);
Real delta_vac_had(kin::Kinematics kin);

// Base functions, for calculating cross-sections directly from leptonic and
// hadronic coefficients.
struct Born {
	Real coeff;
	explicit Born(kin::Kinematics kin);
};

Real born_uu_base(Born b, lep::LepBornUU lep, had::HadUU had);
Real born_ul_base(Born b, lep::LepBornUP lep, had::HadUL had);
Real born_ut1_base(Born b, lep::LepBornUP lep, had::HadUT1 had);
Real born_ut2_base(Born b, lep::LepBornUU lep, had::HadUT2 had);
Real born_lu_base(Born b, lep::LepBornLU lep, had::HadLU had);
Real born_ll_base(Born b, lep::LepBornLP lep, had::HadLL had);
Real born_lt1_base(Born b, lep::LepBornLP lep, had::HadLT1 had);
Real born_lt2_base(Born b, lep::LepBornLU lep, had::HadLT2 had);

struct Amm {
	Real coeff;
	explicit Amm(kin::Kinematics kin);
};

Real amm_uu_base(Amm b, lep::LepAmmUU lep, had::HadUU had);
Real amm_ul_base(Amm b, lep::LepAmmUP lep, had::HadUL had);
Real amm_ut1_base(Amm b, lep::LepAmmUP lep, had::HadUT1 had);
Real amm_ut2_base(Amm b, lep::LepAmmUU lep, had::HadUT2 had);
Real amm_lu_base(Amm b, lep::LepAmmLU lep, had::HadLU had);
Real amm_ll_base(Amm b, lep::LepAmmLP lep, had::HadLL had);
Real amm_lt1_base(Amm b, lep::LepAmmLP lep, had::HadLT1 had);
Real amm_lt2_base(Amm b, lep::LepAmmLU lep, had::HadLT2 had);

struct Rad {
	Real coeff;
	Real R;
	explicit Rad(kin::KinematicsRad kin);
};

Real rad_hard_uu_base(Rad b, lep::LepRadUU lep, had::HadUU shift_had);
math::Vec3 rad_hard_up_base(Rad b, lep::LepRadUX lep, had::HadUP shift_had, math::Transform3 shift_rot);
Real rad_hard_lu_base(Rad b, lep::LepRadLU lep, had::HadLU shift_had);
math::Vec3 rad_hard_lp_base(Rad b, lep::LepRadLX lep, had::HadLP shift_had, math::Transform3 shift_rot);

Real rad_soft_uu_base(Rad b, lep::LepRadUU lep, had::HadUU had, had::HadUU shift_had);
math::Vec3 rad_soft_up_base(Rad b, lep::LepRadUX lep, had::HadUP had, had::HadUP shift_had, math::Transform3 shift_rot);
Real rad_soft_lu_base(Rad b, lep::LepRadLU lep, had::HadLU had, had::HadLU shift_had);
math::Vec3 rad_soft_lp_base(Rad b, lep::LepRadLX lep, had::HadLP had, had::HadLP shift_had, math::Transform3 shift_rot);

Real rad_soft_uu_base_R(Rad b, lep::LepRadUU lep, had::HadUU shift_had);
math::Vec3 rad_soft_up_base_R(Rad b, lep::LepRadUX lep, had::HadUP shift_had, math::Transform3 shift_rot);
Real rad_soft_lu_base_R(Rad b, lep::LepRadLU lep, had::HadLU shift_had);
math::Vec3 rad_soft_lp_base_R(Rad b, lep::LepRadLX lep, had::HadLP shift_had, math::Transform3 shift_rot);

Real rad_soft_uu_base_R0(Rad b, lep::LepRadUU lep, had::HadUU had);
Real rad_soft_ul_base_R0(Rad b, lep::LepRadUP lep, had::HadUL had);
Real rad_soft_ut1_base_R0(Rad b, lep::LepRadUP lep, had::HadUT1 had);
Real rad_soft_ut2_base_R0(Rad b, lep::LepRadUU lep, had::HadUT2 had);
Real rad_soft_lu_base_R0(Rad b, lep::LepRadLU lep, had::HadLU had);
Real rad_soft_ll_base_R0(Rad b, lep::LepRadLP lep, had::HadLL had);
Real rad_soft_lt1_base_R0(Rad b, lep::LepRadLP lep, had::HadLT1 had);
Real rad_soft_lt2_base_R0(Rad b, lep::LepRadLU lep, had::HadLT2 had);

}
}

#endif

