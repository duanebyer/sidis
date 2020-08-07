#ifndef SIDIS_CROSS_SECTION_HPP
#define SIDIS_CROSS_SECTION_HPP

#include "sidis/numeric.hpp"

namespace sidis {

namespace math {
	struct Vec3;
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
}
namespace had {
	struct HadUU;
	struct HadUL;
	struct HadUT1;
	struct HadUT2;
	struct HadLU;
	struct HadLL;
	struct HadLT1;
	struct HadLT2;
}
namespace sf {
	struct SfUU;
	struct SfUL;
	struct SfUT;
	struct SfLU;
	struct SfLL;
	struct SfLT;
	struct SfXU;
	struct SfXL;
	struct SfXT;
	struct SfUX;
	struct SfLX;
	struct Sf;
}

namespace xs {

struct Born {
	Real coeff;
	Real c_1;

	explicit Born(kin::Kinematics kin);
};

struct Amm {
	Real coeff;
	Real c_1;

	explicit Amm(kin::Kinematics kin);
};

Real born(Real lambda_e, math::Vec3 eta, kin::Kinematics kin, sf::Sf sf);
Real delta_vr(kin::Kinematics kin);
Real delta_vac_lep(kin::Kinematics kin);
Real delta_vac_had(kin::Kinematics kin);

Real amm(Real lambda_e, math::Vec3 eta, kin::Kinematics kin, sf::Sf sf);
Real rad(Real lambda_e, math::Vec3 eta, kin::KinematicsRad kin, sf::Sf sf);
Real ex(Real lambda_e, math::Vec3 eta, kin::KinematicsEx kin, sf::Sf sf);

Real born_xu(Real lambda_e, kin::Kinematics kin, sf::SfXU sf);
Real born_xl(Real lambda_e, kin::Kinematics kin, sf::SfXL sf);
Real born_xt(Real lambda_e, Real phi_s, kin::Kinematics kin, sf::SfXT sf);
Real born_xt1(Real lambda_e, kin::Kinematics kin, sf::SfXT sf);
Real born_xt2(Real lambda_e, kin::Kinematics kin, sf::SfXT sf);
Real born_lx(math::Vec3 eta, kin::Kinematics kin, sf::SfLX sf);
Real born_ux(math::Vec3 eta, kin::Kinematics kin, sf::SfUX sf);

Real born_uu(kin::Kinematics kin, sf::SfUU sf);
Real born_ul(kin::Kinematics kin, sf::SfUL sf);
Real born_ut(Real phi_s, kin::Kinematics kin, sf::SfUT sf);
Real born_ut1(kin::Kinematics kin, sf::SfUT sf);
Real born_ut2(kin::Kinematics kin, sf::SfUT sf);
Real born_lu(kin::Kinematics kin, sf::SfLU sf);
Real born_ll(kin::Kinematics kin, sf::SfLL sf);
Real born_lt(Real phi_s, kin::Kinematics kin, sf::SfLT sf);
Real born_lt1(kin::Kinematics kin, sf::SfLT sf);
Real born_lt2(kin::Kinematics kin, sf::SfLT sf);

Real amm_uu(kin::Kinematics kin, sf::SfUU sf);
Real amm_ul(kin::Kinematics kin, sf::SfUL sf);
Real amm_ut(Real phi_s, kin::Kinematics kin, sf::SfUT sf);
Real amm_ut1(kin::Kinematics kin, sf::SfUT sf);
Real amm_ut2(kin::Kinematics kin, sf::SfUT sf);
Real amm_lu(kin::Kinematics kin, sf::SfLU sf);
Real amm_ll(kin::Kinematics kin, sf::SfLL sf);
Real amm_lt(Real phi_s, kin::Kinematics kin, sf::SfLT sf);
Real amm_lt1(kin::Kinematics kin, sf::SfLT sf);
Real amm_lt2(kin::Kinematics kin, sf::SfLT sf);

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

}
}

#endif

