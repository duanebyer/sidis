#include "sidis/leptonic_coeff.hpp"

#include "sidis/kinematics.hpp"
#include "sidis/math.hpp"

using namespace sidis;
using namespace sidis::lep;
using namespace sidis::math;
using namespace sidis::kin;

// Born coefficients. Equation [1.16].
LepBornUU::LepBornUU(kin::Kinematics kin) {
	theta_1 = kin.Q_sq - 2.*sq(kin.m);
	theta_2 = 0.5*(kin.S*kin.X - sq(kin.M)*kin.Q_sq);
	theta_3 = 0.5*(kin.V_1*kin.V_2 - sq(kin.mh)*kin.Q_sq);
	theta_4 = 0.5*(kin.S*kin.V_2 + kin.X*kin.V_1 - kin.z*kin.Q_sq*kin.S_x);
}
LepBornUP::LepBornUP(kin::Kinematics kin) {
	theta_6 = -kin.S_p*kin.vol_phi_h;
	theta_8 = -2.*kin.V_p*kin.vol_phi_h;
}
LepBornLU::LepBornLU(kin::Kinematics kin) {
	theta_5 = (2.*kin.S*kin.vol_phi_h)/kin.lambda_S_sqrt;
}
LepBornLP::LepBornLP(kin::Kinematics kin) {
	theta_7 = kin.S/(4.*kin.lambda_S_sqrt)
		*(kin.lambda_Y*kin.V_p - kin.S_p*kin.S_x*(kin.z*kin.Q_sq + kin.V_m));
	theta_9 = 1./(2.*kin.lambda_S_sqrt)*(
		kin.S*(
			kin.Q_sq*(kin.z*kin.S_x*kin.V_p - sq(kin.mh)*kin.S_p)
			+ kin.V_m*(kin.S*kin.V_2 - kin.X*kin.V_1))
		+ 2.*sq(kin.m)*(
			4.*sq(kin.M)*sq(kin.V_m)
			+ kin.lambda_Y*sq(kin.mh)
			- kin.z*sq(kin.S_x)*(kin.z*kin.Q_sq + 2.*kin.V_m)
			));
}

// AMM coefficients. Equation [1.54].
LepAmmUU::LepAmmUU(kin::Kinematics kin) {
	theta_1 = 6.;
	theta_2 = -kin.lambda_Y/(2.*kin.Q_sq);
	theta_3 = -2.*sq(kin.mh) - 2.*sq(kin.V_m)/kin.Q_sq;
	theta_4 = -2.*kin.S_x*(kin.z + kin.V_m/kin.Q_sq);
}
LepAmmUP::LepAmmUP(kin::Kinematics kin) {
	static_cast<void>(kin);
	theta_6 = 0.;
	theta_8 = 0.;
}
LepAmmLU::LepAmmLU(kin::Kinematics kin) {
	theta_5 = (2.*(2.*kin.S + kin.S_x)*kin.vol_phi_h)
		/(kin.lambda_S_sqrt*kin.Q_sq);
}
LepAmmLP::LepAmmLP(kin::Kinematics kin) {
	theta_7 = (2.*kin.S + kin.S_x)/(4.*kin.lambda_S_sqrt*kin.Q_sq)*(
		kin.S_x*(kin.S*kin.V_2 - kin.X*kin.V_1 - kin.z*kin.S_p*kin.Q_sq)
		+ 4.*sq(kin.M)*kin.Q_sq*kin.V_p);
	theta_9 = 1./(2.*kin.lambda_S_sqrt*kin.Q_sq)*(
		sq(kin.S_x)*(
			4.*sq(kin.m)*(sq(kin.mh) - kin.z*(kin.z*kin.Q_sq + 2.*kin.V_m))
			+ kin.V_1*kin.V_m)
		- 4.*(sq(kin.M)*(kin.Q_sq - 4.*sq(kin.m)) + sq(kin.S))
			*(sq(kin.mh)*kin.Q_sq + sq(kin.V_m))
		+ kin.z*kin.Q_sq*kin.S_x*(
			kin.S_x*(kin.z*kin.Q_sq + kin.V_1 + kin.V_m)
			+ 2.*kin.S*kin.V_p)
		+ 2.*kin.S*kin.S_x*kin.V_m*kin.V_p);
}

