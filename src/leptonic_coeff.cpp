#include "sidis/leptonic_coeff.hpp"

#include "sidis/math.hpp"

using namespace sidis;
using namespace sidis::lep;
using namespace sidis::math;
using namespace sidis::kin;

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

