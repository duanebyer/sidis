#include "sidis/leptonic_coeff.hpp"

#include "sidis/kinematics.hpp"
#include "sidis/extra/math.hpp"

using namespace sidis;
using namespace sidis::lep;
using namespace sidis::math;
using namespace sidis::kin;

// Born coefficients. Equation [1.16].
LepBornBaseUU::LepBornBaseUU(Kinematics const& kin) {
	theta_1 = kin.Q_sq - 2.*sq(kin.m);
	theta_2 = 0.5*(kin.S*kin.X - sq(kin.M)*kin.Q_sq);
	theta_3 = 0.5*(kin.V_1*kin.V_2 - sq(kin.mh)*kin.Q_sq);
	theta_4 = 0.5*(kin.S*kin.V_2 + kin.X*kin.V_1 - kin.z*kin.Q_sq*kin.S_x);
}
LepBornBaseUP::LepBornBaseUP(Kinematics const& kin) {
	theta_6 = -kin.S_p*kin.vol_phi_h;
	theta_8 = -2.*kin.V_p*kin.vol_phi_h;
}
LepBornBaseLU::LepBornBaseLU(Kinematics const& kin) {
	theta_5 = (2.*kin.S*kin.vol_phi_h)/kin.lambda_S_sqrt;
}
LepBornBaseLP::LepBornBaseLP(Kinematics const& kin) {
	theta_7 = kin.S/(4.*kin.lambda_S_sqrt)*(
		kin.lambda_Y*kin.V_p
		- kin.S_p*kin.S_x*(kin.z*kin.Q_sq + kin.V_m));
	theta_9 = 1./(2.*kin.lambda_S_sqrt)*(
		kin.S*(
			kin.Q_sq*(kin.z*kin.S_x*kin.V_p - sq(kin.mh)*kin.S_p)
			+ kin.V_m*(kin.S*kin.V_2 - kin.X*kin.V_1))
		+ 2.*sq(kin.m)*(
			4.*sq(kin.M)*sq(kin.V_m)
			+ kin.lambda_Y*sq(kin.mh)
			- kin.z*sq(kin.S_x)*(kin.z*kin.Q_sq + 2.*kin.V_m)));
}

// AMM coefficients. Equation [1.54].
LepAmmBaseUU::LepAmmBaseUU(Kinematics const& kin) {
	theta_1 = 6.;
	theta_2 = -kin.lambda_Y/(2.*kin.Q_sq);
	theta_3 = -2.*sq(kin.mh) - 2.*sq(kin.V_m)/kin.Q_sq;
	theta_4 = -2.*kin.S_x*(kin.z + kin.V_m/kin.Q_sq);
}
LepAmmBaseUP::LepAmmBaseUP(Kinematics const& kin) {
	static_cast<void>(kin);
	theta_6 = 0.;
	theta_8 = 0.;
}
LepAmmBaseLU::LepAmmBaseLU(Kinematics const& kin) {
	theta_5 = (2.*(2.*kin.S + kin.S_x)*kin.vol_phi_h)
		/(kin.lambda_S_sqrt*kin.Q_sq);
}
LepAmmBaseLP::LepAmmBaseLP(Kinematics const& kin) {
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

// Non-radiative coefficients.
LepNradBaseUU::LepNradBaseUU(Kinematics const& kin) : born(kin), amm(kin) { }
LepNradBaseUP::LepNradBaseUP(Kinematics const& kin) : born(kin), amm(kin) { }
LepNradBaseLU::LepNradBaseLU(Kinematics const& kin) : born(kin), amm(kin) { }
LepNradBaseLP::LepNradBaseLP(Kinematics const& kin) : born(kin), amm(kin) { }
LepNradBaseUU::LepNradBaseUU(LepBornBaseUU const& born, LepAmmBaseUU const& amm) : born(born), amm(amm) { }
LepNradBaseUP::LepNradBaseUP(LepBornBaseUP const& born, LepAmmBaseUP const& amm) : born(born), amm(amm) { }
LepNradBaseLU::LepNradBaseLU(LepBornBaseLU const& born, LepAmmBaseLU const& amm) : born(born), amm(amm) { }
LepNradBaseLP::LepNradBaseLP(LepBornBaseLP const& born, LepAmmBaseLP const& amm) : born(born), amm(amm) { }

// Radiative coefficients. Equations [1.B1] and [1.B2].
LepRadBaseUU::LepRadBaseUU(KinematicsRad const& kin) {
	LepBornBaseUU lep_born(kin.project());
	theta_011 = 4.*kin.F_IR*lep_born.theta_1;
	theta_012 = 4.*kin.tau*kin.F_IR;
	theta_013 = -4. - 2.*kin.F_d*sq(kin.tau);
	theta_021 = 4.*kin.F_IR*lep_born.theta_2;
	theta_022 = 0.5*(
		kin.F_1p*kin.S_x*kin.S_p
		+ 2.*kin.F_2m*sq(kin.m)*kin.S_p
		+ 2.*kin.F_IR*(kin.S_x - 2.*kin.tau*sq(kin.M))
		- kin.F_d*kin.tau*sq(kin.S_p));
	theta_023 = 0.5*(
		kin.F_d*(4.*sq(kin.m) + kin.tau*(2.*kin.tau*sq(kin.M) - kin.S_x))
		- kin.F_1p*kin.S_p
		+ 4.*sq(kin.M));
	theta_031 = 4.*kin.F_IR*lep_born.theta_3;
	theta_032 = 2.*(
		kin.F_IR*(kin.mu*kin.V_m - kin.tau*sq(kin.mh))
		+ kin.V_p*(
			kin.F_2m*kin.mu*sq(kin.m)
			+ kin.F_1p*kin.V_m
			- kin.F_d*kin.tau*kin.V_p));
	theta_033 = kin.F_d*(
			2.*sq(kin.mu)*sq(kin.m)
			+ kin.tau*(kin.tau*sq(kin.mh) - kin.mu*kin.V_m))
		- kin.F_1p*kin.mu*kin.V_p + 2.*sq(kin.mh);
	theta_041 = 4.*kin.F_IR*lep_born.theta_4;
	theta_042 = kin.F_1p*(kin.S*kin.V_1 - kin.X*kin.V_2)
		+ kin.F_2m*sq(kin.m)*(kin.mu*kin.S_p + 2.*kin.V_p)
		- 2.*kin.F_d*kin.tau*kin.S_p*kin.V_p
		+ kin.F_IR*((kin.mu - 2.*kin.tau*kin.z)*kin.S_x + 2.*kin.V_m);
	theta_043 = 0.5*(
		kin.F_d*(
			8.*kin.mu*sq(kin.m)
			+ kin.tau*(
				(2.*kin.tau*kin.z - kin.mu)*kin.S_x
				- 2.*kin.V_m))
		- kin.F_1p*(kin.mu*kin.S_p + 2.*kin.V_p)
		+ 4.*kin.z*kin.S_x);
}
LepRadBaseUP::LepRadBaseUP(KinematicsRad const& kin) {
	LepBornBaseUP lep_born(kin.project());
	theta_061 = 4.*kin.F_IR*lep_born.theta_6;
	theta_062 = 1./(2.*kin.lambda_1)*(
		kin.vol_phi_h*(
			(
				4.*sq(kin.M)*kin.Q_sq*(kin.Q_sq + 4.*sq(kin.m))
				- sq(kin.S_x)*(kin.Q_sq - 4.*sq(kin.m))
				- 8.*kin.Q_sq*kin.S*kin.X)*(
					kin.F_1p*kin.S_x
					+ 2.*kin.F_2m*sq(kin.m)
					- kin.F_d*kin.tau*kin.S_p)
			+ 2.*kin.F_IR*kin.S_p*(
				2.*kin.tau*(
					2.*sq(kin.M)*(kin.Q_sq + 2.*sq(kin.m))
					- kin.S*kin.X)
				- kin.S_x*(kin.Q_sq + 4.*sq(kin.m))))
		+ 2.*kin.S_p*kin.vol_phi_k_R*(
			kin.F_2m*sq(kin.m)*(
				kin.S_x*(kin.z*kin.S_p*kin.Q_sq - kin.S*kin.V_2 + kin.V_1*kin.X)
				- 4.*sq(kin.M)*kin.Q_sq*kin.V_p)
			+ kin.F_IR*(
				(kin.Q_sq + 4.*sq(kin.m))
					*(kin.z*sq(kin.S_x) - 4.*sq(kin.M)*kin.V_m)
				+ kin.S_p*(kin.X*kin.V_1 - kin.S*kin.V_2))));
	theta_063 = 1./(2.*kin.lambda_1)*(
		2.*kin.vol_phi_h*(
			kin.F_1p*(
				2.*kin.Q_sq*(kin.S*kin.X - 2.*sq(kin.M)*kin.Q_sq)
				- kin.tau*kin.S_x*(
					sq(kin.S_x)
					+ 3.*kin.S*kin.X
					- 4.*sq(kin.m)*sq(kin.M))
				- (kin.Q_sq + 2.*sq(kin.m))*sq(kin.S_x))
			+ kin.F_2m*sq(kin.m)*(
				2.*kin.tau*(
					2.*sq(kin.M)*(kin.Q_sq + 2.*sq(kin.m))
					- kin.S*kin.X)
				- kin.S_x*(kin.Q_sq + 4.*sq(kin.m)))
			- kin.F_IR*kin.Q_sq*kin.S_p
			+ kin.F_d*kin.S_p*(
				sq(kin.tau)*(
					sq(kin.S_x)
					+ 2.*kin.S*kin.X
					- 2.*sq(kin.M)*(kin.Q_sq + 4.*sq(kin.m)))
				+ 2.*kin.tau*kin.S_x*(kin.Q_sq + 2.*sq(kin.m))
				- 4.*sq(kin.m)*kin.Q_sq)
			+ kin.S_p*sq(kin.S_x))
		+ kin.vol_phi_k_R*(
			(kin.S_x*kin.F_1p + 2.*sq(kin.m)*kin.F_2m)*(
				(kin.Q_sq + 4.*sq(kin.m))
					*(kin.z*sq(kin.S_x) - 4.*sq(kin.M)*kin.V_m)
				+ kin.S_p*(kin.X*kin.V_1 - kin.S*kin.V_2))
			+ 2.*kin.F_2p*sq(kin.m)*(
				kin.S_x*(kin.z*kin.Q_sq*kin.S_p - kin.S*kin.V_2 + kin.X*kin.V_1)
				- 4.*sq(kin.M)*kin.Q_sq*kin.V_p)
			+ kin.F_d*(
				4.*kin.tau*(
					sq(kin.M)*kin.Q_sq*(
						4.*kin.S*kin.V_m
						+ kin.S_x*(kin.V_2 - kin.V_m))
					+ 2.*kin.S*kin.X*(kin.S*kin.V_2 - kin.X*kin.V_1)
					+ 2.*sq(kin.m)*kin.S_p*(
						4.*sq(kin.M)*kin.V_m
						- kin.z*sq(kin.S_x)))
				+ (3.*kin.tau*kin.S_x + 2.*(kin.Q_sq - 2.*sq(kin.m)))*(
					kin.S*kin.V_2
					- kin.X*kin.V_1
					- kin.z*kin.Q_sq*kin.S_p)*kin.S_x
				+ 8.*(kin.Q_sq - 2.*sq(kin.m))*sq(kin.M)*kin.Q_sq*kin.V_p)));
	theta_064 = 1./(2.*kin.lambda_1)*(
		kin.vol_phi_h*(
			kin.F_1p*(
				(kin.Q_sq + 4.*sq(kin.m))*kin.S_x
				+ 2.*kin.tau*(
					kin.S*kin.X
					- 2.*sq(kin.M)*(kin.Q_sq + 2.*sq(kin.m))))
			+ kin.S_p*(kin.tau*kin.Q_sq*kin.F_d - 2.*kin.S_x))
		+ kin.vol_phi_k_R*(
			kin.F_1p*(
				(kin.Q_sq + 4.*sq(kin.m))
					*(4.*sq(kin.M)*kin.V_m - kin.z*sq(kin.S_x))
				+ kin.S_p*(kin.S*kin.V_2 - kin.X*kin.V_1))
			+ kin.F_d*kin.tau*(
				4.*sq(kin.M)*kin.Q_sq*kin.V_p
				+ kin.S_x*(
					kin.S*kin.V_2
					- kin.X*kin.V_1
					- kin.z*kin.Q_sq*kin.S_p))));
	theta_081 = 4.*kin.F_IR*lep_born.theta_8;
	theta_082 = 1./kin.lambda_1*(
		kin.vol_phi_h*(
			kin.F_1p*(
				kin.Q_sq*kin.S_x*(kin.S_x*kin.V_p - 2.*kin.S*kin.V_2)
				- 2.*kin.V_m*(
					2.*kin.lambda_1
					+ kin.Q_sq*kin.S*kin.S_x))
			- 2.*kin.F_2m*sq(kin.m)*(
				2.*kin.mu*kin.lambda_1
				+ kin.Q_sq*kin.S_p*kin.V_p)
			+ kin.V_p*(
				2.*kin.F_2p*sq(kin.m)*(
					2.*kin.tau*(
						2.*(kin.Q_sq + 2.*sq(kin.m))*sq(kin.M)
						- kin.S*kin.X)
					- (kin.Q_sq + 4.*sq(kin.m))*kin.S_x)
				+ kin.F_d*(
					4.*sq(kin.m)*(
						(3.*kin.Q_sq + 4.*sq(kin.m))*kin.S_x
						+ kin.tau*(
							2.*kin.S*kin.X
							- 4.*(3.*kin.Q_sq + 2.*sq(kin.m))*sq(kin.M)
							- sq(kin.S_x)))
					+ kin.Q_sq*(
						kin.tau*(12.*kin.S*kin.X + sq(kin.S_x))
						+ 2.*kin.Q_sq*(kin.S_x - 6.*kin.tau*sq(kin.M))))))
		+ 2.*kin.V_p*kin.vol_phi_k_R*(
			kin.F_IR*(
				(kin.Q_sq + 4.*sq(kin.m))*(
					kin.z*sq(kin.S_x)
					- 4.*sq(kin.M)*kin.V_m)
				+ kin.S_p*(kin.X*kin.V_1 - kin.S*kin.V_2))
			+ kin.F_2m*sq(kin.m)*(
				kin.S_x*(kin.X*kin.V_1 - kin.S*kin.V_2 + kin.z*kin.S_p*kin.Q_sq)
				- 4.*kin.Q_sq*kin.V_p*sq(kin.M))));
	theta_083 = 1./(2.*kin.lambda_1)*(
		kin.vol_phi_h*(
			kin.F_d*(
				2.*kin.mu*(kin.Q_sq - 2.*sq(kin.m))*kin.Q_sq*kin.S_p
				+ kin.tau*(
					2.*(kin.Q_sq + 8.*sq(kin.m))*kin.S_x*kin.V_p
					+ kin.Q_sq*(
						kin.mu*kin.S_x*kin.S_p
						+ 2.*kin.S*kin.V_1
						- 2.*kin.X*kin.V_2))
				- 2.*sq(kin.tau)*(
					4.*(kin.Q_sq + 4.*sq(kin.m))*kin.V_p*sq(kin.M)
					- kin.S_p*(kin.S*kin.V_1 + kin.X*kin.V_2)))
			+ 2.*kin.F_2m*kin.mu*sq(kin.m)*(
				2.*kin.tau*(
					2.*(kin.Q_sq + 2.*sq(kin.m))*sq(kin.M)
					- kin.S*kin.X)
				- (kin.Q_sq + 4.*sq(kin.m))*kin.S_x)
			+ 4.*kin.S_p*kin.S_x*kin.V_m
			- 2.*kin.F_2p*kin.mu*sq(kin.m)*kin.Q_sq*kin.S_p
			+ kin.F_1p*(
				2.*kin.tau*(
					sq(kin.X)*kin.V_2
					- sq(kin.S)*kin.V_1)
				+ 8.*sq(kin.m)*kin.V_m*(2.*kin.tau*sq(kin.M) - kin.S_x)
				+ kin.Q_sq*(
					4.*kin.mu*(kin.S*kin.X - 2.*kin.Q_sq*sq(kin.M))
					- kin.S_x*(2.*kin.V_m + kin.mu*kin.S_x))))
		+ 2.*kin.vol_phi_k_R*(
			(kin.F_1p*kin.V_m + kin.F_2m*kin.mu*sq(kin.m))*(
				(kin.Q_sq + 4.*sq(kin.m))*(
					kin.z*sq(kin.S_x)
					- 4.*sq(kin.M)*kin.V_m)
				+ kin.S_p*(kin.X*kin.V_1 - kin.S*kin.V_2))
			+ kin.F_2p*kin.mu*sq(kin.m)*(
				kin.S_x*(kin.z*kin.S_p*kin.Q_sq + kin.X*kin.V_1 - kin.S*kin.V_2)
				- 4.*kin.Q_sq*kin.V_p*sq(kin.M))
			+ kin.F_d*(
				kin.mu*(kin.Q_sq - 2.*sq(kin.m))*(
					4.*kin.Q_sq*kin.V_p*sq(kin.M)
					+ kin.S_x*(
						kin.S*kin.V_2
						- kin.X*kin.V_1
						- kin.z*kin.S_p*kin.Q_sq))
				+ kin.tau*(
					kin.S_x*kin.S_p*kin.V_2*(kin.V_1 + kin.V_p)
					+ 2.*kin.V_m*kin.V_p*(
						(kin.S_x - 4.*kin.S)*kin.X
						+ 2.*(3.*kin.Q_sq + 8.*sq(kin.m))*sq(kin.M))
					+ kin.z*kin.S_x*(
						kin.Q_sq*(kin.X*kin.V_2 - kin.S*kin.V_1)
						- (kin.Q_sq + 8.*sq(kin.m))*kin.S_x*kin.V_p)))));
	theta_084 = kin.mu/(2.*kin.lambda_1)*(
		kin.vol_phi_h*(
			kin.F_1p*(
				(kin.Q_sq + 4.*sq(kin.m))*kin.S_x
				+ 2.*kin.tau*(
					kin.S*kin.X
					- 2.*(kin.Q_sq + 2.*sq(kin.m))*sq(kin.M)))
			+ kin.S_p*(kin.F_d*kin.tau*kin.Q_sq - 2.*kin.S_x))
		+ kin.vol_phi_k_R*(
			kin.F_1p*(
				(kin.Q_sq + 4.*sq(kin.m))*(
					4.*sq(kin.M)*kin.V_m
					- kin.z*sq(kin.S_x))
				+ kin.S_p*(kin.S*kin.V_2 - kin.X*kin.V_1))
			+ kin.F_d*kin.tau*(
				4.*kin.Q_sq*kin.V_p*sq(kin.M)
				+ kin.S_x*(
					kin.S*kin.V_2
					- kin.X*kin.V_1
					- kin.z*kin.S_p*kin.Q_sq))));
}
LepRadBaseLU::LepRadBaseLU(KinematicsRad const& kin) {
	LepBornBaseLU lep_born(kin.project());
	theta_051 = 4.*kin.F_IR*lep_born.theta_5;
	theta_052 = kin.S/(kin.lambda_1*kin.lambda_S_sqrt)*(
		kin.vol_phi_h*(
			2.*kin.F_IR*(
				kin.S_x*(kin.Q_sq + 4.*sq(kin.m))
				+ 2.*kin.tau*(
					kin.S*kin.X
					- 2.*sq(kin.M)*(kin.Q_sq + 2.*sq(kin.m))))
			+ kin.Q_sq*(
				kin.S_p*(kin.F_1p*kin.S_x + 2.*kin.F_2m*sq(kin.m))
				- kin.F_d*kin.tau*(4.*kin.S*kin.X + sq(kin.S_x))))
		+ 2.*kin.vol_phi_k_R*(
			kin.F_2m*sq(kin.m)*(
				kin.S_x*(
					kin.S*kin.V_2
					- kin.X*kin.V_1
					- kin.z*kin.Q_sq*kin.S_p)
				+ 4.*sq(kin.M)*kin.Q_sq*kin.V_p)
			+ kin.F_IR*(
				(kin.Q_sq + 4.*sq(kin.m))
					*(4.*sq(kin.M)*kin.V_m - kin.z*sq(kin.S_x))
				+ kin.S_p*(kin.S*kin.V_2 - kin.X*kin.V_1))));
	Real theta_053_hat = (2.*kin.S)/(kin.lambda_1*kin.lambda_S_sqrt)*kin.F_21*(
		kin.vol_phi_k_R*(
			2.*(kin.mu*kin.Q_sq + kin.tau*kin.V_1)
				*(kin.S*kin.X - sq(kin.M)*kin.Q_sq)
			+ (kin.Q_sq + kin.tau*kin.S)
				*(kin.z*kin.Q_sq*kin.S_x - kin.S*kin.V_2 - kin.X*kin.V_1))
		- kin.vol_phi_h*sq(kin.Q_sq + kin.tau*kin.S));
	theta_053 = theta_053_hat + kin.S/(kin.lambda_1*kin.lambda_S_sqrt)*(
		kin.vol_phi_h*(
			8.*kin.F_21*sq(kin.m)*(
				kin.tau*(kin.tau*sq(kin.M) - kin.S_x)
				- kin.Q_sq)
			+ kin.F_1p*(
				kin.Q_sq*(4.*kin.tau*sq(kin.M) + kin.S_p)
				+ 2.*kin.tau*kin.S*kin.S_x)
			+ kin.F_d*kin.tau*(
				4.*sq(kin.m)*(2.*kin.tau*sq(kin.M) - kin.S_x)
				+ kin.Q_sq*(kin.S_x - 4.*kin.S)
				- 2.*kin.tau*sq(kin.S)))
		+ 2.*kin.vol_phi_k_R*(
			2.*kin.F_21*sq(kin.m)*(
				kin.S_x*(
					2.*kin.z*kin.Q_sq
					+ 2.*kin.V_m
					+ (kin.tau*kin.z - kin.mu)*kin.S_x)
				- 4.*sq(kin.M)*(kin.mu*kin.Q_sq + kin.tau*kin.V_m))
			+ kin.F_d*kin.tau*(
				2.*sq(kin.m)*(kin.z*sq(kin.S_x) - 4.*sq(kin.M)*kin.V_m)
				- 2.*sq(kin.M)*kin.Q_sq*kin.V_1
				+ kin.S*(
					kin.z*kin.S_x*kin.Q_sq
					- kin.S*kin.V_2
					+ kin.X*kin.V_1))));
	theta_151 = 0.;
	theta_152 = (2.*sq(kin.m))/(kin.lambda_1*kin.lambda_S_sqrt)*(
		kin.vol_phi_h*(
			2.*kin.F_21*(
				2.*sq(kin.m)*kin.lambda_Y
				+ (kin.Q_sq + kin.tau*kin.S)*(
					2.*sq(kin.M)*kin.Q_sq
					+ kin.S*kin.S_x))
			- kin.F_1p*kin.S_x*kin.lambda_Y
			+ kin.F_d*(
				2.*kin.Q_sq*kin.X*kin.S_x
				+ kin.tau*kin.S_x*(2.*sq(kin.S) - sq(kin.S_p))
				+ 4.*sq(kin.M)*kin.Q_sq*(kin.tau*kin.S - kin.Q_sq)
				- 4.*sq(kin.m)*kin.lambda_Y))
		+ 2.*kin.vol_phi_k_R*(kin.F_d*kin.X - kin.F_21*kin.S)*(
			kin.S_x*(
				kin.z*kin.S_p*kin.Q_sq
				+ kin.X*kin.V_1
				- kin.S*kin.V_2)
			- 4.*sq(kin.M)*kin.Q_sq*kin.V_p));
	theta_153 = (2.*sq(kin.m))/(kin.lambda_1*kin.lambda_S_sqrt)*(
		kin.vol_phi_h*(
			2.*kin.F_21*(
				(kin.Q_sq + 2.*sq(kin.m))*(2.*kin.tau*sq(kin.M) + kin.X)
				- (kin.tau*kin.X + 2.*sq(kin.m))*kin.S)
			- kin.F_1p*kin.lambda_Y
			+ kin.F_d*(
				4.*sq(kin.m)*(kin.S_x - 2.*kin.tau*sq(kin.M))
				+ 2.*kin.Q_sq*kin.S + kin.tau*(sq(kin.S) + sq(kin.X))))
		+ 2.*kin.vol_phi_k_R*(
			kin.F_21*(
				2.*sq(kin.m)*(kin.z*sq(kin.S_x) - 4.*sq(kin.M)*kin.V_m)
				+ 2.*sq(kin.M)*kin.Q_sq*kin.V_2
				+ kin.X*(
					kin.X*kin.V_1
					- kin.S*kin.V_2
					- kin.z*kin.S_x*kin.Q_sq))
			+ kin.F_d*(2.*sq(kin.m)*(4.*sq(kin.M)*kin.V_m - kin.z*sq(kin.S_x))
				+ 2.*sq(kin.M)*kin.Q_sq*kin.V_1
				+ kin.S*(
					kin.S*kin.V_2
					- kin.X*kin.V_1
					- kin.z*kin.Q_sq*kin.S_x))));
}
LepRadBaseLP::LepRadBaseLP(KinematicsRad const& kin) {
	LepBornBaseLP lep_born(kin.project());
	theta_071 = 4.*kin.F_IR*lep_born.theta_7;
	theta_072 = kin.S/(2.*kin.lambda_S_sqrt)*(
		kin.F_1p*kin.Q_sq*(4.*sq(kin.M)*kin.V_m - kin.z*sq(kin.S_x))
		+ kin.F_2m*sq(kin.m)*(
			kin.mu*kin.lambda_Y
			- 2.*kin.S_x*(kin.z*kin.Q_sq + kin.V_m))
		+ kin.F_IR*(
			2.*(4.*kin.tau*sq(kin.M) - kin.S_x)*kin.V_p
			+ (kin.mu - 2.*kin.tau*kin.z)*kin.S_p*kin.S_x
			- 2.*kin.S*kin.V_2
			+ 2.*kin.X*kin.V_1)
		+ kin.F_d*kin.tau*(
			kin.Q_sq*(kin.z*kin.S_x*kin.S_p - 4.*sq(kin.M)*kin.V_p)
			+ kin.S_x*(kin.X*kin.V_1 - kin.S*kin.V_2)));
	theta_073 = kin.S/(4.*kin.lambda_S_sqrt)*(
		kin.F_1p*(
			kin.S_x*(4.*kin.z*kin.Q_sq + 2.*kin.V_m - kin.mu*kin.S_x)
			- 8.*kin.mu*sq(kin.M)*kin.Q_sq)
		+ 2.*kin.F_2m*sq(kin.m)*(
			4.*kin.mu*kin.tau*sq(kin.M)
			+ 2.*kin.V_m
			- (kin.mu + 2.*kin.tau*kin.z)*kin.S_x)
		+ 2.*kin.F_IR*(2.*kin.V_p - kin.mu*kin.S_p)
		+ kin.F_d*kin.tau*(
			4.*(kin.S_x - 2.*kin.tau*sq(kin.M))*kin.V_p
			+ kin.S_p*((2.*kin.tau*kin.z - kin.mu)*kin.S_x - 2.*kin.V_m)));
	theta_074 = kin.S/(4.*kin.lambda_S_sqrt)*(
		kin.F_1p*(
			(kin.mu + 2.*kin.tau*kin.z)*kin.S_x
			- 2.*kin.V_m
			- 4.*kin.mu*kin.tau*sq(kin.M))
		+ kin.F_d*kin.tau*(kin.mu*kin.S_p - 2.*kin.V_p));
	theta_171 = 0.;
	theta_172 = sq(kin.m)/kin.lambda_S_sqrt*(
		kin.F_21*(
			4.*sq(kin.M)*(kin.tau*kin.S*kin.V_m - kin.Q_sq*kin.V_p)
			- sq(kin.S_x)*(kin.tau*kin.z*kin.S + kin.z*kin.Q_sq + kin.V_1)
			+ kin.mu*kin.lambda_Y*kin.S)
		+ kin.F_d*(
			4.*sq(kin.M)*(kin.Q_sq*kin.V_p - kin.tau*kin.X*kin.V_m)
			+ sq(kin.S_x)*(kin.tau*kin.z*kin.X + kin.V_2 - kin.z*kin.Q_sq)
			- kin.mu*kin.lambda_Y*kin.X));
	theta_173 = sq(kin.m)/kin.lambda_S_sqrt*(
		kin.F_21*(
			2.*sq(kin.M)*(
				kin.mu*(kin.Q_sq + kin.tau*kin.S)
				- 2.*kin.tau*kin.V_p)
			+ kin.S_x*(
				(kin.tau*kin.z - 2.*kin.mu)*kin.S
				+ 2.*kin.V_p
				- kin.z*kin.Q_sq)
			+ (kin.mu - kin.tau*kin.z)*sq(kin.S_x))
		+ kin.F_d*(
			2.*sq(kin.M)*(
				kin.mu*(kin.Q_sq - kin.tau*kin.X)
				+ 2.*kin.tau*kin.V_p)
			+ kin.S_x*(
				(2.*kin.mu - kin.tau*kin.z)*kin.S
				- kin.z*kin.Q_sq
				- 2.*kin.V_p)
			- kin.mu*sq(kin.S_x)));
	theta_174 = sq(kin.m)/kin.lambda_S_sqrt*(
		kin.F_21*(
			2.*kin.mu*kin.tau*sq(kin.M)
			+ kin.mu*kin.X
			- kin.tau*kin.z*kin.S_x
			- kin.V_2)
		+ kin.F_d*(
			2.*kin.mu*kin.tau*sq(kin.M)
			+ kin.V_1
			- kin.mu*kin.S
			- kin.tau*kin.z*kin.S_x));
	theta_091 = (2.*kin.S)/kin.lambda_S_sqrt*kin.F_IR*(
		kin.Q_sq*(kin.z*kin.S_x*kin.V_p - sq(kin.mh)*kin.S_p)
		+ kin.V_m*(kin.S*kin.V_2 - kin.X*kin.V_1));
	theta_092 = kin.S/kin.lambda_S_sqrt*(
		kin.F_1p*kin.Q_sq*kin.S_x*(kin.z*kin.V_m - sq(kin.mh))
		+ kin.F_2m*sq(kin.m)*(
			kin.Q_sq*(kin.mu*kin.z*kin.S_x - 2.*sq(kin.mh))
			+ kin.V_m*(kin.mu*kin.S_x - 2.*kin.V_m))
		+ kin.F_d*kin.tau*(
			kin.Q_sq*(sq(kin.mh)*kin.S_p - kin.z*kin.S_x*kin.V_p)
			+ kin.V_m*(kin.X*kin.V_1 - kin.S*kin.V_2))
		+ kin.F_IR*(
			2.*kin.V_m*(2.*kin.mu*kin.S - kin.V_p)
			+ 2.*kin.tau*(kin.z*kin.S_x*kin.V_p - sq(kin.mh)*kin.S_p)
			- kin.mu*(kin.V_1 + kin.V_m)*kin.S_x));
	theta_093 = kin.S/(2.*kin.lambda_S_sqrt)*(
		kin.F_1p*(
			2.*(2.*sq(kin.mh)*kin.Q_sq + sq(kin.V_m))
			- kin.mu*kin.S_x*(2.*kin.z*kin.Q_sq + kin.V_m))
		+ kin.F_2m*sq(kin.m)*(
			kin.mu*((2.*kin.tau*kin.z - kin.mu)*kin.S_x + 2.*kin.V_m)
			- 4.*kin.tau*sq(kin.mh))
		+ kin.F_IR*kin.mu*(2.*kin.V_p - kin.mu*kin.S_p)
		+ kin.F_d*kin.tau*(
			2.*kin.tau*(sq(kin.mh)*kin.S_p - kin.z*kin.S_x*kin.V_p)
			+ kin.mu*kin.S_x*(kin.V_m + kin.V_1)
			+ 2.*kin.V_m*(kin.V_p - 2.*kin.mu*kin.S)));
	theta_094 = kin.S/(4.*kin.lambda_S_sqrt)*(
		kin.F_1p*(
			2.*kin.tau*(2.*sq(kin.mh) - kin.mu*kin.z*kin.S_x)
			+ kin.mu*(kin.mu*kin.S_x - 2.*kin.V_m))
		+ kin.F_d*kin.mu*kin.tau*(kin.mu*kin.S_p - 2.*kin.V_p));
	theta_191 = (4.*sq(kin.m))/kin.lambda_S_sqrt*kin.F_IR*(
		sq(kin.mh)*kin.lambda_Y
		+ 4.*sq(kin.M)*sq(kin.V_m)
		- kin.z*sq(kin.S_x)*(kin.z*kin.Q_sq + 2.*kin.V_m));
	theta_192 = (2.*sq(kin.m))/kin.lambda_S_sqrt*(
		2.*kin.F_2p*sq(kin.m)*(
			2.*sq(kin.mh)*(2.*kin.tau*sq(kin.M) - kin.S_x)
			+ 2.*(kin.z*kin.S_x - 2.*kin.mu*sq(kin.M))*kin.V_m
			+ kin.z*(kin.mu - kin.tau*kin.z)*sq(kin.S_x))
		+ kin.F_21*kin.S_x*(
			kin.z*kin.Q_sq*(kin.mu*kin.S - kin.V_p)
			+ kin.V_m*((kin.mu + kin.tau*kin.z)*kin.S - kin.V_1)
			- sq(kin.mh)*(kin.Q_sq + kin.tau*kin.S))
		+ kin.F_d*(
			kin.S_x*(
				(sq(kin.mh) - kin.z*kin.V_m)*(
					kin.tau*kin.X
					+ 3.*kin.Q_sq
					+ 8.*sq(kin.m))
				+ (kin.V_m + kin.z*kin.Q_sq)*(kin.V_2 - kin.mu*kin.X))
			+ 2.*(kin.Q_sq + 2.*sq(kin.m))*(
				4.*sq(kin.M)*(kin.mu*kin.V_m - kin.tau*sq(kin.mh))
				+ (kin.z*kin.tau - kin.mu)*kin.z*sq(kin.S_x))));
	theta_193 = sq(kin.m)/kin.lambda_S_sqrt*(
		4.*kin.F_2p*sq(kin.m)*(
			sq(kin.mh)
			+ sq(kin.mu)*sq(kin.M)
			- kin.mu*kin.z*kin.S_x)
		+ kin.F_21*(
			2.*sq(kin.mh)*(kin.tau*kin.X - kin.Q_sq)
			+ kin.S_x*(
				kin.mu*(kin.tau*kin.z - kin.mu)*kin.S
				- 2.*kin.tau*kin.z*kin.V_p
				+ kin.mu*kin.z*kin.Q_sq
				+ kin.mu*kin.V_1)
			+ 2.*kin.V_m*(kin.V_2 - kin.mu*kin.X))
		+ kin.F_d*(
			kin.mu*(
				(kin.Q_sq + 2.*sq(kin.m))*(
					5.*kin.z*kin.S_x
					- 4.*kin.mu*sq(kin.M))
				+ (kin.mu - kin.tau*kin.z)*kin.X*kin.S_x)
			- 2.*sq(kin.mh)*(kin.tau*kin.S + 3.*kin.Q_sq + 4.*sq(kin.m))
			+ kin.S_x*(
				2.*kin.tau*kin.z*kin.V_p
				- kin.mu*kin.V_2
				- 2.*kin.mu*kin.z*sq(kin.m))
			+ 2.*kin.V_m*(kin.mu*kin.S - kin.V_1)));
	theta_194 = sq(kin.m)/kin.lambda_S_sqrt*(
		kin.F_21*(
			kin.mu*(kin.tau*kin.z*kin.S_x + kin.mu*kin.X - kin.V_2)
			- 2.*kin.tau*sq(kin.mh))
		+ kin.F_d*(
			kin.mu*(kin.tau*kin.z*kin.S_x + kin.V_1 - kin.mu*kin.S)
			- 2.*kin.tau*sq(kin.mh)));
}

