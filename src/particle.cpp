#include "sidis/particle.hpp"

#include "sidis/extra/exception.hpp"

using namespace sidis;
using namespace sidis::part;

Particles::Particles(
		part::Nucleus target,
		part::Lepton beam,
		part::Hadron hadron,
		Real Mth) :
		target(target),
		beam(beam),
		hadron(hadron),
		M(mass(target)),
		m(mass(beam)),
		mh(mass(hadron)),
		Mth(Mth) {
	if (!(Mth >= M) || !(Mth <= M + mh)) {
		throw MassThresholdOutOfRange(Mth, M, M + mh);
	}
}

