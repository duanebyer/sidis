#ifndef SIDISGEN_GENERATOR_HPP
#define SIDISGEN_GENERATOR_HPP

#include <array>
#include <random>

#include <bubble.hpp>

#include <sidis/sidis.hpp>

#include "utility.hpp"

template<std::size_t D>
using Point = std::array<Double, D>;

class Params;

namespace sidis {
	namespace kin {
		struct Kinematics;
		struct KinematicsRad;
	}
}

// Map the non-radiative cross-section onto unit hypercube for Monte-Carlo.
class NradDensity final {
	sidis::cut::Cut _cut;
	sidis::sf::SfSet const& _sf;
	RcMethod _rc_method;
	Double _soft_threshold;
	sidis::part::Particles _ps;
	Double _S;
	Double _beam_pol;
	sidis::math::Vec3 _target_pol;

public:
	NradDensity(Params& params, sidis::sf::SfSet const& sf);
	// Get the density in the unit hypercube.
	Double eval(Point<6> const& unit_vec, sidis::kin::Kinematics* kin) const noexcept;
	Double eval(Point<6> const& unit_vec) const noexcept;
	// Transform from the unit hypercube into phase space.
	Double transform(Point<6> const& unit_vec, sidis::kin::Kinematics* kin) const noexcept;

	Double operator()(Point<6> const& unit_vec) const noexcept {
		return eval(unit_vec);
	}
};

// Map the radiative cross-section onto unit hypercube for Monte-Carlo.
class RadDensity final {
	sidis::cut::Cut _cut;
	sidis::cut::CutRad _cut_rad;
	sidis::sf::SfSet const& _sf;
	RcMethod _rc_method;
	Double _soft_threshold;
	sidis::part::Particles _ps;
	Double _S;
	Double _beam_pol;
	sidis::math::Vec3 _target_pol;

public:
	RadDensity(Params& params, sidis::sf::SfSet const& sf);
	Double eval(Point<9> const& unit_vec, sidis::kin::KinematicsRad* kin) const noexcept;
	Double eval(Point<9> const& unit_vec) const noexcept;
	Double transform(Point<9> const& unit_vec, sidis::kin::KinematicsRad* kin) const noexcept;

	Double operator()(Point<9> const& unit_vec) const noexcept {
		return eval(unit_vec);
	}
};

// Tagged union for different cross-section types.
struct Density final {
	EventType event_type;
	union Impl {
		NradDensity nrad;
		RadDensity rad;
		Impl() { }
	} density;

	Density(EventType event_type, Params& params, sidis::sf::SfSet const& sf) :
			event_type(event_type) {
		switch (event_type) {
		case EventType::NRAD:
			new (&density.nrad) NradDensity(params, sf);
			break;
		case EventType::RAD:
			new (&density.rad) RadDensity(params, sf);
			break;
		default:
			UNREACHABLE();
		}
	}
};

// Types of probability distributions available.
enum class DistType {
	// Uniform distribution. Default.
	UNIFORM,
	// k-d tree constructed to approximate density.
	FOAM,
	// Variation of FOAM, with different construction process.
	BUBBLE,
	// Stratified VEGAS algorithm.
	VEGAS,
};

// Identifying name for each type of distribution.
inline char const* dist_type_name(DistType type) {
	switch (type) {
	case DistType::UNIFORM:
		return "uniform";
	case DistType::FOAM:
		return "foam";
	case DistType::BUBBLE:
		return "bubble";
	case DistType::VEGAS:
		return "vegas";
	default:
		UNREACHABLE();
	}
}

// Random number engine type used by Monte-Carlo generators.
using RndEngine = std::mt19937_64;

// Monte-Carlo event drawn from the unit hypercube.
template<std::size_t D>
struct UnitEvent final {
	Point<D> vec;
	Double weight;
};

struct UniformParams { };
struct FoamParams { };
using BubbleParams = bubble::CellBuilderParams<Double>;
struct VegasParams { };

// Parameters for building a probability distribution.
struct DistParams final {
	// Tagged union over different distribution types.
	DistType dist_type;
	union Impl {
		UniformParams uniform;
		FoamParams foam;
		BubbleParams bubble;
		VegasParams vegas;
		Impl() { }
	} params;

	DistParams(EventType event_type, Params& params);
};

// Underlying engines.
struct UniformEngine { };
struct FoamEngine { };
template<std::size_t D>
using BubbleEngine = bubble::CellGenerator<D, Double>;
struct VegasEngine { };

// Random number distribution over the unit hypercube of dimension `D`. This
// type wraps several possible underlying engines, as described by the enum
// `DistType`.
template<std::size_t D>
class Dist final {
	// Tagged union over the underlying engine.
	DistType _dist_type;
	// This flag is used to ensure exception safety. If false, then none of the
	// union variants are active. Initialize to false, set to true right after
	// any placement new, set to false right before any destructor. Only
	// destruct if true. That way, even if the constructors/destructors
	// themselves fail, no undefined behaviour can happen.
	bool _dist_valid;
	union Impl {
		UniformEngine uniform;
		FoamEngine foam;
		BubbleEngine<D> bubble;
		VegasEngine vegas;
		Impl() { }
		~Impl() { }
	} _engine;

public:
	// It's annoying to implement these, and for now only the move constructor
	// is needed, so delete the others.
	Dist(Dist<D> const& other) = delete;
	Dist(Dist<D>&& other) noexcept;
	Dist<D>& operator=(Dist<D> const& other) = delete;
	Dist<D>& operator=(Dist<D>&& other) noexcept = delete;
	~Dist();

	// Constructs an empty distribution that provides events evenly throughout
	// the unit hypercube.
	explicit Dist(DistType dist_type);

	// Draws a weighted event from the distribution.
	UnitEvent<D> draw(RndEngine& rnd) const;

	// Average weight of events produced by the distribution.
	Double prime() const;

	// Serialization to binary streams.
	static std::ostream& write(std::ostream& os, Dist<D> const& dist);
	static std::istream& read(std::istream& is, Dist<D>& dist);

	// Constructs a distribution that approximates the provided density.
	template<std::size_t D1, typename F1>
	friend Dist<D1> build_dist_approx(DistParams const& dist_params, F1 const& density);
};

// Generated event from a cross-section.
struct Event final {
	// Tagged union over different event types (non-radiative, radiative, ...).
	EventType event_type;
	Double weight;
	union Impl {
		sidis::kin::Kinematics nrad;
		sidis::kin::KinematicsRad rad;
	} kin;
};

// Generator for producing Monte-Carlo events from a cross-section.
class Generator final {
	EventType _event_type;
	union DensityImpl {
		NradDensity nrad;
		RadDensity rad;
		DensityImpl() { }
	} _density;
	// Need a validity flag for the distribution because it has a non-trivial
	// destructor.
	bool _dist_valid;
	union DistImpl {
		Dist<6> nrad;
		Dist<9> rad;
		DistImpl() { }
		~DistImpl() { }
	} _dist;

public:
	Generator(Generator const& other) = delete;
	Generator(Generator&& other) noexcept;
	Generator& operator=(Generator const& other) = delete;
	Generator& operator=(Generator&& other) noexcept = delete;
	~Generator();

	Generator(Density density);

	EventType event_type() const {
		return _event_type;
	}

	// Builds the underlying distribution to approximate the cross-section.
	void build_dist(DistParams dist_params);

	Event draw(RndEngine& rnd) const;
	Double prime() const;

	// Serialization to binary streams of underlying distribution.
	static std::ostream& write_dist(std::ostream& os, Generator const& gen);
	static std::istream& read_dist(std::istream& is, Generator& gen);
};

#include "generator.ipp"

#endif

