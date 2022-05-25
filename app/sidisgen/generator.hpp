#ifndef SIDISGEN_GENERATOR_HPP
#define SIDISGEN_GENERATOR_HPP

#include <istream>
#include <ostream>

#include <bubble.hpp>

#include <sidis/sidis.hpp>

#include "event_type.hpp"
#include "params.hpp"

using Real = sidis::Real;
using Seed = bubble::Seed;
template<bubble::Dim DIM>
using Point = bubble::Point<DIM, Real>;
using Stats = bubble::Stats<Real>;
using StatsAccum = bubble::StatsAccum<Real>;

// Each of the following functions computes a different cross-section on the
// unit hypercube. The full phase space is rescaled to fit on a hypercube for
// the convenience of the underlying event generator, but a `transform` method
// is provided that takes a unit vector to the full 9-dimensional phase space.
// The `transform` method also returns the Jacobian of the transformation.

class NradDensity {
	Params _params;
	sidis::cut::Cut _cut;
	sidis::part::Particles _ps;
	Real _S;
	sidis::sf::SfSet const& _sf;

public:
	NradDensity(Params const& params, sidis::sf::SfSet const& sf);
	Real operator()(Point<6> const& vec) const noexcept;
	Real transform(Point<6> const& unit_vec, Point<6>* ph_vec) const noexcept;
};

class RadDensity {
	Params _params;
	sidis::cut::Cut _cut;
	sidis::cut::CutRad _cut_rad;
	sidis::part::Particles _ps;
	Real _S;
	sidis::sf::SfSet const& _sf;

public:
	RadDensity(Params const& params, sidis::sf::SfSet const& sf);
	Real operator()(Point<9> const& vec) const noexcept;
	Real transform(Point<9> const& unit_vec, Point<9>* ph_vec) const noexcept;
};

class ExclDensity {
public:
	// For now, the exclusive contribution is not supported.
	// TODO: Throw an exception to ensure this doesn't get called by accident.
	ExclDensity(Params const&, sidis::sf::SfSet const&) { }
	Real operator()(Point<8> const&) const noexcept { return 0.; }
	Real transform(Point<8> const&, Point<8>*) const noexcept { return 0.; };
};

struct BuilderReporters {
	bubble::ExploreProgressReporter<Real>* explore_progress = nullptr;
	bubble::TuneProgressReporter<Real>* tune_progress = nullptr;
};

// Builds a generator and writes it to a stream.
class Builder {
	// Internally, use a tagged union to distinguish the different possible
	// phase space dimensions from each other.
	using NradBuilder = bubble::CellBuilder<6, Real, NradDensity>;
	using RadBuilder = bubble::CellBuilder<9, Real, RadDensity>;
	using ExclBuilder = bubble::CellBuilder<8, Real, ExclDensity>;
	EventType const _type;
	union BuilderImpl {
		NradBuilder nrad;
		RadBuilder rad;
		ExclBuilder excl;
		BuilderImpl() { }
		~BuilderImpl() { }
	} _builder;
	Int_t _seed;

public:

	Builder(
		EventType type,
		BuilderReporters const& reporters,
		Params const& params,
		sidis::sf::SfSet const& sf);
	Builder(Builder const& other) = delete;
	Builder(Builder&& other);
	Builder& operator=(Builder const& other) = delete;
	~Builder();

	EventType type() const {
		return _type;
	}
	Int_t seed() const {
		return _seed;
	}
	// Constructs the generator.
	void explore();
	void tune();
	// Writes the generator to a binary stream.
	void write(std::ostream& os);
	// Calculates the relative variance (with uncertainty).
	Real rel_var(Real* err_out=nullptr) const;
	// Number of cells.
	std::size_t size() const;
};

// Loads a generator from a stream and produces events from it.
class Generator {
	// Uses a tagged union, similar to `Builder`.
	using NradGenerator = bubble::CellGenerator<6, Real, NradDensity>;
	using RadGenerator = bubble::CellGenerator<9, Real, RadDensity>;
	using ExclGenerator = bubble::CellGenerator<8, Real, ExclDensity>;
	union GeneratorImpl {
		NradGenerator nrad;
		RadGenerator rad;
		ExclGenerator excl;
		GeneratorImpl() { }
		~GeneratorImpl() { }
	} _generator;
	EventType _type;
	Int_t _seed;

	// Stores statistics about the produced weights.
	StatsAccum _weights;
	std::size_t _count;
	std::size_t _count_acc;

public:
	Generator(
		EventType type,
		Params const& params,
		sidis::sf::SfSet const& sf,
		std::istream& is);
	Generator(Generator const& other) = delete;
	Generator(Generator&& other);
	Generator& operator=(Generator const& other) = delete;
	~Generator();

	EventType type() const {
		return _type;
	}
	Int_t seed() const {
		return _seed;
	}
	// Produces an event with a certain weight. Sometimes, the generator may
	// produce an event with weight zero. In this case, the event is rejected,
	// and a new event is sampled to replace it. Rejected events are still
	// tracked in the resulting statistics, but never returned from this
	// function.
	Real generate(Real* ph_out, Real* unit_out);
	Real generate(Real* ph_out) {
		Real unit_out[9];
		return generate(ph_out, unit_out);
	}
	// Gets total number of events generated.
	std::size_t count() const {
		return _count;
	}
	// Gets total number of events generated that were accepted.
	std::size_t count_acc() const {
		return _count_acc;
	}
	// Gets acceptance fraction of events.
	Real acceptance() const {
		return static_cast<Real>(_count_acc) / _count;
	}
	// Gets statistics of the provided weights.
	Stats weights() const {
		return _weights.total();
	}
	// Gets statistics of the provided weights that were accepted.
	Stats weights_acc() const {
		Stats stats = _weights.total();
		stats.rescale_count(_count_acc);
		return stats;
	}

	// Gets the overall prime of the generator. Combine with the average weight
	// to get the cross-section.
	Real prime() const;
};

#endif

