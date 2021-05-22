#ifndef SIDISGEN_EVENT_GENERATOR_HPP
#define SIDISGEN_EVENT_GENERATOR_HPP

#include <memory>
#include <vector>

#include <sidis/sidis.hpp>

#include <TFile.h>
#include <TFoam.h>
#include <TFoamIntegrand.h>
#include <TRandom3.h>

#include "event_stats.hpp"
#include "params.hpp"

class EventDensity : public TFoamIntegrand {
public:
	int const dim;
	EventDensity(int dim);

	virtual Double_t Transform(Double_t const* vec_unit, Double_t* vec_out) const = 0;
	virtual Double_t Density(int dim, Double_t* vec_unit) override = 0;
};

class EventGenerator {
private:
	EventType _type;
	std::unique_ptr<TFile> _foam_file;
	std::unique_ptr<EventDensity> _density;
	TRandom3* _random;
	TFoam* _foam;
	EventStats _stats;
	std::vector<Double_t> _vec_unit;

	EventGenerator() :
		_type(EventType::NRAD),
		_foam_file(),
		_density(),
		_foam(nullptr),
		_stats() { }
	EventGenerator(
		EventType type,
		std::unique_ptr<TFile> foam_file,
		std::unique_ptr<EventDensity> density,
		TRandom3* random,
		TFoam* foam);

	static std::unique_ptr<EventDensity> alloc_density(
		EventType type,
		Params const& params,
		sidis::sf::SfSet const& sf);

public:
	static EventGenerator write(
		EventType type,
		Params const& params,
		sidis::sf::SfSet const& sf,
		TRandom3* random);
	static EventGenerator read(
		EventType type,
		Params const& params,
		sidis::sf::SfSet const& sf,
		TRandom3* random);

	Double_t generate(Double_t* vec);

	EventType type() const {
		return _type;
	}
	int dim() const {
		return _density->dim;
	}
	EventStats const& stats() const {
		return _stats;
	}
};

#endif

