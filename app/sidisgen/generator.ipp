#ifndef SIDISGEN_GENERATOR_IPP
#define SIDISGEN_GENERATOR_IPP

#include "generator.hpp"

#include <utility>

template<std::size_t D>
Dist<D>::Dist(DistType dist_type) :
		_dist_type(dist_type),
		_dist_valid(false) {
	switch (_dist_type) {
	case DistType::UNIFORM:
		new (&_engine.uniform) UniformEngine();
		_dist_valid = true;
		break;
	case DistType::FOAM:
		new (&_engine.foam) FoamEngine();
		_dist_valid = true;
		break;
	case DistType::BUBBLE:
		new (&_engine.bubble) BubbleEngine<D>();
		_dist_valid = true;
		break;
	case DistType::VEGAS:
		new (&_engine.vegas) VegasEngine();
		_dist_valid = true;
		break;
	default:
		UNREACHABLE();
	}
}

template<std::size_t D>
Dist<D>::Dist(Dist<D>&& other) noexcept :
		_dist_type(other._dist_type),
		_dist_valid(false) {
	switch (_dist_type) {
	case DistType::UNIFORM:
		new (&_engine.uniform) UniformEngine(std::move(other._engine.uniform));
		_dist_valid = true;
		break;
	case DistType::FOAM:
		new (&_engine.foam) FoamEngine(std::move(other._engine.foam));
		_dist_valid = true;
		break;
	case DistType::BUBBLE:
		new (&_engine.bubble) BubbleEngine<D>(std::move(other._engine.bubble));
		_dist_valid = true;
		break;
	case DistType::VEGAS:
		new (&_engine.vegas) VegasEngine(std::move(other._engine.vegas));
		_dist_valid = true;
		break;
	default:
		UNREACHABLE();
	}
}

template<std::size_t D>
Dist<D>::~Dist() {
	if (_dist_valid) {
		switch (_dist_type) {
		case DistType::UNIFORM:
			_dist_valid = false;
			_engine.uniform.~UniformEngine();
			break;
		case DistType::FOAM:
			_dist_valid = false;
			_engine.foam.~FoamEngine();
			break;
		case DistType::BUBBLE:
			_dist_valid = false;
			_engine.bubble.~BubbleEngine<D>();
			break;
		case DistType::VEGAS:
			_dist_valid = false;
			_engine.vegas.~VegasEngine();
			break;
		default:
			UNREACHABLE();
		}
	}
}

template<std::size_t D> 
UnitEvent<D> Dist<D>::draw(RndEngine& rnd) const {
	UnitEvent<D> event;
	switch (_dist_type) {
	case DistType::UNIFORM:
		{
			std::uniform_real_distribution<Double> dist;
			event.weight = 1.;
			for (std::size_t dim = 0; dim < D; ++dim) {
				event.vec[dim] = dist(rnd);
			}
		}
		break;
	case DistType::FOAM:
		throw "Error";
	case DistType::BUBBLE:
		_engine.bubble.generate(rnd, &event.weight, &event.vec);
		break;
	case DistType::VEGAS:
		throw "Error";
	default:
		UNREACHABLE();
	}
	return event;
}

template<std::size_t D>
Double Dist<D>::prime() const {
	switch (_dist_type) {
	case DistType::UNIFORM:
		return 1.;
	case DistType::FOAM:
		throw "Error";
	case DistType::BUBBLE:
		return _engine.bubble.prime();
	case DistType::VEGAS:
		throw "Error";
	default:
		UNREACHABLE();
	}
}

template<std::size_t D>
std::ostream& Dist<D>::write(std::ostream& os, Dist<D> const& dist) {
	std::size_t dim = D;
	if (!os.write(reinterpret_cast<char const*>(&dist._dist_type), sizeof(DistType))) {
		goto error;
	}
	if (!os.write(reinterpret_cast<char const*>(&dim), sizeof(std::size_t))) {
		goto error;
	}
	switch (dist._dist_type) {
	case DistType::UNIFORM:
		goto error;
	case DistType::FOAM:
		goto error;
	case DistType::BUBBLE:
		return dist._engine.bubble.write(os);
	case DistType::VEGAS:
		goto error;
	default:
		UNREACHABLE();
	}
error:
	os.setstate(std::ios_base::failbit);
	return os;
}

template<std::size_t D>
std::istream& Dist<D>::read(std::istream& is, Dist<D>& dist) {
	DistType dist_type_read;
	std::size_t dim;
	if (!is.read(reinterpret_cast<char*>(&dist_type_read), sizeof(DistType))) {
		goto error;
	}
	if (!is.read(reinterpret_cast<char*>(&dim), sizeof(std::size_t))) {
		goto error;
	}
	if (dim != D) {
		goto error;
	}
	dist.~Dist<D>();
	// Create new engine.
	switch (dist_type_read) {
	case DistType::UNIFORM:
		dist._dist_type = DistType::UNIFORM;
		new (&dist._engine.uniform) UniformEngine();
		dist._dist_valid = true;
		goto error;
	case DistType::FOAM:
		dist._dist_type = DistType::FOAM;
		new (&dist._engine.foam) FoamEngine();
		dist._dist_valid = true;
		goto error;
	case DistType::BUBBLE:
		dist._dist_type = DistType::BUBBLE;
		new (&dist._engine.bubble) BubbleEngine<D>();
		dist._dist_valid = true;
		dist._engine.bubble.read(is);
		break;
	case DistType::VEGAS:
		dist._dist_type = DistType::VEGAS;
		new (&dist._engine.vegas) VegasEngine();
		dist._dist_valid = true;
		goto error;
	default:
		goto error;
	}
	return is;
error:
	is.setstate(std::ios_base::failbit);
	return is;
}

template<std::size_t D, typename F>
inline Dist<D> build_dist_approx(DistParams const& dist_params, F const& density) {
	Dist<D> dist(dist_params.dist_type);
	switch (dist_params.dist_type) {
	case DistType::UNIFORM:
		break;
	case DistType::FOAM:
		throw "Error";
	case DistType::BUBBLE:
		{
			bubble::CellBuilder<D, Double, F> builder(density, nullptr, nullptr);
			builder.par = dist_params.params.bubble;
			builder.explore();
			builder.tune();
			dist._engine.bubble = bubble::make_generator(builder);
		}
		break;
	case DistType::VEGAS:
		throw "Error";
	default:
		UNREACHABLE();
	}
	return dist;
}

#endif

