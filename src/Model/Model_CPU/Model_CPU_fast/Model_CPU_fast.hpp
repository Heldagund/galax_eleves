#ifdef GALAX_MODEL_CPU_FAST

#ifndef MODEL_CPU_FAST_HPP_
#define MODEL_CPU_FAST_HPP_

#include "../Model_CPU.hpp"

class Model_CPU_fast : public Model_CPU
{
public:
    Model_CPU_fast(const Initstate& initstate, Particles& particles);

    virtual ~Model_CPU_fast() = default;

    virtual void step();

private:
	// std::vector<float> posx_ex;
	// std::vector<float> posy_ex;
	// std::vector<float> posz_ex;
	// std::vector<float> mass_ex;
};
#endif // MODEL_CPU_FAST_HPP_

#endif // GALAX_MODEL_CPU_FAST
