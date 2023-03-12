#ifdef GALAX_MODEL_CPU_FMM

#ifndef MODEL_CPU_FMM_HPP_
#define MODEL_CPU_FMM_HPP_

#include "../Model_CPU.hpp"

class Model_CPU_FMM : public Model_CPU
{
public:
    Model_CPU_FMM(const Initstate& initstate, Particles& particles);

    virtual ~Model_CPU_FMM() = default;

    virtual void step();
};
#endif // MODEL_CPU_FMM_HPP_

#endif // GALAX_MODEL_CPU_FMM
