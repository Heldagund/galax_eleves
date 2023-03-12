#ifdef GALAX_MODEL_CPU_FMM

#include <cmath>

#include "Model_CPU_FMM.hpp"

#include <xsimd/xsimd.hpp>
#include <omp.h>

namespace xs = xsimd;
using b_type = xs::batch<float, xs::avx2>;

Model_CPU_FMM
::Model_CPU_FMM(const Initstate& initstate, Particles& particles)
: Model_CPU(initstate, particles)
{
	root = new Cell(R_ROOT);
	r_min = R_ROOT;
}

Model_CPU_FMM::~Model_CPU_FMM()
{
	delete root;
}

void Model_CPU_FMM
::step()
{
	BuildTree();
	MeasureDepth(root);
	std::cout << "r_min = " << r_min <<std::endl;
}

#endif // GALAX_MODEL_CPU_FMM
