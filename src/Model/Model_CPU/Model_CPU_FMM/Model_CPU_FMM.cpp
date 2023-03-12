#ifdef GALAX_MODEL_CPU_FAST

#include <cmath>

#include "Model_CPU_fast.hpp"

#include <xsimd/xsimd.hpp>
#include <omp.h>

namespace xs = xsimd;
using b_type = xs::batch<float, xs::avx2>;

Model_CPU_fast
::Model_CPU_fast(const Initstate& initstate, Particles& particles)
: Model_CPU(initstate, particles)
{
}

void Model_CPU_fast
::step()
{
    std::fill(accelerationsx.begin(), accelerationsx.end(), 0);
    std::fill(accelerationsy.begin(), accelerationsy.end(), 0);
    std::fill(accelerationsz.begin(), accelerationsz.end(), 0);

// OMP  version
    // #pragma omp parallel for // collapse(2)
	// for (int i = 0; i < n_particles; i++)
	// {
	// 	for (int j = 0; j < n_particles; j++)
	// 	{
	// 		if(i != j)
	// 		{
	// 			const float diffx = particles.x[j] - particles.x[i];
	// 			const float diffy = particles.y[j] - particles.y[i];
	// 			const float diffz = particles.z[j] - particles.z[i];

	// 			float dij = diffx * diffx + diffy * diffy + diffz * diffz;

	// 			if (dij < 1.0)
	// 			{
	// 				dij = 10.0;
	// 			}
	// 			else
	// 			{
	// 				dij = std::sqrt(dij);
	// 				dij = 10.0 / (dij * dij * dij);
	// 			}

	// 			accelerationsx[i] += diffx * dij * initstate.masses[j];
	// 			accelerationsy[i] += diffy * dij * initstate.masses[j];
	// 			accelerationsz[i] += diffz * dij * initstate.masses[j];
	// 		}
	// 	}
	// }

    // #pragma omp parallel for
	// for (int i = 0; i < n_particles; i++)
	// {
	// 	velocitiesx[i] += accelerationsx[i] * 2.0f;
	// 	velocitiesy[i] += accelerationsy[i] * 2.0f;
	// 	velocitiesz[i] += accelerationsz[i] * 2.0f;
	// 	particles.x[i] += velocitiesx   [i] * 0.1f;
	// 	particles.y[i] += velocitiesy   [i] * 0.1f;
	// 	particles.z[i] += velocitiesz   [i] * 0.1f;
	// }

// OMP + xsimd version
	const b_type dij_threshold = b_type::broadcast(1.0);
	const b_type dij_min = b_type::broadcast(10.0);
	auto vec_size = n_particles-n_particles%b_type::size;

	#pragma omp parallel for
		for (int i = 0; i < vec_size; i += b_type::size)
		{
			b_type raccx_i = b_type::load_unaligned(&accelerationsx[i]);
			b_type raccy_i = b_type::load_unaligned(&accelerationsy[i]);
			b_type raccz_i = b_type::load_unaligned(&accelerationsz[i]);

			// load registers body i
			const b_type rposx_i = b_type::load_unaligned(&particles.x[i]);
			const b_type rposy_i = b_type::load_unaligned(&particles.y[i]);
			const b_type rposz_i = b_type::load_unaligned(&particles.z[i]);

			for (int j = 0; j < n_particles; j++)
			{
				if(i != j)
				{
					const b_type rposx_j = b_type::load_unaligned(&particles.x[j]);
					const b_type rposy_j = b_type::load_unaligned(&particles.y[j]);
					const b_type rposz_j = b_type::load_unaligned(&particles.z[j]);
					const b_type rmass_j = b_type::load_unaligned(&initstate.masses[j]);

					b_type diffx = rposx_j - rposx_i;
					b_type diffy = rposy_j - rposy_i;
					b_type diffz = rposz_j - rposz_i;
					b_type dij = diffx * diffx + diffy * diffy + diffz * diffz;

					auto mask_identidy = xs::le(dij, b_type::broadcast(0.000001));
					auto mask_threshold = xs::lt(dij, dij_threshold);
					
					dij =  xs::sqrt(dij);
					dij =  dij = 10.0 / (dij * dij * dij);
					dij = xs::select(mask_threshold, dij_min, dij);	

					raccx_i += diffx * dij * rmass_j;
					raccy_i += diffy * dij * rmass_j;
					raccz_i += diffz * dij * rmass_j;
				}
			}

			raccx_i.store_unaligned(&accelerationsx[i]);
			raccy_i.store_unaligned(&accelerationsy[i]);
			raccz_i.store_unaligned(&accelerationsz[i]);
		}

	#pragma omp parallel for // collapse(2)
	for (int i = vec_size; i < n_particles; i++)
	{
		for (int j = 0; j < n_particles; j++)
		{
			if(i != j)
			{
				const float diffx = particles.x[j] - particles.x[i];
				const float diffy = particles.y[j] - particles.y[i];
				const float diffz = particles.z[j] - particles.z[i];

				float dij = diffx * diffx + diffy * diffy + diffz * diffz;

				if (dij < 1.0)
				{
					dij = 10.0;
				}
				else
				{
					dij = std::sqrt(dij);
					dij = 10.0 / (dij * dij * dij);
				}

				accelerationsx[i] += diffx * dij * initstate.masses[j];
				accelerationsy[i] += diffy * dij * initstate.masses[j];
				accelerationsz[i] += diffz * dij * initstate.masses[j];
			}
		}
	}

	#pragma omp parallel for
	for (int i = 0; i < n_particles; i++)
	{
		velocitiesx[i] += accelerationsx[i] * 2.0f;
		velocitiesy[i] += accelerationsy[i] * 2.0f;
		velocitiesz[i] += accelerationsz[i] * 2.0f;
		particles.x[i] += velocitiesx   [i] * 0.1f;
		particles.y[i] += velocitiesy   [i] * 0.1f;
		particles.z[i] += velocitiesz   [i] * 0.1f;
	}


    // #pragma omp parallel for
	// for (int i = 0; i < n_particles; i += b_type::size)
	// {
	// 	b_type rx_i = b_type::load_unaligned(&accelerationsx[i]);
	// 	b_type raccy_i = b_type::load_unaligned(&accelerationsy[i]);
	// 	b_type raccz_i = b_type::load_unaligned(&accelerationsz[i]);

	// }

	// #pragma omp parallel for
	// for (int i = vec_size; i < n_particles; i ++)
	// {
	// 	velocitiesx[i] += accelerationsx[i] * 2.0f;
	// 	velocitiesy[i] += accelerationsy[i] * 2.0f;
	// 	velocitiesz[i] += accelerationsz[i] * 2.0f;
	// 	particles.x[i] += velocitiesx   [i] * 0.1f;
	// 	particles.y[i] += velocitiesy   [i] * 0.1f;
	// 	particles.z[i] += velocitiesz   [i] * 0.1f;
	// }
}

#endif // GALAX_MODEL_CPU_FAST
