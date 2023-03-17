// #ifdef GALAX_MODEL_CPU_FMM

// #include <cmath>

// #include "Model_CPU_FMM.hpp"

// #include <xsimd/xsimd.hpp>
// #include <omp.h>

// namespace xs = xsimd;
// using b_type = xs::batch<float, xs::avx2>;

// Model_CPU_FMM
// ::Model_CPU_FMM(const Initstate& initstate, Particles& particles)
// : Model_CPU(initstate, particles)
// {
// 	root = new Cell(R_ROOT);
// 	r_min = R_ROOT;
// }

// Model_CPU_FMM::~Model_CPU_FMM()
// {
// 	delete root;
// }

// void Model_CPU_FMM
// ::step()
// {
// 	BuildTree();
// 	//MeasureTreeDepth(root);
// 	GetMultipole(root);
// 	//std::cout << "r_min = " << r_min <<std::endl;
// 	ResetTree(root);
// }

// void Model_CPU_FMM::Evaluate(Cell* entry, int index_target)
// {
// 	if(entry->child)
// 	{
// 		for(int octant=0; octant<8; octant++)
// 		{
// 			if(entry->nchild & (1 << octant))
// 			{
// 				Cell* c = entry->child[octant];
// 				float distance = std::sqrt( (particles.x[index_target] - c->x)*(particles.x[index_target] - c->x)
// 											+ (particles.y[index_target] - c->y)*(particles.y[index_target] - c->y)
// 											+ (particles.z[index_target] - c->z)*(particles.z[index_target] - c->z));
// 				if(c.r > THETA*distance)
// 				{
// 					Evaluate(c, index_target)
// 				}
// 				else
// 				{
// 					float dx = particles.x[index_target] - c->x;
// 					float dy = particles.y[index_target] - c->y;
// 					float dz = particles.z[index_target] - c->z;
// 					float r5 = std::powf(distance, 5);
// 					float r7 = std::powf(distance, 7);

					
// 					accelerationsx[i] +=
// 					accelerationsy[i] +=
// 					accelerationsz[i] +=
// 				}
// 			}    
// 		}
// 	}
// 	else
// 	{
// 		for(int i=0; i<entry->nleaf; i++)
// 		{
// 			int index_source = entry->leaf[i];
// 			if(index_source != index_target)
// 			{
// 				const float diffx = particles.x[index_source] - particles.x[index_target];
// 				const float diffy = particles.y[index_source] - particles.y[index_target];
// 				const float diffz = particles.z[index_source] - particles.z[index_target];

// 				float dij = diffx * diffx + diffy * diffy + diffz * diffz;

// 				if (dij < 1.0)
// 				{
// 					dij = 10.0;
// 				}
// 				else
// 				{
// 					dij = std::sqrt(dij);
// 					dij = 10.0 / (dij * dij * dij);
// 				}

// 				accelerationsx[i] += diffx * dij * initstate.masses[j];
// 				accelerationsy[i] += diffy * dij * initstate.masses[j];
// 				accelerationsz[i] += diffz * dij * initstate.masses[j];
// 			}
// 		}
// 	}
// }
// #endif // GALAX_MODEL_CPU_FMM
