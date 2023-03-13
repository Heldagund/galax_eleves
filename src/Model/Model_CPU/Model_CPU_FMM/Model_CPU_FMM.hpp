#ifdef GALAX_MODEL_CPU_FMM

#ifndef MODEL_CPU_FMM_HPP_
#define MODEL_CPU_FMM_HPP_

#include "../Model_CPU.hpp"
#include "../../../Particles.hpp"
#include <bitset>
#include <cassert>
#include <cstring>
#include <iostream>
#include <string>

#define N_CRIT 100
#define R_ROOT 64

class Cell
{
public:
    Cell(float r = 0.0f, float x = 0.0f, float y = 0.0f, float z = 0.0f, Cell* parent = nullptr)
    : r(r), x(x), y(y), z(z), parent(parent), nleaf(0), nchild(0)
    {
        leaf = new int[N_CRIT];
        child = nullptr;
    }
    ~Cell()
    {
        if(leaf)
        {
            delete[] leaf;
        }
        if(child)
        {
            for(int i=0; i<8; i++)
            {
                if(nchild & (1 << i))
                {
                    delete child[i];
                }
            }
            delete[] child;
        }
    }

    void AddChild(int octant)
    {
        if(child == nullptr)
        {
            child = new Cell*[8];
        }
        float r_child = r/2;
        child[octant] = new Cell(r_child, 
                                 x + r_child*((octant & 1) * 2 - 1),
                                 y + r_child*((octant & 2) - 1),
                                 z + r_child*((octant & 4) / 2 - 1), 
                                 this);
        nchild = nchild | (1 << octant);
    }
    
    void SplitCell(Particles& particles)
    {
        for(int i=0; i<nleaf; i++)
        {
            int l = leaf[i];
            int octant = (particles.x[l] > x) + ((particles.y[l] > y) << 1) 
                          + ((particles.z[l] > z) << 2);
            if(!(nchild & (1 << octant)))
            {
                AddChild(octant);
            }
            Cell* c = child[octant];
            c->leaf[c->nleaf] = l;
            c->nleaf += 1;
            if(c->nleaf >= N_CRIT)
            {
                c->SplitCell(particles);
            }
        }
        //leaf is no longer needed
        if(leaf)
        {
            delete[] leaf;
            leaf = nullptr;
        }
    }

public:
    int nleaf;
    int* leaf;
    int nchild;
    Cell** child;
    Cell* parent;
    float x;
    float y;
    float z;
    float r;
    float multipole[10] = {0.0f};
};

class Model_CPU_FMM : public Model_CPU
{
public:
    Model_CPU_FMM(const Initstate& initstate, Particles& particles);

    virtual ~Model_CPU_FMM();

    virtual void step();

    void BuildTree()
    {
        for(int i=0; i<n_particles; i++)
        {
            Cell* curr = root;
            while(curr->child)
            {
                curr->nleaf += 1;
                int octant = (particles.x[i] > curr->x) + ((particles.y[i] > curr->y) << 1)
                              + ((particles.z[i] > curr->z) << 2);
                if(!(curr->nchild & (1 << octant)))
                {
                    curr->AddChild(octant);
                }
                curr = curr->child[octant];
            }
            curr->leaf[curr->nleaf] = i;
            curr->nleaf += 1;
            if(curr->nleaf >= N_CRIT)
            {
                curr->SplitCell(particles);
            }
        }
    }

    void GetMultipole(Cell* entry)
    {
        if(entry->child)
        {
            for(int octant=0; octant<8; octant++)
            {
                if(entry->nchild & (1 << octant))
                {
                    Cell* c = entry->child[octant];
                    GetMultipole(c);

                    float dx = entry->x - c->x;
                    float dy = entry->y - c->y;
                    float dz = entry->z - c->z;
                    entry->multipole[0] += c->multipole[0];
                    entry->multipole[1] += c->multipole[1] + c->multipole[0]*dx;
                    entry->multipole[2] += c->multipole[2] + c->multipole[0]*dy;
                    entry->multipole[3] += c->multipole[3] + c->multipole[0]*dz;
                    entry->multipole[4] += c->multipole[4] + c->multipole[1]*dx + 0.5*c->multipole[0]*dx*dx;
                    entry->multipole[5] += c->multipole[5] + c->multipole[2]*dy + 0.5*c->multipole[0]*dy*dy;
                    entry->multipole[6] += c->multipole[6] + c->multipole[3]*dz + 0.5*c->multipole[0]*dz*dz;
                    entry->multipole[7] += c->multipole[7] + 0.5*dx*c->multipole[2] + 0.5*dy*c->multipole[1] + 0.5*dx*dy*c->multipole[0];
                    entry->multipole[8] += c->multipole[8] + 0.5*dy*c->multipole[3] + 0.5*dz*c->multipole[2] + 0.5*dy*dz*c->multipole[0];
                    entry->multipole[9] += c->multipole[9] + 0.5*dz*c->multipole[1] + 0.5*dx*c->multipole[3] + 0.5*dz*dx*c->multipole[0];
                }
            }
        }
        else
        {
            for(int i=0; i< entry->nleaf; i++)
            {
                int l = entry->leaf[i];
                float dx = entry->x - particles.x[l];
                float dy = entry->y - particles.y[l];
                float dz = entry->z - particles.z[l];
                entry->multipole[0] += initstate.masses[l];
                entry->multipole[1] += initstate.masses[l] * dx;
                entry->multipole[2] += initstate.masses[l] * dy;
                entry->multipole[3] += initstate.masses[l] * dz;
                entry->multipole[4] += initstate.masses[l] * dx*dx/2;
                entry->multipole[5] += initstate.masses[l] * dy*dy/2; 
                entry->multipole[6] += initstate.masses[l] * dz*dz/2;
                entry->multipole[7] += initstate.masses[l] * dx*dy/2;
                entry->multipole[8] += initstate.masses[l] * dy*dz/2;
                entry->multipole[9] += initstate.masses[l] * dz*dx/2;
            }
        }
    }

    void ResetTree(Cell* entry)
    {
        if(entry->child)
        {
            entry->nleaf = 0;
            memset(entry->multipole, 0, 10*sizeof(float));
            for(int octant=0; octant<8; octant++)
            {
                if(entry->nchild & (1 << octant))
                {
                    ResetTree(entry->child[octant]);
                }
            }
        }
        else
        {
            entry->nleaf = 0;
            memset(entry->multipole, 0, 10*sizeof(float));
        }
    }

    void MeasureTreeDepth(Cell* entry)
    {
        if(entry->child)
        {
            for(int octant=0; octant<8; octant++)
            {
                if(entry->nchild & (1 << octant))
                {
                    MeasureTreeDepth(entry->child[octant]);
                }    
            }
        }
        else
        {
            if(entry->r < r_min)
            {
                r_min = entry->r;
            }
        }
    }

    void ShowTree(Cell* entry, std::string prefix="")
    {
        std::cout << "(x=" << entry->x << ", y=" << entry->y << ", z=" << entry->z << ", r=" << entry->r;
        if(entry->child)
        {
            std::bitset<8> nc(entry->nchild);
            std::cout << ", nchild=" << nc.to_string() << ", nleaf=" << entry->nleaf << ")" << std::endl;
            for(int i=0; i<10; i++)
            {
                std::cout << entry->multipole[i] << "  ";
            }
            std::cout <<std::endl;
            prefix += "--";
            for(int octant=0; octant<8; octant++)
            {
                if(entry->nchild & (1 << octant))
                {
                    std::cout << prefix;
                    ShowTree(entry->child[octant], prefix);
                }    
            }
        }
        else
        {
            std::cout << " leaf: ";
            for(int i=0; i< entry->nleaf; i++)
            {
                std::cout << entry->leaf[i] << " ";
            }
            std::cout << ")" << std::endl;

            for(int i=0; i<10; i++)
            {
                std::cout << entry->multipole[i] << "  ";
            }
            std::cout << std::endl;
        }
    }

private:
    Cell* root;
    float r_min;
};

#endif // MODEL_CPU_FMM_HPP_

#endif // GALAX_MODEL_CPU_FMM
