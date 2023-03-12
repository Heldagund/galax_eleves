#ifdef GALAX_MODEL_CPU_FMM

#ifndef MODEL_CPU_FMM_HPP_
#define MODEL_CPU_FMM_HPP_

#include "../Model_CPU.hpp"
#include "../../../Particles.hpp"
#include <bitset>
#include <iostream>
#include <string>

#define N_CRIT 10
#define R_ROOT 64


class Cell
{
public:
    Cell(float r = 0.0, float x = 0.0, float y = 0.0, float z = 0.0, Cell* parent = nullptr)
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
        if(!child)
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
    float multipole[10];
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

    void MeasureDepth(Cell* entry)
    {
        if(entry->child)
        {
            for(int octant=0; octant<8; octant++)
            {
                if(entry->nchild & (1 << octant))
                {
                    MeasureDepth(entry->child[octant]);
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
        std::cout << "(" << entry->x << "," << entry->y << "," << entry->z << "," << entry->r;
        if(entry->child)
        {
            std::bitset<8> nc(entry->nchild);
            std::cout << "," << nc.to_string() << ")" << std::endl;
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
            if(entry->r < r_min)
            {
                r_min = entry->r;
            }
            std::cout << " leaf: ";
            for(int i=0; i< entry->nleaf; i++)
            {
                std::cout << entry->leaf[i] << " ";
            }
            std::cout << ")" << std::endl;
        }
    }

private:
    Cell* root;
    float r_min;
};

#endif // MODEL_CPU_FMM_HPP_

#endif // GALAX_MODEL_CPU_FMM
