#include "Fluid.h"
#include <iostream>
#include <iomanip>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>
#include <ngl/ShaderLib.h>
#include <ngl/VAOFactory.h>
#include <ngl/Util.h>

std::mt19937 Fluid::m_generator;
auto randomPositivezDist=std::uniform_real_distribution<float>(0.0f,1.0f);

////////// Initialisation //////////

void Fluid::initialise(size_t _sizeX, size_t _sizeY, float _gridsize, float _timestep, const ngl::Vec3 &_force, size_t _right, size_t _left, size_t _top, size_t _bottom, bool _type)
{
    // Apply the settings.
    m_type = _type;
    // Split the frame in to cells according to the grid size.
    // Grid size is changed to a nearest value that can divide 1x1 frame evenly.
    m_cell_res = static_cast<size_t>(1/_gridsize);
    // Add 2 to the simulation sizes for solid cells on the edges.
    m_sizeX = _sizeX+2;
    m_sizeY = _sizeY+2;
    m_resolutionX = m_sizeX*m_cell_res;
    m_resolutionY = m_sizeY*m_cell_res;
    m_timestep = _timestep;
    m_force = _force;

    // Initialise vectors.
    m_boundary.resize(m_resolutionX*m_resolutionY, 0);
    m_pressure.resize(m_resolutionX*m_resolutionY, 0.0f);
    m_velocityX.resize((m_resolutionX+1)*m_resolutionY, 0.0f);
    m_velocityY.resize(m_resolutionX*(m_resolutionY+1), 0.0f);
    m_velocityX_diff.reserve((m_resolutionX+1)*m_resolutionY);
    m_velocityY_diff.reserve(m_resolutionX*(m_resolutionY+1));

    // Initialise VAO.
    m_vao = ngl::vaoFactoryCast<ngl::MultiBufferVAO>(ngl::VAOFactory::createVAO(ngl::multiBufferVAO,GL_POINTS));
    m_vao->bind();
        m_vao->setData(ngl::MultiBufferVAO::VertexData(0,0));
        m_vao->setData(ngl::MultiBufferVAO::VertexData(0,0));
    m_vao->unbind();

    // Sides of the initial fluid position.
    // The input values outside the simulation frame are pushed back to the edges of the frame.
    // It also sets the left and bottom values to right and top values respectively when left > right or bottom > top.
    // When right == left or top == bottom, no particle is initialised.
    size_t right = (std::min(_sizeX,_right)+1)*m_cell_res;
    size_t left = (std::max<size_t>(0,std::min(std::min(_sizeX,_right),_left))+1)*m_cell_res;
    size_t top = (std::min(_sizeY,_top)+1)*m_cell_res;
    size_t bottom = (std::max<size_t>(0,std::min(std::min(_sizeY,_top),_bottom))+1)*m_cell_res;

    // Initialise 4 particles in each fluid cell.
    size_t numParticlesInCell = 4;
    m_numParticles = (top-bottom)*(right-left)*numParticlesInCell;
    m_position.reserve(m_numParticles);
    m_velocity.reserve(m_numParticles);
    for(int j=bottom; j<top; ++j)
    {
        for(int i=left; i<right; ++i)
        {
            for(int k=0; k<numParticlesInCell; ++k)
            {
                float x = (randomPositivezDist(m_generator)+i)/m_cell_res;
                float y = (randomPositivezDist(m_generator)+j)/m_cell_res;
                m_position.push_back({x,y,0.0f});
                m_velocity.push_back(0.0f);
            }
            m_boundary[j*m_resolutionX+i] = 1;
        }
    }
    
    // Set solid cells on the edge.
    for(int i=0; i<m_resolutionY; ++i)
    {
        for(int j=0; j<m_cell_res; ++j)
        {
            m_boundary[i*m_resolutionX+j] = -1;
            m_boundary[i*m_resolutionX+m_resolutionX-1-j] = -1;
        }
    }
    for(int i=0; i<m_resolutionX; ++i)
    {
        for(int j=0; j<m_cell_res; ++j)
        {
            m_boundary[j*m_resolutionX+i] = -1;
            m_boundary[m_resolutionX*m_resolutionY-1-i-j*m_resolutionX] = -1;
        }
    }

    // Add the positions of the solid cells to the vector for visualisation.
    m_solid.reserve(m_cell_res*m_resolutionX*2 + m_cell_res*m_resolutionX*2 - m_cell_res*m_cell_res*4);
    for(int j=0; j<m_resolutionY; ++j)
    {
        for(int i=0; i<m_resolutionX; ++i)
        {
            if(m_boundary[j*m_resolutionX+i] == -1)
            {
                m_solid.push_back({(static_cast<float>(i)+0.5f)/m_cell_res,(static_cast<float>(j)+0.5f)/m_cell_res, 0.0f});
            }
        }
    }
}

void Fluid::simulate()
{
    particleToGrid();
    addForce();
    project();
    calculateDifferences();
    gridToParticle();
    advect();
}

////////// Particle to Grid //////////

// This function transfer velocities from particles to grid using equations 8, 9 and 10 in my report.
void Fluid::particleToGrid()
{
    // Reset the vector.
    std::fill(m_velocityX.begin(), m_velocityX.end(), 0.0f);
    std::fill(m_velocityY.begin(), m_velocityY.end(), 0.0f);
    
    // top part of the equation 8 for u (m_velocityX).
    std::vector<float> ut;
    // bottom part of the equation 8 for u (m_velocityX).
    std::vector<float> ub;
    // top part of the equation 8 for v (m_velocityY).
    std::vector<float> vt;
    // bottom part of the equation 8 for v (m_velocityY).
    std::vector<float> vb;

    ut.resize((m_resolutionX+1)*m_resolutionY, 0.0f);
    ub.resize((m_resolutionX+1)*m_resolutionY, 0.0f);
    vt.resize(m_resolutionX*(m_resolutionY+1), 0.0f);
    vb.resize(m_resolutionX*(m_resolutionY+1), 0.0f);

    // Loop through the particles to calculate the ut, ub, vt, vb.
    for(int k=0; k<m_numParticles; ++k)
    {
        // Find which cell the particle is located.
        float x = m_position[k].m_x*m_cell_res;
        float y = m_position[k].m_y*m_cell_res;
        size_t i = static_cast<size_t>(x);
        size_t j = static_cast<size_t>(y);

        // Calculate the top and bottom part of the equation 8 for u(i-1/2,j). 
        ut[j*(m_resolutionX+1)+i] += m_velocity[k].m_x*weight(x-i,y-j);
        ub[j*(m_resolutionX+1)+i] += weight(x-i,y-j);
        // Calculate the top and bottom part of the equation 8 for v(i,j-1/2). 
        vt[j*m_resolutionX+i] += m_velocity[k].m_y*weight(x-i,y-j);
        vb[j*m_resolutionX+i] += weight(x-i,y-j);
        // Calculate the top and bottom part of the equation 8 for u(i+1/2,j). 
        ut[j*(m_resolutionX+1)+i+1] += m_velocity[k].m_x*weight(x-(i+1),y-j);
        ub[j*(m_resolutionX+1)+i+1] += weight(x-(i+1),y-j);
        // Calculate the top and bottom part of the equation 8 for v(i,j+1/2). 
        vt[(j+1)*m_resolutionX+i] += m_velocity[k].m_y*weight(x-i,y-(j+1));
        vb[(j+1)*m_resolutionX+i] += weight(x-i,y-(j+1));
    }

    // Divide the top by bottom to get the velocities on the grid.
    for(int i=0; i<m_velocityX.size(); ++i)
    {
        if(ub[i] != 0)
        {
            m_velocityX[i] += ut[i] / ub[i];
        }
    }
    for(int i=0; i<m_velocityY.size(); ++i)
    {
        if(vb[i] != 0)
        {
            m_velocityY[i] += vt[i] / vb[i];
        }
    }
    // Keep the current velocities to calculate the differences later.
    m_velocityX_old = m_velocityX;
    m_velocityY_old = m_velocityY;
}

// Calculate weight of the particles using the distance between the particle and the cell face (Equation 9).
float Fluid::weight(float _x, float _y)
{
    return hat(_x)*hat(_y);
}

// Hat function used for calculating the weight (Equation 10).
float Fluid::hat(float _r)
{
    if(_r >= 0.0f && _r < 1.0f)
    {
        return 1.0f-_r;
    }
    else if(_r >= -1.0f && _r < 0.0f)
    {
        return 1.0f+_r;
    }
    else
    {
        return 0;
    }
}

////////// Add External Forces //////////

void Fluid::addForce()
{
    // Add forces.
    for(auto &u : m_velocityX)
    {
        u += m_force.m_x * m_timestep;
    }
    for(auto &v : m_velocityY)
    {
        v += m_force.m_y * m_timestep;
    }

    // Set velocities on the solid faces to 0.
    for(int j=0; j<m_resolutionY; ++j)
    {
        for(int i=0; i<m_resolutionX; ++i)
        {
            if(m_boundary[j*m_resolutionX+i] == -1)
            {
                m_velocityX[j*(m_resolutionX+1)+i] = 0.0f;
                m_velocityX[j*(m_resolutionX+1)+i+1] = 0.0f;
                m_velocityY[j*m_resolutionX+i] = 0.0f;
                m_velocityY[(j+1)*m_resolutionX+i] = 0.0f;                  
            }
        }
    }
}

////////// Projection //////////

void Fluid::project()
{ 
    size_t num_fluidcells = 0;
    // vectors to store the indices of the fluid cells.
    std::vector<int> indices;
    indices.resize(m_resolutionX*m_resolutionY,-1);
    for(int j=0; j<m_resolutionY; ++j)
    {
        for(int i=0; i<m_resolutionX; ++i)
        {
            if(m_boundary[j*m_resolutionX+i] == 1)
            {
                indices[j*m_resolutionX+i] = num_fluidcells++;
            }
        }
    }

    // Solve Poisson problem Ap=b (Equation 20)
    Eigen::VectorXf p(num_fluidcells), b(num_fluidcells);
    Eigen::SparseMatrix<float> A(num_fluidcells,num_fluidcells);
    A.reserve(Eigen::VectorXi::Constant(num_fluidcells,5));

    // Calculate coefficients for each fluid cell.
    // Set -1 to fluid neighbours and 4-(number of non-solid neighbours) to the cell.
    for(int j=0; j<m_resolutionY; ++j)
    {
        for(int i=0; i<m_resolutionX; ++i)
        {
            if(m_boundary[j*m_resolutionX+i] == 1)
            {
                size_t num_nonSolidNeighbours = 4;
                if(m_boundary[j*m_resolutionX+i-1] == -1)
                {
                    num_nonSolidNeighbours--;
                }
                else if(m_boundary[j*m_resolutionX+i-1] == 1)
                {
                    A.insert(indices[j*m_resolutionX+i],indices[j*m_resolutionX+i-1]) = -1;
                }

                if(m_boundary[(j-1)*m_resolutionX+i] == -1)
                {
                    num_nonSolidNeighbours--;
                }
                else if(m_boundary[(j-1)*m_resolutionX+i] == 1)
                {
                    A.insert(indices[j*m_resolutionX+i],indices[(j-1)*m_resolutionX+i]) = -1;
                }

                if(m_boundary[j*m_resolutionX+i+1] == -1)
                {
                    num_nonSolidNeighbours--;
                }
                else if(m_boundary[j*m_resolutionX+i+1] == 1)
                {
                    A.insert(indices[j*m_resolutionX+i],indices[j*m_resolutionX+i+1]) = -1;
                }

                if(m_boundary[(j+1)*m_resolutionX+i] == -1)
                {
                    num_nonSolidNeighbours--;
                }
                else if(m_boundary[(j+1)*m_resolutionX+i] == 1)
                {
                    A.insert(indices[j*m_resolutionX+i],indices[(j+1)*m_resolutionX+i]) = -1;
                }
                A.insert(indices[j*m_resolutionX+i],indices[j*m_resolutionX+i]) = num_nonSolidNeighbours;

                // Calculate the divergence of the cell.
                b[indices[j*m_resolutionX+i]] = - (m_velocityX[j*(m_resolutionX+1)+i+1]-m_velocityX[j*(m_resolutionX+1)+i] + m_velocityY[(j+1)*m_resolutionX+i]-m_velocityY[j*m_resolutionX+i]);
            }
        }
    }

    // Compute the Poisson problem using Conjugate Gradient to get pressures.
    Eigen::ConjugateGradient<Eigen::SparseMatrix<float>> cg;
    cg.compute(A);
    p = cg.solve(b);

    // Assign the calculated pressures to the pressure field.
    std::fill(m_pressure.begin(), m_pressure.end(), 0.0f);
    size_t count = 0;
    for(int j=0; j<m_resolutionY; ++j)
    {
        for(int i=0; i<m_resolutionX; ++i)
        {
            if(m_boundary[j*m_resolutionX+i] == 1)
            {
                m_pressure[j*m_resolutionX+i] = p[indices[j*m_resolutionX+i]];
                ++count;
            }
        }
    }

    // Add values to m_pressureViz vector for pressure visualisation.
    m_pressureViz.clear();
    m_pressureViz.reserve(count);
    for(int j=0; j<m_resolutionY; ++j)
    {
        for(int i=0; i<m_resolutionX; ++i)
        {
            if(m_boundary[j*m_resolutionX+i] == 1)
            {
                m_pressureViz.push_back({static_cast<float>(i)/m_cell_res,static_cast<float>(j)/m_cell_res,std::clamp(m_pressure[j*m_resolutionX+i],0.0f,2.55f)});
            }
        }
    }    

    // Update the velocity with the calculated pressures.
    for(int j=0; j<m_resolutionY; ++j)
    {
        for(int i=1; i<m_resolutionX; ++i)
        {
            m_velocityX[j*(m_resolutionX+1)+i] -= m_pressure[j*m_resolutionX+i] - m_pressure[j*m_resolutionX+i-1];
        }
    }
    for(int j=1; j<m_resolutionY; ++j)
    {
        for(int i=0; i<m_resolutionX; ++i)
        {
            m_velocityY[j*m_resolutionX+i] -= m_pressure[j*m_resolutionX+i] - m_pressure[(j-1)*m_resolutionX+i];
        }
    }
}

void Fluid::calculateDifferences()
{
    if(!m_type)
    {
        ////////// FLIP //////////
        // Calculate the velocity changes on the grid.
        for(int i=0; i<m_velocityX.size(); ++i)
        {
            m_velocityX_diff[i] = m_velocityX[i] - m_velocityX_old[i];
        }
        for(int i=0; i<m_velocityY.size(); ++i)
        {
            m_velocityY_diff[i] = m_velocityY[i] - m_velocityY_old[i];
        }
    }
    else
    {
        ////////// PIC //////////
        // Assign the new velocities.
        for(int i=0; i<m_velocityX.size(); ++i)
        {
            m_velocityX_diff[i] = m_velocityX[i];
        }
        for(int i=0; i<m_velocityY.size(); ++i)
        {
            m_velocityY_diff[i] = m_velocityY[i];
        }    
    }

    // Set the velocity on the solid surfaces to 0.
    for(int j=0; j<m_resolutionY; ++j)
    {
        for(int i=0; i<m_resolutionX; ++i)
        {
            if(m_boundary[j*m_resolutionX+i] == -1)
            {
                m_velocityX_diff[j*(m_resolutionX+1)+i] = 0.0f;
                m_velocityX_diff[j*(m_resolutionX+1)+i+1] = 0.0f;
                m_velocityY_diff[j*m_resolutionX+i] = 0.0f;
                m_velocityY_diff[(j+1)*m_resolutionX+i] = 0.0f;                  
            }
        }
    }
}

////////// Grid to Particle //////////

// Update the particle velocities using the velocity fields on the grid.
void Fluid::gridToParticle()
{
    for(int k=0; k<m_numParticles; ++k)
    {
        if(!m_type)
        {
            ////////// FLIP //////////
            // Add the interpolated velocity changes to particle velocities.
            m_velocity[k].m_x += gridToParticleX(m_position[k]);
            m_velocity[k].m_y += gridToParticleY(m_position[k]);            
        }
        else
        {
            ////////// PIC //////////
            // Assign the interpolated velocities to particle velocities.
            m_velocity[k].m_x = gridToParticleX(m_position[k]);
            m_velocity[k].m_y = gridToParticleY(m_position[k]);
        }

    }    
}

// Interpolate the velocities from grid to particles.
float Fluid::gridToParticleX(const ngl::Vec3 &_p)
{
    // Find which cell the particle is located.
    float x = _p.m_x*m_cell_res;
    float y = _p.m_y*m_cell_res;
    size_t i = static_cast<size_t>(x);
    size_t j = static_cast<size_t>(y);

    // Interpolate the velocities according to the equation 21 in my report.
    return (m_velocityX_diff[j*(m_resolutionX+1)+i]*weight(x-i,y-j) + m_velocityX_diff[j*(m_resolutionX+1)+(i+1)]*weight(x-(i+1),y-j))
                    / (weight(x-i,y-j) + weight(x-(i+1),y-j));
}

float Fluid::gridToParticleY(const ngl::Vec3 &_p)
{
    // Find which cell the particle is located.
    float x = _p.m_x*m_cell_res;
    float y = _p.m_y*m_cell_res;
    size_t i = static_cast<size_t>(x);
    size_t j = static_cast<size_t>(y);

    // Interpolate the velocities according to the equation 21 in my report.
    return (m_velocityY_diff[j*m_resolutionX+i]*weight(x-i,y-j) + m_velocityY_diff[(j+1)*m_resolutionX+i]*weight(x-i,y-(j+1)))
                    / (weight(x-i,y-j) + weight(x-i,y-(j+1)));
}

////////// Advection //////////

void Fluid::advect()
{
    // Reset the boundary field.
    std::replace (m_boundary.begin(), m_boundary.end(), 1, 0);
    for(int k=0; k<m_numParticles; ++k)
    {
        // Update the particles positions using Euler forward method.
        m_position[k].m_x += m_timestep * m_velocity[k].m_x;
        m_position[k].m_y += m_timestep * m_velocity[k].m_y;

        // Push the partilces outside the frame to inside.
        m_position[k].m_x = std::max( (1.0f+0.3f/m_cell_res), std::min( m_sizeX-(1.0f+0.3f/m_cell_res), m_position[k].m_x ) );
        m_position[k].m_y = std::max( (1.0f+0.3f/m_cell_res), std::min( m_sizeY-(1.0f+0.3f/m_cell_res), m_position[k].m_y ) );

        // Mark the cell as a fluid cell.
        size_t i = static_cast<size_t>(m_position[k].m_x*m_cell_res);
        size_t j = static_cast<size_t>(m_position[k].m_y*m_cell_res);
        m_boundary[j*m_resolutionX+i] = 1;
    }
}

////////// Render //////////

void Fluid::render(size_t _w, size_t _h, bool _flu, bool _vel, bool _pre)
{
    const auto ColourShader = "ColourShader";
    const auto GridViz = "GridViz";
    const auto PressureViz = "PressureViz";
    const auto SolidShader = "SolidShader";

    float width = static_cast<float>(m_sizeX);
    float height = static_cast<float>(m_sizeY);

    // Use orthographic projection for 2D simulation.
    auto view = ngl::lookAt({width/2.0f,height/2.0f,10}, {width/2.0f,height/2.0f,0}, {0,1,0});
    auto project = ngl::ortho(-width/2.0f, width/2.0f, -height/2.0f, height/2.0f, 0.1f, 50.0f);

    // Visualise solid cells.
    ngl::ShaderLib::use(SolidShader);
    ngl::ShaderLib::setUniform("MVP",project*view);
    ngl::ShaderLib::setUniform("size",ngl::Vec2(m_resolutionX,m_resolutionY));

    m_vao->bind();
        m_vao->setData(0,ngl::MultiBufferVAO::VertexData(m_solid.size()*sizeof(ngl::Vec3),m_solid[0].m_x));
        m_vao->setVertexAttributePointer(0,3,GL_FLOAT,0,0);
        m_vao->setNumIndices(m_solid.size());
        m_vao->draw();
    m_vao->unbind();   
    
    // Visualise fluid particles when fluid flag is set.
    if(_flu)
    {
        ngl::ShaderLib::use(ColourShader);
        ngl::ShaderLib::setUniform("MVP",project*view);
        glPointSize(3);

        m_vao->bind();
            m_vao->setData(0,ngl::MultiBufferVAO::VertexData(m_numParticles*sizeof(ngl::Vec3),m_position[0].m_x));
            m_vao->setVertexAttributePointer(0,3,GL_FLOAT,0,0);          
            m_vao->setNumIndices(m_numParticles);
            m_vao->draw();
        m_vao->unbind();
    }

    // Visualise velocity vectors when velocity flag is set.
    if(_vel)
    {
        ngl::ShaderLib::use(GridViz);
        ngl::ShaderLib::setUniform("thickness",2.0f);
        ngl::ShaderLib::setUniform("thickness2",1.0f);
        ngl::ShaderLib::setUniform("viewportSize",ngl::Vec2(_w,_h));
        ngl::ShaderLib::setUniform("MVP",project*view);

        m_vao->bind();
            m_vao->setData(0,ngl::MultiBufferVAO::VertexData(m_numParticles*sizeof(ngl::Vec3),m_position[0].m_x));
            m_vao->setVertexAttributePointer(0,3,GL_FLOAT,0,0);
            m_vao->setData(1,ngl::MultiBufferVAO::VertexData(m_numParticles*sizeof(ngl::Vec3),m_velocity[0].m_x));
            m_vao->setVertexAttributePointer(1,3,GL_FLOAT,0,0);              
            m_vao->setNumIndices(m_numParticles);
            m_vao->draw();
        m_vao->unbind();        
    }

    // Visualise pressure when pressure flag is set.
    if(_pre)
    {
        ngl::ShaderLib::use(PressureViz);
        glPointSize(10);
        ngl::ShaderLib::setUniform("MVP",project*view);

        m_vao->bind();
            m_vao->setData(0,ngl::MultiBufferVAO::VertexData(m_pressureViz.size()*sizeof(ngl::Vec3),m_pressureViz[0].m_x));
            m_vao->setVertexAttributePointer(0,3,GL_FLOAT,0,0);
            m_vao->setNumIndices(m_pressureViz.size());
            m_vao->draw();
        m_vao->unbind();        
    }
}