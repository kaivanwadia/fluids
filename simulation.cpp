#include "simulation.h"
#include <QGLWidget>
#include <QPainter>
#include "simparameters.h"
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Geometry>

const double PI = 3.1415926535898;

using namespace Eigen;
using namespace std;

Simulation::Simulation(const SimParameters &params) : params_(params), time_(0)
{
    fluid_ = new Fluid();


}

void Simulation::render()
{
    double baseradius = 0.02;
    double pulsefactor = 0.1;
    double pulsespeed = 50.0;
    double baselinewidth = 0.5;
    int numcirclewedges = 20;

    renderLock_.lock();
    {
        fluid_->render();
        // Mass Spring System Render
        // Rendering Rods
        for(vector<Rod>::iterator it = rods_.begin(); it != rods_.end(); ++it)
        {
            glColor3f(0.75, 0.0, 0.75);
            Vector2d sourcepos = particles_[it->p1].pos;
            Vector2d destpos   = particles_[it->p2].pos;

            glLineWidth(8);

            if (it->isHinge)
            {
                glColor3f(1.0, 1.0, 0.0);
                glLineWidth(4);
            }

            glBegin(GL_LINES);
            glVertex2f(sourcepos[0], sourcepos[1]);
            glVertex2f(destpos[0], destpos[1]);
            glEnd();
        }
        //Rendering Springs
        for(vector<Spring>::iterator it = springs_.begin(); it != springs_.end(); ++it)
        {
            glColor3f(0.0, 0.0, 1.0);
            Vector2d sourcepos = particles_[it->p1].pos;
            Vector2d destpos   = particles_[it->p2].pos;

            double dist = (sourcepos-destpos).norm();

            glLineWidth(baselinewidth/dist);

            if (it->unsnappable)
            {
                glColor3f(0.752941, 0.752941, 0.752941);
                glLineWidth(14);
            }

            glBegin(GL_LINES);
            glVertex2f(sourcepos[0], sourcepos[1]);
            glVertex2f(destpos[0], destpos[1]);
            glEnd();
        }
        //Rendering Particles
        for(vector<Particle>::iterator it = particles_.begin(); it != particles_.end(); ++it)
        {
            double radius = baseradius*sqrt(it->mass);
            double pulse = pulsefactor*sin(pulsespeed*time_);
            radius *= (1.0 + pulse);

            if(it->fixed)
            {
                radius = baseradius;
                glColor3f(1.0,0,0);
            }
            int swap = 1;
            glBegin(GL_TRIANGLE_FAN);
            {
                glVertex2f(it->pos[0], it->pos[1]);
                for(int i=0; i<=numcirclewedges; i++)
                {
                    if (i%5==0 && !it->fixed)
                    {
                        swap =swap*-1;
                        if(swap>0)
                        {
                            glColor3f(0.9+0.2*pulse, 0.4+0.2*pulse, 0);

                        }
                        else
                        {
                            glColor3f(0, 0.9+0.2*pulse, 0);

                        }
                    }
                    glVertex2f(it->pos[0] + radius * cos(2*PI*i/numcirclewedges),
                               it->pos[1] + radius * sin(2*PI*i/numcirclewedges));
                }
            }
            glEnd();
        }
    }
    renderLock_.unlock();
}

void Simulation::takeSimulationStep()
{
    fluidSimulationStep();
    massSpringSimulationStep();

    time_ += params_.timeStep;
}

void Simulation::fluidSimulationStep()
{
    addSource(fluid_->vx, fluid_->vxOld);
    addSource(fluid_->vy, fluid_->vyOld);
    // Velocity Code

    this->swap(fluid_->vx, fluid_->vxOld);
    this->swap(fluid_->vy, fluid_->vyOld);
//    cout<<"Velocity X : \n"<<fluid_->vx<<endl;
//    cout<<"Velocity Y : \n"<<fluid_->vy<<endl;

    advection(1, fluid_->vx, fluid_->vxOld, fluid_->vxOld, fluid_->vyOld );
    advection(2, fluid_->vy, fluid_->vyOld, fluid_->vxOld, fluid_->vyOld );

    project(fluid_->vx, fluid_->vy, fluid_->vxOld, fluid_->vyOld);

    this->swap(fluid_->vx, fluid_->vxOld);
    this->swap(fluid_->vy, fluid_->vyOld);

    diffuse(1, fluid_->vx, fluid_->vxOld, params_.viscosityFluid);
    diffuse(2, fluid_->vy, fluid_->vyOld, params_.viscosityFluid);

    project(fluid_->vx, fluid_->vy, fluid_->vxOld, fluid_->vyOld);

    //Density Code
    addSource(fluid_->fluidDensity, fluid_->fluidDensityOld);
    this->swap(fluid_->fluidDensity, fluid_->fluidDensityOld);


    diffuse(0, fluid_->fluidDensity, fluid_->fluidDensityOld, params_.diffusionConstant);
    this->swap(fluid_->fluidDensity, fluid_->fluidDensityOld);

    advection(0, fluid_->fluidDensity, fluid_->fluidDensityOld, fluid_->vx, fluid_->vy);
    fluid_->fluidDensityOld.setZero();
    fluid_->vxOld.setZero();
    fluid_->vyOld.setZero();

    cout << "fluid density: " << fluid_->getTotalDensity() << endl;
}

void Simulation::addSource(MatrixXd &d, MatrixXd &dOld)
{

    d += params_.timeStep * dOld;

}

void Simulation::diffuse(int boundary, MatrixXd &d, MatrixXd &dOld, double diffFactor)
{
    double a = params_.timeStep * diffFactor * fluid_->n * fluid_->n;

    for(int k = 0; k < 20; k++)
    {
        for(int i  =1; i <= fluid_->n; i++)
        {
            for(int j = 1; j <= fluid_->n; j++)
            {
                d.coeffRef(i,j) = (dOld.coeff(i,j) + a * (d.coeff(i-1,j) + d.coeff(i+1,j) +
                                                          d.coeff(i,j-1) + d.coeff(i,j+1)))/(1+4*a);
            }
        }
        setBoundry(boundary, d);
    }
}

void Simulation::project(Eigen::MatrixXd &x, Eigen::MatrixXd &y, Eigen::MatrixXd &xOld, Eigen::MatrixXd &yOld)
{
    for(int i = 1; i <= fluid_->n; i++)
    {
        for(int j = 1; j <= fluid_->n; j++)
        {
            yOld.coeffRef(i,j) = (x.coeff(i+1, j) - x.coeff(i-1, j)
                                    + y.coeff(i,j+1) - y.coeff(i,j-1)) * -0.5f / fluid_->n;

        }
    }

    xOld.setZero();

    setBoundry(0,yOld);
    setBoundry(0,xOld);

    for(int k = 0; k < 20; k++)
    {
        for(int i  =1; i <= fluid_->n; i++)
        {
            for(int j = 1; j <= fluid_->n; j++)
            {
                xOld.coeffRef(i,j) =  ( xOld.coeff(i-1, j) + xOld.coeff(i+1, j)
                                         + xOld.coeff(i,j-1) + xOld.coeff(i,j+1)
                                         + yOld.coeff(i,j)) / 4 ;
            }
        }
        setBoundry(0, xOld);
    }

    for(int i = 1; i <= fluid_->n; i++)
    {
        for(int j = 1; j <= fluid_->n; j++)
        {
            x.coeffRef(i,j) -= 0.5f * fluid_->n * (xOld.coeff(i+1, j) - xOld.coeff(i-1, j));
            y.coeffRef(i,j) -= 0.5f * fluid_->n * (xOld.coeff(i, j+1) - xOld.coeff(i, j-1));
        }
    }
    setBoundry(1,x);
    setBoundry(2,y);
}

void Simulation::setBoundry(int b, MatrixXd& m)
{
    for(int i = 1; i <= fluid_->n; i++)
    {
        m.coeffRef(0,i) =               (b==1) ? -m.coeff(1,i) : m.coeff(1,i);
        m.coeffRef(fluid_->n, i) =    (b==1) ? -m.coeff(fluid_->n-1,i) : m.coeff(fluid_->n-1,i);
        m.coeffRef(i,0) =               (b==2) ? -m.coeff(i,1) : m.coeff(i,1);
        m.coeffRef(i, fluid_->n) =    (b==2) ? -m.coeff(i, fluid_->n-1) : m.coeff(i,fluid_->n-1);
    }
    m.coeffRef(0,0) =                   0.5f *(m.coeff(1,0) + m.coeff(0,1));
    m.coeffRef(0,fluid_->n) =           0.5f *(m.coeff(1,fluid_->n) + m.coeff(0, fluid_->n-1));
    m.coeffRef(fluid_->n,0) =           0.5f *(m.coeff(fluid_->n-1,0) + m.coeff(fluid_->n,1));
    m.coeffRef(fluid_->n,fluid_->n) =   0.5f *(m.coeff(fluid_->n-1,fluid_->n) + m.coeff(fluid_->n,fluid_->n-1));

}

void Simulation::swap(MatrixXd &left, MatrixXd &right)
{
    MatrixXd tmp(left);
    left = right;
    right = tmp;
}

void Simulation::linearSolver(int b, Eigen::MatrixXd &x, Eigen::MatrixXd &xOld, double a, double c)
{
    for(int k = 0; k < 20; k++)
    {
        for(int i  =1; i <= fluid_->n; i++)
        {
            for(int j = 1; j <= fluid_->n; j++)
            {
                x.coeffRef(i,j) = (a * ( x.coeff(i-1, j) + x.coeff(i+1, j)
                                         + x.coeff(i,j-1) + x.coeff(i,j+1))
                                         + xOld.coeff(i,j)) / c ;
            }
        }
        setBoundry(b, x);
    }
}

void Simulation::advection(int boundry, MatrixXd &d, MatrixXd &dOld, MatrixXd &xOld, MatrixXd &yOld)
{
    double dt0 = params_.timeStep * fluid_->n;
    int n = fluid_->n;

    for(int i = 1; i <= n; i++)
    {
        for(int j = 1; j <= n; j++)
        {
            double x = i - dt0 * xOld.coeff(i,j);
            double y = j - dt0 * yOld.coeff(i,j);

            if(x > n + 0.5f)
            {
                x = n + 0.5f;
            }
            if(x < 0.5f)
            {
                x = 0.5f;
            }

            int i0 = (int) x;
            int i1 = i0 + 1;

            if(y > n + 0.5f)
            {
                y = n + 0.5f;
            }
            if(y < 0.5f)
            {
                y = 0.5f;
            }

            int j0 = (int) y;
            int j1 = j0 +1;

            double s1 = x - i0;
            double s0 = 1 - s1;
            double t1 = y - j0;
            double t0 = 1 - t1;

            d.coeffRef(i,j) = s0 * (t0 * dOld.coeff(i0, j0) + t1 * dOld.coeff(i0,j1))
                    + s1 * (t0 * dOld.coeff(i1,j0) + t1 * dOld.coeff(i1, j1));
        }
    }
    setBoundry(boundry, d);
}

void Simulation::addDensity(double x, double y)
{
    int i = floor((x+1)/fluid_->sizeOfVoxel);
    int j = -1* floor((y-1)/fluid_->sizeOfVoxel);

    if(i >= 0 && i <fluid_->n && j >= 0 && j < fluid_->n)
    {
        fluid_->fluidDensityOld.coeffRef(i,j) += params_.densityMagnitude;
    }
    for(int r = 1; r <= params_.densityRadius; r++)
    {
        if(i-r >= 0 && i+r <fluid_->n && j-r >= 0 && j+r < fluid_->n)
        {
            fluid_->fluidDensityOld.coeffRef(i-r,j-r) += params_.densityMagnitude;
            fluid_->fluidDensityOld.coeffRef(i-r,j) += params_.densityMagnitude;
            fluid_->fluidDensityOld.coeffRef(i-r,j+r) += params_.densityMagnitude;
            fluid_->fluidDensityOld.coeffRef(i,j+r) += params_.densityMagnitude;
            fluid_->fluidDensityOld.coeffRef(i+r,j+r) += params_.densityMagnitude;
            fluid_->fluidDensityOld.coeffRef(i+r,j) += params_.densityMagnitude;
            fluid_->fluidDensityOld.coeffRef(i+r,j-r) += params_.densityMagnitude;
            fluid_->fluidDensityOld.coeffRef(i,j-r) += params_.densityMagnitude;
        }
    }
}

void Simulation::addVelocity(double x, double y, double velX, double velY)
{
    int i = floor((x+1)/fluid_->sizeOfVoxel);
    int j = -1* floor((y-1)/fluid_->sizeOfVoxel);
    if(i >= 0 && i <fluid_->n && j >= 0 && j < fluid_->n)
    {
        fluid_->vxOld.coeffRef(i,j) += velX;
        fluid_->vyOld.coeffRef(i,j) += velY;
    }
    for(int r = 1; r <= params_.velocityRadius; r++)
    {
        if(i-r >= 0 && i+r <fluid_->n && j-r >= 0 && j+r < fluid_->n)
        {
//            cout<<"Here : "<<velX<<" : "<<velY<<endl;
            fluid_->vxOld.coeffRef(i-r,j-r) += velX;
            fluid_->vxOld.coeffRef(i-r,j) += velX;
            fluid_->vxOld.coeffRef(i-r,j+r) += velX;
            fluid_->vxOld.coeffRef(i,j+r) += velX;
            fluid_->vxOld.coeffRef(i+r,j+r) += velX;
            fluid_->vxOld.coeffRef(i+r,j) += velX;
            fluid_->vxOld.coeffRef(i+r,j-r) += velX;
            fluid_->vxOld.coeffRef(i,j-r) += velX;

            fluid_->vyOld.coeffRef(i,j) += velY;
            fluid_->vyOld.coeffRef(i-r,j-r) += velY;
            fluid_->vyOld.coeffRef(i-r,j) += velY;
            fluid_->vyOld.coeffRef(i-r,j+r) += velY;
            fluid_->vyOld.coeffRef(i,j+r) += velY;
            fluid_->vyOld.coeffRef(i+r,j+r) += velY;
            fluid_->vyOld.coeffRef(i+r,j) += velY;
            fluid_->vyOld.coeffRef(i+r,j-r) += velY;
            fluid_->vyOld.coeffRef(i,j-r) += velY;
        }
    }
}

void Simulation::clearScene()
{
    renderLock_.lock();
    {
        fluid_->zeroEverything();
        particles_.clear();
        springs_.clear();
        flexibleRodHinges_.clear();
        ropeHinges_.clear();
        rods_.clear();
    }
    renderLock_.unlock();
}

void Simulation::printSimParameters() {
    cout<<"\n\nSIM PARAMETERS:"<<endl;
    cout<<"Sim running:"<<params_.simRunning<<endl;
    cout<<"Timestep:"<<params_.timeStep<<endl;
    cout<<"Newton tol:"<<params_.NewtonTolerance<<endl;
    cout<<"Newton Iters:"<<params_.NewtonMaxIters<<endl;
    cout<<"ClickMode: "<<params_.clickMode<<endl;
    cout<<"densityMagnitude:"<<params_.densityMagnitude<<endl;
    cout<<"densityRadius:"<<params_.densityRadius<<endl;
    cout<<"velocityMagnitude:"<<params_.velocityMagnitude<<endl;
    cout<<"velocityRadius:"<<params_.velocityRadius<<endl;
    cout<<"diffusionConstant:"<<params_.diffusionConstant<<endl;
    cout<<"viscosityConstant:"<<params_.viscosityFluid<<endl;

    cout<<"Connector: "<<params_.connector<<endl;

    cout<<"springStiffness:"<<params_.springStiffness<<endl;
    cout<<"maxSpringStrain: "<<params_.maxSpringStrain<<endl;
    cout<<"Damping stiffness:"<<params_.dampingStiffness<<endl;

    cout<<"particleMass:"<<params_.particleMass<<endl;
    cout<<"particleFixed: "<<params_.particleFixed<<endl;
    cout<<"maxSpringDist:"<<params_.maxSpringDist<<endl;

    cout<<"rodDensity:"<<params_.rodDensity<<endl;
    cout<<"rodStretchStiffness:"<<params_.rodStretchStiffness<<endl;
    cout<<"rodBendingStiffness:"<<params_.rodBendingStiffness<<endl;
    cout<<"rodSegments:"<<params_.rodSegments<<endl;

    cout<<"ropeDensity:"<<params_.ropeDensity<<endl;
    cout<<"ropeBend:"<<params_.ropeBend<<endl;
    cout<<"ropeSegments:"<<params_.ropeSegments<<endl;
}

// Particle Stuff
void Simulation::massSpringSimulationStep()
{
    VectorXd q, qprev, v;
    buildConfiguration(q, qprev, v);
    numericalIntegration(q, qprev, v);
    unbuildConfiguration(q, v);
    pruneOverstrainedSprings();
    removeOutsideParticles();
}

void Simulation::numericalIntegration(VectorXd &q, VectorXd &qprev, VectorXd &v)
{
    VectorXd F;
    SparseMatrix<double> H;
    SparseMatrix<double> Minv;

    computeMassInverse(Minv);

    VectorXd oldq = q;
    q += params_.timeStep*v;
    computeForceAndHessian(q, oldq, F, H, v);
    v += params_.timeStep*Minv*F;
    computeStepProject(q, oldq, v);
}

void Simulation::computeLagrangeMultipliers(const VectorXd &qVV, const VectorXd &F, VectorXd &v)
{
    VectorXd lamGuess(rods_.size());
    lamGuess.setZero();
    SparseMatrix<double> massInverseMatrix(qVV.rows(), qVV.rows());
    massInverseMatrix.setZero();
    computeMassInverse(massInverseMatrix);
    VectorXd fOfx(rods_.size());
    SparseMatrix<double> gradF(rods_.size(), rods_.size());

    VectorXd c = qVV + params_.timeStep*v + params_.timeStep * params_.timeStep * massInverseMatrix * F;
    SparseMatrix<double> gradGTransposeofqVV = computeGradGTranspose(qVV);
    for(int newtonIt =0 ; newtonIt< params_.NewtonMaxIters; newtonIt++)
    {
        VectorXd qInside(qVV.rows());
        qInside.setZero();
        fOfx.setZero();
        gradF.setZero();
        // Compute F of Lambda(i+1)
        qInside = c + (params_.timeStep * params_.timeStep) * massInverseMatrix * gradGTransposeofqVV * lamGuess;

        for(int i=0; i<rods_.size(); i++)
        {
            int p1pos = rods_[i].p1*2;
            int p2pos = rods_[i].p2*2;
            Vector2d p1 = qInside.segment<2>(p1pos);
            Vector2d p2 = qInside.segment<2>(p2pos);
            fOfx[i] += (p2-p1).squaredNorm() - rods_[i].restlen*rods_[i].restlen;
        }
        if(fOfx.norm() < params_.NewtonTolerance)
        {
            break;
        }
        //Compute Gradient of F(lambda i + 1)
        gradF = computeGradGTranspose(qInside).transpose() * params_.timeStep * params_.timeStep * massInverseMatrix * gradGTransposeofqVV;
        gradF.makeCompressed();
        SparseQR<SparseMatrix<double>, COLAMDOrdering<int> > solver;
        solver.compute(gradF);
        VectorXd deltaLamda = solver.solve(-fOfx);
        lamGuess += deltaLamda;
    }
    v = v + params_.timeStep * massInverseMatrix * F + params_.timeStep * massInverseMatrix * gradGTransposeofqVV * lamGuess;
}

Eigen::SparseMatrix<double> Simulation::computeGradGTranspose(const VectorXd &q)
{
    SparseMatrix<double> gradGTranspose(q.rows(), rods_.size());
    gradGTranspose.setZero();

    for(int i=0; i<rods_.size(); i++)
    {
        int p1index = rods_[i].p1*2;
        int p2index = rods_[i].p2*2;
        Vector2d p1 = q.segment<2>(p1index);
        Vector2d p2 = q.segment<2>(p2index);

        Vector2d gradient1 = (p2 - p1)*-2;
        Vector2d gradient2 = (p1 - p2)*-2;

        gradGTranspose.coeffRef(p1index, i) += gradient1[0];
        gradGTranspose.coeffRef(p1index+1, i) += gradient1[1];

        gradGTranspose.coeffRef(p2index, i) += gradient2[0];
        gradGTranspose.coeffRef(p2index+1, i) += gradient2[1];
    }
    return gradGTranspose;
}

void Simulation::computeStepProject(VectorXd &q, VectorXd &oldq, VectorXd &v)
{
    VectorXd qGuess = q;
    VectorXd lamGuess(rods_.size());
    lamGuess.setZero();
    SparseMatrix<double> massInverseMatrix(q.rows(), q.rows());
    massInverseMatrix.setZero();
    computeMassInverse(massInverseMatrix);
    VectorXd forceX(qGuess.rows() + rods_.size());

    for(int i = 0; i < params_.NewtonMaxIters; i++)
    {
        forceX.setZero();
        for(int j = 0; j < rods_.size(); j++)
        {
            int p1pos = rods_[j].p1*2;
            int p2pos = rods_[j].p2*2;
            Vector2d p1 = qGuess.segment<2>(2*rods_[j].p1);
            Vector2d p2 = qGuess.segment<2>(2*rods_[j].p2);
            forceX.segment<2>(rods_[j].p1*2) += lamGuess[j] * massInverseMatrix.coeff(p1pos, p1pos) * (p2 - p1) * 2;
            forceX.segment<2>(rods_[j].p2*2) += lamGuess[j] * massInverseMatrix.coeff(p2pos, p2pos) * (p1 - p2) * 2;
        }

        VectorXd qDifference(qGuess.rows());
        qDifference.setZero();
        qDifference = q - qGuess;
        for(int j =0; j < qDifference.rows(); j++)
        {
            forceX[j] += qDifference[j];
        }

        for(int j = 0; j < rods_.size(); j++)
        {
            Vector2d dstpos(qGuess[rods_[j].p1*2], qGuess[rods_[j].p1*2+1]);
            Vector2d srcpos(qGuess[rods_[j].p2*2], qGuess[rods_[j].p2*2+1]);
            double dist = (dstpos-srcpos).norm();
            forceX[j+qGuess.rows()] += dist*dist - rods_[j].restlen*rods_[j].restlen;
        }
        if(forceX.norm() < params_.NewtonTolerance)
        {
            break;
        }

        // GRADIENT CALCULATION

        SparseMatrix<double> forceGradient(qGuess.rows()+rods_.size(), qGuess.rows()+rods_.size());
        forceGradient.setZero();

        // Calculating Top Left of Force Gradient Matrix
        vector< Triplet<double> > topLeftTriplets;
        topLeftTriplets.reserve(rods_.size()*8);
        double dfxidxi, dfxidxj, dfyidyi, dfyidyj, dfxjdxi, dfxjdxj, dfyjdyi, dfyjdyj;
        for(int j=0; j<rods_.size(); j++)
        {
            int p1pos = rods_[j].p1*2;
            int p2pos = rods_[j].p2*2;
            dfxidxi = -2*lamGuess[j]*massInverseMatrix.coeff(p1pos, p1pos);
            dfxidxj = -dfxidxi;
            dfyidyi = dfxidxi;
            dfyidyj = -dfyidyi;

            dfxjdxi = 2*lamGuess[j]*massInverseMatrix.coeff(p2pos, p2pos);
            dfxjdxj = -dfxjdxi;
            dfyjdyi = dfxjdxi;
            dfyjdyj = -dfyjdyi;
            topLeftTriplets.push_back(Triplet<double>(p1pos, p1pos, dfxidxi));
            topLeftTriplets.push_back(Triplet<double>(p1pos, p2pos, dfxidxj));
            topLeftTriplets.push_back(Triplet<double>(p1pos + 1, p1pos + 1, dfyidyi));
            topLeftTriplets.push_back(Triplet<double>(p1pos + 1, p2pos + 1, dfyidyj));

            topLeftTriplets.push_back(Triplet<double>(p2pos, p1pos, dfxjdxi));
            topLeftTriplets.push_back(Triplet<double>(p2pos, p2pos, dfxjdxj));
            topLeftTriplets.push_back(Triplet<double>(p2pos + 1, p1pos + 1, dfyjdyi));
            topLeftTriplets.push_back(Triplet<double>(p2pos + 1, p2pos + 1, dfyjdyj));
        }
        SparseMatrix<double> topLeftMatrix(qGuess.rows(), qGuess.rows());
        topLeftMatrix.setZero();
        topLeftMatrix.setFromTriplets(topLeftTriplets.begin(), topLeftTriplets.end());
        SparseMatrix<double> identity(qGuess.rows(), qGuess.rows());
        identity.setIdentity();
        topLeftMatrix += identity;
        for(int j = 0; j < topLeftMatrix.rows(); j++)
        {
            for(int k = 0; k < topLeftMatrix.cols(); k++)
            {
                forceGradient.coeffRef(j,k) = topLeftMatrix.coeff(j,k);
            }
        }
        // Calculating Top Right and Bottom Left of Force Gradient Matrix
        for(int j = 0; j < rods_.size(); j++)
        {
            int n = j + qGuess.rows();
            int p1pos = rods_[j].p1 *2;
            int p2pos = rods_[j].p2 *2;
            Vector2d p1Grad = (qGuess.segment<2>(p2pos) - qGuess.segment<2>(p1pos)) * 2;
            Vector2d p2Grad = (qGuess.segment<2>(p1pos) - qGuess.segment<2>(p2pos)) * 2;
            forceGradient.coeffRef(n,p1pos) += p1Grad[0];
            forceGradient.coeffRef(n,p1pos + 1) += p1Grad[1];
            forceGradient.coeffRef(n,p2pos) += p2Grad[0];
            forceGradient.coeffRef(n,p2pos + 1) += p2Grad[1];
            forceGradient.coeffRef(p1pos,n) += p1Grad[0] * massInverseMatrix.coeff(p1pos,p1pos);
            forceGradient.coeffRef(p1pos + 1,n) += p1Grad[1] * massInverseMatrix.coeff(p1pos,p1pos);
            forceGradient.coeffRef(p2pos,n) += p2Grad[0] * massInverseMatrix.coeff(p2pos,p2pos);
            forceGradient.coeffRef(p2pos + 1,n) += p2Grad[1] * massInverseMatrix.coeff(p2pos,p2pos);
        }
        forceGradient.makeCompressed();
        SparseQR<SparseMatrix<double>, COLAMDOrdering<int> > solver;
        solver.compute(forceGradient);
        VectorXd deltaguess = solver.solve(-forceX);
        for(int j = 0; j < qGuess.rows(); j++)
        {
            qGuess[j] -= deltaguess[j];
        }
        int n = qGuess.rows();
        for(int j = 0; j < lamGuess.rows(); j++)
        {
            lamGuess[j] += deltaguess[j+n];
        }
    }
    v = v + (qGuess - q)/params_.timeStep;
    q = qGuess;
}

void Simulation::buildConfiguration(VectorXd &q, VectorXd &qprev, VectorXd &v)
{
    int ndofs = 2*particles_.size();
    q.resize(ndofs);
    qprev.resize(ndofs);
    v.resize(ndofs);

    for(int i=0; i<(int)particles_.size(); i++)
    {
        q.segment<2>(2*i) = particles_[i].pos;
        qprev.segment<2>(2*i) = particles_[i].prevpos;
        v.segment<2>(2*i) = particles_[i].vel;
    }
}

void Simulation::unbuildConfiguration(const VectorXd &q, const VectorXd &v)
{
    int ndofs = q.size();
    assert(ndofs == int(2*particles_.size()));
    for(int i=0; i<ndofs/2; i++)
    {
        particles_[i].prevpos = particles_[i].pos;
        particles_[i].pos = q.segment<2>(2*i);
        particles_[i].vel = v.segment<2>(2*i);
    }
}

void Simulation::computeForceAndHessian(const VectorXd &q, const VectorXd &qprev, Eigen::VectorXd &F, SparseMatrix<double> &H, Eigen::VectorXd &v)
{
    F.resize(q.size());
    F.setZero();
    H.resize(q.size(), q.size());
    H.setZero();

    vector<Tr> Hcoeffs;
    processSpringForce(q, F, Hcoeffs);
    processDampingForce(q, qprev, F, Hcoeffs);
    processElasticBendingForce(q, F);
    processFluidForce(q, F, v);
    H.setFromTriplets(Hcoeffs.begin(), Hcoeffs.end());
}

void Simulation::processFluidForce(const VectorXd &q, VectorXd &F, const VectorXd &v)
{
    int nparticles = (int)particles_.size();
    for(int i=0; i<nparticles; i++){
        if(!particles_[i].fixed)
        {
            int fluid_i = floor((q[2*i]+1)/fluid_->sizeOfVoxel);
            int fluid_j = -1* floor((q[2*i+1]-1)/fluid_->sizeOfVoxel);
            Vector2d relVel(-(fluid_->vx.coeff(fluid_i, fluid_j) - v[2*i]), -((-fluid_->vy.coeff(fluid_i, fluid_j)) - v[2*i+1]));
            F.coeffRef(2*i) += -6.0*PI*params_.viscosityFluid*(relVel[0])*particles_[i].mass;
            F.coeffRef(2*i+1) += -6.0*PI*params_.viscosityFluid*(relVel[1])*particles_[i].mass;
        }
    }
}

void Simulation::processElasticBendingForce(const VectorXd &q, VectorXd &F)
{
    Vector3d pi,pj,pk,zUnit;
    zUnit.setZero();
    zUnit[2] = 1;
    int s1Id, s2Id, piIndex, pjIndex, pkIndex;

    for (int i=0; i<flexibleRodHinges_.size(); i++)
    {
        pi.setZero();
        pj.setZero();
        pk.setZero();
        s1Id = flexibleRodHinges_[i].s1;
        s2Id = flexibleRodHinges_[i].s2;

        if (springs_[s1Id].p1 == springs_[s2Id].p1)
        {
            pj.segment<2>(0) = q.segment<2>(springs_[s1Id].p1*2);
            pi.segment<2>(0) = q.segment<2>(springs_[s1Id].p2*2);
            pk.segment<2>(0) = q.segment<2>(springs_[s2Id].p2*2);
            pjIndex = springs_[s1Id].p1;
            piIndex = springs_[s1Id].p2;
            pkIndex = springs_[s2Id].p2;
        }
        else if (springs_[s1Id].p1 == springs_[s2Id].p2)
        {
            pj.segment<2>(0) = q.segment<2>(springs_[s1Id].p1*2);
            pi.segment<2>(0) = q.segment<2>(springs_[s1Id].p2*2);
            pk.segment<2>(0) = q.segment<2>(springs_[s2Id].p1*2);
            pjIndex = springs_[s1Id].p1;
            piIndex = springs_[s1Id].p2;
            pkIndex = springs_[s2Id].p1;
        }
        else if (springs_[s1Id].p2 == springs_[s2Id].p1)
        {
            pj.segment<2>(0) = q.segment<2>(springs_[s1Id].p2*2);
            pi.segment<2>(0) = q.segment<2>(springs_[s1Id].p1*2);
            pk.segment<2>(0) = q.segment<2>(springs_[s2Id].p2*2);
            pjIndex = springs_[s1Id].p2;
            piIndex = springs_[s1Id].p1;
            pkIndex = springs_[s2Id].p2;
        }
        else if (springs_[s1Id].p2 == springs_[s2Id].p2)
        {
            pj.segment<2>(0) = q.segment<2>(springs_[s1Id].p2*2);
            pi.segment<2>(0) = q.segment<2>(springs_[s1Id].p1*2);
            pk.segment<2>(0) = q.segment<2>(springs_[s2Id].p1*2);
            pjIndex = springs_[s1Id].p2;
            piIndex = springs_[s1Id].p1;
            pkIndex = springs_[s2Id].p1;
        }
        double y = ((pj - pi).cross(pk - pj)).dot(zUnit);
        double x = ((pj - pi).norm() * (pk - pj).norm()) + ((pj - pi).dot(pk - pj));
        double theta = 2 * atan2(y, x);
        Vector3d Fi = (flexibleRodHinges_[i].stiffness * theta * (pj - pi).cross(zUnit)) / (pj - pi).squaredNorm();
        Vector3d Fk = (flexibleRodHinges_[i].stiffness * theta * (pk - pj).cross(zUnit)) / (pk - pj).squaredNorm();
        F.segment<2>(piIndex * 2) += Fi.segment<2>(0);
        F.segment<2>(pkIndex * 2) += Fk.segment<2>(0);
        F.segment<2>(pjIndex * 2) += (-Fi-Fk).segment<2>(0);
    }
    for (int i=0; i<ropeHinges_.size(); i++)
    {
        pi.setZero();
        pj.setZero();
        pk.setZero();
        s1Id = ropeHinges_[i].s1;
        s2Id = ropeHinges_[i].s2;

        if (rods_[s1Id].p1 == rods_[s2Id].p1)
        {
            pj.segment<2>(0) = q.segment<2>(rods_[s1Id].p1*2);
            pi.segment<2>(0) = q.segment<2>(rods_[s1Id].p2*2);
            pk.segment<2>(0) = q.segment<2>(rods_[s2Id].p2*2);
            pjIndex = rods_[s1Id].p1;
            piIndex = rods_[s1Id].p2;
            pkIndex = rods_[s2Id].p2;
        }
        else if (rods_[s1Id].p1 == rods_[s2Id].p2)
        {
            pj.segment<2>(0) = q.segment<2>(rods_[s1Id].p1*2);
            pi.segment<2>(0) = q.segment<2>(rods_[s1Id].p2*2);
            pk.segment<2>(0) = q.segment<2>(rods_[s2Id].p1*2);
            pjIndex = rods_[s1Id].p1;
            piIndex = rods_[s1Id].p2;
            pkIndex = rods_[s2Id].p1;
        }
        else if (rods_[s1Id].p2 == rods_[s2Id].p1)
        {
            pj.segment<2>(0) = q.segment<2>(rods_[s1Id].p2*2);
            pi.segment<2>(0) = q.segment<2>(rods_[s1Id].p1*2);
            pk.segment<2>(0) = q.segment<2>(rods_[s2Id].p2*2);
            pjIndex = rods_[s1Id].p2;
            piIndex = rods_[s1Id].p1;
            pkIndex = rods_[s2Id].p2;
        }
        else if (rods_[s1Id].p2 == rods_[s2Id].p2)
        {
            pj.segment<2>(0) = q.segment<2>(rods_[s1Id].p2*2);
            pi.segment<2>(0) = q.segment<2>(rods_[s1Id].p1*2);
            pk.segment<2>(0) = q.segment<2>(rods_[s2Id].p1*2);
            pjIndex = rods_[s1Id].p2;
            piIndex = rods_[s1Id].p1;
            pkIndex = rods_[s2Id].p1;
        }
        double y = ((pj - pi).cross(pk - pj)).dot(zUnit);
        double x = ((pj - pi).norm() * (pk - pj).norm()) + ((pj - pi).dot(pk - pj));
        double theta = 2 * atan2(y, x);
        Vector3d Fi = (ropeHinges_[i].stiffness * theta * (pj - pi).cross(zUnit)) / (pj - pi).squaredNorm();
        Vector3d Fk = (ropeHinges_[i].stiffness * theta * (pk - pj).cross(zUnit)) / (pk - pj).squaredNorm();
        F.segment<2>(piIndex * 2) += Fi.segment<2>(0);
        F.segment<2>(pkIndex * 2) += Fk.segment<2>(0);
        F.segment<2>(pjIndex * 2) += (-Fi-Fk).segment<2>(0);
    }
}

void Simulation::processSpringForce(const VectorXd &q, VectorXd &F, std::vector<Tr> &H)
{
    int nsprings = (int)springs_.size();

    for(int i=0; i<nsprings; i++)
    {
        Vector2d p1 = q.segment<2>(2*springs_[i].p1);
        Vector2d p2 = q.segment<2>(2*springs_[i].p2);
        double dist = (p2-p1).norm();
        Vector2d localF = springs_[i].stiffness*(dist-springs_[i].restlen)/dist * (p2-p1);
        F.segment<2>(2*springs_[i].p1) += localF;
        F.segment<2>(2*springs_[i].p2) -= localF;

        Matrix2d I;
        I << 1, 0, 0, 1;
        Matrix2d localH = springs_[i].stiffness * (1.0 - springs_[i].restlen/dist)*I;
        localH += springs_[i].stiffness*springs_[i].restlen*(p2-p1)*(p2-p1).transpose()/dist/dist/dist;

        for(int j=0; j<2; j++)
            for(int k=0; k<2;k++)
            {
                H.push_back(Tr(2*springs_[i].p1+j, 2*springs_[i].p1+k, localH.coeff(j,k)));
                H.push_back(Tr(2*springs_[i].p2+j, 2*springs_[i].p2+k, localH.coeff(j,k)));
                H.push_back(Tr(2*springs_[i].p1+j, 2*springs_[i].p2+k, -localH.coeff(j,k)));
                H.push_back(Tr(2*springs_[i].p2+j, 2*springs_[i].p1+k, -localH.coeff(j,k)));
            }
    }
}

void Simulation::processDampingForce(const VectorXd &q, const VectorXd &qprev, VectorXd &F, std::vector<Tr> &H)
{
    int nsprings = (int)springs_.size();

    for(int i=0; i<nsprings; i++)
    {
        if (springs_[i].unsnappable)
        {
            continue;
        }
        Vector2d p1 = q.segment<2>(2*springs_[i].p1);
        Vector2d p2 = q.segment<2>(2*springs_[i].p2);
        Vector2d p1prev = qprev.segment<2>(2*springs_[i].p1);
        Vector2d p2prev = qprev.segment<2>(2*springs_[i].p2);

        Vector2d relvel = (p2 - p2prev)/params_.timeStep - (p1 - p1prev)/params_.timeStep;
        Vector2d localF = params_.dampingStiffness*relvel;
        F.segment<2>(2*springs_[i].p1) += localF;
        F.segment<2>(2*springs_[i].p2) -= localF;

        Matrix2d I;
        I << 1, 0, 0, 1;
        Matrix2d localH = params_.dampingStiffness*I/params_.timeStep;

        for(int j=0; j<2; j++)
            for(int k=0; k<2;k++)
            {
                H.push_back(Tr(2*springs_[i].p1+j, 2*springs_[i].p1+k, localH.coeff(j,k)));
                H.push_back(Tr(2*springs_[i].p2+j, 2*springs_[i].p2+k, localH.coeff(j,k)));
                H.push_back(Tr(2*springs_[i].p1+j, 2*springs_[i].p2+k, -localH.coeff(j,k)));
                H.push_back(Tr(2*springs_[i].p2+j, 2*springs_[i].p1+k, -localH.coeff(j,k)));
            }
    }
}

void Simulation::addParticle(double x, double y)
{
    addParticle(x, y, params_.particleFixed);
}

void Simulation::addParticle(double x, double y, bool particleFixed)
{
    renderLock_.lock();
    {
        Vector2d newParticlePos(x,y);
        double particleMass = params_.particleMass;
        int newParticleIndex = particles_.size();
        if (params_.connector == params_.CT_FLEXIBLE_ROD || params_.connector == params_.CT_ROPE)
        {
            particles_.push_back(Particle(newParticlePos, particleMass, particleFixed, false));
        }
        for(int i=0; i<(int)particles_.size(); i++)
        {
            if (i == newParticleIndex)
            {
                continue;
            }
            Vector2d pos = particles_[i].pos;
            double dist = (pos-newParticlePos).norm();
            if(!particles_[i].inert && dist <= params_.maxSpringDist)
            {
                if (params_.connector == params_.CT_RIGID_ROD)
                {
                    rods_.push_back(Rod(particles_.size(), i, dist));
                }
                else if (params_.connector == params_.CT_SPRING)
                {
                    springs_.push_back(Spring(particles_.size(), i, params_.springStiffness/dist, dist));
                }
                else if (params_.connector == params_.CT_FLEXIBLE_ROD)
                {
                    int rodSegs = params_.rodSegments;
                    if (rodSegs <= 1)
                    {
                        rodSegs = 2;
                    }
                    double springLength = dist/rodSegs;
                    double springMass = params_.rodDensity * springLength;
                    double springStiffness = params_.rodStretchStiffness/springLength;
                    double hingeStiffness = 0;
                    Vector2d unitVector = (newParticlePos - pos)/dist;
                    Vector2d distanceToMove = unitVector * (dist/rodSegs);
                    Vector2d newInertParticlePos = (distanceToMove * 1) + pos;
                    springs_.push_back(Spring(particles_.size(), i, springStiffness, springLength, springMass, true));
                    particles_.push_back(Particle(newInertParticlePos, springMass, false, true));
                    int j;
                    for (j=2; j<=rodSegs - 1; j++)
                    {
                        newInertParticlePos = (distanceToMove * j) + pos;
                        springs_.push_back(Spring(particles_.size(), particles_.size() - 1, springStiffness, springLength, springMass, true));
                        particles_.push_back(Particle(newInertParticlePos, springMass, false, true));
                        hingeStiffness = (params_.rodBendingStiffness * 2)/(springs_[springs_.size() - 1].restlen + springs_[springs_.size() - 2].restlen);
                        flexibleRodHinges_.push_back(FlexibleRodHinge(springs_.size() - 1, springs_.size() - 2, hingeStiffness));
                    }
                    springs_.push_back(Spring(newParticleIndex, particles_.size() - 1, springStiffness, springLength, springMass, true));
                    hingeStiffness = (params_.rodBendingStiffness * 2)/(springs_[springs_.size() - 1].restlen + springs_[springs_.size() - 2].restlen);
                    flexibleRodHinges_.push_back(FlexibleRodHinge(springs_.size() - 1, springs_.size() - 2, hingeStiffness));
                    particleMass += springMass/2;
                }
                else if (params_.connector == params_.CT_ROPE)
                {
                    int rodSegs = params_.ropeSegments;
                    if (rodSegs <= 1)
                    {
                        rodSegs = 2;
                    }
                    double rodLength = dist/rodSegs;
                    double rodMass = params_.rodDensity * rodLength;
                    double hingeStiffness = 0;
                    Vector2d unitVector = (newParticlePos - pos)/dist;
                    Vector2d distanceToMove = unitVector * (dist/rodSegs);
                    Vector2d newInertParticlePos = (distanceToMove * 1) + pos;
                    rods_.push_back(Rod(particles_.size(), i, rodLength, rodMass, true));
                    particles_.push_back(Particle(newInertParticlePos, rodMass, false, true));
                    int j;
                    for (j=2; j<=rodSegs - 1; j++)
                    {
                        newInertParticlePos = (distanceToMove * j) + pos;
                        rods_.push_back(Rod(particles_.size(), particles_.size() - 1, rodLength, rodMass, true));
                        particles_.push_back(Particle(newInertParticlePos, rodMass, false, true));
                        hingeStiffness = (params_.rodBendingStiffness * 2)/(rods_[rods_.size() - 1].restlen + rods_[rods_.size() - 2].restlen);
                        ropeHinges_.push_back(RopeHinge(rods_.size() - 1, rods_.size() - 2, hingeStiffness));
                    }
                    rods_.push_back(Rod(newParticleIndex, particles_.size() - 1, rodLength, rodMass, true));
                    hingeStiffness = (params_.rodBendingStiffness * 2)/(rods_[rods_.size() - 1].restlen + rods_[rods_.size() - 2].restlen);
                    ropeHinges_.push_back(RopeHinge(rods_.size() - 1, rods_.size() - 2, hingeStiffness));
                    particleMass += rodMass/2;
                }
            }
        }
        if(particleFixed)
        {
            particleMass = std::numeric_limits<double>::infinity();
        }
        if (params_.connector != params_.CT_FLEXIBLE_ROD && params_.connector != params_.CT_ROPE)
        {
            particles_.push_back(Particle(newParticlePos, particleMass, particleFixed, false));
        }
        else
        {
            particles_[newParticleIndex].mass = particleMass;
        }
    }
    renderLock_.unlock();
}

void Simulation::computeMassInverse(Eigen::SparseMatrix<double> &Minv)
{
    int ndofs = 2*int(particles_.size());

    Minv.resize(ndofs, ndofs);
    Minv.setZero();

    vector<Tr> Minvcoeffs;
    for(int i=0; i<ndofs/2; i++)
    {
        Minvcoeffs.push_back(Tr(2*i,   2*i,   1.0/particles_[i].mass));
        Minvcoeffs.push_back(Tr(2*i+1, 2*i+1, 1.0/particles_[i].mass));
    }

    Minv.setFromTriplets(Minvcoeffs.begin(), Minvcoeffs.end());
}

void Simulation::pruneOverstrainedSprings()
{
    set<int> springstodelete;

    vector<Spring> newsprings;
    vector<int> remainingspringmap;
    int nsprings = springs_.size();
    for(int i=0; i<nsprings; i++)
    {
        Vector2d p1 = particles_[springs_[i].p1].pos;
        Vector2d p2 = particles_[springs_[i].p2].pos;
        double d = (p2-p1).norm();

        double strain = (d - springs_[i].restlen)/springs_[i].restlen;
        if(!springs_[i].unsnappable && strain > params_.maxSpringStrain)
        {
            springstodelete.insert(i);
        }
    }

    renderLock_.lock();
    {
        if(!springstodelete.empty())
        {
            for(int i=0; i<(int)springs_.size(); i++)
            {
                if(springstodelete.count(i) == 0)
                {
                    remainingspringmap.push_back(newsprings.size());
                    newsprings.push_back(springs_[i]);
                }
                else
                {
                    remainingspringmap.push_back(-1);
                }
            }
        }
        if(!springstodelete.empty())
        {
            springs_ = newsprings;
            for(vector<FlexibleRodHinge>::iterator hinge = flexibleRodHinges_.begin(); hinge != flexibleRodHinges_.end(); ++hinge)
            {
                hinge->s1 = remainingspringmap[hinge->s1];
                hinge->s2 = remainingspringmap[hinge->s2];
            }
        }
    }
    renderLock_.unlock();
}

void Simulation::removeOutsideParticles()
{
    set<int> particlestodelete;
    set<int> springstodelete;
    set<int> springHingesToDelete;
    set<int> ropeHingesToDelete;
    set<int> rodsToDelete;

    vector<Particle> newparticles;
    vector<Spring> newsprings;
    vector<FlexibleRodHinge> newSpringHinges;
    vector<RopeHinge> newRopeHinges;
    vector<Rod> newrods;
    vector<int> remainingparticlemap;
    vector<int> remainingspringmap;
    vector<int> remainingrodmap;
    for(int i=0; i<(int)particles_.size(); i++)
    {
        Vector2d particlePos = particles_[i].pos;

        if(fabs(particlePos[0]) > 2 || fabs(particlePos[1]) > 2)
        {
            particlestodelete.insert(i);
            break;
        }
    }
    if(!particlestodelete.empty())
    {
        for(int i=0; i<(int)springs_.size(); i++)
        {
            if(particlestodelete.count(springs_[i].p1) || particlestodelete.count(springs_[i].p2))
            {
                springstodelete.insert(i);
                for(int j=0; j<(int)flexibleRodHinges_.size(); j++)
                {
                    if (springstodelete.count(flexibleRodHinges_[j].s1) || springstodelete.count(flexibleRodHinges_[j].s2))
                    {
                        springHingesToDelete.insert(j);
                    }
                }
            }

        }
        for(int i=0; i<(int)rods_.size(); i++)
        {
            if(particlestodelete.count(rods_[i].p1) || particlestodelete.count(rods_[i].p2))
            {
                rodsToDelete.insert(i);
                for(int j=0; j<(int)ropeHinges_.size(); j++)
                {
                    if (rodsToDelete.count(ropeHinges_[j].s1) || rodsToDelete.count(ropeHinges_[j].s2))
                    {
                        ropeHingesToDelete.insert(j);
                    }
                }
            }
        }
        for(int i=0; i<(int)particles_.size(); i++)
        {
            if(particlestodelete.count(i) == 0)
            {
                remainingparticlemap.push_back(newparticles.size());
                newparticles.push_back(particles_[i]);
            }
            else
            {
                remainingparticlemap.push_back(-1);
            }
        }
    }
    if(!springstodelete.empty())
    {
        for(int i=0; i<(int)springs_.size(); i++)
        {
            if(springstodelete.count(i) == 0)
            {
                remainingspringmap.push_back(newsprings.size());
                newsprings.push_back(springs_[i]);
            }
            else
            {
                remainingspringmap.push_back(-1);
            }
        }
    }
    if(!rodsToDelete.empty())
    {
        for(int i=0; i<(int)rods_.size(); i++)
        {
            if(rodsToDelete.count(i) == 0)
            {
                remainingrodmap.push_back(newrods.size());
                newrods.push_back(rods_[i]);
            }
            else
            {
                remainingrodmap.push_back(-1);
            }
        }
    }
    if(!springHingesToDelete.empty())
    {
        for(int i=0; i<(int)flexibleRodHinges_.size(); i++)
        {
            if(springHingesToDelete.count(i) == 0)
            {
                newSpringHinges.push_back(flexibleRodHinges_[i]);
            }
        }
    }
    if(!ropeHingesToDelete.empty())
    {
        for(int i=0; i<(int)ropeHinges_.size(); i++)
        {
            if(ropeHingesToDelete.count(i) == 0)
            {
                newRopeHinges.push_back(ropeHinges_[i]);
            }
        }
    }
    if(!springstodelete.empty() || !particlestodelete.empty() || !rodsToDelete.empty() || !springHingesToDelete.empty() || !ropeHingesToDelete.empty())
    {
        renderLock_.lock();
        {
            if(!ropeHingesToDelete.empty())
                ropeHinges_ = newRopeHinges;
            if(!springHingesToDelete.empty())
                flexibleRodHinges_ = newSpringHinges;
            if(!rodsToDelete.empty())
            {
                rods_ = newrods;
                for(vector<RopeHinge>::iterator hinge = ropeHinges_.begin(); hinge != ropeHinges_.end(); ++hinge)
                {
                    hinge->s1 = remainingrodmap[hinge->s1];
                    hinge->s2 = remainingrodmap[hinge->s2];
                }
            }
            if(!springstodelete.empty())
            {
                springs_ = newsprings;
                for(vector<FlexibleRodHinge>::iterator hinge = flexibleRodHinges_.begin(); hinge != flexibleRodHinges_.end(); ++hinge)
                {
                    hinge->s1 = remainingspringmap[hinge->s1];
                    hinge->s2 = remainingspringmap[hinge->s2];
                }
            }
            if(!particlestodelete.empty())
            {
                particles_ = newparticles;
                for(vector<Spring>::iterator it = springs_.begin(); it != springs_.end(); ++it)
                {
                    it->p1 = remainingparticlemap[it->p1];
                    it->p2 = remainingparticlemap[it->p2];
                }
                for(vector<Rod>::iterator it = rods_.begin(); it != rods_.end(); ++it)
                {
                    it->p1 = remainingparticlemap[it->p1];
                    it->p2 = remainingparticlemap[it->p2];
                }
            }
        }
        renderLock_.unlock();
    }
}

double Simulation::ptSegmentDist(const Vector2d &p, const Vector2d &q1, const Vector2d &q2)
{
    double t = (p-q1).dot(q2-q1) / (q2-q1).dot(q2-q1);
    double linedistsq = (q1 + t*(q2-q1) - p).squaredNorm();
    double q1dist = (p-q1).squaredNorm();
    double q2dist = (p-q2).squaredNorm();
    double mindistsq = min(linedistsq, min(q1dist, q2dist));
    return sqrt(mindistsq);
}
