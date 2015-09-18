#ifndef SIMULATION_H
#define SIMULATION_H

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>
#include <set>
#include <QMutex>
#include "fluid.h"
#include "simparameters.h"

typedef Eigen::Triplet<double> Tr;

struct Particle
{
public:
    Particle(Eigen::Vector2d pos, double mass, bool isFixed, bool inert) : pos(pos), mass(mass), fixed(isFixed), inert(inert)
    {
        vel.setZero();
        prevpos = pos;
    }
    Particle(Eigen::Vector2d pos, double mass, bool isFixed) : pos(pos), mass(mass), fixed(isFixed)
    {
        vel.setZero();
        prevpos = pos;
    }

    Eigen::Vector2d pos;
    Eigen::Vector2d prevpos;
    Eigen::Vector2d vel;
    double mass;
    bool fixed;
    bool inert;
};

struct Connector
{
public:
    Connector(int p1, int p2, double restlen, double mass) : p1(p1), p2(p2), restlen(restlen), mass(mass) {}
    Connector(int p1, int p2, double restlen) : p1(p1), p2(p2), restlen(restlen) {
        mass = 0;
    }
    int p1;
    int p2;
    double restlen;
    double mass;
};

struct Rod: public Connector
{
public:
    Rod(int p1, int p2, double restlen) : Connector(p1, p2, restlen) {
        isHinge = false;
    }
    Rod(int p1, int p2, double restlen, double mass, bool isHinge) : Connector(p1, p2, restlen, mass), isHinge(isHinge) {
    }
    bool isHinge;
};

struct FlexibleRodHinge
{
public:
    FlexibleRodHinge(int s1, int s2, double stiffness) : s1(s1), s2(s2), stiffness(stiffness) {
    }

    int s1;
    int s2;
    double stiffness;
};

struct RopeHinge
{
public:
    RopeHinge(int s1, int s2, double stiffness) : s1(s1), s2(s2), stiffness(stiffness) {
    }

    int s1;
    int s2;
    double stiffness;
};

struct Spring: public Connector
{
public:
    Spring(int p1, int p2, double stiffness, double restlen) : Connector(p1, p2, restlen), stiffness(stiffness) {
        unsnappable = false;
    }
    Spring(int p1, int p2, double stiffness, double restlen, double mass, bool unsnappable) : Connector(p1, p2, restlen, mass), stiffness(stiffness), unsnappable(unsnappable) {}
    bool unsnappable;
    double stiffness;
};

class Simulation
{
public:
    Simulation(const SimParameters &params);

    void takeSimulationStep();
    void render();
    void clearScene();
    void fluidSimulationStep();
    void advection(int boundary, Eigen::MatrixXd &d, Eigen::MatrixXd &dOld, Eigen::MatrixXd &uOld, Eigen::MatrixXd &vOld);
    void diffuse(int boundary, Eigen::MatrixXd &d, Eigen::MatrixXd &dOld, double diffFactor);
    void linearSolver(int b, Eigen::MatrixXd &x, Eigen::MatrixXd &xOld, double a, double c);
    void project(Eigen::MatrixXd &x, Eigen::MatrixXd &y, Eigen::MatrixXd &xOld, Eigen::MatrixXd &yOld);
    void swap(Eigen::MatrixXd &left, Eigen::MatrixXd &right);
    void addSource(Eigen::MatrixXd &d, Eigen::MatrixXd& dOld);
    void addVelocity(double x, double y, double velX, double velY);
    void addDensity(double x, double y);
    void setBoundry(int b, Eigen::MatrixXd& m);

    void massSpringSimulationStep();

    // Particle Stuff
    void addParticle(double x, double y);

private:
    const SimParameters &params_;
    QMutex renderLock_;

    double time_;
    Fluid *fluid_;

    std::vector<Particle> particles_;
    std::vector<Spring> springs_;
    std::vector<FlexibleRodHinge> flexibleRodHinges_;
    std::vector<RopeHinge> ropeHinges_;
    std::vector<Rod> rods_;

    void addParticle(double x, double y, bool particleFixed);
    void buildConfiguration(Eigen::VectorXd &q, Eigen::VectorXd &qprev, Eigen::VectorXd &v);
    void unbuildConfiguration(const Eigen::VectorXd &q, const Eigen::VectorXd &v);
    void computeForceAndHessian(const Eigen::VectorXd &q, const Eigen::VectorXd &qprev, Eigen::VectorXd &F, Eigen::SparseMatrix<double> &H, Eigen::VectorXd &v);
    void processElasticBendingForce(const Eigen::VectorXd &q, Eigen::VectorXd &F);
    void processSpringForce(const Eigen::VectorXd &q, Eigen::VectorXd &F, std::vector<Tr> &H);
    void processDampingForce(const Eigen::VectorXd &q, const Eigen::VectorXd &qprev, Eigen::VectorXd &F, std::vector<Tr> &H);
    void processFloorForce(const Eigen::VectorXd &q, const Eigen::VectorXd &qprev, Eigen::VectorXd &F, std::vector<Tr> &H);
    void processFluidForce(const Eigen::VectorXd &q,  Eigen::VectorXd &F, const Eigen::VectorXd &v);
    void computeMassInverse(Eigen::SparseMatrix<double> &Minv);
    Eigen::SparseMatrix<double> computeGradGTranspose(const Eigen::VectorXd &q);

    void numericalIntegration(Eigen::VectorXd &q, Eigen::VectorXd &qprev, Eigen::VectorXd &v);
    void computeStepProject(Eigen::VectorXd &q, Eigen::VectorXd &oldq, Eigen::VectorXd &v);
    void computeLagrangeMultipliers(const Eigen::VectorXd &q, const Eigen::VectorXd &F, Eigen::VectorXd &v);

    void pruneOverstrainedSprings();
    void removeOutsideParticles();
    double ptSegmentDist(const Eigen::Vector2d &p, const Eigen::Vector2d &q1, const Eigen::Vector2d &q2);
    void printSimParameters();

};

#endif // SIMULATION_H
