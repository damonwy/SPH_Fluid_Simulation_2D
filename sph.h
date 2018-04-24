#pragma once
#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include <vector>

class Particle;
enum Scheme {EL, LF, SIXTH};

class SPH {

public:
	SPH(int n, float l, float u, float h, float _epsilon, float _viscosity_coeff, float _den_bar, float _pres_bar, float _upsilon, float _dt);
	~SPH();

	// fixed values
	Eigen::Vector2f LLC;
	Eigen::Vector2f URC;

	// computed each time step
	Eigen::Vector2f llc;
	Eigen::Vector2f urc;
	int nx;
	int ny;
	float dx;
	float dy;

	// user input
	float epsilon;
	float viscosity_coeff;
	float upsilon;
	float den_bar;
	float pres_bar;
	float h;
	float dt;
	float wallsticky;
	Eigen::Vector2f grav;
	Scheme update_scheme;

	std::vector<std::shared_ptr<Particle> > particles;
	std::vector<std::vector<size_t> > occupancy_volume; // list of nearby particles = occupancy_volume[i]
	
	void updateFluid();

private:
	void sixth(float dt);
	void leapfrog(float _dt);
	void euler(float dt);

	// apply to all particles
	void updateDensity();
	void updateForces();
	void checkBoundary(float wallsticky);
	void updateOVindices();

	// apply to one particle
	void findNeighbors(std::shared_ptr<Particle> pa, std::vector<size_t> *neighbor_indices);
	int getGridIndex(std::shared_ptr<Particle> pa, int &ix, int &iy);

	// apply between two particles
	float computeSV(std::shared_ptr<Particle> pa, std::shared_ptr<Particle> pb);
	float computeSP(std::shared_ptr<Particle> pa, std::shared_ptr<Particle> pb);
	float getWeight(std::shared_ptr<Particle> pa, std::shared_ptr<Particle> pb);
	Eigen::Vector2f getGradWeight(std::shared_ptr<Particle> pa, std::shared_ptr<Particle> pb); 

	// helper
	float randFloat(float l, float h);
	void updateOVboundary(Eigen::Vector2f _llc, Eigen::Vector2f _urc);
};