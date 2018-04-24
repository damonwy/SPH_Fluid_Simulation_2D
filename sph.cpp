#include "sph.h"
#include "Particle.h"

#include <iostream>
#include <cmath>
#include <memory>

using namespace std;
using namespace Eigen;
const float  PI_F = 3.14159265358979f;

SPH::SPH(int n, 
	float l, 
	float u, 
	float h, 
	float _epsilon, 
	float _viscosity_coeff, 
	float _den_bar, 
	float _pres_bar, 
	float _upsilon, 
	float _dt){
	
	this->update_scheme = EL;
	this->epsilon = _epsilon;
	this->viscosity_coeff = _viscosity_coeff;
	this->den_bar = _den_bar;
	this->upsilon = _upsilon;
	this->pres_bar = _pres_bar;
	this->llc << l, l;
	this->urc << u, u;
	this->h = h;
	this->dt = _dt;

	this->grav << 0.0f, 98.0f;
	this->wallsticky  = 0.5f;

	Vector2f v;
	v.setZero();

	Vector2f x;
	x.setZero();

	for(int i = 0; i < n; ++i) {
		x << randFloat(l + 1.0f, u - 1.0f), randFloat(l + 1.0f, u - 1.0f);
		auto p = make_shared<Particle>(x, v, h);

		// give particle random color
		for(int j = 0; j < 3; j++){
			p->color(j) = randFloat( 0.05, 0.95);
		}
		particles.push_back(p);
	}

	updateOVboundary(llc,  urc);
	updateOVindices();
}

float SPH::getWeight(shared_ptr<Particle> pa, shared_ptr<Particle> pb) {
	float w;
	float r = (pa->x - pb->x).norm();
	float h = pa->h;

	if (r / h < 1.0) {
		w = pow((1.0 - r / h), 3) * 10.0f / (pow(h, 2) * PI_F);
	}
	else {
		w = 0.0f;
	}
	return w;
}

Vector2f SPH::getGradWeight(shared_ptr<Particle> pa, shared_ptr<Particle> pb) {
	Vector2f gw;
	gw.setZero();

	float r = (pa->x - pb->x).norm();
	float h = pa->h;

	if (r / h < 1.0f) {
		gw = - pow((1.0f - r / h), 2)* 30.0f / (PI_F * pow(h, 3)) * (pa->x - pb->x).normalized();
	}
	return gw;
}

void SPH::updateDensity(){
	for(int i = 0; i < (int)particles.size(); i++){

		vector<size_t> neighbor_indices;
		auto pa = particles[i];

		// neighbor_indices has all the indices of nearby particles
		findNeighbors(pa, &neighbor_indices);
		pa->den = 0.0f;

		for(int j = 0; j < neighbor_indices.size(); j++){
			int ib = neighbor_indices[j];
			auto pb = particles[ib];
			pa->den += pb->m * getWeight(pa, pb);
		}

		// compute pressure using Tait equation where
		// pres_bar: strength  upsilon: power  den_bar: base density
		pa->pres = pres_bar * (pow(pa->den/den_bar, upsilon) - 1.0f); 
	}
}

void SPH::updateOVboundary(Vector2f _llc, Vector2f _urc){

	nx = int(((_urc - _llc)(0) / h) + 1);
	ny = int(((_urc - _llc)(1) / h) + 1);
	dx = (_urc - _llc)(0) / (float)(nx - 1);
	dy = (_urc - _llc)(1) / (float)(ny - 1);
	occupancy_volume.clear();
	occupancy_volume.resize(nx * ny);
}

void SPH::updateOVindices(){
	for(int i = 0; i < (int)particles.size(); i++){
		auto pa = particles[i];
		int ix, iy;
		int id = getGridIndex(pa, ix, iy);
		occupancy_volume[id].push_back(i);
	}
}

float SPH::computeSV(shared_ptr<Particle> pa, shared_ptr<Particle> pb){
	Vector2f rab = pa->x - pb->x;
	Vector2f vab = pa->v - pb->v;

	float vr = rab.dot(vab);
	float SV;
	float gamma = viscosity_coeff * pa->h / (pa->den + pb->den);
	if(vr < 0.0f){
		SV = -gamma * vr / (pow(rab.norm(), 2) + epsilon * pow(pa->h, 2));
	}else{
		SV = 0.0f;
	}
	return SV;
}

float SPH::computeSP(shared_ptr<Particle> pa, shared_ptr<Particle> pb){
	float SP;
	SP = pa->pres / pow(pa->den, 2) + pb->pres / pow(pb->den, 2);
	return SP;
}

void SPH::updateForces(){

	for(int i = 0; i < (int)particles.size(); i++){
		Vector2f force;
		force.setZero();
		
		Vector2f gw;	// gradient weight
		float SP, SV; 	// parameters for force_pressure, force_viscosity

		vector<size_t> neighbor_indices;
		auto pa = particles[i];
		findNeighbors(pa, &neighbor_indices);

		// add forces from neighbors
		for(int j = 0; j < neighbor_indices.size(); j++){

			int ib = neighbor_indices[j];
			auto pb = particles[ib];
			gw = getGradWeight(pa, pb);
			SP = computeSP(pa, pb);
			SV = computeSV(pa, pb);
			force -= pb->m * (SP + SV) * gw;
		}
		force *= 1e-1;

		// add gravity
		force += grav;
		pa->f = force; 
	}
}

void SPH::findNeighbors(shared_ptr<Particle> pa, vector<size_t> *neighbor_indices){

	// For particle i want a list of particle indices for particles with 2h
    int grid_id;
    int ix, iy;
    grid_id = getGridIndex(pa, ix, iy);

    int l, r, u, d;
    l = ix - 1; 
    r = ix + 1;
    u = iy + 1;
    d = iy - 1;

    vector<size_t> nearby_grids;

    int index;

    //nearby_grids[i] stores the indices of grids around grid i
    /*
	-------------
		|	|
	-------------
		|i  |
	-------------
		|	|
	-------------
    */

    // first check the four corners
    if (l >= 0 && u < ny ){ 
    	index = l + nx * u;
    	nearby_grids.push_back((size_t) index); 
    }

	if (l >= 0 && d >= 0 ) {
    	index = l + nx * d;
     	nearby_grids.push_back((size_t) index); 
 	}

	if (r < nx && d >= 0 ) { 
    	index = r + nx * d;
    	nearby_grids.push_back((size_t) index); 
    }

    if (r < nx && u <  ny ) { 
    	index = r + nx * u;
    	nearby_grids.push_back((size_t) index); 
    }

    // then check four middle grids
    if ( u <  ny ) { 
    	index = ix + nx * u;
    	nearby_grids.push_back((size_t) index);
   	}

 	if (r < nx) { 
    	index = r + nx * iy;
    	nearby_grids.push_back((size_t) index); 
    }

    if (l >= 0) { 
    	index = l + nx * iy;
    	nearby_grids.push_back((size_t) index); 
    }

    if ( d >= 0 ) { 
    	index = ix + nx * d;
    	nearby_grids.push_back((size_t) index); 
    }
  
	nearby_grids.push_back((size_t) grid_id);

    for(int i=0; i < nearby_grids.size(); i++) {

    	int n_pts = occupancy_volume[nearby_grids[i]].size();

        for(int j=0; j < n_pts; j++) {

            neighbor_indices->push_back(occupancy_volume[nearby_grids[i]][j]);
        }
    }
}

void SPH::checkBoundary(float wallsticky){

	for(int i = 0; i < (int)particles.size(); i++){
  		auto pa = particles[i];
  		
  		// modify the velocity and position of particles when out of boundary
  		if(pa->x(0) < LLC(0)){
			pa->x(0) = LLC(0) + (LLC(0) - pa->x(0)) * wallsticky; 
			pa->v(0) = -pa->v(0) * wallsticky;
      	}
  		if(pa->x(0) > URC(0)){
  			pa->x(0) = URC(0) + (URC(0) - pa->x(0)) * wallsticky; 
			pa->v(0) = -pa->v(0) * wallsticky;
  		}
  		if(pa->x(1) < LLC(1)){
			pa->x(1) = LLC(1) + (LLC(1) - pa->x(1)) * wallsticky; 
			pa->v(1) = -pa->v(1) * wallsticky;
      	}
  		if(pa->x(1) > URC(1)){
  			pa->x(1) = URC(1) + (URC(1) - pa->x(1)) * wallsticky; 
			pa->v(1) = -pa->v(1) * wallsticky;
  		}

  		// update the occupancy volume boundary
  		if(pa->x(0) < llc(0)){
      		llc(0) = pa->x(0);
      	}
      	if(pa->x(0) > urc(0)){
      		urc(0) = pa->x(0);
      	}
      	if(pa->x(1) < llc(1)){
      		llc(1) = pa->x(1);
      	}
      	if(pa->x(1) > urc(1)){
      		urc(1) = pa->x(1);
      	}

    }
}

int SPH::getGridIndex(shared_ptr<Particle> pa, int &ix, int &iy){

	// compute the index of a particle given an occupancy volume grid
	Vector2f pos = pa->x;

	ix = int((pos-llc)(0)/(dx));
	iy = int((pos-llc)(1)/(dy));
	int index = ix + nx * iy;
	return index;
}

void SPH::updateFluid(){

	if(update_scheme == LF){
		leapfrog(dt);
	}

	if(update_scheme == SIXTH){
		sixth(dt);
	}

	if(update_scheme == EL){
		euler(dt);
	}
}

void SPH::euler(float dt){

	updateDensity();
    updateForces();

    for(int i = 0; i < (int)particles.size(); i++){
      		auto pa = particles[i];
      		pa->v += pa->f * dt;
      		pa->x += pa->v * dt;
    }

    checkBoundary(wallsticky);
    updateOVboundary(llc, urc);
    updateOVindices();
}

void SPH::leapfrog(float _dt){
		// non-implicit method

		updateDensity();
    	updateForces();

	 	_dt = _dt / 2.0f;

      	for(int i = 0; i < (int)particles.size(); i++){
      		auto pa = particles[i];
      		pa->x += pa->v * _dt; // half step
      	}

      	checkBoundary(wallsticky);
      	updateOVboundary(llc, urc);
    	updateOVindices();

		updateDensity();
    	updateForces();

    	for(int i = 0; i < (int)particles.size(); i++){
      		auto pa = particles[i];
      		pa->v += pa->f * _dt * 2.0f; // full step
      		pa->x += pa->v * _dt;		// half step
      	}

		checkBoundary(wallsticky);
      	updateOVboundary(llc, urc);
    	updateOVindices();
}

void SPH::sixth(float dt) {
  float a, b;
  a = 1.0f / (4.0f - pow(4.0f, 1.0f / 3.0f));
  b = 1.0f - 4.0f * a;

  for(int i = 0; i < 5; i++){
  	  leapfrog(a * dt);
  }
}

float SPH::randFloat(float l, float h){

	float r = rand() / (float)RAND_MAX;
	return (1.0f - r) * l + r * h;
}