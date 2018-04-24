#include "Particle.h"
using namespace Eigen;

Particle::Particle(Vector2f _x, Vector2f _v, float _h){
	this->m = 1.0f;
	this->color << 1.0f, 0.0f, 0.0f;
	this->v = _v;
	this->x = _x;
	this->h = _h;
	this->den = 0.0f;

}

Particle:: ~Particle(){
	
}