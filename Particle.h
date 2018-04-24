#pragma once
#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>

class Particle {
	public:
		Particle(Eigen::Vector2f _x, Eigen::Vector2f _v, float _h);
		~Particle();

		Eigen::Vector2f x;
		Eigen::Vector2f v;
		Eigen::Vector2f f;
		Eigen::Vector3f color;
		float pres;
		float den;
		float m;
		float h;
};