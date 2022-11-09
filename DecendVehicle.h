#pragma once

#include <math.h>

#include "SFML/Graphics.hpp"
#include "Planet.h"
#include "excelExporter.h"


#define PI 3.1415926535897932384
#define H_IGNITION_INITIAL 0
#define G_EARTH 9.81

class DecendVehicle
{
public:

	ExcelExporter excelExporter;
	float totalMass;
	float massSC;
	float cX = 0;
	float cY = 0;
	float midSurfase = 0;
	double phi = 0;
	float alpha = 0;
	float alphaDescend = -1.0 / 180.0 * PI;
	double dt = 0.01;
	double time = 0;
	double overLoad;
	double maxOverLoad = 0;
	bool bcolibration = false;
	bool landed = false;


	sf::Vector2f velocity;
	sf::Vector2f position;
	sf::Vector2f relVelocity;
	sf::Vector2f windForce;

	//virtual void initialize(Planet& planet);

	//virtual void control();

	//virtual void dynamic(Planet& planet);

	//virtual double cxCoefficient(double v, double angleA, double h, Planet& planet);

	//virtual double cyCoefficient(double v, double angleA, double h, Planet& planet, double cxC);

	//virtual void calculateOptimalMass(Planet& planet);

	double dot(sf::Vector2f vec1, sf::Vector2f vec2) {
		double dott = vec1.x * vec2.x + vec1.y * vec2.y;
		return dott;
	}

	double vecLength(sf::Vector2f vec) {
		return sqrt(vec.x * vec.x + vec.y * vec.y);
	}

	sf::Vector2f normalize(sf::Vector2f vec) {
		sf::Vector2f newVec;
		newVec.x = vec.x / vecLength(vec);
		newVec.y = vec.y / vecLength(vec);
		return(newVec);
	}

	sf::Vector2f turnVector(sf::Vector2f vec, double a) {

		sf::Vector2f newVec;
		newVec.x = vec.x * cos(a) - vec.y * sin(a);
		newVec.y = vec.x * sin(a) + vec.y * cos(a);
		return newVec;
	}



};

