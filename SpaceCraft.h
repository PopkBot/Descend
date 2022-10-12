#pragma once
//Импорт -начало-

#include <math.h>

#include "SFML/Graphics.hpp"
#include "excelExporter.h"
#include "Planet.h"

//Импорт -конец-

//Объявление констант -начало-

#define SC_PARAMS_COUNT 8
#define SC_DIRECTORY "Space_Crafts\\"
#define PI 3.1415926535897932384

//Объявление констант -конец-
class SpaceCraft
{
public:
	//Объявление переменных -начало-

	ExcelExporter excelExporter;

	double massSC;
	double cX;
	double cY;
	double midSurfase;
	double phi;
	double dt = 0.05;
	double time = 0;

	sf::Vector2f velocity;
	sf::Vector2f position;
	//sf::Vector2f pressureCenter;
	



	//Объявление переменных -конец-


	//Методы -начало-
	SpaceCraft(string fileName) {

		double params[SC_PARAMS_COUNT];
		excelExporter.extractDataFromFile(SC_DIRECTORY, fileName, params, SC_PARAMS_COUNT);
		
		massSC=params[0];
		cX = params[1];
		cY = params[2];
		midSurfase = params[3];
		phi = params[5]*PI/180;
		velocity.x = params[4]*sin(phi);
		velocity.y = -params[4]*cos(phi);
		position.x = params[6];
		position.y = params[7];
		//pressureCenter.x = params[8];
		//pressureCenter.y = params[9];

		printf("\nSpaceCraft parametrs:\nmassSC %2f\ncX %3f\ncY %3f\nVx %4f\nVy %4f\nphi %3f\nX %3f\nY %3f\n",
			massSC, cX, cY, velocity.x, velocity.y, phi*180/3.1415, position.x, position.y);

	}

	void dynamic(Planet &planet) {

		
		sf::Vector2f windForceLocal;
		sf::Vector2f windForceGlobal;
		sf::Vector2f acselGlobal;
		sf::Vector2f gravityForce;
		double qForce = cX * (pow(velocity.x - planet.windSpeed, 2) + velocity.y * velocity.y) * planet.getDensityByHeight(position.y) / 2 * midSurfase;
		double yForce = cY * (pow(velocity.x - planet.windSpeed, 2) + velocity.y * velocity.y) * planet.getDensityByHeight(position.y) / 2 * midSurfase;

		double alpha;
		if (planet.windSpeed < 0)
			alpha =acos( dot(normalize( velocity),normalize( sf::Vector2f{ velocity.x -(float) planet.windSpeed,velocity.y })));
		else if (planet.windSpeed > 0)
			alpha = -acos(dot(normalize(velocity), normalize(sf::Vector2f{ velocity.x - (float)planet.windSpeed,velocity.y })));
		else
			alpha = 0;

		windForceLocal.x = qForce * cos(alpha) - yForce * sin(alpha);
		windForceLocal.y = qForce * sin(alpha) + yForce * cos(alpha);

		//windForceLocal = -normalize(sf::Vector2f{ velocity.x - (float)planet.windSpeed,velocity.y })*(float)qForce
		//				+(float)yForce*normalize(turnVector(sf::Vector2f{ velocity.x - (float)planet.windSpeed,velocity.y }, 0.5 * PI));
		
		windForceGlobal.y = windForceLocal.x * cos(phi) + windForceLocal.y * sin(phi);
		windForceGlobal.x = -windForceLocal.x * sin(phi) + windForceLocal.y * cos(phi);

		//windForceGlobal = turnVector(windForceLocal,- PI / 2 + phi);

		gravityForce.y =- planet.gravPar * massSC / pow(planet.radius+position.y, 2);
		gravityForce.x = 0;
		
		acselGlobal = (gravityForce + windForceGlobal) / (float) massSC;

		velocity =velocity+ acselGlobal *(float) dt;

		position = position + velocity * (float) dt;

		if (velocity.x > 0)
			phi = acos(dot(normalize(velocity), sf::Vector2f{ 0,-1 }));
		else if (velocity.x < 0)
			phi = -acos(dot(normalize(velocity), sf::Vector2f{ 0,-1 }));
		else
			phi = 0;
		time += dt;

		if (position.y < 0)
			position.y = 0;

		printf("t %2f\t", time);
	}

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

	sf::Vector2f turnVector(sf::Vector2f vec, double a){
	
		sf::Vector2f newVec;
		newVec.x = vec.x * cos(a) - vec.y * sin(a);
		newVec.y = vec.x * sin(a) + vec.y * cos(a);
		return newVec;
	}









	//Методы -конец-
};

