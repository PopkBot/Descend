#pragma once

#include <math.h>

#include "SFML/Graphics.hpp"
#include "Planet.h"



#define PI 3.1415926535897932384
#define H_IGNITION_INITIAL 0
#define G_EARTH 9.81

class DescendVehicle
{
public:

	
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

	double machArray[FILE_RESOLUTION];
	double angleAttackArray[FILE_RESOLUTION];
	double cXMachArray[FILE_RESOLUTION][FILE_RESOLUTION];
	double adRatioAngleArray[FILE_RESOLUTION][FILE_RESOLUTION];

	virtual void initialize(Planet& planet) {}

	virtual void control(Planet& planet) {}

	virtual void dynamic(Planet& planet) {}

	virtual void calculateOptimalMass(Planet& planet) {}

	virtual void printStats() {}

	virtual void printFinalStats() {}

	void printMatrixParams() {

		printf("\n\nCx\nMach No\t");
		for (int i = 0; i < FILE_RESOLUTION; i++) {
			printf("%-4.2f ", machArray[i]);
		}
		printf("\nAoA");
		for (int i = 0; i < FILE_RESOLUTION; i++) {
			printf("\n%-4.2f\t", angleAttackArray[i] / PI * 180);
			for (int j = 0; j < FILE_RESOLUTION; j++) {
				printf("%-4.2f ", cXMachArray[i][j]);

			}
		}
		printf("\n\nL/D\nMach No\t");
		for (int i = 0; i < FILE_RESOLUTION; i++) {
			printf("%-4.2f ", machArray[i]);
		}
		printf("\nAoA");
		for (int i = 0; i < FILE_RESOLUTION; i++) {
			printf("\n%-4.2f\t", angleAttackArray[i] / PI * 180);
			for (int j = 0; j < FILE_RESOLUTION; j++) {
				printf("%-4.2f ", adRatioAngleArray[i][j]);

			}
		}


	}

	double cxCoefficient(double v, double angleA, double h, Planet& planet) {

		double cxCoef;
		int cxAIndex = FILE_RESOLUTION - 1;
		int cxMIndex = FILE_RESOLUTION - 1;
		double machNo = v / planet.getSonicSpeed(h);

		bool cxAFound = false;
		bool cxMFound = false;

		for (int i = 0; i < FILE_RESOLUTION; i++) {

			if (angleA <= angleAttackArray[i] && !cxAFound) {
				cxAIndex = i;
				cxAFound = true;
			}
			if (machNo <= machArray[i] && !cxMFound) {
				cxMIndex = i;
				cxMFound = true;
			}

		}
		if (cxMIndex > 0) {
			//cxCoef = cX * (cXMachArray[cxAIndex][cxMIndex] - cXMachArray[cxAIndex][cxMIndex - 1]) / (machArray[cxMIndex] - machArray[cxMIndex - 1]) * (machNo - machArray[cxMIndex - 1]) + cXMachArray[cxAIndex][cxMIndex - 1];
			cxCoef = cX * cXMachArray[cxAIndex][cxMIndex];
			return(cxCoef);
		}
		else {
			cxCoef = cX * cXMachArray[cxAIndex][cxMIndex];
		}


		return(cxCoef);

	}

	double cyCoefficient(double v, double angleA, double h, Planet& planet, double cxC) {

		double LDCoef;
		int cxAIndex = FILE_RESOLUTION - 1;
		int cxMIndex = FILE_RESOLUTION - 1;
		double machNo = v / planet.getSonicSpeed(h);

		bool cxAFound = false;
		bool cxMFound = false;

		for (int i = 0; i < FILE_RESOLUTION; i++) {

			if (angleA <= angleAttackArray[i] && !cxAFound) {
				cxAIndex = i;
				cxAFound = true;
			}
			if (machNo <= machArray[i] && !cxMFound) {
				cxMIndex = i;
				cxMFound = true;
			}

		}
		if (cxMIndex > 0) {
			//LDCoef = (adRatioAngleArray[cxAIndex][cxMIndex] - adRatioAngleArray[cxAIndex][cxMIndex - 1]) / (machArray[cxMIndex] - machArray[cxMIndex - 1]) * (machNo - machArray[cxMIndex - 1]) + adRatioAngleArray[cxAIndex][cxMIndex - 1];
			LDCoef = adRatioAngleArray[cxAIndex][cxMIndex];
		}
		else {
			LDCoef = adRatioAngleArray[cxAIndex][cxMIndex];
		}
		//printf("\n cyccc %.5f\tcxC %0.2f\tld %.2f\tang %.2f\n", cxC * ld,cxC,ld,angleA);
		return (LDCoef * cxC);


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

	sf::Vector2f turnVector(sf::Vector2f vec, double a) {

		sf::Vector2f newVec;
		newVec.x = vec.x * cos(a) - vec.y * sin(a);
		newVec.y = vec.x * sin(a) + vec.y * cos(a);
		return newVec;
	}



};

