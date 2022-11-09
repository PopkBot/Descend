#pragma once


#include "DecendVehicle.h"
#include <math.h>


#define PARACHUTE_PARAMS_COUNT 13
#define PARACHUTE_MATRIX_PARAPMS 5

#define H_EPS 0.1
#define MASS_EPS 1
#define V_EPS 0.1
#define V_DECEND 1
#define V_LANDING 6

#define LANDING_OVERLOAD 6
#define THRUST_SPEC_IMP 2300


#define PARACHUTE_DIRECTORY "Space_Crafts\\"
#define PARACHUTE_MATRIX_PARAMS_DIR "Space_Crafts\\SC_matrix_parametrs\\"

enum ParashuteType {
	DRAG = 0,
	MAIN = 1,
	SPARE = 2
};


class ParachuteDecend : public DecendVehicle 
{

public:

	double params[PARACHUTE_PARAMS_COUNT];

	double machArray[FILE_RESOLUTION];
	double cXMachArray[FILE_RESOLUTION];
	double angleAttackArray[FILE_RESOLUTION];
	double cXAtackArray[FILE_RESOLUTION];
	double adRatioAngleArray[FILE_RESOLUTION];

	double parachuteDensity;
	double parashuteDEployHight=0;
	double parSysMass = 0;
	double hz = 100;
	double g0;

	float fuelMass = 0;
	float mfuel = 0;
	float thrust=0;
	float ignitionHight = 0;
	float IgnitionVel = 0;
	float thrustMass = 0;
	
	float maxQH = 0;
	float cP[3];
	float parashuteDensity[3];
	double parashuteSurface[3] = {1,1,1};
	double maxQForParashute[3] = { 24060,			// 0.66 * 72900 / 2.0
								1000,				// 0.8*2500/2.0
								1125 };				// 0.9 * 2500 / 2.0

	bool parashuteDeployed = false;
	bool mainParashuteAvailable = true;
	bool onlyDragPar = false;
	bool thrustIsOn = false;
	bool isParashute =true;

	
	

	ParashuteType parashuteType;

	ParachuteDecend(string fileName, string matrixParFileName) {

		double dataArray[PARACHUTE_MATRIX_PARAPMS][FILE_RESOLUTION];

		excelExporter.extractMatrixFromFile(PARACHUTE_MATRIX_PARAMS_DIR, matrixParFileName, dataArray, FILE_RESOLUTION, PARACHUTE_MATRIX_PARAPMS);

		for (int i = 0; i < FILE_RESOLUTION; i++) {
			machArray[i] = dataArray[0][i];
			cXMachArray[i] = dataArray[1][i];
			angleAttackArray[i] = dataArray[2][i];
			cXAtackArray[i] = dataArray[3][i];
			adRatioAngleArray[i] = dataArray[4][i];
		}

		excelExporter.extractDataFromFile(PARACHUTE_DIRECTORY, fileName, params, PARACHUTE_PARAMS_COUNT);

		massSC = params[0];						//Масса КА
		totalMass = massSC + parSysMass;
		cX = params[1];							//сХ
		midSurfase = params[2];					//Площадь Миделя
		phi = params[4] * PI / 180;				//Начальный угол входа 
		velocity.x = params[3] * sin(phi);		//Vx
		velocity.y = -params[3] * cos(phi);		//Vy
		relVelocity = velocity;
		position.x = params[5];					//X
		position.y = params[6];					//Y
		parashuteDensity[ParashuteType::DRAG] = params[7];			//Плотность ткани тормозного парашюта
		parashuteDensity[ParashuteType::MAIN] = params[8];			//Плотность ткани основного парашюта
		parashuteDensity[ParashuteType::SPARE] = params[9];			//Плотность ткани запасного парашюта
		cP[ParashuteType::DRAG] =  params[10];	//сП тормозного парашюта
		cP[ParashuteType::MAIN] = params[11];	//сП основного парашюта
		cP[ParashuteType::SPARE] = params[12];	//сП запасного парашюта

		parashuteDeployed = false;
		mainParashuteAvailable = true;


	}

	void initialize(Planet& planet) {


		massSC = params[0];						//Масса КА
		totalMass = massSC;
		cX = params[1];							//сХ
		midSurfase = params[2];					//Площадь Миделя
		phi = params[4] * PI / 180;				//Начальный угол входа 
		velocity.x = params[3] * sin(phi);		//Vx
		velocity.y = -params[3] * cos(phi);		//Vy
		relVelocity = velocity;
		position.x = params[5];					//X
		position.y = params[6];					//Y
		parashuteDensity[ParashuteType::DRAG] = params[7];			//Плотность ткани тормозного парашюта
		parashuteDensity[ParashuteType::MAIN] = params[8];			//Плотность ткани основного парашюта
		parashuteDensity[ParashuteType::SPARE] = params[9];			//Плотность ткани запасного парашюта
		cP[ParashuteType::DRAG] = params[10];	//сП тормозного парашюта
		cP[ParashuteType::MAIN] = params[11];	//сП основного парашюта
		cP[ParashuteType::SPARE] = params[12];	//сП запасного парашюта

		parashuteDeployed = false;
		mainParashuteAvailable = true;

		dt = 0.01;
		time = 0;
		maxOverLoad = 0;
		landed = false;
		mfuel = fuelMass;
		g0 = (planet.gravPar / planet.radius) / planet.radius;
		thrustIsOn = false;
		overLoad = 6;
		if (!bcolibration) {
			totalMass = massSC + parSysMass + thrustMass;
		}
	}

	void control(double density) {

		double currentQ = density * pow(vecLength(velocity), 2) / 2.0;

		parashuteDeployed = false;
		//printf("\nflag00");
		if (position.y > parashuteDEployHight) {
			alpha = alphaDescend;
		}
		else {
			alpha = 0;
		}
		if (position.y < maxQH && isParashute) {
			if (position.y < parashuteDEployHight && maxQForParashute[ParashuteType::DRAG] >= currentQ) {
				parashuteDeployed = true;
				parashuteType = ParashuteType::DRAG;
				//printf(" current q %10.2f maxQ %10.2f", currentQ, maxQForParashute[ParashuteType::DRAG]);
			}
			if (position.y < parashuteDEployHight && maxQForParashute[ParashuteType::MAIN] >= currentQ && !onlyDragPar) {
				parashuteDeployed = true;
				parashuteType = ParashuteType::MAIN;
				//printf("\nflag01 current q %10.2f maxQ %10.2f", currentQ, maxQForParashute[ParashuteType::MAIN]);
			}
			if (!mainParashuteAvailable && position.y < parashuteDEployHight && maxQForParashute[ParashuteType::SPARE] >= currentQ) {
				parashuteDeployed = true;
				parashuteType = ParashuteType::SPARE;
				//printf("\nflag02");
			}
			if (parashuteDeployed && parashuteType == ParashuteType::MAIN && maxQForParashute[ParashuteType::MAIN] < currentQ) {
				mainParashuteAvailable = false;
				//printf("\nflag03");
			}
		}
		if (position.y <= ignitionHight) {

			thrustIsOn = true;

		}

		//printf("\t depl %2d", parashuteDeployed);


	}

	void dynamic(Planet& planet) {


		sf::Vector2f acselGlobal;
		sf::Vector2f gravityForce;
		sf::Vector2f centripicalAcs = sf::Vector2f{ 0,0 };
		float windSpeed = planet.getWindSpeed(position.y);

		float g0 = (planet.gravPar / planet.radius) / planet.radius;


		float dAlpha = windSpeed * cos(phi) / vecLength(velocity);


		relVelocity = velocity - sf::Vector2f{ windSpeed,0 };


		float cxCoef = cxCoefficient(vecLength(relVelocity), alpha + dAlpha, position.y, planet);
		float cyCoef = cyCoefficient(vecLength(relVelocity), alpha + dAlpha, position.y, planet, cxCoef);

		float qForce = cxCoef * pow(vecLength(relVelocity), 2) * planet.getDensityByHeight(position.y) / 2 * midSurfase;
		float yForce = cyCoef * pow(vecLength(relVelocity), 2) * planet.getDensityByHeight(position.y) / 2 * midSurfase;


		

		windForce = -normalize(relVelocity) * qForce + turnVector(normalize(relVelocity), PI / 2) * yForce;

		

		gravityForce.y = -planet.gravPar * totalMass / pow(planet.radius + position.y, 2);
		gravityForce.x = 0;


		
		centripicalAcs.y = pow(velocity.x, 2) / (planet.radius + position.y);

		sf::Vector2f dragF = parashuteDragForce(planet.getDensityByHeight(position.y), relVelocity);

		acselGlobal = (gravityForce + windForce + dragF + getThrustForce(g0)) / totalMass + centripicalAcs;


		//printf("\nh %5.1f acs %10.2f v %5.2f gF %10.2f wF %10.2f dF %10.2f depl %2d thr %-10.2f",position.y,acselGlobal.y,vecLength(velocity),
		//	gravityForce.y,windForce.y, dragF.y,parashuteDeployed, getThrustForce(g0).y);


		
		overLoad = vecLength(acselGlobal) / G_EARTH;
		if (overLoad > maxOverLoad) {
			maxOverLoad = overLoad;
		}

		velocity = velocity + acselGlobal * (float)dt;

		position = position + velocity * (float)dt;

	

		if (velocity.x > 0)
			phi = acos(dot(normalize(velocity), sf::Vector2f{ 0,-1 }));
		else if (velocity.x < 0)
			phi = -acos(dot(normalize(velocity), sf::Vector2f{ 0,-1 }));
		else
			phi = 0;


		if (abs(position.y - parashuteDEployHight) < 100 || position.y < (hz * 1.1)) {
			dt = 0.0001;
		}
		else {
			dt = 0.01;
		}


		time += dt;

		/*
		if (position.y < H_EPS && bcolibration) {

		}
		else if (position.y < H_EPS / 100 && vecLength(velocity) < V_EPS * 10 && !bcolibration) {

			position.y = 0;
			velocity.x = 0;
			velocity.y = 0;
			phi = 0;
		}
		else if (position.y < -H_EPS) {
			position.y = 0;

			landed = true;

		}
		*/

		//printf("t %2f\t", time);
	
	}

	double cxCoefficient(double v, double angleA, double h, Planet& planet) {

		double cxCoef;
		double cxA = 1;
		double cxM = 1;
		double machNo;
		if (planet.getSonicSpeed(h) > 0) {
			machNo = v / planet.getSonicSpeed(h);
		}
		else {
			machNo = 0;
		}

		bool cxAFound = false;
		bool cxMFound = false;

		for (int i = 0; i < FILE_RESOLUTION; i++) {

			if (angleA < angleAttackArray[i] && !cxAFound) {
				cxA = (cXAtackArray[i] - cXAtackArray[i - 1]) / (angleAttackArray[i] - angleAttackArray[i - 1]) * (angleA - angleAttackArray[i - 1]) + cXAtackArray[i - 1];
				cxAFound = true;
			}
			if (machNo < machArray[i] && !cxMFound) {
				cxM = (cXMachArray[i] - cXMachArray[i - 1]) / (machArray[i] - machArray[i - 1]) * (machNo - machArray[i - 1]) + cXMachArray[i - 1];
				cxMFound = true;
			}

		}

		if (!cxAFound) {
			cxA = cXAtackArray[FILE_RESOLUTION - 1];
		}
		if (!cxMFound) {
			cxM = cXMachArray[FILE_RESOLUTION - 1];
		}
		cxCoef = cX * cxA * cxM;
		return(cxCoef);

	}

	double cyCoefficient(double v, double angleA, double h, Planet& planet, double cxC) {

		double ld = 0;
		bool ldFound = false;

		for (int i = 0; i < FILE_RESOLUTION; i++) {

			if (angleA <= angleAttackArray[i] && !ldFound) {
				ld = (adRatioAngleArray[i] - adRatioAngleArray[i - 1]) / (angleAttackArray[i] - angleAttackArray[i - 1]) * (angleA - angleAttackArray[i - 1]) + adRatioAngleArray[i - 1];
				ldFound = true;

			}
		}
		if (!ldFound) {
			ld = adRatioAngleArray[FILE_RESOLUTION - 1];
		}
		//printf("\n cyccc %.5f\tcxC %0.2f\tld %.2f\tang %.2f\n", cxC * ld,cxC,ld,angleA);
		return (cxC * ld);


	}


	void calculateOptimalMass(Planet& planet) {

		double g0 = (planet.gravPar / planet.radius) / planet.radius;

		//while (true) {

			initialize(planet);
			double maxQ = 0;

			//double mTPx;
			double mTP;
			double mOP;
			double mTPxSys;
			double mTPSys;
			double mOpSys;
			double mTPxSpareSys;
			double mOpSpareSys;

			float minMass;
			float freeFallVelocity;
			float optimalParVel=0;
			
			isParashute = true;
		

			
			while ( position.y > H_EPS ) {

				control(planet.getDensityByHeight(position.y));
				dynamic(planet);
				if (planet.getDensityByHeight(position.y) * pow(vecLength(velocity), 2) / 2.0 > maxQ) {

					maxQ = planet.getDensityByHeight(position.y) * pow(vecLength(velocity), 2) / 2.0;
					maxQH = position.y;
					//printf("\nmaxQ = %10.2f maxQ H = %10.2f", maxQ, maxQH);
				}

			}

			freeFallVelocity = abs(velocity.y);
			
			printf("\nmaxQ = %10.2f maxQ H = %10.2f freeFallVel %-10.2f", maxQ, maxQH, freeFallVelocity);

			initialize(planet);
			parashuteDEployHight = hz;

			//double overLoad = 2;

			bcolibration = true;

			float cxCoef = cxCoefficient(vecLength(sf::Vector2f{ 0,-V_LANDING }), 0,0, planet);
			float cyCoef = cyCoefficient(vecLength(sf::Vector2f{ 0,-V_LANDING }), 0,0, planet, cxCoef);

			float qForce = cxCoef * pow(vecLength(sf::Vector2f{ 0,-V_LANDING }), 2) * planet.getDensityByHeight(0) / 2 * midSurfase;
			float yForce = cyCoef * pow(vecLength(sf::Vector2f{ 0,-V_LANDING }), 2) * planet.getDensityByHeight(0) / 2 * midSurfase;

			windForce = sf::Vector2f{0,0};
			

		

			calculateParashute(windForce.y, V_LANDING, planet, massSC);
			totalMass += parSysMass;
			minMass = totalMass;

			printf("\nOnly Par parSysMass %5.2f", parSysMass);
			
			for (float parDecVel = V_LANDING; parDecVel <= freeFallVelocity; parDecVel += 0.1) {

				initialize(planet);


				cxCoef = cxCoefficient(parDecVel, 0, 0, planet);
				cyCoef = cyCoefficient(parDecVel, 0, 0, planet, cxCoef);

				qForce = cxCoef * pow(parDecVel, 2) * planet.getDensityByHeight(0) / 2 * midSurfase;
				yForce = cyCoef * pow(parDecVel, 2) * planet.getDensityByHeight(0) / 2 * midSurfase;

				windForce = -normalize(sf::Vector2f{ 0,-parDecVel }) * qForce + turnVector(normalize(sf::Vector2f{ 0,-parDecVel }), PI / 2) * yForce;

				float thrMassBuf = getThrusterMass(parDecVel - V_LANDING, totalMass);
				totalMass += getThrusterMass(parDecVel - V_LANDING, totalMass);

	
			

				calculateParashute(windForce.y, parDecVel, planet, massSC);
				totalMass += parSysMass;

				totalMass -=  thrMassBuf;
				totalMass += getThrusterMass(parDecVel - V_LANDING, totalMass);

				if (totalMass < minMass) {
					minMass = totalMass;
					optimalParVel = parDecVel;
				}
				//printf("\nparVel %-5.2f totalMass %-10.2f parSysMass %-10.2f thrMass %-5.2f", parDecVel, totalMass,parSysMass, getThrusterMass(parDecVel - V_LANDING, totalMass));

			}

			initialize(planet);


			
			
			totalMass += getThrusterMass(freeFallVelocity - V_LANDING, totalMass);

			if (totalMass < minMass) {
				initialize(planet);
				thrustMass = getThrusterMass(freeFallVelocity - V_LANDING, totalMass);
				minMass = totalMass;
				optimalParVel = freeFallVelocity;
				isParashute = false;
				parSysMass = 0;
				//printf("\nflag0");
			}
			else {
				initialize(planet);
				float thrMassBuf = getThrusterMass(optimalParVel - V_LANDING, totalMass);
				totalMass += getThrusterMass(optimalParVel - V_LANDING, totalMass);
				cxCoef = cxCoefficient(optimalParVel, 0, 0, planet);
				cyCoef = cyCoefficient(optimalParVel, 0, 0, planet, cxCoef);
				qForce = cxCoef * pow(optimalParVel, 2) * planet.getDensityByHeight(0) / 2 * midSurfase;
				yForce = cyCoef * pow(optimalParVel, 2) * planet.getDensityByHeight(0) / 2 * midSurfase;
				windForce = -normalize(sf::Vector2f{ 0,-optimalParVel }) * qForce + turnVector(normalize(sf::Vector2f{ 0,-optimalParVel }), PI / 2) * yForce;
				calculateParashute(windForce.y, optimalParVel, planet, totalMass);
				totalMass += parSysMass;

				totalMass -= thrMassBuf;
				thrustMass = getThrusterMass(optimalParVel - V_LANDING, totalMass);
				totalMass += getThrusterMass(optimalParVel - V_LANDING, totalMass);
				isParashute = true;
				//printf("\nflag1");
			}

			printf("\nMin parVel %-5.2f totalMass %-10.2f parSysMass %-10.2f thrMass %-5.2f", optimalParVel, totalMass, parSysMass, getThrusterMass(optimalParVel - V_LANDING, totalMass));


		
			printf("\n only thrust: thrust mass %5.2f", getThrusterMass(freeFallVelocity, massSC));
			//totalMass = minMass;
			
			

			thrust = totalMass * g0 * LANDING_OVERLOAD;
			fuelMass = (optimalParVel-V_LANDING) * totalMass / THRUST_SPEC_IMP * LANDING_OVERLOAD / sqrtf(LANDING_OVERLOAD - 1) * atan((1 - V_LANDING / optimalParVel) / (sqrtf(LANDING_OVERLOAD - 1) * (V_LANDING / optimalParVel + 1)));
			mfuel = fuelMass;
			IgnitionVel = optimalParVel;
			ignitionHight = getIgnitionHight(g0, optimalParVel, V_LANDING);
			//findOptimalDescendAngle(planet);
			//printf("\nmTPx %10.2f mTP %10.2f mOP %10.2f\n", mTPx / parashuteDensity[ParashuteType::DRAG], mTP / parashuteDensity[ParashuteType::DRAG], mOP / parashuteDensity[ParashuteType::MAIN]);
			bcolibration = false;
			float fuelBuff = fuelMass;
			fuelMass = 20*fuelMass;
			mfuel = fuelMass;
			printf("\noptParVel %-5.2f totalMass %-10.2f parSysMass %-5.2f trhustMass %-5.2f ignH %-5.2f", optimalParVel, totalMass, parSysMass,
				getThrusterMass(optimalParVel - V_LANDING, totalMass), getIgnitionHight(g0, IgnitionVel, V_LANDING));

			int drawCount = 0;

			while (true)
			{
				initialize(planet);
				//while (abs(velocity.y)>0.1) {
				while (!(position.y < (H_EPS * 2) && abs(velocity.y) - V_LANDING < V_EPS)) {

					control(planet.getDensityByHeight(position.y));
					dynamic(planet);
					//printf("hIgn %.4f\t\tY %.4f\t\tvy %.5f fuel mass %-10.2f\n", ignitionHight, position.y, velocity.y, mfuel);
					if (drawCount >= 20) {
						
						drawCount = 0;
						//printf("hIgn %.4f\t\tY %.4f\t\tvy %.5f fuel mass %-10.2f\n", ignitionHight, position.y, velocity.y, mfuel);


					}
					drawCount++;
				}


				if (0< (position.y) < H_EPS) {
					//if(abs(position.y/10) > 1){
					ignitionHight -= position.y / abs(position.y)*0.01;
					printf("hIgn %.4f\t\tY %.4f\t\tvy %.5f tot mass %-10.2f\n", ignitionHight, position.y, velocity.y,totalMass);
				}
				else {
					//printf("\nFFFFFFFFhIgn %.4f\t\tY %.4f\t\tvy %.5f tot mass %-10.2f\n", ignitionHight, position.y, velocity.y, totalMass);
					break;
				}
			}
			initialize(planet);
			while (!(position.y < (H_EPS * 2) && abs(velocity.y) - V_LANDING < V_EPS)) {

				control(planet.getDensityByHeight(position.y));
				dynamic(planet);
				//printf("hIgn %.4f\t\tY %.4f\t\tvy %.5f tot mass %-10.2f\n", ignitionHight, position.y, velocity.y, totalMass);
			}
			fuelMass -= mfuel;
			printf("hIgn %.4f\t\tY %.4f\t\tvy %.5f tot mass %-10.2f\n", ignitionHight, position.y, velocity.y, totalMass);





























			printf("\noptParVel %-5.2f totalMass %-10.2f parSysMass %-5.2f trhustMass %-5.2f ignH %-5.2f", optimalParVel,totalMass, parSysMass,
						getThrusterMass(optimalParVel - V_LANDING, totalMass), getIgnitionHight(g0, IgnitionVel, V_LANDING));
			/*
			while (true)
			{
				initialize(planet);
				
				
				while (!(position.y < hz && abs(velocity.y) < V_DECEND) && position.y > H_EPS) {

					control(planet.getDensityByHeight(position.y));
					dynamic(planet);

				}


				if (abs(velocity.y) > V_DECEND) {
					//if(abs(position.y/10) > 1){
					parashuteDEployHight += 1;
					//parashuteSurface[ParashuteType::DRAG] += 1;
					printf("hParDep %.4f\t\tsurf %.4f\t\tvy %.5f\tdragForce %10.2f\tdens %5.2f\tparDep %d\n", parashuteDEployHight, parashuteSurface[ParashuteType::DRAG], velocity.y,vecLength( parashuteDragForce(planet.getDensityByHeight(0),velocity)), planet.getDensityByHeight(0),parashuteDeployed);
				}
				else {
					//printf("FFFFFFFFhIgn %.4f\t\tY %.4f\t\tvy %.5f\n", hIgnition, position.y, velocity.y);
					bcolibration = false;
					printf("hParDep %.4f\t\tsurf %.4f\t\tvy %.5f\tdragForce %10.2f\tdens %5.2f\tparDep %d\n", parashuteDEployHight, parashuteSurface[ParashuteType::DRAG], velocity.y, vecLength(parashuteDragForce(planet.getDensityByHeight(0), velocity)), planet.getDensityByHeight(0), parashuteDeployed);

					//printf("\nflag 002\n");
					break;
				}
			}
			*/


		//}

	}

	void calculateParashute(float dragForce,float vel,Planet& planet,float totalMass1) {

		double mTPx;
		double mTP;
		double mOP;
		double mTPxSys;
		double mTPSys;
		double mOpSys;
		double mTPxSpareSys;
		double mOpSpareSys;

		mTPx = parashuteDensity[ParashuteType::DRAG] * (2 * (totalMass1 * g0-dragForce)) / (cP[ParashuteType::DRAG] * planet.getDensityByHeight(hz) * vel * vel);
		mTP = parashuteDensity[ParashuteType::DRAG] * (2 * (totalMass1 * g0 - dragForce)) / (cP[ParashuteType::DRAG] * maxQForParashute[ParashuteType::MAIN]);
		mOP = parashuteDensity[ParashuteType::MAIN] * (2 * (totalMass1 * g0 - dragForce)) / (cP[ParashuteType::MAIN] * planet.getDensityByHeight(0) * vel * vel);


		mTPxSys = 0.0018 * totalMass1 + mTPx * 1.09;
		mTPSys = 0.0018 * totalMass1 + mTP * 1.09;
		mOpSys = 12.52 / 10000.0 * overLoad * pow(totalMass1, 1.5) / vel + mOP * 1.05;

		mTPxSpareSys = mTPxSys * 0.64;
		mOpSpareSys = (mTPSys + mOpSys) * 0.64;

		
			if (mTPxSys < mTPSys + mOpSys) {
				onlyDragPar = true;
				parSysMass = mTPxSys * 1.64;
				//totalMass += parSysMass;
				parashuteSurface[ParashuteType::DRAG] = (2 * (totalMass1 * g0 - dragForce)) / (cP[ParashuteType::DRAG] * planet.getDensityByHeight(hz) * vel * vel);
				parashuteSurface[ParashuteType::MAIN] = 0;

				parashuteSurface[ParashuteType::SPARE] = parashuteSurface[ParashuteType::DRAG] * 0.64;
				maxQForParashute[ParashuteType::SPARE] = maxQForParashute[ParashuteType::DRAG];
				cP[ParashuteType::SPARE] = cP[ParashuteType::DRAG];
				//printf("\n flag000\n");
			}
			else {
				onlyDragPar = false;
				parSysMass = (mTPSys + mOpSys) * 1.64;
				//totalMass += parSysMass;
				parashuteSurface[ParashuteType::DRAG] = (2 * (totalMass1 * g0 - dragForce)) / (cP[ParashuteType::DRAG] * maxQForParashute[ParashuteType::MAIN]);
				parashuteSurface[ParashuteType::MAIN] = (2 * (totalMass1 * g0 - dragForce)) / (cP[ParashuteType::MAIN] * planet.getDensityByHeight(0) * vel * vel);

				parashuteSurface[ParashuteType::SPARE] = parashuteSurface[ParashuteType::MAIN] * 0.64;
				maxQForParashute[ParashuteType::SPARE] = maxQForParashute[ParashuteType::MAIN];
				cP[ParashuteType::SPARE] = cP[ParashuteType::MAIN];
				//printf("\n flag001\n");
			}
			//isParashute = true;
		


	}

	void findOptimalDescendAngle(Planet& planet) {

		double descendAngleMinOL;
		double lowestMaxOL = 1000000;
		double angle = 0;
		double hIgnBuff = parashuteDEployHight;
		parashuteDEployHight = 0;

		if (planet.isAtmosphere) {

			for (int desAng = -10; desAng <= (int)(angleAttackArray[FILE_RESOLUTION - 1] * 180 / PI); desAng += 1) {

				alphaDescend = desAng / 180.0 * PI;
				initialize(planet);
				while (!(vecLength(velocity) < 100) && position.y > 0)
				{
					control(planet.getDensityByHeight(position.y));
					dynamic(planet);
					//printf("y %.2f\n", position.y);

				}
				//printf("desAng %-5f  maxOL %-5.2f\n", desAng*1.0,maxOverLoad);
				if (maxOverLoad < lowestMaxOL) {
					lowestMaxOL = maxOverLoad;
					descendAngleMinOL = alphaDescend;
				}
			}
		}
		else {
			alphaDescend = 0;
		}
		parashuteDEployHight = hIgnBuff;

		//printf("Descent calculation complete descent angle = %-5.2f maxOL = %-5.2f\n", alphaDescend * 180 / PI, lowestMaxOL);



	}

	sf::Vector2f parashuteDragForce(double density , sf::Vector2f relVelocity ) {

		float vel = vecLength(relVelocity);
		sf::Vector2f dragForce=-normalize(relVelocity);

		if (parashuteDeployed == true) {

			if (parashuteType == ParashuteType::DRAG) {

				dragForce *= (float)(cP[ParashuteType::DRAG] * density * vel * vel * parashuteSurface[ParashuteType::DRAG] / 2);
				//printf("\nflag0 cP %-5.2f den %-5.2f vel %-5.2f surf %-5.2f",
					//cP[ParashuteType::DRAG], density, vel, parashuteSurface[ParashuteType::DRAG]);
			}
			else if (parashuteType == ParashuteType::MAIN) {

				dragForce *= (float)(cP[ParashuteType::MAIN] * density * vel * vel * parashuteSurface[ParashuteType::MAIN] / 2);
				//printf("\nflag1");
			}
			else if (parashuteType == ParashuteType::SPARE) {

				dragForce *= (float)(cP[ParashuteType::SPARE] * density * vel * vel * parashuteSurface[ParashuteType::SPARE] / 2);
				//printf("\nflag2");
			}
		}
		else {
			//printf("\nfock\n");
			return sf::Vector2f{ 0,0 };
		}

		return dragForce;

	}

	sf::Vector2f getThrustForce(float g0) {

		float n = LANDING_OVERLOAD;
		float specImp = THRUST_SPEC_IMP;
		
		

		if (mfuel > 0 && thrustIsOn) {
			mfuel -= thrust / specImp * dt;
			totalMass-= thrust / specImp * dt;
			return (-normalize(velocity) * thrust);
		}
		else {
			return sf::Vector2f{ 0,0 };
		}

	}

	float getIgnitionHight(double g0 , float v0 , float vk) {

		float underLog = LANDING_OVERLOAD / (LANDING_OVERLOAD - 1 + vk * vk / (v0 * v0));
		float h = v0 * v0 / (2.0 * g0) * log(underLog);
		return h;

	}

	float getThrusterMass(float v0, float mass) {
		if (v0 != 0) {
			return 1.19 / 1000.0 * v0 * mass + 7.8;
		}
		else {
			return (0);
		}
	}


};

