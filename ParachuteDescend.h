#pragma once

#include "excelExporter.h"
#include "DescendVehicle.h"
#include <math.h>


#define PARACHUTE_PARAMS_COUNT 13
#define PARACHUTE_MATRIX_PARAPMS 2

#define H_EPS 0.1
#define MASS_EPS 1
#define V_EPS 0.5
#define V_DECEND 10
#define V_LANDING 0.1
#define G0 9.81
#define MAX_AOA 60

#define LANDING_OVERLOAD 6
#define THRUST_SPEC_IMP 2300


#define PARACHUTE_DIRECTORY "Space_Crafts\\"
#define PARACHUTE_MATRIX_PARAMS_DIR "Space_Crafts\\SC_matrix_parametrs\\"

enum ParashuteType {
	DRAG = 0,
	MAIN = 1,
	SPARE = 2
};


class ParachuteDescend : public DescendVehicle 
{

public:

	ExcelExporter excelExporter;
	double params[PARACHUTE_PARAMS_COUNT];

	

	double parachuteDensity;
	double parashuteDEployHight=0;
	double parSysMass = 0;
	double hz = 100;
	double g0;
	double currentQ = 0;

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
	double maxQForParashute[3] = { 25060,			// 0.66 * 72900 / 2.0
								3000,				// 0.8*2500/2.0
								4000 };				// 0.9 * 2500 / 2.0

	bool parashuteDeployed = false;
	bool mainParashuteAvailable = true;
	bool onlyDragPar = false;
	bool thrustIsOn = false;
	bool isParashute =true;
	bool isDragCut = false;

	
	

	ParashuteType parashuteType;

	ParachuteDescend(string fileName, string matrixParFileName) {

		printf("\nPar vehicle\n");
		double dataArray[PARACHUTE_MATRIX_PARAPMS][FILE_RESOLUTION];

		ExcelExporter:: extractMatrixFromFile_ForSpaceCraft(PARACHUTE_MATRIX_PARAMS_DIR, matrixParFileName, dataArray, PARACHUTE_MATRIX_PARAPMS, cXMachArray, adRatioAngleArray);

		for (int i = 0; i < FILE_RESOLUTION; i++) {
			machArray[i] = dataArray[0][i];
			angleAttackArray[i] = dataArray[1][i];
		}

		ExcelExporter::extractDataFromFile(PARACHUTE_DIRECTORY, fileName, params, PARACHUTE_PARAMS_COUNT);

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
		printMatrixParams();

	}

	

	void initialize(Planet& planet)
		override
	{


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
		isDragCut = false;
		overLoad = 6;
		if (!bcolibration) {
			totalMass = massSC + parSysMass + thrustMass;
		}
	}

	void printStats()
		override 
	{
		printf("t %-10.3f  y %-10.5f  v %-10.2f  phi %-5.2f  mt %-10.3f\tq %-10.3f deployed %-5d parType %2d\n",
			time,
			position.y,
			vecLength(velocity),
			phi * 180 / PI,
			mfuel,
			currentQ,
			parashuteDeployed,
			parashuteType);
	}

	void printFinalStats()
		override
	{
		printf("y %-10.2f  v %-10.2f  phi %-5.2f  mt %-10.3f\n",
			position.y,
			vecLength(velocity),
			phi * 180 / PI,
			0);

		printf("x %2f\ty %2f\tv %3f\tVx %2f\tVy %4f\tphi %2.2f\t fuel mass %2.2f\n", position.x, position.y,
			vecLength(velocity), velocity.x,
			velocity.y, phi * 180 / PI, mfuel);
		printf("max Overload %-5.2f", maxOverLoad);
	}

	void control(Planet& planet)
		override {

		double currentQ = planet.getDensityByHeight(position.y) * pow(vecLength(velocity), 2) / 2.0;

		parashuteDeployed = false;
		//printf("\nflag00");
		if (position.y > parashuteDEployHight) {
			alpha = alphaDescend;
		}
		else {
			alpha = 0;
		}
		if (position.y < maxQH && isParashute) {
			if (/*position.y < parashuteDEployHight &&*/ maxQForParashute[ParashuteType::DRAG] >= currentQ) {
				parashuteDeployed = true;
				parashuteType = ParashuteType::DRAG;
				//printf(" \ny %-5.2f current q %10.2f maxQ %10.2f",position.y, currentQ, maxQForParashute[ParashuteType::DRAG]);
			}
			if (position.y < parashuteDEployHight && maxQForParashute[ParashuteType::MAIN] >= currentQ && !onlyDragPar) {
				parashuteDeployed = true;
				parashuteType = ParashuteType::MAIN;
				if (!isDragCut) {
					isDragCut = true;
					totalMass -= parashuteSurface[ParashuteType::DRAG] * parashuteDensity[ParashuteType::DRAG];
				}
				//printf("\nflag01 current q %10.2f maxQ %10.2f", currentQ, maxQForParashute[ParashuteType::MAIN]);
			}
			//if (!mainParashuteAvailable && position.y < parashuteDEployHight && maxQForParashute[ParashuteType::SPARE] >= currentQ) {
			//	parashuteDeployed = true;
			//	parashuteType = ParashuteType::SPARE;
				//printf("\nflag02");
			//}
			//if (parashuteDeployed && parashuteType == ParashuteType::MAIN && maxQForParashute[ParashuteType::MAIN] < currentQ) {
			//	mainParashuteAvailable = false;
				//printf("\nflag03");
			//}
		}
		if (position.y <= ignitionHight) {

			thrustIsOn = true;

		}

		//printf("\t depl %2d", parashuteDeployed);


	}

	void dynamic(Planet& planet)override
	{


		sf::Vector2f acselGlobal;
		sf::Vector2f gravityForce;
		sf::Vector2f centripicalAcs = sf::Vector2f{ 0,0 };
		float windSpeed = planet.getWindSpeed(position.y);

		float g0 = (planet.gravPar / planet.radius) / planet.radius;


		float dAlpha = windSpeed * cos(phi) / vecLength(velocity);

		currentQ = pow(vecLength(velocity), 2) * planet.getDensityByHeight(position.y) * 0.5;
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

		if (position.y < H_EPS/10.0  ) {	//&& abs(velocity.y) - V_LANDING < V_EPS
			landed = true;
		}

		//printf("t %2f\t", time);
	
	}


	



	void calculateOptimalMass(Planet& planet)
		override
	{

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

				control(planet);
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
			
			if (V_DECEND == 0) {

				for (float parDecVel = V_LANDING; parDecVel <= freeFallVelocity; parDecVel += 0.1) {

					initialize(planet);


					cxCoef = cxCoefficient(parDecVel, 0, 0, planet);
					cyCoef = cyCoefficient(parDecVel, 0, 0, planet, cxCoef);

					qForce = cxCoef * pow(parDecVel, 2) * planet.getDensityByHeight(0) / 2 * midSurfase;
					yForce = cyCoef * pow(parDecVel, 2) * planet.getDensityByHeight(0) / 2 * midSurfase;

					windForce = -normalize(sf::Vector2f{ 0,-parDecVel }) * qForce + turnVector(normalize(sf::Vector2f{ 0,-parDecVel }), PI / 2) * yForce;

					float thrMassBuf = getThrusterMass(parDecVel - V_LANDING, massSC);
					totalMass += getThrusterMass(parDecVel - V_LANDING, massSC);




					calculateParashute(windForce.y, parDecVel, planet, totalMass);
					totalMass += parSysMass;

					totalMass -= thrMassBuf;
					totalMass += getThrusterMass(parDecVel - V_LANDING, totalMass);

					if (totalMass < minMass) {
						minMass = totalMass;
						optimalParVel = parDecVel;
					}
					printf("\nparVel %-5.2f totalMass %-10.2f parSysMass %-10.2f thrMass %-5.2f", parDecVel, totalMass, parSysMass, getThrusterMass(parDecVel - V_LANDING, totalMass));

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
			}
			else {
				initialize(planet);

				minMass = totalMass;
				optimalParVel = V_DECEND;

				float thrMassBuf = getThrusterMass(V_DECEND - V_LANDING, massSC);
				totalMass += getThrusterMass(V_DECEND - V_LANDING, massSC);
				cxCoef = cxCoefficient(V_DECEND, 0, 0, planet);
				cyCoef = cyCoefficient(V_DECEND, 0, 0, planet, cxCoef);
				qForce = cxCoef * pow(V_DECEND, 2) * planet.getDensityByHeight(0) / 2 * midSurfase;
				yForce = cyCoef * pow(V_DECEND, 2) * planet.getDensityByHeight(0) / 2 * midSurfase;
				windForce = -normalize(sf::Vector2f{ 0,-V_DECEND }) * qForce + turnVector(normalize(sf::Vector2f{ 0,-V_DECEND }), PI / 2) * yForce;
				calculateParashute(windForce.y, V_DECEND, planet, totalMass);
				totalMass += parSysMass;
				totalMass -= thrMassBuf;
				thrustMass = getThrusterMass(V_DECEND - V_LANDING, totalMass);
				totalMass += getThrusterMass(V_DECEND - V_LANDING, totalMass);
				isParashute = true;




			}



			printf("\nMin parVel %-5.2f totalMass %-10.2f parSysMass %-10.2f thrMass %-5.2f", optimalParVel, totalMass, parSysMass, getThrusterMass(optimalParVel - V_LANDING, totalMass));


		
			printf("\n only thrust: thrust mass %5.2f", getThrusterMass(freeFallVelocity, massSC));
			//totalMass = minMass;
			
			

			thrust = totalMass * G0 * LANDING_OVERLOAD;

			//fuelMass = (optimalParVel-V_LANDING) * totalMass / THRUST_SPEC_IMP * LANDING_OVERLOAD / sqrtf(LANDING_OVERLOAD - 1) * atan((1 - V_LANDING / optimalParVel) / (sqrtf(LANDING_OVERLOAD - 1) * (V_LANDING / optimalParVel + 1)));
			fuelMass = (exp((optimalParVel - V_LANDING) / THRUST_SPEC_IMP) - 1) * totalMass;
			mfuel = fuelMass;
			IgnitionVel = optimalParVel;
			ignitionHight = 0;// getIgnitionHight(g0, optimalParVel, V_LANDING);
			findOptimalDescendAngle(planet);
			//printf("\nmTPx %10.2f mTP %10.2f mOP %10.2f\n", mTPx / parashuteDensity[ParashuteType::DRAG], mTP / parashuteDensity[ParashuteType::DRAG], mOP / parashuteDensity[ParashuteType::MAIN]);
			bcolibration = false;
			

			int drawCount = 0;

			









			printf("\nhIgn %.4f\t\tY %.4f\t\tvy %.5f tot mass %-10.2f\n", ignitionHight, position.y, velocity.y, totalMass);

			float vyPrev = freeFallVelocity;

			while (true)
			{
				initialize(planet);


				while (!(position.y < hz && abs(velocity.y) > optimalParVel) && position.y > H_EPS) {

					control(planet);
					dynamic(planet);
					//printf("\ny %-10.2f q %-10.2f", position.y, planet.getDensityByHeight(position.y)* pow(vecLength(velocity), 2) / 2.0);
				}


				if (abs(velocity.y) > optimalParVel*1.2 ) {
					//if(abs(position.y/10) > 1){
					//if (vyPrev < velocity.y) {
						//parashuteSurface[ParashuteType::MAIN] += 1;
					//}
					//else {
						parashuteDEployHight += 1;
					//}
					
					//parashuteSurface[ParashuteType::DRAG] += 1;
					printf("hParDep %.2f\tq %.2f\tvy %.2f\tdragForce %10.2f\tdens %5.2f\tparDepType %d parashuteDeployed %d onlyDrag %d\n",
						parashuteDEployHight, planet.getDensityByHeight(position.y)* pow(vecLength(velocity), 2) / 2.0,
						vecLength(velocity), vecLength(parashuteDragForce(planet.getDensityByHeight(position.y), velocity)),
						planet.getDensityByHeight(position.y), parashuteType, parashuteDeployed, onlyDragPar);
				}
				else {
					//printf("FFFFFFFFhIgn %.4f\t\tY %.4f\t\tvy %.5f\n", hIgnition, position.y, velocity.y);
					bcolibration = false;
					//printf("Fin hParDep %.4f\t\tsurf %.4f\t\tvy %.5f\tdragForce %10.2f\tdens %5.2f\tparDep %d\n", parashuteDEployHight, parashuteSurface[ParashuteType::DRAG], velocity.y, vecLength(parashuteDragForce(planet.getDensityByHeight(0), velocity)), planet.getDensityByHeight(0), parashuteDeployed);

					//printf("\nflag 002\n");
					break;
				}
				vyPrev = velocity.y;
			}
			//float fuelBuff = fuelMass;
			//fuelMass = 2 * fuelMass;
			//mfuel = fuelMass*2;
			

			while (true)
			{
				initialize(planet);
				//while (abs(velocity.y)>0.1) {
				//while (!(position.y < (H_EPS * 2) && abs(velocity.y) - V_LANDING < V_EPS)) {
				while(mfuel>0){
					control(planet);
					dynamic(planet);
					//printf("hIgn %.4f\t\tY %.4f\t\tvy %.5f fuel mass %-10.2f\n", ignitionHight, position.y, velocity.y, mfuel);
					if (drawCount >= 20) {

						drawCount = 0;
						//printf("hIgn %.4f\t\tY %.4f\t\tvy %.5f fuel mass %-10.2f\n", ignitionHight, position.y, velocity.y, mfuel);


					}
					drawCount++;
				}

				if (vecLength(velocity) - V_LANDING > V_EPS) {
					float dm = (exp((vecLength(velocity) - V_LANDING) / THRUST_SPEC_IMP) - 1) * massSC * 0.5;
					fuelMass +=dm+0.1 ;
					thrustMass = getThrusterMass(V_DECEND-V_LANDING, totalMass + fuelMass);
					printf("f mass %-5.2f thr mass %-5.2f\t", fuelMass, thrustMass);
				}
				if (0 < (position.y) < H_EPS) {
					//if(abs(position.y/10) > 1){
					ignitionHight -= position.y / abs(position.y) * 0.1 + position.y * 0.1;
					printf("hIgn %.4f\t\tY %.4f\t\tvy %.5f tot mass %-10.2f\n", ignitionHight, position.y, velocity.y, totalMass);
				}
				else {
					initialize(planet);
					printf("\nFFFFFFFFhIgn %.4f\t\tY %.4f\t\tvy %.5f tot mass %-10.2f\n", ignitionHight, position.y, velocity.y, totalMass);
					break;
				}
			}
			initialize(planet);
			while (!(position.y < (H_EPS * 2) && abs(velocity.y) - V_LANDING < V_EPS  || position.y < 0)) {

				control(planet);
				dynamic(planet);
				//printf("hIgn %.4f\t\tY %.4f\t\tvy %.5f tot mass %-10.2f\n", ignitionHight, position.y, velocity.y, totalMass);
			}
			//fuelMass -= mfuel;
			initialize(planet);
			printf("\noptParVel %-5.2f totalMass %-10.2f parSysMass %-5.2f trhustMass %-5.2f ignH %-5.2f", optimalParVel, totalMass, parSysMass,
				getThrusterMass(optimalParVel - V_LANDING, totalMass), getIgnitionHight(g0, IgnitionVel, V_LANDING));



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

		mTPx = parashuteDensity[ParashuteType::DRAG] * (2 * (totalMass1 * g0 - dragForce)) / (cP[ParashuteType::DRAG] * planet.getDensityByHeight(hz) * vel * vel);
		mTP = parashuteDensity[ParashuteType::DRAG] * (  (totalMass1 * g0 - dragForce)) / (cP[ParashuteType::DRAG] * maxQForParashute[ParashuteType::MAIN]);
		mOP = parashuteDensity[ParashuteType::MAIN] * (2 * (totalMass1 * g0 - dragForce)) / (cP[ParashuteType::MAIN] * planet.getDensityByHeight(hz) * vel * vel);


		mTPxSys = 0.0018 * totalMass1 + mTPx * 1.09;
		mTPSys = 0.0018 * totalMass1 + mTP * 1.09;
		mOpSys = 12.52 / 10000.0 * overLoad * pow(totalMass1, 1.5) / vel + mOP * 1.05;

		mTPxSpareSys = mTPxSys * 0.77;
		mOpSpareSys = (mTPSys + mOpSys) * 0.77;
		//printf("\t\tSTPx %10.2f STP %10.2f SOP %10.2f", mTPx / parashuteDensity[ParashuteType::DRAG], mTP / parashuteDensity[ParashuteType::DRAG], mOP / parashuteDensity[ParashuteType::MAIN]);
		//printf("\t\tmTPx %10.2f mTP %10.2f mOP %10.2f", mTPx, mTP , mOP );
			if (mTPxSys < mTPSys + mOpSys) {
				onlyDragPar = true;
				parSysMass = mTPxSys * 1.77;
				//totalMass += parSysMass;
				parashuteSurface[ParashuteType::DRAG] = (2 * (totalMass1 * g0 - dragForce)) / (cP[ParashuteType::DRAG] * planet.getDensityByHeight(hz) * vel * vel);
				parashuteSurface[ParashuteType::MAIN] = 0;

				parashuteSurface[ParashuteType::SPARE] = parashuteSurface[ParashuteType::DRAG] * 0.77;
				maxQForParashute[ParashuteType::SPARE] = maxQForParashute[ParashuteType::DRAG];
				cP[ParashuteType::SPARE] = cP[ParashuteType::DRAG];
				//printf(" flag000 parSys %-5.2f", mTPxSys * 1.64);
			}
			else {
				onlyDragPar = false;
				parSysMass = (mTPSys + mOpSys) * 1.77;
				//totalMass += parSysMass;
				parashuteSurface[ParashuteType::DRAG] = (2 * (totalMass1 * g0 - dragForce)) / (cP[ParashuteType::DRAG] * maxQForParashute[ParashuteType::MAIN]);
				parashuteSurface[ParashuteType::MAIN] = (2 * (totalMass1 * g0 - dragForce)) / (cP[ParashuteType::MAIN] * planet.getDensityByHeight(hz) * vel * vel);

				parashuteSurface[ParashuteType::SPARE] = parashuteSurface[ParashuteType::MAIN] * 0.77;
				maxQForParashute[ParashuteType::SPARE] = maxQForParashute[ParashuteType::MAIN];
				cP[ParashuteType::SPARE] = cP[ParashuteType::MAIN];
				//printf(" flag001 parSys %-5.2f", (mTPSys + mOpSys) * 1.64);
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

			for (int desAng = -10; desAng <=(int) MAX_AOA; desAng += 1) {

				alphaDescend = desAng / 180.0 * PI;
				initialize(planet);
				while (!(vecLength(velocity) < 100) && position.y > 0)
				{
					control(planet);
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
				//	cP[ParashuteType::DRAG], density, vel, parashuteSurface[ParashuteType::DRAG]);
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

