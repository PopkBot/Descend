#pragma once
//Импорт -начало-

#include <math.h>

#include "SFML/Graphics.hpp"
#include "excelExporter.h"
#include "Planet.h"
#include "DescendVehicle.h"

//Импорт -конец-

//Объявление констант -начало-

#define SC_LRE_PARAMS_COUNT 11
#define SC_LRE_MATRIX_PARAPMS 2
#define SC_DIRECTORY "Space_Crafts\\"
#define SC_LRE_MATRIX_PARAMS_DIR "Space_Crafts\\SC_matrix_parametrs\\"
#define PI 3.1415926535897932384
#define H_EPS 0.1
#define MASS_EPS 3
#define V_EPS 0.5
#define H_IGNITION_INITIAL 0
#define G_EARTH 9.81
#define MAX_AOA 60

//Объявление констант -конец-


class SC_LRE : public DescendVehicle
{

public:
	//Объявление переменных -начало-

	ExcelExporter excelExporter;

	float totalMass;
	float massSC;
	double mt;
	double mtDebt;
	double fuelMass;
	float mdu;
	float specImpuls=2800;
	float massCoefsPropelSys[3];
	//double atoCoeff;
	float massCoefsEng[3];
	double gammaCoef;
	float oxCoef;
	float overLoadCoef=0;
	float cX=0;
	float cY=0;
	float midSurfase=0;
	
	float alpha=0;
	float gamma=0;
	float alphaDescend =  -1.0 / 180.0 * PI;
	float engSweepAngle=10/180*PI;
	float thrust=0;
	double dt = 0.01;
	double time = 0;
	double hz=20;
	double lowThrustCoef=0;
	double overLoad;
	double maxOverLoad=0;


	float engOn = 0.0f;
	float hIgnition = 0;

	bool bcolibration = false;
	//bool landed = false;

	
	
	sf::Vector2f thrustVector = sf::Vector2f{ 0,0 };


	double machArray[FILE_RESOLUTION];
	double angleAttackArray[FILE_RESOLUTION];
	double cXMachArray[FILE_RESOLUTION][FILE_RESOLUTION];
	double adRatioAngleArray[FILE_RESOLUTION][FILE_RESOLUTION];

	
	

	double params[SC_LRE_PARAMS_COUNT];

	//sf::Vector2f pressureCenter;




	//Объявление переменных -конец-


	//Методы -начало-

	SC_LRE(string fileName, string matrixParFileName) {

		printf("\LRE vehicle\n");
		double dataArray[SC_LRE_MATRIX_PARAPMS][FILE_RESOLUTION];

		ExcelExporter::extractMatrixFromFile_ForSpaceCraft(SC_LRE_MATRIX_PARAMS_DIR, matrixParFileName, dataArray, SC_LRE_MATRIX_PARAPMS, cXMachArray, adRatioAngleArray);

		for (int i = 0; i < FILE_RESOLUTION; i++) {
			machArray[i] = dataArray[0][i];
			angleAttackArray[i] = dataArray[1][i];
			
		}
		printMatrixParams();
		
		ExcelExporter::extractDataFromFile(SC_DIRECTORY, fileName, params, SC_LRE_PARAMS_COUNT);

		massSC = params[0];						//Масса КА
		totalMass = massSC;
		cX = params[1];							//сХ
		specImpuls = params[2];							//сУ
		midSurfase = params[3];					//Площадь Миделя
		phi = params[5] * PI / 180;				//Начальный угол входа 
		velocity.x = params[4] * sin(phi);		//Vx
		velocity.y = -params[4] * cos(phi);		//Vy
		relVelocity = velocity;
		position.x = params[6];					//X
		position.y = params[7];					//Y
		lowThrustCoef = params[9];				//Коэффициент минимальной тяги
		overLoadCoef = params[10];				//Коэффициент перегрузки
		//engSweepAngle = params[11];				//Предельный угол отклонения тяги
		if (params[8] == 1) {					//АТ НДМГ
			massCoefsPropelSys[0] = 0.043;
			massCoefsPropelSys[1] = 1.6;
			massCoefsPropelSys[2] = -0.11;

			massCoefsEng[0] = 0.021;
			massCoefsEng[1] = 4.1;
			massCoefsEng[2] = -0.019;
		}
		else if (params[8] == 2) {				//О2 Керосин
			massCoefsPropelSys[0] = 0.054;
			massCoefsPropelSys[1] = 1.6;
			massCoefsPropelSys[2] = -0.12;

			massCoefsEng[0] = 0.029;
			massCoefsEng[1] = 4.1;
			massCoefsEng[2] = -0.02;
		}
		else if (params[8] == 3) {				//О2 Н2
			massCoefsPropelSys[0] = 0.105;
			massCoefsPropelSys[1] = 1.6;
			massCoefsPropelSys[2] = -0.28;

			massCoefsEng[0] = 0.49;
			massCoefsEng[1] = 4.1;
			massCoefsEng[2] = -0.023;
		}
		else {									//О2 СПГ
			massCoefsPropelSys[0] = 0.062;
			massCoefsPropelSys[1] = 1.6;
			massCoefsPropelSys[2] = -0.15;

			massCoefsEng[0] = 0.032;
			massCoefsEng[1] = 4.1;
			massCoefsEng[2] = -0.021;
		}

		

		//pressureCenter.x = params[8];
		//pressureCenter.y = params[9];

		printf("\nSpaceCraft LRE parametrs:\nmassSC %2f\ncX %3f\ncY %3f\nVx %4f\nVy %4f\nphi %3f\nX %3f\nY %3f\nlowThustCoeff %2f\nfuel %2f\n",
			massSC, cX, cY, velocity.x, velocity.y, phi * 180 / 3.1415, position.x, position.y,lowThrustCoef,params[8]);

	}

	void initialize(Planet &planet)override {

		
		//totalMass = massSC;
		midSurfase = params[3];					//Площадь Миделя
		phi = params[5] * PI / 180;				//Начальный угол входа 
		velocity.x = params[4] * sin(phi);		//Vx
		velocity.y = -params[4] * cos(phi);		//Vy
		relVelocity = velocity;
		position.x = params[6];					//X
		position.y = params[7];					//Y
		lowThrustCoef = params[9];				//Коэффициент минимальной тяги
		overLoadCoef = params[10];				//Коэффициент перегрузки
		dt = 0.01;
		time = 0;
		engOn = 0;
		maxOverLoad = 0;
		mt = fuelMass;
		mtDebt = 0;

		double g0 = (planet.gravPar / planet.radius) / planet.radius;
		totalMass = findTotalMass(fuelMass, planet, MASS_EPS);
		thrust = totalMass * G_EARTH * overLoadCoef;
		landed = false;

		//printf("total mass %2f\tg0 %.2f\tthrust %.2f\n", totalMass, g0, thrust);
	}

	void printStats()override {
		printf("t %-10.3f  y %-10.5f  v %-10.2f  phi %-5.2f  mt %-10.3f\n",
			time,
			position.y,
			vecLength(velocity),
			phi * 180 / PI,
			mt);
			
	}

	void printFinalStats()
		override {

		printf("y %-10.2f  v %-10.2f  phi %-5.2f  mt %-10.3f\n",
			position.y,
			vecLength(velocity),
			phi * 180 / PI,
			mt);

		printf("x %2f\ty %2f\tv %3f\tVx %2f\tVy %4f\tphi %2.2f\t fuel mass %2.2f\n", position.x, position.y,
			vecLength(velocity), velocity.x,
			velocity.y, phi * 180 / PI, fuelMass);
		printf("max Overload %-5.2f", maxOverLoad);

	}

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

	void dynamic(Planet& planet)override {

		
		sf::Vector2f acselGlobal;
		sf::Vector2f gravityForce;
		sf::Vector2f centripicalAcs = sf::Vector2f{ 0,0 };
		float windSpeed = planet.getWindSpeed(position.y);
		

		float dAlpha = windSpeed * cos(phi) / vecLength(velocity);


		relVelocity = velocity - sf::Vector2f{ windSpeed,0 };


		float cxCoef = cxCoefficient(vecLength(relVelocity), alpha+ dAlpha, position.y, planet);
		float cyCoef = cyCoefficient(vecLength(relVelocity), alpha+ dAlpha, position.y, planet, cxCoef);

		float qForce = cxCoef * pow(vecLength(relVelocity), 2)* planet.getDensityByHeight(position.y) / 2 * midSurfase;
		float yForce = cyCoef * pow(vecLength(relVelocity), 2)* planet.getDensityByHeight(position.y) / 2 * midSurfase;

		
		//windForceLocal.x = qForce * cos(alpha) - yForce * sin(alpha);
		//windForceLocal.y = qForce * sin(alpha) + yForce * cos(alpha);

		//windForceGlobal.y = windForceLocal.x * cos(phi) + windForceLocal.y * sin(phi);
		//windForceGlobal.x = -windForceLocal.x * sin(phi) + windForceLocal.y * cos(phi);

		windForce= -normalize(relVelocity) * qForce + turnVector(normalize(relVelocity), PI / 2) * yForce;
		
		//windForceLocal = -normalize(sf::Vector2f{ velocity.x - (float)planet.windSpeed,velocity.y })*(float)qForce
		//				+(float)yForce*normalize(turnVector(sf::Vector2f{ velocity.x - (float)planet.windSpeed,velocity.y }, 0.5 * PI));

		

		//printf("\n WF local %.2f\tWF global %.2f", vecLength(windForceLocal), vecLength(windForceGlobal));
		//printf("\nproj WF local %.2f\tproj WF global %.2f\talpha %.2f", dot(windForceLocal, sf::Vector2f{0,1}), dot(windForceGlobal, sf::Vector2f{ 0,1 }),alpha);

		//windForceGlobal = turnVector(windForceLocal,- PI / 2 + phi);

		gravityForce.y = -planet.gravPar * totalMass / pow(planet.radius + position.y, 2);
		gravityForce.x = 0;
		
		
		//thrustForceGlobal = -normalize(velocity) *(float) thrustForce(position.y,hz,thrust,thrust*lowThrustCoef) *engOn;

		centripicalAcs.y = pow(velocity.x, 2) / (planet.radius + position.y);

		acselGlobal = (gravityForce + windForce  + thrustVector*engOn) / totalMass+centripicalAcs;
		//printf("\nWForce %.5f\tacs %.5f\tthru %.5f\n", vecLength(windForce),vecLength(acselGlobal),vecLength(thrustVector));
		overLoad = vecLength(acselGlobal) / G_EARTH;
		if (overLoad > maxOverLoad) {
			maxOverLoad = overLoad;
		}

		velocity = velocity + acselGlobal * (float)dt;

		position = position + velocity * (float)dt;
		
	//	position.y += planet.radius * (1 - 0.5 * pow(vecLength(velocity) / (planet.radius + position.y), 2)) - 
	//					(position.y+ planet.radius) * (1 - 0.5 * pow(vecLength(velocity) / (planet.radius + position.y), 2));

		if (velocity.x > 0)
			phi = acos(dot(normalize(velocity), sf::Vector2f{ 0,-1 }));
		else if (velocity.x < 0)
			phi = -acos(dot(normalize(velocity), sf::Vector2f{ 0,-1 }));
		else
			phi = 0;


		if (abs(position.y - hIgnition) < 500 || position.y < (hz * 2)) {
			dt = 0.0001;
		}
		else {
			dt = 0.01;
		}


		time += dt;
		/*
		if (position.y < H_EPS && bcolibration) {
			
		}
		else if (position.y < H_EPS/100 && vecLength(velocity) < V_EPS * 10 && !bcolibration) {
			
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
		if (position.y < H_EPS && !bcolibration ) {	//&& abs(velocity.y) < V_EPS
			landed = true;
			totalMass = findTotalMass(fuelMass, planet, MASS_EPS);
		}

		//printf("t %2f\t", time);
	}

	void control(Planet& planet)override {

		

		if (position.y <= hIgnition) {
			engOn = 1;
			if (velocity.y > 0) {
				engOn = 0;
				thrustVector *=0.0f;
			}
			if (engOn == 1) {
				
				if (abs(windForce.x / thrustForce(position.y, hz, thrust, thrust * lowThrustCoef)) <= 1) {
					gamma = asin(windForce.x / thrustForce(position.y, hz, thrust, thrust * lowThrustCoef));
					thrustVector = -normalize(turnVector(velocity, gamma)) * (float)thrustForce(position.y, hz, thrust, thrust * lowThrustCoef);
					
					
				}
				else {
					gamma = 0;
					thrustVector = -normalize(turnVector(velocity, gamma)) * (float)thrustForce(position.y, hz, thrust, thrust * lowThrustCoef);
					//thrustVector =normalize( sf::Vector2f(-windForce.x,0) )* (float)thrustForce(position.y, hz, thrust, thrust * lowThrustCoef);
				}
				

				
			//	thrustVector = -normalize(turnVector(velocity, gamma)) * (float)thrustForce(position.y, hz, thrust, thrust * lowThrustCoef);
				//float angle = acos(dot(normalize(thrustVector), normalize(-velocity)));
				if (abs(gamma) <= engSweepAngle) {
					alpha = 0;
				}
				else {
					alpha = (abs(gamma) - engSweepAngle) * gamma / abs(gamma);
				}
			}
		}
		if (position.y > hIgnition) {
			
			alpha = alphaDescend;
			
			
		}



		if (engOn == 1) {
			//fuelBurn();
			double dm = thrustForce(position.y, hz, thrust, thrust * lowThrustCoef) / specImpuls*dt;
			mt -= dm;
			if (mt <= 0) {
				if (bcolibration) {
					mt = 0;
					mtDebt += dm;
					dm = 0;
					
				}
				else {
					mt = 0;
					dm = 0;
					engOn = 0;
				}

			}
			totalMass -= dm;
		}
		
		//printf("\nalpha %.5f\tad %.5f", alpha, alphaDescend);


	}

	void fuelBurn() {

		double dm= thrustForce(position.y, hz, thrust, thrust * lowThrustCoef)/specImpuls;

		mt -= dm;
		totalMass -= dm;

	}

	double findPropelMass(float vv,float imp,double eps) {

		double ee = exp(-vv / imp);
		double mm = 0;
		double delta;
		delta = mm * (ee - (1 - ee) * aCoef(mm)) -
			(1 - ee) * (massSC + mdu);

		while (delta > eps) {
			mm += 0.1;
			delta = mm * (ee - (1 - ee) * (massCoefsPropelSys[0] * (1 + massCoefsPropelSys[2] * exp(massCoefsPropelSys[3] * mm / 1000)))) -
				(1 - ee) * (massSC + mdu);
		}
		return mm;

	}

	double findTotalMass(double propellantMass,Planet &planet,double eps) {
		
		double totMass= massSC;
		double deltaM = 0.1;
		double g0= (planet.gravPar / planet.radius) / planet.radius;

		double delta = massSC + propellantMass * (aCoef(propellantMass) + 1) + totMass*(gammaCoefSC(thrust) * overLoadCoef - 1);

		while (delta > eps) {
			totMass += deltaM;
			delta = massSC + propellantMass * (aCoef(propellantMass) + 1) + totMass * (gammaCoefSC(totMass*G_EARTH * overLoadCoef) * overLoadCoef - 1);
			//printf("delta %-10.2f\ttotM %-10.2f\tfuelSys %-10.2f\tpropS %-10.2f\n", delta, totMass, propellantMass * (aCoef(propellantMass) + 1), gammaCoefSC(totMass * g0 * overLoadCoef) * overLoadCoef*totMass);
		}
		if (propellantMass == 0) {
			totMass = massSC;
		}

		return totMass;
	}

	void calculateOptimalMass(Planet& planet)override {

		
		double h;
		double mtt0;
		double mtt;
		double deltaMT;
		double v0;
		double v;
		double dv=0;
		double phi0;
		double h0;
		double eps=0.1;
		double epsH = 50;
		double dh=0;
		double dphi=-0.0001;
		double delta;
		double g0;
		
		double dV = V_EPS * 10;
		bcolibration = true;
		//dt = 0.001;

		g0= (planet.gravPar / planet.radius) / planet.radius;
		
		fuelMass = 0;
		//totalMass = massSC;// +thrust / 1000 * gammaCoefSC(thrust) / g0;
		//g0 = (planet.gravPar / planet.radius) / planet.radius;
		//totalMass = findTotalMass(fuelMass, planet, MASS_EPS);
		//thrust = totalMass * g0 * overLoadCoef;

		
		hIgnition = 0;// params[7] / 2;
		initialize(planet);
		printf("total mass %2f\tg0 %.2f\tthrust %.2f\n", totalMass, g0, thrust);
			

		/*

			//while (position.y>h0)
			//{
			//	dynamic(planet);
			//}

			
			//v = vecLength(velocity);
			//phi0 = phi;
			//printf("flag0 h0 %.2f\tphi0 %.5f\tv %.2f\th %.5f\n", h0, phi0 / PI * 180, v, h);
			//h = h0;
			//if (phi0 <= eps / 100) {
			//	phi0 = eps / 100;
			//}
			
			while(abs(v)>eps){


				//dv = ((liftForce(planet, h, v, phi0) * cos(phi0) - thrustForce(h, hz, thrust, thrust * lowThrustCoef) - dragForce(planet, h, v, phi0) * sin(phi0)) / totalMass *
				//	v / (-g0 * sin(phi0) + liftForce(planet, h, v, phi0) / totalMass) - v * cos(phi0)) / sin(phi0) * dphi;

				dv = ((liftForce(planet, h, v, phi0) * cos(phi0) - thrustForce(h, hz, thrust, thrust * lowThrustCoef)*sin(phi0) - dragForce(planet, h, v, phi0) * sin(phi0)) / totalMass *
					v / (-g0 * sin(phi0) + liftForce(planet, h, v, phi0) / totalMass) - v * cos(phi0)) / sin(phi0) * dphi;

				dh = -dv * totalMass * v * cos(phi0) / (totalMass * g0 * cos(phi0) - dragForce(planet, h, v, phi0) - thrustForce(h, hz, thrust, thrust * lowThrustCoef));

				//dh = -(thrustForce(h, hz, thrust, thrust * lowThrustCoef) * cos(phi0) + dragForce(planet, h, v, phi0) * cos(phi0) - totalMass * g0) /
				//	(totalMass * g0 * cos(phi0) - dragForce(planet, h, v, phi0) - thrustForce(h, hz, thrust, thrust * lowThrustCoef) + liftForce(planet, h, v, phi0) * sin(phi0)) * dv;

				//printf("flag1 h0 %.2f\tphi0 %.5f\th %.2f\tv %.2f\tdv %.2f\tdh %.2f\toLoad %.2f\n",h0, phi0 / PI * 180, h, v, dv,dh,thrustForce(h, hz, thrust, thrust * lowThrustCoef)/ (totalMass * g0));
				//printf("drag %.2f\tlift %.2f\tthrust %.2f\tmg %.2f\n", dragForce(planet, h, v, phi0), liftForce(planet, h, v, phi0), thrustForce(h, hz, thrust, thrust * lowThrustCoef), totalMass * g0);
				/*
				dv = -1000;
				delta = (v * cos(phi0) + sin(phi0) * abs(dv / dphi)) * (g0 * sin(phi0) - liftForce(planet, h, v, phi0) / totalMass) / v -
					(liftForce(planet, h, v, phi0) * cos(phi0) - thrustForce(h, hz, thrust, thrust * lowThrustCoef) * sin(phi0) - dragForce(planet, h, v, phi0) * sin(phi0)) / totalMass;


				
				while (abs(delta) > eps) {

					dv += 0.1;
					delta = (v * cos(phi0) + sin(phi0) * (dv / dphi)) * (g0 * sin(phi0) - liftForce(planet, h, v, phi0) / totalMass) / v -
						(liftForce(planet, h, v, phi0) * cos(phi0) - thrustForce(h, hz, thrust, thrust * lowThrustCoef) * sin(phi0) - dragForce(planet, h, v, phi0) * sin(phi0)) / totalMass;

					printf("flag1 delta %.5f\tphi0 %.2f\th %.2f\tv %.2f\tdv %.1f\n",delta,phi0,h,v,dv);
				}
				

				

				dh = -1000;

				delta = (thrustForce(h, hz, thrust, thrust * lowThrustCoef) * cos(phi0) + dragForce(planet, h, v, phi0) * cos(phi0) - totalMass * g0) /
					(totalMass * g0 * cos(phi0) - dragForce(planet, h, v, phi0) - thrustForce(h, hz, thrust, thrust * lowThrustCoef) + liftForce(planet, h, v, phi0) * sin(phi0)) * dv - dh;
				while (abs(delta) > eps) {

					dh += 0.1;
					delta = (thrustForce(h, hz, thrust, thrust * lowThrustCoef) * cos(phi0) + dragForce(planet, h, v, phi0) * cos(phi0) - totalMass * g0) /
						(totalMass * g0 * cos(phi0) - dragForce(planet, h, v, phi0) - thrustForce(h, hz, thrust, thrust * lowThrustCoef) + liftForce(planet, h, v, phi0) * sin(phi0)) * dv - dh;
					//printf("flag2 delta %.5f\tphi0 %.2f\th %.2f\tv %.2f\tdv %.1f\tdh %.2f\n", delta, phi0, h, v, dv,dh);
					//printf("flag2 delta %.5f\tdh %.5f\n", delta,dh);
				}
			

			
				h += dh;
				v += dv;
				phi0 += dphi;
				if (phi0 <= eps / 100) {
					phi0 = eps/100;

				}
				

			}

			if (abs(h) > epsH) {
				h0 -= (h) / 100;
			
			}
			


		}
		printf("flag0 h0 %.2f\tphi0 %.5f\tv %.2f\th %.5f\n", h0, phi0 / PI * 180, v, h);
		
		
		//printf("\n\nIgnition height caclulation completed hIgnition = %.3f\n\n", hIgnition);

	
		engOn = 0;
		massSC = params[0];						//Масса КА
		phi = params[5] * PI / 180;				//Начальный угол входа 
		velocity.x = params[4] * sin(phi);		//Vx
		velocity.y = -params[4] * cos(phi);		//Vy
		position.x = params[6];					//X
		position.y = params[7];

		*/



		

			double dmt=0;

			while (true) {


				fuelMass += dmt;
				hIgnition = 0;
				initialize(planet);



				findOptimalDescendAngle(planet);
				findOptimalIgnitionHeight(planet);

				if (mtDebt > 0) {
					dmt = mtDebt;
				}
				else {
					dmt = -mt;
				}
				printf("mt = %-10.3f  mtD = %-10.3f  dmt = %-10.3f  fM = %-10.3f\n", mt, mtDebt, dmt, fuelMass);

				if (abs(dmt) < MASS_EPS) {
					break;
				}

			}

			printf("mt = %-10.3f  mtD = %-10.3f  dmt = %-10.3f  fM = %-10.3f\n", mt, mtDebt, dmt, fuelMass);




		

		findOptimalDescendAngle(planet);
		
		findOptimalIgnitionHeight(planet);

		
		
		initialize(planet);
		bcolibration = false;
		
		printf("\n\nIgnition height caclulation completed hIgnition = %.3f\n\n", hIgnition);
		printf("\nCalculated parametrs:\nhIgnition = %-10.2f\nDescent Angle of attack = %-5.2f\n\n",
			hIgnition,alphaDescend/PI*180);
		




	}

	void findOptimalDescendAngle(Planet &planet) {

		double descendAngleMinOL;
		double lowestMaxOL = 1000000;
		double angle = 0;
		double hIgnBuff = hIgnition;
		hIgnition = 0;

		if (planet.isAtmosphere) {

			for (int desAng = -10; desAng <= (int)MAX_AOA; desAng += 1) {

				alphaDescend = desAng / 180.0 * PI;
				initialize(planet);
				while (!(vecLength(velocity) < 100) && position.y>0)
				{
					control(planet);
					dynamic(planet);
					//printf("y %.2f\n", position.y);

				}
				//printf("decAng %-5f  maxOL %-5.2f\n", desAng*1.0,maxOverLoad);
				if (maxOverLoad < lowestMaxOL) {
					lowestMaxOL = maxOverLoad;
					descendAngleMinOL = alphaDescend;
				}
			}
		}
		else {
			alphaDescend = 0;
		}
		hIgnition = hIgnBuff;

		//printf("Descent calculation complete descent angle = %-5.2f maxOL = %-5.2f\n", alphaDescend * 180 / PI, lowestMaxOL);

		

	}

	void findOptimalIgnitionHeight(Planet &planet) {

		int drawCount = 0;
		while (true)
		{
			initialize(planet);
			//while (abs(velocity.y)>0.1) {
			while (!(position.y < (H_EPS * 2) && abs(velocity.y) < V_EPS)) {

				control(planet);
				dynamic(planet);
				

				if (drawCount >= 50) {

					drawCount = 0;
					//printf("hIgn %.4f\t\tY %.4f\t\tvy %.5f cx %-5.2f  cy %-5.2f\n", hIgnition, position.y, velocity.y, cxCoefficient(1, alpha, position.y, planet),
					//	cyCoefficient(vecLength(relVelocity), alpha, position.y, planet, cxCoefficient(vecLength(relVelocity), alpha, position.y, planet)));


				}
				drawCount++;
			}

			
			if (abs(position.y)  > H_EPS) {
				//if(abs(position.y/10) > 1){
				hIgnition -= position.y / 10;
				//printf("hIgn %.4f\t\tY %.4f\t\tvy %.5f\n", hIgnition, position.y, velocity.y);
			}
			else {
				//printf("FFFFFFFFhIgn %.4f\t\tY %.4f\t\tvy %.5f\n", hIgnition, position.y, velocity.y);
				break;
			}
			
			/*
			if (abs(position.y )< H_EPS ) {//position.y < H_EPS && position.y >= -H_EPS/10
				printf("FFFFFFFFhIgn %.4f\t\tY %.4f\t\tvy %.5f\n", hIgnition, position.y, velocity.y);
				break;
			}
			else {
				if (abs(position.y / 10) < 0.01) {
					hIgnition -= position.y / abs(position.y) * 0.01;
				}
				else {
					hIgnition -= position.y / 10 ;
				}
				printf("hIgn %.4f\t\tY %.4f\t\tvy %.5f\n", hIgnition, position.y, velocity.y);
			}
			*/


		}

	}

	double thrustForce(double h, double hzz, double p0, double pk) {

		if (h > hzz) {
			return p0;
		}
		else if (h < 0) {
			return p0;
		}
		else {
			return ((p0 - pk) / hzz * h + pk);
			//return pk;
		}

	}

	double cxCoefficient(double v,double angleA,double h,Planet &planet) {

		double cxCoef;
		int cxAIndex= FILE_RESOLUTION - 1;
		int cxMIndex= FILE_RESOLUTION - 1;
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

	double cyCoefficient(double v, double angleA, double h, Planet &planet,double cxC) {

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
			LDCoef =  adRatioAngleArray[cxAIndex][cxMIndex];
		}
		//printf("\n cyccc %.5f\tcxC %0.2f\tld %.2f\tang %.2f\n", cxC * ld,cxC,ld,angleA);
		return (LDCoef* cxC);


	}









	double gammaCoefSC(float tt) {

		return (massCoefsEng[0] * (1 + massCoefsEng[1] * exp(massCoefsEng[2] * tt / 1000)));
	}

	double aCoef(double mm) {
		return(massCoefsPropelSys[0] * (1 + massCoefsPropelSys[1] * exp(massCoefsPropelSys[2] * mm / 1000)));
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


	//Методы -конец-





};

