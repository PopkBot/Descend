#pragma once
/*
Venus atmosphere	http://marsmeta.narod.ru/venera.html



*/





//Импорт -начало-

#include "excelExporter.h"
//Импорт -конец-

//Объявление констант -начало-


#define PARAMETRS_COUNT 3
#define ATMOSPHERE_PARAMETRS_COUNT 3
#define PLANETS_DIRECTORY "Planets\\"
#define MATRIX_DIRECTORY "Planets\\Planets_matrix_parametrs\\"
#define ATMOSPHERE_FILE_NAME "AtmoSphere"
#define SONIC_SPEED_FILE "MachNo"

//Объявление констант -конец-

class Planet
{
public:

	//Объявление переменных -начало-

	ExcelExporter excelExporter;
	double atmosphereHight=0;
	double gravPar = 0;
	double radius=1;
	double heightStep[FILE_RESOLUTION];
	double densityStep[FILE_RESOLUTION];
	double sonicSpeedStep[FILE_RESOLUTION];
	float windSpeed = 0;

	bool isAtmosphere = false;

	//Объявление переменных -конец-


	//Методы -начало-

	Planet(string fileName) {

		double dataArray[ATMOSPHERE_PARAMETRS_COUNT][FILE_RESOLUTION];
		excelExporter.extractMatrixFromFile(MATRIX_DIRECTORY,fileName+ ATMOSPHERE_FILE_NAME,dataArray, FILE_RESOLUTION, ATMOSPHERE_PARAMETRS_COUNT);
		
		for (int i = 0; i < FILE_RESOLUTION; i++) {
			heightStep[i] = dataArray[0][i];
			densityStep[i] = dataArray[1][i];
			sonicSpeedStep[i] = dataArray[2][i];
		}
		
		printAtmospereParametrs();
		
		
		atmosphereHight = heightStep[FILE_RESOLUTION - 1];
		double params[PARAMETRS_COUNT];
		excelExporter.extractDataFromFile(PLANETS_DIRECTORY,fileName, params, PARAMETRS_COUNT);
		
		gravPar = params[0]*pow(10,9);
		radius = params[1]*1000;
		if (params[2] > 0) {
			isAtmosphere = true;
		}
		else {
			isAtmosphere = false;
		}

		printf("\nPlanet parametrs:\ngravPar %5f\nradius %2f\n",
			gravPar, radius);

	}

	void printAtmospereParametrs() {
		for (int i = 0; i < FILE_RESOLUTION; i++) {
			
			printf("HS %-20.2f  DS %-20.2f  SSS %-20.2f\n", heightStep[i], densityStep[i], sonicSpeedStep[i]);
		}
	}

	double getDensityByHeight(double height) {


		if(height>=heightStep[0])
			for (int i = 1; i < FILE_RESOLUTION; i++) {
				if (height < heightStep[i]) {
					return (	(densityStep[i]-densityStep[i-1])/(heightStep[i] - heightStep[i - 1])*(height-heightStep[i-1])+densityStep[i-1]	);
				}
			}
		return(0);

	}

	float getWindSpeed(float height) {


		if (height > 10) {
			return windSpeed;
		}
		else {
			return 0;
		}

	}

	float getSonicSpeed(float height) {

		if (height >= heightStep[0]) {
			for (int i = 1; i < FILE_RESOLUTION; i++) {
				if (height < heightStep[i]) {
					return ((sonicSpeedStep[i] - sonicSpeedStep[i - 1]) / (heightStep[i] - heightStep[i - 1]) * (height - heightStep[i - 1]) + sonicSpeedStep[i - 1]);
				}
			}
			return(heightStep[FILE_RESOLUTION-1]);
		}
		else {
			return(heightStep[0]);
		}
		return(heightStep[0]);

	}

	//Методы -конец-



};

