#pragma once
#include "DescendVehicle.h"
#include "Planet.h"
#include "excelExporter.h"



#define GLIDE_PARAMS_COUNT 13
#define GLIDE_MATRIX_PARAPMS 2

#define H_EPS 0.1
#define MASS_EPS 1
#define V_EPS 0.5
#define V_DECEND 1
#define V_LANDING 0

#define GLIDE_DIRECTORY "Space_Crafts\\"
#define GLIDE_MATRIX_PARAMS_DIR "Space_Crafts\\SC_matrix_parametrs\\"








class GlideDescend :
    public DescendVehicle
{
public:

	ExcelExporter excelExporter;

	double params[GLIDE_PARAMS_COUNT];


	GlideDescend(string fileName, string matrixParFileName) {

		printf("\nGlide vehicle\n");
		double dataArray[GLIDE_PARAMS_COUNT][FILE_RESOLUTION];

		ExcelExporter::extractMatrixFromFile_ForSpaceCraft(GLIDE_MATRIX_PARAMS_DIR, matrixParFileName, dataArray, GLIDE_MATRIX_PARAPMS, cXMachArray, adRatioAngleArray);

		for (int i = 0; i < FILE_RESOLUTION; i++) {
			machArray[i] = dataArray[0][i];
			angleAttackArray[i] = dataArray[1][i];
			
		}

		ExcelExporter::extractDataFromFile(GLIDE_DIRECTORY, fileName, params, GLIDE_PARAMS_COUNT);

		massSC = params[0];						//����� ��
		totalMass = massSC;
		cX = params[1];							//��
		midSurfase = params[2];					//������� ������
		phi = params[4] * PI / 180;				//��������� ���� ����� 
		velocity.x = params[3] * sin(phi);		//Vx
		velocity.y = -params[3] * cos(phi);		//Vy
		relVelocity = velocity;
		position.x = params[5];					//X
		position.y = params[6];					//Y

		//parashuteDensity[ParashuteType::DRAG] = params[7];			//��������� ����� ���������� ��������
		//parashuteDensity[ParashuteType::MAIN] = params[8];			//��������� ����� ��������� ��������
		//parashuteDensity[ParashuteType::SPARE] = params[9];			//��������� ����� ��������� ��������
		//cP[ParashuteType::DRAG] = params[10];	//�� ���������� ��������
		//cP[ParashuteType::MAIN] = params[11];	//�� ��������� ��������
		//cP[ParashuteType::SPARE] = params[12];	//�� ��������� ��������

		//parashuteDeployed = false;
		//mainParashuteAvailable = true;
		//printMatrixParams();

	}

	void initialize(Planet& planet) {}

	void control(Planet& planet) {}

	void dynamic(Planet& planet) {}

	void calculateOptimalMass(Planet& planet) {}

	void printStats() {}

	void printFinalStats() {}



};

