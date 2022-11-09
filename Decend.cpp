// Decend.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

//Импорт -начало-
#include <iostream>
#include <thread>
#include <sstream>


#include "Graph.h"
//#include "excelExporter.h"
#include "Planet.h"
#include "SpaceCraft.h"
#include "SC_LRE.h"
#include "ParachuteDecend.h"

//Импорт -конецо-

//Объявление констант -начало-

#define PLANET "Venus"
#define SPACE_CRAFT "SpaceCraft"
#define SPACE_CRAFT_LRE "Soyuz_LRE"
#define SPACE_CRAFT_PAR "Soyuz_Parashute"

#define DRAW_RATE 20




//Объявление констант -конец-



//Объявление переменных -начало-

Planet planet(PLANET);
SpaceCraft spaceCraft(SPACE_CRAFT);
SC_LRE sc_LRE(SPACE_CRAFT_LRE, SPACE_CRAFT_LRE);
ParachuteDecend sc_PAR(SPACE_CRAFT_PAR, SPACE_CRAFT_LRE);
Graph graph;
ExcelExporter excelExporter;

double simSpeed = 0;


//Объявление переменных -конец-



// Методы -начало-


void graphs() {
	graph.initialize();
}

void drawSpaceCraft(int lineNum,SpaceCraft &sc) {

	char buffer[50];

	std::snprintf(buffer, 50, "%.2f", sc.position.y);
	graph.lines[lineNum].text.setString(buffer);

	std::snprintf(buffer, 50, "%.2f", sc.position.x);
	graph.lines[lineNum].textX.setString(buffer);

	graph.drawline(sc.position.x, sc.position.y, 0);

	graph.drawIcon(lineNum, sc.phi);

	graph.lines[lineNum].bdraw = true;

}

void drawSpaceCraft(int lineNum, SC_LRE &sc) {

	char buffer[50];

	std::snprintf(buffer, 50, "%.2f", sc.position.y);
	graph.lines[lineNum].text.setString(buffer);

	std::snprintf(buffer, 50, "%.2f", sc.position.x);
	graph.lines[lineNum].textX.setString(buffer);

	graph.drawline(sc.position.x, sc.position.y, 0);

	graph.drawIcon(lineNum, sc.phi);

	graph.lines[lineNum].bdraw = true;

}

void drawSpaceCraft(int lineNum, ParachuteDecend& sc) {

	char buffer[50];

	std::snprintf(buffer, 50, "%.2f", sc.position.y);
	graph.lines[lineNum].text.setString(buffer);

	std::snprintf(buffer, 50, "%.2f", sc.position.x);
	graph.lines[lineNum].textX.setString(buffer);

	graph.drawline(sc.position.x, sc.position.y, 0);

	graph.drawIcon(lineNum, sc.phi);

	graph.lines[lineNum].bdraw = true;

}

void printMenu() {
	cout <<"\n\n1 - Определить оптимальную массу\n"
		<< "2 - Запустить симмуляцию\n"
		<< "0 - Выход\n";
}

void menuPos1() {
	//sc_LRE.calculateOptimalMass(planet);
	sc_PAR.calculateOptimalMass(planet);
}

void menuPos2() {
	/*

	int delayCount=0;
	int drawCount = 0;
	graph.lines[0].vertex.clear();
	sc_LRE.sc_Initialize(planet);


	while (!(sc_LRE.position.y < H_EPS * 2 && abs(sc_LRE.velocity.y) < V_EPS) && ! sc_LRE.landed) {

		sc_LRE.control();
		sc_LRE.dynamic(planet);
		
		if (drawCount >= DRAW_RATE) {
			drawSpaceCraft(0, sc_LRE);
			drawCount = 0;
			printf("t %-10.3f  y %-10.5f  v %-10.2f  phi %-5.2f  mt %-10.3f\tq %-10.3f\n",
				sc_LRE.time,
				sc_LRE.position.y,
				sc_LRE.vecLength(sc_LRE.velocity),
				sc_LRE.phi * 180 / PI,
				sc_LRE.mt,
				pow(sc_LRE.vecLength(sc_LRE.velocity),2)*planet.getDensityByHeight(sc_LRE.position.y)*0.5);
			

		}
		drawCount++;
		
	}
	drawSpaceCraft(0, sc_LRE);
	printf("y %-10.2f  v %-10.2f  phi %-5.2f  mt %-10.3f\n",
		sc_LRE.position.y,
		sc_LRE.vecLength(sc_LRE.velocity),
		sc_LRE.phi * 180 / PI,
		sc_LRE.mt);

	printf("x %2f\ty %2f\tv %3f\tVx %2f\tVy %4f\tphi %2.2f\t fuel mass %2.2f\n", sc_LRE.position.x, sc_LRE.position.y,
		sc_LRE.vecLength(sc_LRE.velocity), sc_LRE.velocity.x,
		sc_LRE.velocity.y, sc_LRE.phi * 180 / PI,sc_LRE.fuelMass);
	printf("max Overload %-5.2f", sc_LRE.maxOverLoad);
	*/

	int delayCount = 0;
	int drawCount = 0;
	graph.lines[0].vertex.clear();
	sc_PAR.initialize(planet);


	while (!(sc_PAR.position.y < H_EPS * 2 && abs(sc_PAR.velocity.y) - V_LANDING < V_EPS ) && !sc_PAR.landed) {

		sc_PAR.control(planet.getDensityByHeight(sc_PAR.position.y));
		sc_PAR.dynamic(planet);

		if (drawCount >= DRAW_RATE) {
			drawSpaceCraft(0, sc_PAR);
			drawCount = 0;
			printf("t %-10.3f  y %-10.5f  v %-10.2f  phi %-5.2f  mt %-10.3f\tq %-10.3f deployed %-5d parType %2d\n",
				sc_PAR.time,
				sc_PAR.position.y,
				sc_PAR.vecLength(sc_PAR.velocity),
				sc_PAR.phi * 180 / PI,
				sc_PAR.mfuel,
				pow(sc_PAR.vecLength(sc_PAR.velocity), 2) * planet.getDensityByHeight(sc_PAR.position.y) * 0.5,
				sc_PAR.parashuteDeployed,
				sc_PAR.parashuteType);


		}
		drawCount++;

	}
	drawSpaceCraft(0, sc_PAR);
	printf("y %-10.2f  v %-10.2f  phi %-5.2f  mt %-10.3f\n",
		sc_PAR.position.y,
		sc_PAR.vecLength(sc_PAR.velocity),
		sc_PAR.phi * 180 / PI,
		0);

	printf("x %2f\ty %2f\tv %3f\tVx %2f\tVy %4f\tphi %2.2f\t fuel mass %2.2f\n", sc_PAR.position.x, sc_PAR.position.y,
		sc_PAR.vecLength(sc_PAR.velocity), sc_PAR.velocity.x,
		sc_PAR.velocity.y, sc_PAR.phi * 180 / PI, sc_PAR.mfuel);
	printf("max Overload %-5.2f", sc_PAR.maxOverLoad);
}








// Методы -конец-

int main()
{

	setlocale(LC_ALL, "Russian");
	
	string anyThing;
	string input;
	int command;
	bool runLoop = true;
	
   
	graph.initializeAxis();
	thread thread1(graphs);
	thread1.detach();
	

	
	drawSpaceCraft(0, spaceCraft);
	graph.initializeIcon(0, 3, 5, sf::Color{ 255,0,130 });
	graph.scale_y_max();
	graph.qx = graph.qy;
	cout << "Graph_init_complete" << "\n";

	

	while (runLoop) {

		printMenu();
		
		while (!(cin >> command)) {
			cin.clear();
			cin.ignore(numeric_limits<streamsize>::max(), '\n');
			cout << "Invalid input\n";
		}

		switch (command) {
		case 1: 
			menuPos1();
			break;
		case 2:
			menuPos2();
			break;
		case 0:
			runLoop = false;
			break;
		default:
			break;

		}



	}

	

	

	return 0;

}

