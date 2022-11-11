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
#include "DecendVehicle.h"

//Импорт -конецо-

//Объявление констант -начало-

#define PLANET "Venus"
#define SPACE_CRAFT "SpaceCraft"
#define SPACE_CRAFT_LRE "Soyuz_LRE"
#define SPACE_CRAFT_PAR "Soyuz_Parashute"

#define DRAW_RATE 100




//Объявление констант -конец-



//Объявление переменных -начало-

Planet planet(PLANET);
//SpaceCraft spaceCraft(SPACE_CRAFT);


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

void drawSpaceCraft(int lineNum, DecendVehicle* vehicle) {

	char buffer[50];

	std::snprintf(buffer, 50, "%.2f", vehicle->position.y);
	graph.lines[lineNum].text.setString(buffer);

	std::snprintf(buffer, 50, "%.2f", vehicle->position.x);
	graph.lines[lineNum].textX.setString(buffer);

	graph.drawline(vehicle->position.x, vehicle->position.y, lineNum);

	graph.drawIcon(lineNum, vehicle->phi + vehicle->alpha);

	graph.lines[lineNum].bdraw = true;

}

void drawSpaceCraft( DecendVehicle* vehicle[]) {

	char buffer[50];
	for (int i = 0; i < lineCount; i++) {
		std::snprintf(buffer, 50, "%.2f", vehicle[i]->position.y);
		graph.lines[i].text.setString(buffer);

		std::snprintf(buffer, 50, "%.2f", vehicle[i]->position.x);
		graph.lines[i].textX.setString(buffer);

		graph.drawline(vehicle[i]->position.x, vehicle[i]->position.y, i);

		graph.drawIcon(i, vehicle[i]->phi + vehicle[i]->alpha);

		graph.lines[i].bdraw = true;
	}
}

void printMenu() {
	cout <<"\n\n1 - Определить оптимальную массу\n"
		<< "2 - Запустить симмуляцию\n"
		<< "0 - Выход\n";
}

void menuPos1() {
	
	
}

void menuPos2() {
	
}

void simulateVehicle(DecendVehicle* vehicle,int lineNum) {

	int delayCount = 0;
	int drawCount = 0;
	graph.lines[lineNum].vertex.clear();
	vehicle->initialize(planet);


	while (!vehicle->landed) {

		vehicle->control(planet);
		vehicle->dynamic(planet);

		if (drawCount >= DRAW_RATE) {
			drawSpaceCraft(lineNum, vehicle);
			drawCount = 0;
			

			vehicle->printStats();


		}
		drawCount++;

	}
	drawSpaceCraft(lineNum, vehicle);
	vehicle->printFinalStats();

}

void calculateOptimal(DecendVehicle* vehicle) {

	vehicle->calculateOptimalMass(planet);
}

void initializeVehicleIcon(DecendVehicle* vehicle,int lineNum, int shapeEdges, float radius, sf::Color iconColor, sf::Color lineColor) {

	
	graph.initializeIcon(lineNum, shapeEdges, radius, iconColor,lineColor);
	graph.scale_y_max();
	graph.qx = graph.qy;
	//drawSpaceCraft(lineNum, vehicle);
}







// Методы -конец-

int main()
{

	

	SC_LRE sc_LRE(SPACE_CRAFT_LRE, SPACE_CRAFT_LRE);
	ParachuteDecend sc_PAR(SPACE_CRAFT_PAR, SPACE_CRAFT_LRE);

	DecendVehicle* p[lineCount] = { &sc_LRE,&sc_PAR };

	setlocale(LC_ALL, "Russian");
	
	string anyThing;
	string input;
	int command;
	bool runLoop = true;
	
   
	graph.initializeAxis();
	
	
	initializeVehicleIcon(&sc_PAR, 1, 5, 5, sf::Color::Red, sf::Color::Red);
	
	initializeVehicleIcon(&sc_LRE, 0, 3, 5, sf::Color::Yellow, sf::Color::Yellow);
	thread thread1(graphs);
	thread1.detach();
	
	cout << "Graph_init_complete" << "\n";
	//p = &sc_PAR;
	//simulateVehicle(&sc_PAR);
	//simulateVehicle(&sc_LRE);
	
	while (runLoop) {

		printMenu();


		
		while (!(cin >> command)) {
			cin.clear();
			cin.ignore(numeric_limits<streamsize>::max(), '\n');
			cout << "Invalid input\n";
		}

		switch (command) {
		case 1: 
			calculateOptimal(&sc_LRE);
			calculateOptimal(&sc_PAR);
			break;
		case 2:
			//graph.clearLines();
			
			simulateVehicle(&sc_PAR,1);
			simulateVehicle(&sc_LRE, 0);
			drawSpaceCraft(p);
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

