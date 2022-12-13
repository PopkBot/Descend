﻿// Descend.cpp : This file contains the 'main' function. Program execution begins and ends there.
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
#include "ParachuteDescend.h"
#include "DescendVehicle.h"

//Импорт -конецо-

//Объявление констант -начало-

#define PLANET "Earth"
#define SPACE_CRAFT "SpaceCraft"
#define SPACE_CRAFT_LRE "First_Stage_LRE"
#define SPACE_CRAFT_LRE_2 "First_Stage_LRE_2"
#define SPACE_CRAFT_PAR "First_Stage_PAR"

#define DRAW_RATE 50




//Объявление констант -конец-



//Объявление переменных -начало-

Planet planet(PLANET);



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

void drawSpaceCraft(int lineNum, DescendVehicle* vehicle) {

	char buffer[50];

	std::snprintf(buffer, 50, "%.2f", vehicle->position.y);
	graph.lines[lineNum].text.setString(buffer);

	std::snprintf(buffer, 50, "%.2f", vehicle->position.x);
	graph.lines[lineNum].textX.setString(buffer);

	
	graph.drawline(vehicle->position.x, vehicle->position.y, lineNum);
	graph.drawIcon(lineNum, vehicle->phi + vehicle->alpha);
	

	graph.lines[lineNum].bdraw = true;

}

void drawSpaceCraft( DescendVehicle* vehicle[]) {

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

void drawSpaceCraft(int lineNum, SpaceCraft* vehicle) {

	char buffer[50];

	std::snprintf(buffer, 50, "%.2f", vehicle->position.y);
	graph.lines[lineNum].text.setString(buffer);

	std::snprintf(buffer, 50, "%.2f", vehicle->position.x);
	graph.lines[lineNum].textX.setString(buffer);


	graph.drawline(vehicle->position.x, vehicle->position.y, lineNum);
	graph.drawIcon(lineNum, vehicle->phi);


	graph.lines[lineNum].bdraw = true;

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

void simulateVehicle(DescendVehicle* vehicle,int lineNum) {

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

void simulateVehicle(SpaceCraft* vehicle, int lineNum) {

	int delayCount = 0;
	int drawCount = 0;
	graph.lines[lineNum].vertex.clear();
	vehicle->initialize();


	while (!vehicle->landed) {

		
		vehicle->dynamic(planet);

		if (drawCount >= DRAW_RATE) {
			drawSpaceCraft(lineNum, vehicle);
			drawCount = 0;


			vehicle->printStats();


		}
		drawCount++;

	}
	drawSpaceCraft(lineNum, vehicle);
	

}

void calculateOptimal(DescendVehicle* vehicle) {

	vehicle->calculateOptimalMass(planet);
}

void initializeVehicleIcon(DescendVehicle* vehicle,int lineNum, int shapeEdges, float radius, sf::Color iconColor, sf::Color lineColor) {

	
	graph.initializeIcon(lineNum, shapeEdges, radius, iconColor,lineColor);
	graph.scale_y_max();
	graph.qx = graph.qy;
	//drawSpaceCraft(lineNum, vehicle);
}

void initializeVehicleIcon( int lineNum, int shapeEdges, float radius, sf::Color iconColor, sf::Color lineColor) {


	graph.initializeIcon(lineNum, shapeEdges, radius, iconColor, lineColor);
	graph.scale_y_max();
	graph.qx = graph.qy;
	//drawSpaceCraft(lineNum, vehicle);
}







// Методы -конец-

int main()
{

	

	SC_LRE sc_LRE(SPACE_CRAFT_LRE, SPACE_CRAFT_LRE);
	SC_LRE sc_LRE_2(SPACE_CRAFT_LRE_2, SPACE_CRAFT_LRE);
	ParachuteDescend sc_PAR(SPACE_CRAFT_PAR, SPACE_CRAFT_LRE);

	SpaceCraft spaceCraft(SPACE_CRAFT);

	//DescendVehicle* p[lineCount] = { &sc_LRE,&sc_PAR,&sc_LRE_2 };

	setlocale(LC_ALL, "Russian");
	
	string anyThing;
	string input;
	int command;
	bool runLoop = true;
	
	
   
	graph.initializeAxis();
	
	
	//initializeVehicleIcon(&sc_PAR, 1, 5, 5, sf::Color::Red, sf::Color::Red);
	
	//initializeVehicleIcon(&sc_LRE, 0, 3, 5, sf::Color::Yellow, sf::Color::Yellow);

	//initializeVehicleIcon(&sc_LRE_2, 2, 3, 5, sf::Color::Green, sf::Color::Green);
	
	initializeVehicleIcon( 0,6, 5, sf::Color::Blue, sf::Color::Blue);

	
	thread thread1(graphs);
	thread1.detach();
	
	
	cout << "Graph_init_complete" << "\n";
	//simulateVehicle(&sc_LRE_2, 2);
	/*


	//p = &sc_PAR;
	//simulateVehicle(&sc_PAR);
	//simulateVehicle(&sc_LRE);

//	calculateOptimal(&sc_PAR);
	calculateOptimal(&sc_LRE_2);
	calculateOptimal(&sc_LRE);
	

	simulateVehicle(&sc_LRE_2, 2);
//	simulateVehicle(&sc_PAR, 1);
	simulateVehicle(&sc_LRE, 0);

	drawSpaceCraft(p);
	printf("\n\nLRE vel %-5.2f  LRE_2 vel %-5.2f  Par vel %-5.2f", sc_LRE.vecLength(sc_LRE.velocity), sc_LRE_2.vecLength(sc_LRE_2.velocity), sc_PAR.vecLength(sc_PAR.velocity));
	sc_PAR.initialize(planet);
	sc_LRE.initialize(planet);
	sc_LRE_2.initialize(planet);
	printf("\n\nLRE totalMass %-10.2f  LRE_2 totalMass %-10.2f  Par totalMass %-10.2f", sc_LRE.totalMass, sc_LRE_2.totalMass, sc_PAR.totalMass);
	*/
	
	while (runLoop) {

		printMenu();


		
		while (!(cin >> command)) {
			cin.clear();
			cin.ignore(numeric_limits<streamsize>::max(), '\n');
			cout << "Invalid input\n";
		}

		switch (command) {
		case 1: 
			
			calculateOptimal(&sc_PAR);
			calculateOptimal(&sc_LRE);
			calculateOptimal(&sc_LRE_2);
			break;
		case 2:
			
			//graph.clearLines();
			//drawSpaceCraft(p);
			
			simulateVehicle(&sc_PAR,1);
			simulateVehicle(&sc_LRE, 0);
			simulateVehicle(&sc_LRE_2, 2);
			//drawSpaceCraft(p);
			
			printf("\n\nLRE vel %-5.2f  LRE_2 vel %-5.2f  Par vel %-5.2f", sc_LRE.vecLength(sc_LRE.velocity), sc_LRE_2.vecLength(sc_LRE_2.velocity), sc_PAR.vecLength(sc_PAR.velocity));
			printf("\n\nLRE x %-10.2f  LRE_2 x %-10.2f  Par x %-10.2f", sc_LRE.position.x, sc_LRE_2.position.x, sc_PAR.position.x);
			sc_PAR.initialize(planet);
			sc_LRE.initialize(planet);
			sc_LRE_2.initialize(planet);
			
			printf("\n\nLRE totalMass %-10.2f  LRE_2 totalMass %-10.2f  Par totalMass %-10.2f", sc_LRE.totalMass, sc_LRE_2.totalMass, sc_PAR.totalMass);
			
			
			
			break;

		case 3:
			spaceCraft.initialize();
			simulateVehicle(&spaceCraft, 0);
			break;
		
		case 0:
			runLoop = false;
			break;
		default:
			break;

		}



	}

	printf("\n\n EXIT");

	

	return 0;

}

