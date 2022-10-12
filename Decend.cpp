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

//Импорт -конецо-

//Объявление констант -начало-

#define PLANET "Venus"
#define SPACE_CRAFT "SpaceCraft"
#define SPACE_CRAFT_LRE "Soyuz_LRE"

#define DRAW_RATE 10


//Объявление констант -конец-



//Объявление переменных -начало-

Planet planet(PLANET);
SpaceCraft spaceCraft(SPACE_CRAFT);
SC_LRE sc_LRE(SPACE_CRAFT_LRE);
Graph graph;
ExcelExporter excelExporter;


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

// Методы -конец-

int main()
{

	
	
	string anyThing;
	string input;
	
   
	graph.initializeAxis();
	thread thread1(graphs);
	thread1.detach();
	

	
	drawSpaceCraft(0, spaceCraft);
	graph.initializeIcon(0, 3, 5, sf::Color{ 255,0,130 });
	graph.scale_y_max();
	graph.qx = graph.qy;
	cout << "Graph_init_complete" << "\n"<<"Type any key to start\n";
	cin >> anyThing;

	sc_LRE.calculateOptimalMass(planet);

	cout  << "Type any key to continue\n";
	cin >> anyThing;

	int drawCount=0;

	while (true) {
		
		sc_LRE.sc_Initialize(planet);


		while (!(sc_LRE.position.y < H_EPS * 2 && abs(sc_LRE.velocity.y) < V_EPS)) {//(sc_LRE.position.y<H_EPS*2 && abs(sc_LRE.velocity.y)<0.1)

			sc_LRE.control();
			sc_LRE.dynamic(planet);
			//printf("x %2f\ty %2f\tv %3f\tphi %3f\tVx %2f\tVy %4f\n", spaceCraft.position.x, spaceCraft.position.y,
			//	spaceCraft.vecLength(spaceCraft.velocity), spaceCraft.phi*180/3.1415,spaceCraft.velocity.x,
			//	spaceCraft.velocity.y);
			if (drawCount >= DRAW_RATE) {
				drawSpaceCraft(0, sc_LRE);
				drawCount = 0;
				printf("t %-10.3f  y %-10.5f  v %-10.2f  phi %-5.2f  mt %-10.3f\n",
					sc_LRE.time,
					sc_LRE.position.y,
					sc_LRE.vecLength(sc_LRE.velocity),
					sc_LRE.phi * 180 / PI,
					sc_LRE.mt);
				//printf("t %.2f\tfoverLOAD %0.2f\tmaxOL %.2f\n",sc_LRE.time, sc_LRE.overLoad,sc_LRE.maxOverLoad);
				//printf("alpha %0.3f\n", sc_LRE.alpha);

			}
			drawCount++;

			//printf("thrustX %.2f\tthrustY %.2f\n", sc_LRE.thrustVector.x * sc_LRE.engOn, sc_LRE.thrustVector.y * sc_LRE.engOn);
		}
		drawSpaceCraft(0, sc_LRE);
		printf("y %-10.2f  v %-10.2f  phi %-5.2f  mt %-10.3f\n",
			sc_LRE.position.y,
			sc_LRE.vecLength(sc_LRE.velocity),
			sc_LRE.phi * 180 / PI,
			sc_LRE.mt);

		printf("x %2f\ty %2f\tv %3f\tVx %2f\tVy %4f\tphi %3f\n", sc_LRE.position.x, sc_LRE.position.y,
			sc_LRE.vecLength(sc_LRE.velocity), sc_LRE.velocity.x,
			sc_LRE.velocity.y, sc_LRE.phi * 180 / PI);
		printf("max Overload %-5.2f", sc_LRE.maxOverLoad);

		printf("\nType any key to restart = == = == = == = Type 0 to exit\n");
		cin >> input;
		if (input=="0") {
			printf("\nExit\n");
			break;
		}
		else {
			printf("\n Restart\n");
			graph.lines[0].vertex.clear();
		}

	}

	/*
	while (spaceCraft.position.y > 0) {
		spaceCraft.dynamic(planet);
		//printf("x %2f\ty %2f\tv %3f\tphi %3f\tVx %2f\tVy %4f\n", spaceCraft.position.x, spaceCraft.position.y,
		//	spaceCraft.vecLength(spaceCraft.velocity), spaceCraft.phi*180/3.1415,spaceCraft.velocity.x,
		//	spaceCraft.velocity.y);
		drawSpaceCraft(0, spaceCraft);
		printf("x %2f\ty %2f\tv %3f\tVx %2f\tVy %4f\tphi %3f\n",spaceCraft.position.x,  spaceCraft.position.y,
			spaceCraft.vecLength(spaceCraft.velocity),spaceCraft.velocity.x,
			spaceCraft.velocity.y,spaceCraft.phi*180/PI);
	}

	*/

	cout << "\n\nPrint anything\n";
	cin >> anyThing;

	return 0;

}

