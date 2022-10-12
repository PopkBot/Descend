#pragma once
//Импорт -начало-

#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>


using namespace std;
//using std::filesystem::directory_iterator;
//namespace fs = std::filesystem;
//Импорт -конец-

//Объявление констант -начало-

#define FILE_RESOLUTION 51


//Объявление констант -конец-








class ExcelExporter
{
public:


//Объявление переменных -начало-

	ifstream file;



//Объявление переменных -конец-


// Методы -начало-

	void extractMatrixFromFile(string dir, string fileName, double dataArray[][FILE_RESOLUTION], int dataLength, int dataCount) {

		file.open(dir + fileName + ".csv");
		string dataString;
		string::size_type stringPosition;
		file >> dataString;
		
		for (int i = 0; i < dataLength; i++) {
			file >> dataString;
			for (int j = 0; j < dataCount; j++) {
				stringPosition = dataString.find(';');
				//cout << dataString.substr(0, stringPosition)<<"\n";
				
				dataArray[j][i] = stod(dataString.substr(0, stringPosition));
				dataString.erase(0, stringPosition + 1);

			}

		}

		file.close();


	}

	void showFileData(string dir,string fileName,int fileLength) {

		string data;
		string::size_type stringPosition;
		double dataDouble;
		file.open(dir + fileName + ".csv");
		for(int i=0;i< fileLength;i++){
			file >> data;
			stringPosition = data.find(';');
			 data=data.substr(stringPosition+1);
			 dataDouble = stod(data);
			cout << dataDouble << endl;
		}

	}

	void showFileData(string dir,string fileName) {

		string data;
		string::size_type stringPosition;
		double dataDouble;
		file.open(dir + fileName + ".csv");
		while(file.eof()) {
			file >> data;
			stringPosition = data.find(';');
			data = data.substr(stringPosition + 1);
			dataDouble = stod(data) + 0.001;
			cout << dataDouble << endl;
		}

	}

	void extractDataFromFile(string dir, string fileName, double data[], int dataLength) {

		file.open(dir + fileName + ".csv");
		string dataString;
		string::size_type stringPosition;

		for (int i = 0; i < dataLength; i++) {
			file >> dataString;
			stringPosition = dataString.find(';');
			dataString = dataString.substr(stringPosition + 1);
			data[i] = stod(dataString);
			//cout << dataString << "\n";
		}
		file.close();

	}

// Методы -конец-



};

