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

	static void extractMatrixFromFile(string dir, string fileName, double dataArray[][FILE_RESOLUTION], int dataLength, int dataCount) {

		ifstream file;
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

	static	void extractMatrixFromFile_ForSpaceCraft(string dir, string fileName, double dataArray[][FILE_RESOLUTION], int dataCount,double matrixData[FILE_RESOLUTION][FILE_RESOLUTION]) {

		ifstream file;
		file.open(dir + fileName + ".csv");
		string dataString;
		string::size_type stringPosition;
		file >> dataString;

		for (int i = 0; i < FILE_RESOLUTION; i++) {
			file >> dataString;
			for (int j = 0; j < dataCount; j++) {
				stringPosition = dataString.find(';');
				//cout << dataString.substr(0, stringPosition)<<"\n";

				dataArray[j][i] = stod(dataString.substr(0, stringPosition));
				dataString.erase(0, stringPosition + 1);

			}
			for (int j = 0; j < FILE_RESOLUTION; j++) {
				stringPosition = dataString.find(';');
				//cout << dataString.substr(0, stringPosition)<<"\n";

				matrixData[i][j] = stod(dataString.substr(0, stringPosition));
				dataString.erase(0, stringPosition + 1);

			}

		}

		file.close();


	}

	static void extractMatrixFromFile_ForSpaceCraft(string dir, string fileName, double dataArray[][FILE_RESOLUTION], int dataCount, double matrixData1[FILE_RESOLUTION][FILE_RESOLUTION], double matrixData2[FILE_RESOLUTION][FILE_RESOLUTION]) {
		
		ifstream file;
		file.open(dir + fileName + ".csv");
		string dataString;
		string::size_type stringPosition;
		file >> dataString;

		for (int i = 0; i < FILE_RESOLUTION; i++) {
			file >> dataString;
			for (int j = 0; j < dataCount; j++) {
				stringPosition = dataString.find(';');
				//cout << dataString.substr(0, stringPosition)<<"\n";

				dataArray[j][i] = stod(dataString.substr(0, stringPosition));
				dataString.erase(0, stringPosition + 1);

			}
			for (int j = 0; j < FILE_RESOLUTION; j++) {
				stringPosition = dataString.find(';');
				//cout << dataString.substr(0, stringPosition)<<"\n";

				matrixData1[i][j] = stod(dataString.substr(0, stringPosition));
				dataString.erase(0, stringPosition + 1);

			}
			for (int j = 0; j < FILE_RESOLUTION; j++) {
				stringPosition = dataString.find(';');
				//cout << dataString.substr(0, stringPosition)<<"\n";

				matrixData2[i][j] = stod(dataString.substr(0, stringPosition));
				dataString.erase(0, stringPosition + 1);

			}

		}

		file.close();


	}

	static void showFileData(string dir,string fileName,int fileLength) {
		ifstream file;
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

	static void showFileData(string dir,string fileName) {
		ifstream file;
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

	static void extractDataFromFile(string dir, string fileName, double data[], int dataLength) {
		ifstream file;
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

