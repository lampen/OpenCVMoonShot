//Includes
#include<opencv2\opencv.hpp>﻿
#include<opencv2\core\core.hpp>
#include<opencv2\highgui\highgui.hpp>
#include<opencv2\imgproc\imgproc.hpp>
#include<iostream>
#include<iomanip>
#include<sstream>
#include<string>
#include<opencv\highgui.h>
#include<opencv\cv.h>
#include<vector>
#include<cstdlib>
#include<stdio.h>
#include<stdlib.h>
#include<fstream>

#include "SunMoon.h"
#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>
#include <limits>
#include <string>
#include <iostream>
#include <time.h>
#include <chrono>
#include <ctime>

#pragma warning(disable : 4996)

//namespace
using namespace cv;
using namespace std;

//MoonLocation
using namespace Astronomy;
struct tm m_datumzeit;
int m_Zone;
double latitude = 55.782604799999994; //Output Datahub
double longitude = 12.528025600000001;  //Output Datahub
double deltaT = 65.0;

//pos
double posX = 4;  //Output Datahub
double posY = 5;  //Output Datahub

//moonPos
double moonPosAz;  //Input Datahub
double moonPosAlt;  //Intput Datahub
double lookingAtAz = 6; //Output Datahub
double lookingAtAl = 7; //Output Datahub

//Batteri
double bat1Voltage = 10; //Output Datahub
double bat2Voltage = 11; //Output Datahub

//clock
struct tm gmt;

//txt file
std::ofstream outfile;
std::string logLineFirst = "Clock \t \t \t OurPosX \t \t OurPosY \t \t LookingAtAz \t \t LookingAtAl \t \t MoonPosAz \t \t MoonPosAl \t \t Batteri1 \t \t Batteri2\n";
int uncertainty = 2;

const std::string currentDateTime() {
	time_t     now = time(0);
	struct tm  tstruct;
	char       buf[80];
	tstruct = *localtime(&now);
	strftime(buf, sizeof(buf), "%d-%X", &tstruct);
	return buf;
}

void logFile() {
	std::string logLine = "";
	outfile.open("LogFile.txt", std::ios_base::app);
	std::string TimeStamp = currentDateTime();
	logLine += (TimeStamp) + "\t\t";
	logLine += std::to_string(posX) + "\t\t" + std::to_string(posY) + "\t\t";
	logLine += std::to_string(lookingAtAz) + "\t\t" + std::to_string(lookingAtAl) + "\t\t";
	logLine += std::to_string(moonPosAz) + "\t\t" + std::to_string(moonPosAlt) + "\t\t";
	logLine += std::to_string(bat1Voltage) + "\t\t" + std::to_string(bat2Voltage) + "\n";
	outfile << logLine;
	//std::cout << logLine;
	outfile.close();
}

void Time() {
	time_t now = time(NULL);
#ifdef _MSC_VER
	memset(&gmt, 0, sizeof(struct tm));
	memset(&m_datumzeit, 0, sizeof(struct tm));
	localtime_s(&m_datumzeit, &now);
	gmtime_s(&gmt, &now);
#else
	m_datumzeit = *localtime(&now);
	struct tm gmt = *gmtime(&now);
#endif
	time_t utc = mktime(&gmt);
	if (gmt.tm_isdst == 1) utc -= 3600;
	m_Zone = (int)(difftime(now, utc) / 3600.0);
	
}

void MoonLocation() {
	Astronomy::SunMoon::SunMoon(latitude, longitude, &m_datumzeit, deltaT, m_Zone);

	&Astronomy::SunMoon::GetMoonAlt;
	&Astronomy::SunMoon::GetMoonAz;
	//std::cout << "MoonAlt:  "<< m_MoonAlt << "		Moonz:  " << m_MoonAz << "\n";


}

int main() {
	outfile.open("LogFile.txt", std::ios_base::app);
	outfile << logLineFirst;
	outfile.close();

	VideoCapture cap(0); // open the default camera
	if (!cap.isOpened()) {  // check if we succeeded
		return -1;
	}

	while (1) {
		//checking for position comment out
		/*if (lookingAtAz <= moonPosAz - uncertainty && lookingAtAz >= moonPosAz + uncertainty && lookingAtAl <= moonPosAl - uncertainty && lookingAtAl >= moonPosAl + uncertainty) {
			std::cout << "Az" << "		" << lookingAtAz << "		" << moonPosAz <<"\n";
			std::cout << "Al" << "		" << lookingAtAl << "		" << moonPosAl <<"\n";
			_sleep(1000);
		}*/
		

		//Time
		Time();

		//MoonLocation
		MoonLocation();

		//picture
		Mat imgImport;

		cap.read(imgImport);

		imshow("img", imgImport);

		if (imgImport.empty()) {
			std::cout << "error: image did'nt load ";
			return -2;
		}

		std::string saveNameMoon = "../Testing with moonPos/MoonFolder/";
		cv::String Path2(saveNameMoon);
		vector<cv::String> fileNameSaveMoon;
		cv::glob(Path2, fileNameSaveMoon, true);
		cv::Mat imgGrayScale;
		double meanM;

		//const
		float moonThreshold = 20;
		int printedLine = 0;

		//Convert to grayscale
		cv::cvtColor(imgImport, imgGrayScale, CV_BGR2GRAY);

		//prints one line with:
		if (printedLine == 0) {
			//std::cout << "Mean" << "		" << "Over threshold\n";
			printedLine = 1;
		}

		//mean
		Scalar average = mean(imgGrayScale);

		//Scalar->string
		string AVG = std::to_string(average[0]);


		//prints AVG
		String width = "	";
		//std::cout << AVG << width;

		//Is the picture within the threshold
		String TorF;

		meanM = average[0];
		if (moonThreshold <= meanM) {
			TorF = "True";
		}
		else {
			TorF = "false";
		}
		//std::cout << TorF + "\n";
		cv::imwrite(saveNameMoon + AVG + ".jpg", imgImport);


		//running logfile
		logFile();		
	}
	// hold windows open
	system("pause");
	return(0);
}
