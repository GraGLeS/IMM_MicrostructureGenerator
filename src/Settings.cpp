/*
 * Settings.cpp
 *
 *  Created on: Apr 25, 2016
 *      Author: Jonathan Nierengarten
 */
//*********************************************************************************************************************************************************************************************
#include "Settings.h"
#include "rapidxml.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <cstring>
#include <stdexcept>
#include <string>

using namespace std;
using namespace rapidxml;

//*********************************************************************************************************************************************************************************************

//Initialization of the setting variables *****************************************************************************************************************************************************

unsigned int Settings::NumberOfGrains = 0;
unsigned int Settings::NumberOfSubgrains = 0;
unsigned int Settings::NumberOfGridpoints = 0;
unsigned long Settings::NumberOfPointsPerSubGrain = 0;
unsigned long Settings::MaximumNumberOfThreads = 1;
double Settings::Settings::SubgrainOriScatter = 0.0;
double Settings::StoredElasticEnergyMax = 0.0;
double Settings::StoredElasticEnergyMin = 0.0;
double Settings::StoredElasticScatterGrain = 0.0;
double Settings::StoredElasticScatterSubgrain = 0.0;


bool Settings::ExecuteInParallel = false;
bool Settings::PlotIPF2DSection = false;
double Settings::PlotWindowXMin = 0.0;
double Settings::PlotWindowXMax = 1.0;
double Settings::PlotWindowYMin = 0.0;
double Settings::PlotWindowYMax = 1.0;

//bool Settings::UseOrientationSpace = false;

E_CRYSTAL_STRUCTURE Settings::CrystalStructure = E_DEFAULT_STRUCTURE;
E_MICROGENMODE Settings::MicroGenMode = E_DEFAULT_GENERATOR;
E_TEXTURE Settings::TextureGEN = E_DEFAULT_TEXTURE;
E_PLOT_DIMENSION Settings::PlotDimension = E_DEFAULT_DIMENSION;

string Settings::ReadFromFilename;
string Settings::AdditionalFilename;

//Definition of the needed functions **********************************************************************************************************************************************************

void Settings::readXML(string filename) {

	//Find the wished .xml file ***********************************************************************************************************************************************************

	if (0 == filename.compare(""))
		filename = string("parameters.xml");
	ifstream file(filename);
	if (file.fail()) {
		throw runtime_error(string("Unable to locate file ") + filename);
	}
	stringstream contents;
	contents << file.rdbuf();
	string xmlDocument(contents.str());

	xml_document<> tree;
	tree.parse<0>(&xmlDocument[0]);

	xml_node<>* rootNode = tree.first_node();
	if (0 != strcmp(rootNode->name(), "Parameters")) {
		throw runtime_error("undefined parameters file!");
	}

	//Read all the required parameters ****************************************************************************************************************************************************

	if (0 != rootNode->first_node("NumberOfGrains")) {
		NumberOfGrains = std::stoul(
				rootNode->first_node("NumberOfGrains")->value());
	}

	if (0 != rootNode->first_node("NumberOfSubgrains")) {
		NumberOfSubgrains = std::stoul(
				rootNode->first_node("NumberOfSubgrains")->value());
	}
	//*****************************************************************************************************************************

	if (0 != rootNode->first_node("NumberOfPointsPerSubGrain")) {
		NumberOfPointsPerSubGrain = std::stoul(
				rootNode->first_node("NumberOfPointsPerSubGrain")->value());
		if ( NumberOfPointsPerSubGrain <= MINIMUM_DISCRETIZATION )
			NumberOfPointsPerSubGrain = MINIMUM_DISCRETIZATION;
	}
	if (0 != rootNode->first_node("StoredElasticEnergyMin")) {
		StoredElasticEnergyMin = std::stod(
				rootNode->first_node("StoredElasticEnergyMin")->value());
		if ( StoredElasticEnergyMin <= 0.0 ) 
			StoredElasticEnergyMin = 0.0;
	}
	if (0 != rootNode->first_node("StoredElasticEnergyMax")) {
		StoredElasticEnergyMax = std::stod(
				rootNode->first_node("StoredElasticEnergyMax")->value());
		if ( StoredElasticEnergyMax <= 0.0 ) 
			StoredElasticEnergyMax = 0.0;
		if ( StoredElasticEnergyMax <= StoredElasticEnergyMin ) 
			StoredElasticEnergyMax = StoredElasticEnergyMin;
	}
	if (0 != rootNode->first_node("StoredElasticScatterGrain")) {
		StoredElasticScatterGrain = std::stod(
				rootNode->first_node("StoredElasticScatterGrain")->value());
		if ( StoredElasticScatterGrain <= 0.0 )
			StoredElasticScatterGrain = 0.0;
	}
	if (0 != rootNode->first_node("StoredElasticScatterSubgrain")) {
		StoredElasticScatterSubgrain = std::stod(
				rootNode->first_node("StoredElasticScatterSubgrain")->value());
		if ( StoredElasticScatterSubgrain <= 0.0 ) 
			StoredElasticScatterSubgrain = 0.0;
	}

	if (0 != rootNode->first_node("SubgrainOriScatter")) {
		SubgrainOriScatter = std::stod(
				rootNode->first_node("SubgrainOriScatter")->value());
		if ( SubgrainOriScatter <= 0.5 )
			SubgrainOriScatter = 0.5;
	}
	//*****************************************************************************************************************************

	if (0 != rootNode->first_node("ExecuteInParallel")) {
		ExecuteInParallel = (bool) std::stoul(
				rootNode->first_node("ExecuteInParallel")->value());
	}

	//if (0 != rootNode->first_node("UseOrientationSpace")) {
	//	UseOrientationSpace = (bool) std::stoul(
	//			rootNode->first_node("UseOrientationSpace")->value());
	//}

	//*****************************************************************************************************************************

	if (0 != rootNode->first_node("CrystalStructure")) {
		CrystalStructure = (E_CRYSTAL_STRUCTURE) std::stoi(
				rootNode->first_node("CrystalStructure")->value());

		if (CrystalStructure >= E_DEFAULT_STRUCTURE)
			CrystalStructure = E_DEFAULT_STRUCTURE;
	}

	if (0 != rootNode->first_node("MicroGenMode")) {
		MicroGenMode = (E_MICROGENMODE) std::stoi(
				rootNode->first_node("MicroGenMode")->value());

		if (MicroGenMode >= E_DEFAULT_GENERATOR)
			MicroGenMode = E_DEFAULT_GENERATOR;
	}

	if (0 != rootNode->first_node("TextureGEN")) {
		TextureGEN = (E_TEXTURE) std::stoi(
				rootNode->first_node("TextureGEN")->value());

		if (TextureGEN > E_USE_PREFERENCEORI)
			TextureGEN = E_DEFAULT_TEXTURE;
	}

	if (0 != rootNode->first_node("PlotDimension")) {
		PlotDimension = (E_PLOT_DIMENSION) std::stoi(
				rootNode->first_node("PlotDimension")->value());

		if (PlotDimension >= E_DEFAULT_DIMENSION)
			PlotDimension = E_DEFAULT_DIMENSION;
	}

	if (0 != rootNode->first_node("PlotIPF2DSection")) {
		PlotIPF2DSection = (bool) std::stoul(
				rootNode->first_node("PlotIPF2DSection")->value());
	}
	if (0 != rootNode->first_node("PlotWindowXMin")) {
			PlotWindowXMin = std::stod(
					rootNode->first_node("PlotWindowXMin")->value());
	}
	if (0 != rootNode->first_node("PlotWindowXMax")) {
				PlotWindowXMax = std::stod(
						rootNode->first_node("PlotWindowXMax")->value());
	}
	if (0 != rootNode->first_node("PlotWindowYMin")) {
				PlotWindowYMin = std::stod(
						rootNode->first_node("PlotWindowYMin")->value());
	}
	if (0 != rootNode->first_node("PlotWindowYMax")) {
				PlotWindowYMax = std::stod(
						rootNode->first_node("PlotWindowYMax")->value());
	}
	if ( PlotWindowXMin <= 0.0 ) PlotWindowXMin = 0.0;
	if ( PlotWindowXMax >= 1.0 ) PlotWindowXMax = 1.0;
	if ( PlotWindowYMin <= 0.0 ) PlotWindowYMin = 0.0;
	if ( PlotWindowYMax >= 1.0 ) PlotWindowYMax = 1.0;

	//*****************************************************************************************************************************

	if (0 != rootNode->first_node("ReadFromFilename")) {
		ReadFromFilename = rootNode->first_node("ReadFromFilename")->value();
	}

	if (0 != rootNode->first_node("AdditionalFilename")) {
		AdditionalFilename =
				rootNode->first_node("AdditionalFilename")->value();
	}

}

