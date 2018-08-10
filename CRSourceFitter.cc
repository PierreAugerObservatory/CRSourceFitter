// ******************************************************************************
// * License and Disclaimer                                                     *
// *                                                                            *
// * Copyright 2011 Simone Riggi																			          *
// *																																	          *
// * This file is part of CRSourceFitter																        *
// * CRSourceFitter is free software: you can redistribute it and/or modify it  *
// * under the terms of the GNU General Public License as published by          *
// * the Free Software Foundation, either * version 3 of the License,           *
// * or (at your option) any later version.                                     *
// * CRSourceFitter is distributed in the hope that it will be useful, but 			*
// * WITHOUT ANY WARRANTY; without even the implied warranty of                 * 
// * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                       *
// * See the GNU General Public License for more details. You should            * 
// * have received a copy of the GNU General Public License along with          * 
// * CRSourceFitter. If not, see http://www.gnu.org/licenses/.                  *
// ******************************************************************************
/**
* @mainpage CRSourceFitter tool

* <img src="SourceCompositionFitter_screenshot.jpg" alt="Screenshot">
* <b>Description of the project </b> <br>
* This project is developed to perform a combined fit of spectrum and Xmax data. <br>
The project settings can be specified in a configuration file which is parsed by the ConfigParser class.<br>
The MCReader and DataReader classes are provided to read and get data from ROOT files containing the data (Xmax, spectrum) and MC information (both reconstructed and generated). A reader is provided also to get data from a ROOT file containing the simulated propagation data. Its primary function is to build the transfer matrix used to propagate the spectrum from the source to the Earth, and to calculate spectra observed @ Earth. <br>
A global model with a set of parameters is assumed and fitted to spectrum+Xmax data. The fitted parameters are: injection index gamma, maximum acceleration energy Emax @ source, flux normalization, mass fractions. Eventually the shift of Xmax distribution can be included as fit parameter. The likelihood fit of Xmax, spectrum, Xmax+spectrum data is done by the SourceFitter class. Eventually the EarthFitter class perform the Xmax fit @ Earth level <br>
As utilities the Drawer class plus some other ROOT macros are provided to draw the obtained fit results and other operations.
* 
* @author S. Riggi <br> Dept. of Physics and Astronomy, University of Catania <br> INFN, Section of Catania  <br> Italy
*/

/**
* @file CRSourceFitter.cc
* @class main
* @brief main class used to handle the entire project
*
* ...
* @author S. Riggi
* @date 22/04/2010
*/

#include "ConfigParser.h"
#include "SourceFitter.h"
#include "DrawFitResults.h"
#include "PropagationMCReader.h"
#include "AnalysisConsts.h"
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TCut.h>
#include <TEventList.h>
#include <TMath.h>
#include <TPad.h>
#include <TVirtualPad.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <TDirectory.h>
#include <TStyle.h>

#include<TPaveText.h>
#include<TVirtualFitter.h>
#include<TObjArray.h>
#include<TMatrixD.h>
#include<TColor.h>

#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>
#include <getopt.h>

using namespace std;


//FileNames
std::string DataFileName;
std::string MCFileName;
std::string MCPropFileName;
std::string OutputFileName;
std::string configFileName= "";
std::string inputFileName= "";


//Source params
int Nmass;
double InjIndexAtSource;
double EmaxAtSource;
double RedShiftEvolution;
double FluxNorm;
int SourceEnergyCutoffMode;
double SourceEnergyCutoffShape;
bool EmaxADependent;
std::vector<double> FractAtSource;
std::vector<double> FractParamsAtSource;

//Fit options
bool FitGalSpectrum;
bool FitExtragalFractions;
bool FitGalFractions;
bool FitSourceCutoffShape;
bool CombinedFit;
bool StoreFitResults;
bool FixMassFractions;
bool FixGamma;
bool FixEmax;
bool FixCutoffShape;
bool WeightLikelihood;
double SpectrumLikelihoodWeight;
double XmaxLikelihoodWeight;

void DoSourceFit();
void Draw();

void Usage(char* exeName){
	cout<<"=========== USAGE ==========="<<endl;
	cout<<"Usage: "<<exeName<<" [options]"<<endl;
	cout<<endl;
	cout<<"Options:"<<endl;
  cout<<"-h, --help \t Show help message and exit"<<endl;
	cout<<"-c, --config=[CONFIG FILENAME] \t Configuration file name with source fit run options"<<endl;
	cout<<"-i, --input=[INPUT FILENAME] \t Input file name (.root) containing fit info to be drawn"<<endl;
	cout<<"-d, --draw \t Draw fitted results instead of running source fit"<<endl;
	cout<<"=============================="<<endl;
}

static const struct option options_tab[] = {
  /* name, has_arg, &flag, val */
  { "help", no_argument, 0, 'h' },
	{ "input", required_argument, 0, 'i' },
	{ "config", required_argument, 0, 'c' },
	{ "draw", no_argument, 0, 'd' },
  {(char*)0, (int)0, (int*)0, (int)0}
};


int main(int argc, char **argv) 
{
	//====================================================
	//==         PARSE ARGS
	//=====================================================
	//## Check args
	if(argc<2){
		cout<<endl;
		cerr<< "ERROR: Incorrect number of arguments...see program usage!"<<endl;
		Usage(argv[0]);		
		exit(1);
	}

	//## Get command args
	int c = 0;
  int option_index = 0;
	bool runFit= true;

	while((c = getopt_long(argc, argv, "hi:c:d",options_tab, &option_index)) != -1) {
    
    switch (c) {
			case 0 : 
			{
				break;
			}
			case 'h':
			{
      	Usage(argv[0]);	
				exit(0);
			}
    	case 'i':	
			{
				inputFileName= std::string(optarg);	
				break;	
			}
			case 'c':	
			{
				configFileName= std::string(optarg);	
				break;	
			}
			case 'd':	
			{
				runFit= false;
				break;	
			}
    	default:
			{
      	Usage(argv[0]);	
				exit(0);
			}
    }//close switch
	}//close while
	
	//====================================================
	//==         CHECK ARGS
	//=====================================================
	CRSourceFitter_ns::ConfigParser* parser= nullptr;
	if(runFit){
		//Check config file name given
		if(configFileName==""){
			cerr<<"ERROR: Empty configuration filename given!"<<endl;
			exit(1);
		}

		//Parsing info
		parser= new CRSourceFitter_ns::ConfigParser(configFileName.c_str());
		if(!parser){
			cerr<<"ERROR: Failed to create parser or to read config file "<<configFileName<<"!"<<endl;
			exit(1);
		}
		parser->SetVerbosity(1);
		parser->ReadConfig();
	}//close if
	else{
		if(inputFileName==""){
			cerr<<"ERROR: Empty input filename given!"<<endl;
			exit(1);
		}
	}	


	//====================================================
	//==         RUN FIT/DRAW
	//=====================================================
	if(runFit){
		cout<<"************************************"<<endl;
		cout<<"** SOURCE FITTING                ***"<<endl;
		cout<<"************************************"<<endl;
		DoSourceFit();
	}
	else{
		cout<<"************************************"<<endl;
		cout<<"** DRAWING                       ***"<<endl;
		cout<<"************************************"<<endl;
		Draw();
	}

  return 0;
  
}//close main



void DoSourceFit(){

  //Create a CompositionFitter class
  CRSourceFitter_ns::SourceFitter theSourceFitter;
	theSourceFitter.SetVerbosity(2);
  theSourceFitter.DoSourceFit();	
  
}//close DoSourceFit()


void Draw(){

	//Create a CompositionFitter class
  CRSourceFitter_ns::DrawFitResults theDrawer;
  
  //Set some class variables
  theDrawer.SetNumberOfComponents(Nmass);
  theDrawer.SetDataFileName(inputFileName.c_str());
	theDrawer.SetOutputFileName("Plot.root");	
	theDrawer.MakePlots();

}//close Draw()


