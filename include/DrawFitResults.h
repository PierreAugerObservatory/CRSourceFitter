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
* @file DrawFitResults.h
* @class DrawFitResults
* @brief Draw the results of the fit (fitted data and MC, elongation rate, rms, chi square, ...)
* 
* @author S. Riggi
* @date 22/04/2010
*/

#ifndef _DRAWFITRESULTS_H_
#define _DRAWFITRESULTS_H_

#include <vector>
#include <string>
#include <TH1D.h>
#include <TFile.h>
#include <TTree.h>
#include <TMatrixD.h>
#include <TF1.h>
#include <TStyle.h>
#include <TApplication.h>

namespace CRSourceFitter_ns {

class DrawFitResults : public TObject {

	public:

		/** 
		\brief Class constructor
 		*/
		DrawFitResults();
		/** 
		\brief Class destructor
 		*/
		virtual ~DrawFitResults();

	public:
		
		/** 
		\brief Draw all
 		*/
		void MakePlots();

		/**
		* \brief Set the number of mass components to be drawn
		*/
		void SetNumberOfComponents(int comp) { fNmass = comp;}
		
		/**
		* \brief Set the input filename containing data info to be drawn
		*/
		void SetDataFileName(const char * fileName) { fDataFileName = fileName;}
		/**
		* \brief Set the output filename containing plots of the fit results
		*/
		void SetOutputFileName(const char * fileName) { fOutputFileName = fileName;}

		/**
		* \brief Set the Xmax MC shift
		*/
  	void SetXmaxShift(double val) {fXmaxShift= val;}
	

	private:
  
  
		/**
		* \brief Inititialize data structures
		*/
		void Init();
		/**
		* \brief Draw the Xmax fit results (distributions, chisquare,...)
		*/
		void DrawXmaxFit();
		/**
		* \brief Draw the fitted mass fractions
		*/
		void DrawMassFractions();
		/**
		* \brief Draw the spectrum fit results
		*/
		void DrawSpectrumFit();
		/**
		* \brief Draw the fitted elongation rate
		*/
		void DrawElongationRate();
		/**
		* \brief Draw the fitted Xmax RMS
		*/
		void DrawRMS();

		/**
		* \brief Function to draw the ER reference given by the models 
		*/
		TF1* defineElong(double* p, int color, int style);
		/**
		* \brief Function to draw the RMS reference given by the models
		*/
		TF1* defineRMSref(double* p, int color, int style);

	
	private:

  	int fNmass; 
		const char * fDataFileName;
		const char * fOutputFileName;

		TFile* fInputFile;
		TFile* fPlotFile;
		TApplication* fApplication;
		TStyle* myStyle;
		double fXmaxShift;

	ClassDef(DrawFitResults,1)
	
};

}//close namespace

#endif

