#ifndef __IMAGEEDITOR__
#define __IMAGEEDITOR__

#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <string>
#include <vector>

#include "TApplication.h"
#include <TCanvas.h>
#include <TH2F.h>
#include <TStyle.h>
#include <TF1.h>

using namespace std;

struct Image;

//IMAGE READ, MANIPULATION AND PGM

//Read a pgm from a file and stores the information on M and N
void ReadImage(string filenam, Image& image);
void ReadImage(string filename, vector<vector<int>>& M, int& N);
//Uses the M and N values to create a Image struct
Image MakeImage(vector<vector<int>>& M, int& N);
//Produces a pgm file with inverted colors
Image InvertImage(Image image);
void WriteImage(string filename, const vector<vector<int>>& M, int& N);
//sums two images values, producing a new one
Image ImageSum(const Image& img1, const Image& img2);
//add a segment in the matrix of a Image
void addSegment(Image& image, vector<int> start, vector<int> end, int color);
//add a segment in the matrix of a Image (more specifications)
void addSegment(Image& image, vector<int> start, vector<int> end, int color, int size, string cluster_orientation);

//AUXILIAR

// makes a new image with the settings of another, doesnt copy image matrix
Image CopySettingsImage(Image& image);
//produces a pgm file with image info
void ProduceImage(string filename, Image& image);
//verify if two images have the same settings
bool SameSettings(Image image1, Image image2);
//creates a white image
void CreateEmpty(Image& image, int ncols, int nrows, int N);

//NOISE

//adds noise to a image
void AddNoise(Image& image, string noise_kind, double noise_par);

//generator of random bool value considering a probability
bool GenerateRandomBool(double p);
//random centered gaussian value
double Gaussian(double stdv);
//???
int GetRandomExtreme(int N);
//cartesian to matrix coordinates

//CALCULATIONS

//calculates absolute frequencies
vector<int> GetColourFreq(Image& image);
//calculates relative frequencies
vector<double> GetColourRelFreq(Image& image);
//calculates the media of the colors in a image
double MedColor(Image image);
//calculates the variance of the image values
double Variance(Image image);

vector<int> cart_to_matrix(int nrows, vector<int> c);

//HOUGH

//filters an image with default median filter
void Mediana_Quadrado(Image& image);
//threshholds an image and leaves the edges with a tolerance
void Threshholding(Image& image, double tolerance);
//applies hough algorithm to find lines and returns the p,theta parameter histogram
TH2D* HoughTransform(Image& image);
TH2D* hough_to_img(TH2D* hough_space,TH2D* image);


//PLOTS

//makes a 2D root histogram using a image info
TH2D* GetImageHistogram(Image& image);
int MaxCalc(vector<int> vec);
void CounterVector(vector<int> vec, int val);
void CreateFileHistogram(string title, string filename, vector<int> size, TH2D* hist2D); // JUST TO CHECK

void MakeHistogram(std::vector<int> info, std::string filename);

void MakeHistogram(std::vector<double> info, std::string filename);

void MakeHistogram2D(std::vector<std::vector<int>> info, std::string filename);


struct Image {
    string cod; // codification (not using)
    int nrows, ncols; //numb of rows and columns
    int N; // max color value (white)
    vector<vector<int>> C; // matrix of pixel colors
};

#endif

// FUNCTIONS NOT USED:
    // void PrintColorCoordinates(string stored, int color);
    // void PlotAbsFreq(string stored, string filename);
    // void PlotRelFreq(string stored, string filename);
    // void MakeHistogram(vector<int> info, string filename);
    // void MakeHistogram(vector<double> info, string filename);
    // void MakePointsPlot(vector<int> info, string filename);
    // void MakeHistogram2D(vector<vector<int>> info, string filename);
    // void Image2DHist(string stored, string filename);