#include "ImageProcessor.h"
#include <TSystem.h>
//main
int main(){
    //Criar matrix para armazenar os piexeis da imagem
    vector<vector<int>> M;
    int N = 0;
    
    ReadImage("Trab01_imagem_N.ascii.pgm", M, N);
    Image image = MakeImage(M, N);
    AddNoise(image, "G",  30);

    TH2D *NoiseGaussiano = GetImageHistogram(image);
    TCanvas *c = new TCanvas("c", "NoiseGaussaino", 200,10,600,600);
    NoiseGaussiano->SetMinimum(-0.1);
    gStyle->SetPalette(kGreyScale);
    NoiseGaussiano->Draw("COLZ");
    c->Update();
    gSystem->ProcessEvents();
    c->WaitPrimitive();
    c->SaveAs("../Trab01_Hough_imagem_1.png");
    
    ///////////////////

    Image n = MakeImage(M, N);
    TH2D *HoughHist = HoughTransform(n);

    TCanvas c3("c3", "", 600, 600);
    HoughHist->SetMinimum(-0.1);
    gStyle->SetPalette(kBird);
    HoughHist->Draw("surf1z");
    c3.Update();
    c3.SaveAs("../Trab01_Hough_imagem_2.png");
    c3.WaitPrimitive();

    TH2D *imagem = GetImageHistogram(image);
    TH2D *HoughLines = hough_to_img(HoughHist,imagem);

}