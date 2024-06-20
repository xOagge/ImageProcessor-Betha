#include "ImageProcessor.h"


Image CopySettingsImage(Image& image){
    Image result;
    result.nrows = image.nrows;
    result.ncols = image.ncols;
    result.N = image.N;
    result.C.resize(result.ncols, vector<int>(result.nrows, 0));
    return result;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ProduceImage(string filename, Image& image) {
    // Check if the image exists
    ofstream outFile("../" + filename, ios::binary);
    if (!outFile) {
        cerr << "Error: Unable to open file '" << filename << "for writing." << endl;
        return;
    }

    // Write the PGM file header
    outFile << "P2" << endl;
    outFile << image.ncols << " " << image.nrows << endl;
    outFile << image.N << endl;

    // Write pixel data to the file
    for (int row = 0; row < image.nrows; row++) {
        for (int col = 0; col < image.ncols; col++) {
            outFile << image.C[col][row] << " ";
        }
        outFile << endl; // Newline after each row
    }

    // Close the file
    outFile.close();
}

bool SameSettings(Image image1, Image image2){
    if(image1.nrows == image2.nrows && image1.ncols == image2.ncols && image1.N == image2.N){
        return true;
    }
    return false;
}

void CreateEmpty(Image& image, int ncols, int nrows, int N){
    image.nrows = nrows;
    image.ncols = ncols;
    image.N = N;
    image.C.resize(image.ncols, vector<int>(image.nrows, N));
}

bool GenerateRandomBool(double p) {
    if (p <= 0 || p > 100) {
        cerr << "Error: Probability must be above 0 and equal or below 100." << endl;
        return false;
    }

    // Seed the random number generator
    random_device rd;
    mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()

    // Generate a random integer between 0 and 99
    uniform_int_distribution<> dis(0, 100);
    int randNum = dis(gen);

    // Return true with probability p%
    return randNum <= p;
}

int GetRandomExtreme(int N) {
    // Seed the random number generator
    random_device rd;
    mt19937 gen(rd());

    // Generate a random boolean value
    uniform_int_distribution<> dis(0, 1);
    bool value = dis(gen);

    // Return 45 if the random boolean value is false (0), otherwise return 88
    return value ? 0 : N;
}

double Gaussian(double stdv) {
    // Create a random number generator
    random_device rd;
    mt19937 gen(rd());
    normal_distribution<double> distribution(0.0, stdv);

    // Generate a random value from the normal distribution
    return distribution(gen);
}

vector<int> cart_to_matrix(int nrows, vector<int> c){
    // takes a point in cartesian (x,y) and returns the corresponding point in the image matrix (y,nrows-x)
    vector<int> m;
    m.push_back(nrows - c[1]);
    m.push_back(c[0]);
    return m;
}

// HOUGH STEPS
//apply squared median with dimension one for threshholding
void Mediana_Quadrado(Image& image) {
    Image result = CopySettingsImage(image);
    result.C=image.C;

    int r = 1;
    //iteração
    for (int i = 0; i < result.ncols; i++) {
        for (int j = 0; j < result.nrows; j++) {
            vector<int> ord; //vetor para recolher os valores à volta

            for (int ni = -r; ni <= r; ni++) { 
                for (int nj = -r; nj <= r; nj++) { 
                    if ((i+ni) >= 0 && (j+nj) >=0 && (i+ni) < image.ncols && (j+nj) < image.nrows) {
                        ord.push_back(image.C[i+ni][j+nj]); //adicionar os valores
                    } 
                } 
            }
            //ordenar e recolher o elemnto central
            sort(ord.begin(),ord.end());
            if (ord.size()%2) {result.C[i][j] = ord[ord.size()/2]; } //verificar se é par ou impar
            else {result.C[i][j] = (ord[ord.size()/2]+ord[ord.size()/2-1])/2; }
        }
    }
    //guardar a imagem ou é alterado?
    image.C=result.C;
}

void Threshholding(Image& image, double threshhold){
    //verificação?
    Image result = CopySettingsImage(image);
    result.C=image.C;
    // Iterate over each pixel in the image
    for (int col = 1; col < image.ncols-1; col++) {
        for (int row = 1; row < image.nrows-1; row++) {
            double x = abs(image.C[col-1][row]-image.C[col+1][row])/2;
            double y = abs(image.C[col][row-1]-image.C[col][row+1])/2;
            double gradient = sqrt(x*x+y*y);
            if(gradient>threshhold){
                result.C[col][row]=image.N;;
            }
        }
    }
    image.C=result.C;
}

TH2D* HoughTransform(Image& image){
    TH2D *imagem = GetImageHistogram(image);

    //contabilizar num histograma de parâmteros
    double a_min=-M_PI/2;
    double a_max=M_PI/2;
    double p_min=0;
    double p_max=sqrt((imagem->GetNbinsX()*imagem->GetNbinsX())+(imagem->GetNbinsY()*imagem->GetNbinsY()));
    int as=100;
    int ps=100;
    double a_step = (a_max-a_min)/as;
    double p_step = p_max/ps;

    //criar o histograma que vai ser retornado
    TH2D *hough_hist = new TH2D("hough_hist","hough_hist",as,a_min,a_max,ps,0,p_max);
    double p=p_min+p_step/2;
    while(p<p_max) {
        double a=a_min+a_step/2;
        while (a<a_max) {
            int counter = 0;
            for (int binX = 1; binX <= imagem->GetNbinsX(); ++binX) {
                for (int binY = 1; binY <= imagem->GetNbinsY(); ++binY) {
                    if (imagem->GetBinContent(binX, binY) == 0) {
                        std::function<double(double)> function = [binX, binY](double t) -> double {return binX * cos(t) + binY * sin(t);};
                        double value = function(a);
                        if(p - p_step/2 < value && value < p + p_step/2){counter+=1;}
                    }
                }
            }
            hough_hist->Fill(a, p, counter);
            a+=a_step;
        }
        p+=p_step;
    }
    return hough_hist;
}



TH2D* hough_to_img(TH2D* hough_space, TH2D* image){
    //localizar os spots do histograma (assumem-se linhas os maximos locais)
    vector<vector<double>> yellow_spots;

    //determinar o bin de maior valor
    int max=0;
    for (int binX = 1; binX <= hough_space->GetNbinsX(); ++binX) {
        for (int binY = 1; binY <= hough_space->GetNbinsY(); ++binY) {
            double count = hough_space->GetBinContent(binX, binY);
            if (count>max) {
                max=count;
            }
        }
    }
    //maximos locais são aprox valores próximo em relação ao máximo
    for (int binX = 1; binX <= hough_space->GetNbinsX(); ++binX) {
        for (int binY = 1; binY <= hough_space->GetNbinsY(); ++binY) {
            double count = hough_space->GetBinContent(binX, binY);
            if (count>0.4*max) {
                vector<double> spot = {binX,binY};
                yellow_spots.push_back(spot);
            }
        }
    }

    //desenhar o hist e adicionar segmentos
    TCanvas c3("c3", "", 600, 600);
    image->SetMinimum(-0.1);
    gStyle->SetPalette(kGreyScale);
    image->Draw("colz");

    for (int i=0; i<yellow_spots.size();i++){
        double a_step=M_PI/hough_space->GetNbinsX();
        double p_step=sqrt((image->GetNbinsX()*image->GetNbinsX())+(image->GetNbinsX()*image->GetNbinsX()))/hough_space->GetNbinsY();
        double a=yellow_spots[i][0]*a_step - M_PI/2 - a_step/2;
        double p=yellow_spots[i][1]*p_step - p_step/2;
        double m = -cos(a)/sin(a);
        double b = p/sin(a);
        TF1 *line = new TF1("line", "[0]*x + [1]", image->GetXaxis()->GetXmin(), image->GetXaxis()->GetXmax());
        line->SetParameters(m, b);
        line->SetLineWidth(4);
        line->SetLineColor(8);
        line->Draw("same");
        cout<<a<<","<<p<<endl;
        cout<<m<<","<<b<<endl;
    }

    c3.Update();
    c3.SaveAs("../Trab01_Hough_imagem_3.png");
    c3.WaitPrimitive();

    return image;

}


int MaxCalc(vector<int> vec){
    int result = 0;
    for (int i = 0; i < vec.size(); i++)
        if (vec[i] > result)
            result = vec[i];
    return result;
}

void CounterVector(vector<int> vec, int val) {
    for (int i = 0; i < vec.size(); i++)
        if (val == vec[i])
            cout << i << " ";
    cout << endl;        
}

void CreateFileHistogram(string title, string filename, vector<int> size, TH2D* hist2D) {
    string destination = "../" + filename;
    TCanvas* c1 = new TCanvas("c1", "2D Histogram", size[0], size[1]);
    // Draw the 2D histogram
    hist2D->Draw("COLZ");
    // Update the canvas to display the histogram
    c1->Update();
    // Save the canvas as an image file
    c1->SaveAs(destination.c_str());
    // Wait for user input to exit the application
    c1->WaitPrimitive();
    // Clean up
    delete c1;
}

void MakeHistogram(std::vector<int> info, std::string filename){
    std::string destination = "../" + filename;
    // Create a ROOT canvas
    TCanvas *c1 = new TCanvas("c1", "Histogram", 800, 600);
    // Create a histogram
    TH1F *hist = new TH1F("hist", filename.c_str(), info.size(), 0, info.size());
    // Fill the histogram with the info
    for (int i = 0; i < info.size(); ++i) {
        hist->SetBinContent(i+1, info[i]);
    }
    // Set histogram style
    hist->SetFillColor(38);
    // Draw the histogram
    hist->Draw();
    // Update the canvas to display the histogram
    c1->Update();
    c1->SaveAs(destination.c_str());
    // Wait for user input to exit the application
    c1->WaitPrimitive();
    delete hist;
    delete c1;
}

void MakeHistogram(std::vector<double> info, std::string filename){
    std::string destination = "../" + filename;
    // Create a ROOT canvas
    TCanvas *c1 = new TCanvas("c1", "Histogram", 800, 600);
    // Create a histogram
    TH1F *hist = new TH1F("hist", filename.c_str(), info.size(), 0, info.size());
    // Fill the histogram with the info
    for (int i = 0; i < info.size(); ++i) {
        hist->SetBinContent(i+1, info[i]);
    }
    // Set histogram style
    hist->SetFillColor(38);
    // Draw the histogram
    hist->Draw();
    // Update the canvas to display the histogram
    c1->Update();
    c1->SaveAs(destination.c_str());
    // Wait for user input to exit the application
    c1->WaitPrimitive();
    delete hist;
    delete c1;
}

void MakeHistogram2D(std::vector<std::vector<int>> info, std::string filename){
    std::string destination = "../" + filename;
    
    // Determine the dimensions of the histogram
    int nx = info.size(); // Number of bins along x-axis
    int ny = info[0].size(); // Number of bins along y-axis

    // Create a ROOT canvas
    TCanvas *c1 = new TCanvas("c1", "2D Histogram", nx, ny);
    
    // Set grayscale color palette
    gStyle->SetPalette(kGreyScale);
    
    // Create a 2D histogram
    TH2F *hist2D = new TH2F("hist2D", filename.c_str(), nx, 0, nx, ny, 0, ny);
    hist2D->SetMinimum(-0.1);
    
    // Fill the histogram with the info
    for (int col = 0; col < nx; col++) {
        for (int row = 0; row < ny; row++) {
            hist2D->Fill(col+0.5, ny-row-0.5, info[col][row]);
        }
    }

    int total = nx*ny;
    int count = 0;
    for (int col = 0; col < nx; col++) {
        for (int row = 0; row < ny; row++) {
            if(info[col][row] == hist2D->GetBinContent(col+0.5, ny-row-0.5)){count += 1;}
        }
    }
    
    std::cout << "TOTAL EQUALITY: " << (double)count/(double)total << std::endl;


    // Draw the 2D histogram
    hist2D->Draw("COLZ");
    
    // Update the canvas to display the histogram
    c1->Update();
    
    // Save the canvas as an image file
    c1->SaveAs(destination.c_str());
    
    // Wait for user input to exit the application
    c1->WaitPrimitive();
    
    // Clean up
    delete hist2D;
    delete c1;
}

// functions
void ReadImage(string filename, Image& img){
    ifstream inFile("../" + filename);

    //Make sure file is found
    if (!inFile) {
        cerr << "Error: Unable to open file for reading." << endl;
        return;
    }

    // Read PGM header
    inFile >> img.cod; // Read codification (P2)
    inFile >> img.ncols >> img.nrows; // Read number of columns and rows
    inFile >> img.N; // Read max color value (white)

    // Resize image matrix to match dimensions
    img.C.resize(img.ncols, vector<int>(img.nrows, 0));

    // Read pixel values
    for (int row = 0; row < img.nrows; row++) 
        for (int col = 0; col < img.ncols; col++) 
            inFile >> img.C[col][row];

    inFile.close();
}

void ReadImage(string filename, vector<vector<int>>& M, int& N){
    Image image;
    ifstream inFile("../" + filename);

    //Make sure file is found
    if (!inFile) {
        cerr << "Error: Unable to open file for reading." << endl;
        return;
    }

    // Read PGM header
    inFile >> image.cod; // Read codification (P2)
    inFile >> image.ncols >> image.nrows; // Read number of columns and rows
    inFile >> image.N; // Read max color value (white)

    N = image.N;

    // Resize image matrix to match dimensions
    image.C.resize(image.ncols, vector<int>(image.nrows, 0));
    M.resize(image.ncols, vector<int>(image.nrows, 0));

    // Read pixel values
    for (int row = 0; row < image.nrows; row++) {
        for (int col = 0; col < image.ncols; col++) {
            inFile >> image.C[col][row];
            //cout<< image.C[col][row] << std::endl;
            M[col][row] = image.C[col][row];
        }
    }
    inFile.close();
}

Image MakeImage(vector<vector<int>>& M, int& N){
    Image image;
    image.ncols = M.size();
    image.nrows = M[0].size();
    image.N = N;
    image.C = M;
    return image;
}

Image InvertImage(Image image){
    Image inverted = CopySettingsImage(image);

    for (int col = 0; col < image.ncols; col++) {
        for (int row = 0; row < image.nrows; row++) {
            inverted.C[col][row] = image.N - image.C[col][row];
        }
    }
    return inverted;
}

void WriteImage(string filename, const vector<vector<int>>& M, int& N){

    // Check if the image exists
    ofstream outFile("../" + filename, ios::binary);
    if (!outFile) {
        cerr << "Error: Unable to open file '" << filename << "for writing." << endl;
        return;
    }

    // Write the PGM file header
    outFile << "P2" << endl;
    outFile << M.size() << " " << M[0].size() << endl;
    outFile << N << endl;

    // Write pixel data to the file
    for (int row = 0; row < M[0].size(); row++){
        for (int col = 0; col < M.size(); col++){
            outFile << M[col][row] << " ";
        }
        outFile << endl; // Newline after each row
    }

    // Close the file
    outFile.close();
}

Image ImageSum(const Image& img1, const Image& img2) {
    //check if image to invert exists
    if(!SameSettings(img1,img2)){
        cout << SameSettings(img1,img2) << endl;
        cout << "Images have different settings" << endl;
        return img1; // returns the first image
    }

    Image image;
    image.ncols = img1.C.size();
    image.nrows = img1.C[0].size();
    image.N = img1.N;
    image.C.resize(img1.ncols, vector<int>(img1.nrows, img1.N));

    for (int col = 0; col < image.ncols; col++) {
        for (int row = 0; row < image.nrows; row++) {
            image.C[col][row] = img1.C[col][row] + img2.C[col][row];
        }
    }
    return image;
}

void addSegment(Image& image, vector<int> start, vector<int> end, int color, int size, string cluster_orientation) {
     //verify that the points are valid (x,y) points
    if (start.size() != 2 || end.size() != 2) {
        cout << "Error: Invalid start or end points." << endl;
        return;
    }

    //////////////// X IN THE MATRIX CORRESPONDS TO -Y CARTESIAN /////////////////////////////
    //////////////// Y IN THE MATRIX CORRESPONDES TO X CARTESIAN /////////////////////////////
    ////////////////  CORRESPONDES TO A ROTATION OF -90 DEGREES  /////////////////////////////
    // AN AUXILIARY FUNCTION WAS CREATED TO TAKE A POINT IN CARTESIAN AND DRAW IT IN MATRIX //

    //store the points
    int xi = start[0];
    int yi = start[1];
    int xf = end[0];
    int yf = end[1];

    //swap them to start from left to right
    if(xi>xf){
        xi = end[0];
        yi = end[1];
        xf = start[0];
        yf = start[1];
    }

    //calculate slope
    double slope=0;
    if(xf-xi!=0){slope= static_cast<double> (yf-yi)/(xf-xi);}

    double x = xi;
    double y = yi;

    //if the slope is >1, intersect with horizontal lines
    if(abs(slope)<=1){
        while (x <= xf) {
            for(int z=x-size/2;z<=x+size/2;z++){
                vector<int> c = {z,int(y)};
                vector<int> m = cart_to_matrix(image.nrows, c);
                image.C[m[1]][m[0]] = color;
            }
            x++;
            y=(x-xi)*slope+yi;
        }
    }
    //if the slope is <1, intersect with vertical lines
    if(abs(slope)>1){
        while (x <= xf) {
            for(int z=x-size/2;z<=x+size/2;z++){
                vector<int> c = {z,int(y)};
                vector<int> m = cart_to_matrix(image.nrows, c);
                image.C[m[1]][m[0]] = color;
            }
            x+=1/abs(slope);
            y+=slope/abs(slope);
        }
    }
}

void addSegment(Image& image, vector<int> start, vector<int> end, int color) {
    addSegment(image, start, end, color, 0, "H");
}

void AddNoise(Image& image, string noise_kind, double noise_par) {
    int color;

    if (noise_kind == "SP") { //iterate image to add SaltAndPepper
        for (int col = 0; col < image.ncols; col++) {
            for (int row = 0; row < image.nrows; row++) {
                if(GenerateRandomBool(noise_par)){image.C[col][row] = GetRandomExtreme(image.N);}
            }
        }
    }

    else if (noise_kind == "G") {
        for (int col = 0; col < image.ncols; col++) {
            for (int row = 0; row < image.nrows; row++) {
                int color = image.C[col][row] + (int)Gaussian(noise_par);
                if(color<0){image.C[col][row] = 0;}
                if(color>image.N){image.C[col][row] = image.N;}
                else{image.C[col][row] = color;}
            }
        }
    }

    else {
        cout << "Invalid noise kind" << endl;
        return;
    }
}

vector<int> GetColourFreq(Image& image){
    //create vector to store frequencies
    vector<int> AbsFrequencies;
    AbsFrequencies.resize(image.N+1, 0);

    //iterate image
    for (int col = 0; col < image.ncols; col++){
        for (int row = 0; row < image.nrows; row++) {
            AbsFrequencies[image.C[col][row]] += 1;
        }
    }
    return AbsFrequencies;
}

vector<double> GetColourRelFreq(Image& image){
    //get absolute frequencies and total number of pixels
    vector<int> AbsFrequencies = GetColourFreq(image);
    int Npixels = image.ncols * image.nrows;
    //to store relative frequencies
    vector<double> RelFrequencies;
    RelFrequencies.resize(image.N+1, 0);
    //iterate absolute frequencies to get relative frequencies
    for(int col = 0; col<image.N+1; col++){
        RelFrequencies[col] = (double)AbsFrequencies[col] / (double)Npixels;
    }
    return RelFrequencies;
}

double MedColor(Image image){
    int Sum = 0;
    int Npixels = image.ncols * image.nrows;
    for (int col = 0; col < image.ncols; col++){
        for (int row = 0; row < image.nrows; row++) {
            Sum += image.C[col][row];
        }
    }
    return ((double)Sum / (double)Npixels);
}

double Variance(Image image){
    double diff;
    double Media = MedColor(image);
    int Npixels = image.ncols * image.nrows;
    double Sum = 0;
    for (int col = 0; col < image.ncols; col++){
        for (int row = 0; row < image.nrows; row++) {
            diff = (double)image.C[col][row] - Media;
            Sum += diff*diff;
        }
    }
    return (Sum/(double)(Npixels-1));
}

TH2D* GetImageHistogram(Image& image) {
    string title = "2D Histogram";
    
    // Determine the dimensions of the histogram
    int nx = image.C.size(); // Number of bins along x-axis
    int ny = image.C[0].size(); // Number of bins along y-axis

    TH2D* hist2D = new TH2D("hist2D", title.c_str(), nx, 0, nx, ny, 0, ny);

    // Fill the histogram with the info
    for (int col = 0; col < nx; col++)
        for (int row = 0; row < ny; row++)
            hist2D->Fill(col+0.5, ny-row+0.5, image.C[col][row]);

    int total = nx*ny;
    int count = 0;
    for (int col = 0; col < nx; col++) {
        for (int row = 0; row < ny; row++) {
            if(image.C[col][row] == hist2D->GetBinContent(row+1, ny-col+1)){count += 1;}
        }
    }
    return hist2D;
}