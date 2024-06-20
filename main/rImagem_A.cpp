/* COMPUTAIONAL PHYSICS - PROJECT 1
AUTHORS:
 - 106631, BHAVIN GULAB
 - 106643, JOSÉ GONÇALO MACHADO
 - 106881, YE JINGHAO
*/

#include "ImageProcessor.h"

// main
int main() {
    // Creates the images to store the information: codification, columns, lines, pixels and maximum white color
    Image moldura, letterA, molduraInverted, letterAInverted, invertedMolduraLetterA;

    // A.a. Read Image and set other images to be worked with
    ReadImage("Trab01_imagem_moldura.ascii.pgm", moldura);

    // A.b. Letter A
    CreateEmpty(letterA, 256 , 256, 255);
    vector<int> P1 = {68,28};
    vector<int> P2 = {128,228};
    vector<int> P3 = {188,28};
    vector<int> P4 = {98,128};
    vector<int> P5 = {158,128};
    addSegment(letterA, P1, P2, 0, 3, "H");
    addSegment(letterA, P3, P2, 0, 3, "H");
    addSegment(letterA, P1, P2, 0, 3, "H");
    addSegment(letterA, P4, P5, 0, 3, "H");
    
    // A.c. Inversion and creation of images
    molduraInverted = InvertImage(moldura);
    letterAInverted = InvertImage(letterA);
    WriteImage("Trab01_imagem_moldura_invertida.ascii.pgm", molduraInverted.C, molduraInverted.N);
    WriteImage("Trab01_imagem_A_invertido.ascii.pgm", letterAInverted.C, letterAInverted.N);

    // A.d. Sums the inverted 
    invertedMolduraLetterA = ImageSum(molduraInverted, letterAInverted);

    // A.e. Creates the image of sum of the inverted images
    WriteImage("Trab01_imagem_somada.ascii.pgm", invertedMolduraLetterA.C, invertedMolduraLetterA.N);

    // A.f Adds noise to the sum of the inverted image and creates a new image with noise
    AddNoise(invertedMolduraLetterA, "SP", 30);
    WriteImage("Trab01_imagem_somada_ruido.ascii.pgm", invertedMolduraLetterA.C, invertedMolduraLetterA.N);
    
    return 0;
}