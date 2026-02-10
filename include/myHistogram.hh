#ifndef myHistogram_h
#define myHistogram_h 1

#include <string>
#include <fstream>

// Written by Grant Berland, modified by Julia Claxton
// Header-only histogramming class to record energy deposition per altitude bin

class myHistogram
{
public:

  myHistogram(); // Constructor
  ~myHistogram(); // Destructor

  void AddCountToBin(unsigned int, double); // Fills histogram with data
  void WriteHistogramToFile(std::string); // Writes histogram to path provided

private:  
  // Array initialized to zeros (with fixed resolution)
  double histogramArray[1000] = {};
  double binEdges[101] = {};
};

// Inline constructor and destructor methods
inline myHistogram::myHistogram():histogramArray(){}
inline myHistogram::~myHistogram(){}

inline void myHistogram::AddCountToBin(unsigned int binAddress, double amountToAdd)
{
  histogramArray[binAddress] += amountToAdd;
}

inline void myHistogram::WriteHistogramToFile(std::string filename)
{
  // Open file
  std::ofstream outputFile;
  outputFile.open(filename, std::ios_base::out); // Open in write mode to overwrite any previous results

  // Write data
  outputFile << "altitude_km,energy_deposition_kev" << "\n"; // Header
  for(unsigned int i=0; i<1000; i++)
  {
    outputFile << i << "," << histogramArray[i] << "\n";
  }
  // Close file
  outputFile.close();
}

#endif