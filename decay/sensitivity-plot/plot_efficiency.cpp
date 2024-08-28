#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "TCanvas.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TApplication.h"

int main(int argc, char** argv) {
    // Check if the correct number of arguments are provided
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <path_to_total_efficiency_output.txt>" << std::endl;
        return 1;
    }

    // Open the input file
    std::ifstream inputFile(argv[1]);
    if (!inputFile.is_open()) {
        std::cerr << "Error opening file: " << argv[1] << std::endl;
        return 1;
    }

    std::vector<double> masses;
    std::vector<double> efficiencies;
    double mass, efficiency;

    // Read the data from the file
    while (inputFile >> mass >> efficiency) {
        masses.push_back(mass);
        efficiencies.push_back(efficiency);
    }

    inputFile.close();

    // Create a TApplication object
    TApplication app("app", &argc, argv);



    // Create a canvas
    TCanvas* c1 = new TCanvas("c1", "J/psi Efficiency vs Mass", 800, 600);

    // Set logarithmic scales for both x and y axes
    c1->SetLogx();
    c1->SetLogy();

    // Create a TGraph from the vectors
    int nPoints = masses.size();
    TGraph* graph = new TGraph(nPoints, &masses[0], &efficiencies[0]);

    // Set graph titles and labels
    graph->SetTitle("J/psi Efficiency vs Mass;Mass (GeV/c^2);Efficiency");
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(0.8);
    graph->SetMarkerColor(kBlue);
    graph->SetLineColor(kBlue);

    // Draw the graph
    graph->Draw("APL");

    // Set the minimum and maximum for the y-axis
    double yMin = 1e-1;  // Replace with your desired minimum y-axis value
    double yMax = 1;    // Replace with your desired maximum y-axis value
    graph->GetYaxis()->SetRangeUser(yMin, yMax);

    // Update the canvas
    c1->Update();

    // Save the plot as a PNG file
    c1->SaveAs("efficiency_vs_mass.png");

    // Run the application
    app.Run();

    return 0;
}