#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "TCanvas.h"
#include "TGraph.h"
#include "TMultiGraph.h"  // Include TMultiGraph
#include "TAxis.h"
#include "TApplication.h"
#include "TLegend.h"

int main(int argc, char** argv) {
    // Check if the correct number of arguments are provided
    if (argc != 8) {
        std::cerr << "Usage: " << argv[0] << " <path_to_total_efficiency_output1.txt> <path_to_total_efficiency_output2.txt> <path_to_total_efficiency_output3.txt> <path_to_total_efficiency_output4.txt> <path_to_total_efficiency_output5.txt> <path_to_total_efficiency_output6.txt> <path_to_total_efficiency_output7.txt>" << std::endl;
        return 1;
    }

    // Vectors and file reading for each input file
    std::vector<double> masses1, efficiencies1;
    double mass1, efficiency1;
    std::ifstream inputFile1(argv[1]);
    if (!inputFile1.is_open()) {
        std::cerr << "Error opening file: " << argv[1] << std::endl;
        return 1;
    }
    while (inputFile1 >> mass1 >> efficiency1) {
        masses1.push_back(mass1);
        efficiencies1.push_back(efficiency1);
    }
    inputFile1.close();

    std::vector<double> masses2, efficiencies2;
    double mass2, efficiency2;
    std::ifstream inputFile2(argv[2]);
    if (!inputFile2.is_open()) {
        std::cerr << "Error opening file: " << argv[2] << std::endl;
        return 1;
    }
    while (inputFile2 >> mass2 >> efficiency2) {
        masses2.push_back(mass2);
        efficiencies2.push_back(efficiency2);
    }
    inputFile2.close();

    std::vector<double> masses3, efficiencies3;
    double mass3, efficiency3;
    std::ifstream inputFile3(argv[3]);
    if (!inputFile3.is_open()) {
        std::cerr << "Error opening file: " << argv[3] << std::endl;
        return 1;
    }
    while (inputFile3 >> mass3 >> efficiency3) {
        masses3.push_back(mass3);
        efficiencies3.push_back(efficiency3);
    }
    inputFile3.close();

    std::vector<double> masses4, efficiencies4;
    double mass4, efficiency4;
    std::ifstream inputFile4(argv[4]);
    if (!inputFile4.is_open()) {
        std::cerr << "Error opening file: " << argv[4] << std::endl;
        return 1;
    }
    while (inputFile4 >> mass4 >> efficiency4) {
        masses4.push_back(mass4);
        efficiencies4.push_back(efficiency4);
    }
    inputFile4.close();

    std::vector<double> masses5, efficiencies5;
    double mass5, efficiency5;
    std::ifstream inputFile5(argv[5]);
    if (!inputFile5.is_open()) {
        std::cerr << "Error opening file: " << argv[5] << std::endl;
        return 1;
    }
    while (inputFile5 >> mass5 >> efficiency5) {
        masses5.push_back(mass5);
        efficiencies5.push_back(efficiency5);
    }
    inputFile5.close();

    std::vector<double> masses6, efficiencies6;
    double mass6, efficiency6;
    std::ifstream inputFile6(argv[6]);
    if (!inputFile6.is_open()) {
        std::cerr << "Error opening file: " << argv[6] << std::endl;
        return 1;
    }
    while (inputFile6 >> mass6 >> efficiency6) {
        masses6.push_back(mass6);
        efficiencies6.push_back(efficiency6);
    }
    inputFile6.close();

    // Add the seventh input file (Drell-Yan)
    std::vector<double> masses7, efficiencies7;
    double mass7, efficiency7;
    std::ifstream inputFile7(argv[7]);
    if (!inputFile7.is_open()) {
        std::cerr << "Error opening file: " << argv[7] << std::endl;
        return 1;
    }
    while (inputFile7 >> mass7 >> efficiency7) {
        masses7.push_back(mass7);
        efficiencies7.push_back(efficiency7);
    }
    inputFile7.close();

    // Create a TApplication object
    TApplication app("app", &argc, argv);

    // Create a canvas
    TCanvas* c1 = new TCanvas("c1", "Efficiency vs Mass", 800, 600);

    // Set logarithmic scales for both x and y axes
    c1->SetLogx();
    c1->SetLogy();

    // Create TGraphs from the vectors
    TGraph* graph1 = new TGraph(masses1.size(), &masses1[0], &efficiencies1[0]);
    TGraph* graph2 = new TGraph(masses2.size(), &masses2[0], &efficiencies2[0]);
    TGraph* graph3 = new TGraph(masses3.size(), &masses3[0], &efficiencies3[0]);
    TGraph* graph4 = new TGraph(masses4.size(), &masses4[0], &efficiencies4[0]);
    TGraph* graph5 = new TGraph(masses5.size(), &masses5[0], &efficiencies5[0]);
    TGraph* graph6 = new TGraph(masses6.size(), &masses6[0], &efficiencies6[0]);
    TGraph* graph7 = new TGraph(masses7.size(), &masses7[0], &efficiencies7[0]);  // For Drell-Yan

    // Create a TMultiGraph to hold all the TGraphs
    TMultiGraph* mg = new TMultiGraph();
    mg->Add(graph1);
    mg->Add(graph2);
    mg->Add(graph3);
    mg->Add(graph4);
    mg->Add(graph5);
    mg->Add(graph6);
    mg->Add(graph7);

    // Draw the multigraph
    mg->Draw("AL");

    // Set graph styles
    graph1->SetTitle("Efficiency vs Mass;Mass (GeV/c^2);Efficiency");
    graph1->SetLineColor(kBlue);
    graph1->SetLineWidth(2);

    graph2->SetLineColor(kRed);
    graph2->SetLineWidth(2);

    graph3->SetLineColor(kGreen);
    graph3->SetLineWidth(2);

    graph4->SetLineColor(kMagenta);
    graph4->SetLineWidth(2);

    graph5->SetLineColor(kCyan);
    graph5->SetLineWidth(2);

    graph6->SetLineColor(kOrange);
    graph6->SetLineWidth(2);

    graph7->SetLineColor(kBlack);
    graph7->SetLineWidth(2);

    // Set the x-axis and y-axis range for the multigraph
    double xMin = 1e-2;
    double xMax = 1e1;
    double yMin = 0.005;
    double yMax = 1;

    mg->GetXaxis()->SetRangeUser(xMin, xMax);  // Set x-axis range for the multigraph
    mg->GetYaxis()->SetRangeUser(yMin, yMax);  // Set y-axis range for the multigraph

    // Create a legend
    TLegend* legend = new TLegend(0.7, 0.155, 0.85, 0.45);
    legend->AddEntry(graph1, "J/psi", "l");
    legend->AddEntry(graph2, "Pion", "l");
    legend->AddEntry(graph3, "Rho", "l");
    legend->AddEntry(graph4, "Omega", "l");
    legend->AddEntry(graph5, "Phi", "l");
    legend->AddEntry(graph6, "Eta", "l");
    legend->AddEntry(graph7, "Drell-Yan", "l");
    legend->Draw();

    // Update the canvas
    c1->Update();

    // Save the plot as a PNG file
    c1->SaveAs("efficiency_vs_mass.png");

    // Run the application
    app.Run();

    return 0;
}
