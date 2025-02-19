#include <fstream>

void getAgeo() {
    int nMass = 10;
    
    ofstream ageotxt("ageo_dy.txt");
    if(!ageotxt.is_open()) {
        std::cout << "Failed to make output file" << std::endl;
        return -1;
    }

    std::vector< vector<Double_t> > ageoRecords;


    for (int i = 1; i <=nMass; i++) {
        TFile *input = new TFile( Form("Events/mass_%d/unweighted_events.root", i));

	     TTree *tree = (TTree *) input->Get("LHEF");
         Int_t entries = tree->GetEntries();

	     Double_t dumpLength = 620;
	     Double_t px[4], py[4], pz[4], m[4], e[4];
	     Int_t id[4];

	     tree->SetBranchAddress("Particle.PID", &id);
	     tree->SetBranchAddress("Particle.Px", &px);
	     tree->SetBranchAddress("Particle.Py", &py);
	     tree->SetBranchAddress("Particle.Pz", &pz);
	     tree->SetBranchAddress("Particle.E", &e);
	     tree->SetBranchAddress("Particle.M", &m);

	     Int_t counter_on = 0;
	     Int_t counter_50 = 0;
	     Int_t counter_100 = 0;
         Double_t mass;
	     for (int i = 0; i < entries; i++) {
	     	tree->GetEntry(i);

	     	for (int j = 2; j < 4; j++) {
                mass = m[j];
	     		Double_t bx = px[j]/e[j];
	     		Double_t by = py[j]/e[j];

	     		Double_t x = bx * dumpLength;
	     		Double_t y = by * dumpLength;
	     		

	     		if (abs(x) < 0.5 && abs(y) < 0.5) 
	     			counter_on++;
	     		if (abs(x+0.5) < 0.5 && abs(y) < 0.5) 
	     			counter_50++;
	     		if (abs(x+1.0) < 0.5 && abs(y) < 0.5) 
	     			counter_100++;
	        }
	    }	

        vector<Double_t> ageo; 
        ageo.push_back(mass);
        ageo.push_back(counter_on * 1.0 / entries);
        ageo.push_back(counter_50 * 1.0 / entries);
        ageo.push_back(counter_100 * 1.0 / entries);
        ageoRecords.push_back(ageo);

        input->Close();
    }

    for (const auto &row : ageoRecords) {
        for (Double_t value : row)
            ageotxt << value << " ";
        ageotxt << std::endl;
    }

    ageotxt.close();
    return 0;
}
