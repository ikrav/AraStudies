void singleton(int entry, char *inpath, char *outpath) {

	// load reconstruction library
	gSystem->Load("$ARA_UTIL_INSTALL_DIR/lib/libAraEvent.so");

	// load source file into a tree
	TFile *infile = new TFile(inpath);
	TTree *intree = (TTree*) infile->Get("eventTree");

	// associate pointer with source data
	RawAtriStationEvent *event = 0;
	intree->SetBranchAddress("event", &event);

	// create new tree associated with the same event pointer
	TFile *outfile = new TFile(outpath, "RECREATE");
	TTree *outtree = new TTree("eventTree", "");
	outtree->Branch("event", &event);

	// retrieve singleton and store to output tree
	intree->GetEntry(entry);
	outtree->Fill();

	// write output tree to disk
	outfile->cd();
	outfile->Write();
	outfile->Close();
	infile->Close();
}
