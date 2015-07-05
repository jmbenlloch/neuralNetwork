void root_hist(){
	ifstream testFile("testResults.txt");
	ifstream trainFile("trainResults.txt");

	double test[1437];
	double train[3354];
	double read;

	TH1D *histTest = new TH1D("test","x-x0 (Test set)",100,-60,60); 
	TH1D *histTrain = new TH1D("train","x-x0 (Train set)",100,-60,60); 

	for(unsigned int i=0; i<1437;i++){
		testFile >> read;
		cout << read << std::endl;
		test[i] = read;
		histTest->Fill(read);
	}
	for(unsigned int i=0; i<3354;i++){
		trainFile >> read;
		train[i] = read;
		histTrain->Fill(read);
	}

}
