#include <string>
using namespace std;


void LineParser(ifstream&, vector<string>&);						//Function to turn a multiple lines file into a vector of strings 
void Parser(ifstream& file_to_parse, vector<vector<string> >& );	//Function to turn a vector of lines coming from a text in a vector
							 										//of vectors containing the fields. spreadsheet should be empty
void histProducer(string fileName, string histName, string outputFileName){

	bool verbose = false;

	ifstream inputFile(fileName);
	vector< vector<string> > spreadsheet;

	Parser(inputFile, spreadsheet);
	
	const int points = spreadsheet.size();
	float binXLow[points+1];
	float binCenter[points];
	float crossSection[points];
	
	if(verbose) cout << "Found " << points << " points" << endl;

	for(int i = 0; i <= points; i++){
		if (verbose) cout << "Line " << i << "/" << points <<endl;
		if (i == points){
			binXLow[i] = binXLow[i-1] +2*(binCenter[i-1] - binXLow[i-1]); 
			break;	
		}
		binXLow[i] = stof(spreadsheet.at(i).at(0));
		binCenter[i] = stof(spreadsheet.at(i).at(1));
		crossSection[i] = stof(spreadsheet.at(i).at(2));
		if(verbose ) cout << binXLow[i] << "\t" << binCenter[i] << "\t" << crossSection[i] << endl;
	}
// && TMath::IsNaN(crossSection[i])
	TH1F *hist = new TH1F(histName.c_str(),histName.c_str(),points,binXLow);
	for(int i = 0; i < points; i++) hist->Fill(binCenter[i],crossSection[i]);

	hist->SaveAs(outputFileName.c_str());
}

void LineParser(ifstream& file_to_parse, vector<string>& lines){
  if (!file_to_parse){
    printf("ERROR in LineParser > file not found.\n");
    exit(1);
  }
  while (true){
    if (file_to_parse.eof())break;
    string line;
    getline(file_to_parse,line);
    if (line[0] != '#' && line.size() > 0){
      lines.push_back(line);
    }
  }
}

void Parser(ifstream& file_to_parse,vector<vector<string> >& spreadsheet){
  vector<string> lines;
  LineParser(file_to_parse,lines);
  for (vector<string>::iterator line=lines.begin();line!=lines.end();++line){
    string field_content;
    istringstream linestream(*line);
    if (!linestream.eof()){
      spreadsheet.push_back(vector<string>());
      while(1){
	getline(linestream,field_content,' ');
	spreadsheet.back().push_back(field_content);
	if (linestream.eof()){break;}
      }
    }
  }
}