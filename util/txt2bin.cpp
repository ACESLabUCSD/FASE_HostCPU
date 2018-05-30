#include <boost/format.hpp>
#include <boost/program_options.hpp>
#include "garbled_circuit/garbled_circuit.h"

#include "util/log.h"
#include "util/common.h"
#include "util/util.h"

namespace po = boost::program_options;

#ifdef HW_ACLRTR
using std::ifstream;
using std::ofstream;
#endif

#define SIZE 512

int main(int argc, char* argv[]) {
	LogInitial(argc, argv);
	
	string text_file, bin_file;
	bool bin2text = false;
	uint64_t num_block = SIZE;
	
	boost::format fmter("Binary to text or text to binary conversion");
	po::options_description desc(fmter.str());
	desc.add_options()  //
	("help,h", "produce help message")  //
	("text,t", po::value<string>(&text_file)->default_value(string(TINYGARBLE_SOURCE_DIR) + "/hw_aclrtr/Labels.txt"), "text file")  //
	("bin,b", po::value<string>(&bin_file)->default_value(string(TINYGARBLE_SOURCE_DIR) + "/hw_aclrtr/Labels.bin"), "binary file")  //
	("num_block,n", po::value<uint64_t>(&num_block)->default_value(SIZE), "number of blocks")  //
	("bin2text,r", "binary to text conversion, the default is text to binary");

	po::variables_map vm;
	try {
		po::parsed_options parsed = po::command_line_parser(argc, argv).options(desc).allow_unregistered().run();
		po::store(parsed, vm);
		if (vm.count("help")) {
			std::cout << desc << endl;
			return SUCCESS;
	}
		po::notify(vm);
	}catch (po::error& e) {
		cout << "ERROR: " << e.what() << endl << endl;
		cout << desc << endl;
		return FAILURE;
	}
  
	if (vm.count("bin2text")){
		LOG(INFO) << "Converting binary to text";
		bin2text = true;
	}
	else{
		LOG(INFO) << "Converting text to binary" << endl;
	}
	
	ifstream fin, fcheck;
	ofstream fout;
	
	string row;
	block *labels = nullptr, *labels_check = nullptr;
	CHECK_ALLOC((labels) = new block[num_block]);
	CHECK_ALLOC((labels_check) = new block[num_block]);
	
	if(!bin2text){		
		fin.open(text_file.c_str(), std::ios::in);
		if (!fin.good()) {
			LOG(ERROR) << "file not found:" << text_file << endl;
			return -1;
		}			
		uint64_t text_rd_start_time = RDTSC;
		for (uint i = 0; i < num_block; i++) {
			fin >> row;
			Str2Block(row, &labels[i]);
		}
		uint64_t text_rd_time =  RDTSC - text_rd_start_time;		
		fin.close();

		fout.open(bin_file.c_str(), std::ios::out|std::ios::binary);
		if (!fout.good()) {
			LOG(ERROR) << "file could not be created:" << bin_file << endl;
			return -1;
		}			
		uint64_t bin_wt_start_time = RDTSC;
		fout.write((char*)labels, sizeof(block)*num_block); 
		uint64_t bin_wt_time =  RDTSC - bin_wt_start_time;		
		fout.close();
		
		fcheck.open(bin_file.c_str(), std::ios::in|std::ios::binary);
		if (!fcheck.good()) {
			cout << "file not found:" << bin_file << endl;
			return -1;
		}		
		uint64_t bin_rd_start_time = RDTSC;
		fcheck.read((char*)labels_check, sizeof(block)*num_block); 
		uint64_t bin_rd_time =  RDTSC - bin_rd_start_time;			
		fcheck.close();	

		bool err = false;
		for (uint i = 0; i < num_block; i++) {
			if(!CmpBlock(labels[i], labels_check[i])){
				err = true;
				LOG(ERROR) << "mismatch at label " << i << endl;
			}
#ifdef HW_ACLRTR_PRINT
			printBlock(labels_check[i]);
#endif
		}
		if(!err){
			LOG(INFO) << "conversion successful" << endl;
		}
		
		LOG(INFO) << "text read time:\t" << text_rd_time << endl;
		LOG(INFO) << "binary write time:\t" << bin_wt_time << endl;
		LOG(INFO) << "binary read time:\t" << bin_rd_time << endl;	
		LOG(INFO) << "read time improvement:" << (double)text_rd_time/(double)bin_rd_time << endl;
	}
	else{
		fin.open(bin_file.c_str(), std::ios::in|std::ios::binary);
		if (!fin.good()) {
			cout << "file not found:" << bin_file << endl;
			return -1;
		}		
		uint64_t bin_rd_start_time = RDTSC;
		fin.read((char*)labels, sizeof(block)*num_block); 
		uint64_t bin_rd_time =  RDTSC - bin_rd_start_time;			
		fin.close();

		fout.open(text_file.c_str(), std::ios::out);
		if (!fout.good()) {
			LOG(ERROR) << "file could not be created:" << text_file << endl;
			return -1;
		}			
		uint64_t text_wt_start_time = RDTSC;
		for (uint i = 0; i < num_block; i++) {
			printBlock(labels[i], fout);
		}
		uint64_t text_wt_time =  RDTSC - text_wt_start_time;		
		fout.close();	
		
		LOG(INFO) << "binary read time:\t" << bin_rd_time << endl;
		LOG(INFO) << "text write time:\t" << text_wt_time << endl;
		LOG(INFO) << "conversion successful" << endl;		
	}

	return 0;
}
