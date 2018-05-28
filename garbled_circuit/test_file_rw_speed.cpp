#include <boost/format.hpp>
#include "garbled_circuit/garbled_circuit.h"

#include "scd/scd.h"
#include "scd/scd_evaluator.h"
#include "util/log.h"
#include "crypto/aes.h"
#include "crypto/BN.h"
#include "crypto/OT.h"
#include "crypto/OT_extension.h"
#include "garbled_circuit/garbled_circuit_low_mem.h"
#include "tcpip/tcpip.h"
#include "util/common.h"
#include "util/util.h"

#ifdef HW_ACLRTR
using std::ifstream;
using std::ofstream;
#endif

#define SIZE 512

int main(){
	
	string label_file;
	ifstream flin;
	string label = "";
	label_file = string(TINYGARBLE_SOURCE_DIR) + "/hw_aclrtr/Labels.txt";
	flin.open(label_file.c_str(), std::ios::in);
	if (!flin.good()) {
		cout << "file not found:" << label_file << endl;
		return -1;
	}	
	
	block* init_labels = nullptr;
	CHECK_ALLOC((init_labels) = new block[SIZE]);
	
	uint64_t text_rd_start_time = RDTSC;
	for (uint i = 0; i < SIZE; i++) {
		flin >> label;
		Str2Block(label, &init_labels[i]);
    }
	uint64_t text_rd_time =  RDTSC - text_rd_start_time;
	
	cout << text_rd_time << endl;	
	flin.close();

	ofstream of1;
	string ofile = string(TINYGARBLE_SOURCE_DIR) + "/hw_aclrtr/ofile";
	of1.open(ofile.c_str(), std::ios::out|std::ios::binary);
	if (!of1.good()) {
		cout << "file could not be created:" << ofile << endl;
		return -1;
	}	
	
	uint64_t bin_wt_start_time = RDTSC;
	//for (uint i = 0; i < SIZE; i++) {
		of1.write((char*)init_labels, sizeof(block)*SIZE); 
    //}
	uint64_t bin_wt_time =  RDTSC - bin_wt_start_time;
	of1.close();
	
	
	ifstream if1;
	string ifile = string(TINYGARBLE_SOURCE_DIR) + "/hw_aclrtr/ofile";
	if1.open(ifile.c_str(), std::ios::in|std::ios::binary);
	if (!if1.good()) {
		cout << "file was not found:" << ifile << endl;
		return -1;
	}	
	
	uint64_t bin_rd_start_time = RDTSC;
	//for (uint i = 0; i < SIZE; i++) {
		if1.read((char*)init_labels, sizeof(block)*SIZE); 
    //}
	uint64_t bin_rd_time =  RDTSC - bin_rd_start_time;	
	
	if1.close();
	
	for (uint i = 0; i < SIZE; i++) {
		printBlock(init_labels[i]);
    }
	
	cout << text_rd_time << endl;
	cout << bin_wt_time << endl;
	cout << bin_rd_time << endl;
	
	cout << "read time improvement:" << (double)text_rd_time/(double)bin_rd_time << endl;

	return 0;
}
