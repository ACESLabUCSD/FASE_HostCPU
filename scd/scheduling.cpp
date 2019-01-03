/*
 This file is part of TinyGarble.

 TinyGarble is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 TinyGarble is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with TinyGarble.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "scd/scheduling.h"

#include <string>
#include <cstring>
#include <algorithm>
#include "scd/parse_netlist.h"
#include "util/common.h"
#include "util/log.h"
#include "garbled_circuit/garbled_circuit.h"

using namespace std;

///////////////////////////////////
enum Mark {
  UnMarked = 0,
  TempMarked = 1,  // Temporary mark
  PerMarked = 2   // Permanently mark
};

int TopologicalSortVisit(const ReadCircuit &read_circuit, vector<Mark>* marks,
                         int64_t current_unmark_index,
                         vector<int64_t>* sorted_list, vector<int64_t>* loop) {

  int64_t init_input_size = read_circuit.get_init_input_size();
  int64_t init_input_dff_size = init_input_size + read_circuit.dff_size;

  if (current_unmark_index < 0) {  // CONSTANT or IV 2nd input (-1) are sorted.
    return SUCCESS;
  } else if (marks->at(current_unmark_index) == Mark::TempMarked) {
    LOG(ERROR) << "There is a loop in the circuit." << endl;
    loop->push_back(current_unmark_index);
    return FAILURE;
  } else if (marks->at(current_unmark_index) == Mark::UnMarked) {
    marks->at(current_unmark_index) = Mark::TempMarked;

    if (TopologicalSortVisit(
        read_circuit,
        marks,
        read_circuit.gate_list[current_unmark_index - init_input_dff_size].input[0],
        sorted_list, loop) == FAILURE) {
      loop->push_back(current_unmark_index);
      return FAILURE;
    }
    if (TopologicalSortVisit(
        read_circuit,
        marks,
        read_circuit.gate_list[current_unmark_index - init_input_dff_size].input[1],
        sorted_list, loop) == FAILURE) {
      loop->push_back(current_unmark_index);
      return FAILURE;
    }

    marks->at(current_unmark_index) = Mark::PerMarked;
    sorted_list->push_back(current_unmark_index);
  }

  return SUCCESS;
}

int TopologicalSort(const ReadCircuit &read_circuit,
                    vector<int64_t>* sorted_list,
                    const ReadCircuitString& read_circuit_string) {

  int64_t init_input_size = read_circuit.get_init_input_size();
  int64_t init_input_dff_size = init_input_size + read_circuit.dff_size;

  sorted_list->clear();
  vector<Mark> marks(init_input_dff_size + read_circuit.gate_size,
                     Mark::UnMarked);

  // inputs are already sorted
  for (int64_t i = 0; i < init_input_dff_size; i++) {
    marks[i] = Mark::PerMarked;
    sorted_list->push_back(i);
  }

  while (true) {
    int64_t unmark_index = -1;
    // Everything is sorted.
    if (sorted_list->size() == init_input_dff_size + read_circuit.gate_size) {
      break;
    }
    for (int64_t i = 0; i < (int64_t) read_circuit.gate_size; i++) {
      //TODO: use a list to store unmarked.
      if (marks[init_input_dff_size + i] == Mark::UnMarked) {
        unmark_index = init_input_dff_size + i;
        break;
      }
    }
    if (unmark_index != -1) {  // There is an unmarked gate.
      vector<int64_t> loop;  // for detecting loop
      if (TopologicalSortVisit(read_circuit, &marks, unmark_index, sorted_list,
                               &loop) == FAILURE) {
        string loop_id_str = "";
        string loop_name_str = "";
        for (int64_t i = (int64_t) loop.size() - 1; i > 0; i--) {
          loop_name_str += read_circuit_string.gate_list_string[loop[i]
              - init_input_dff_size].output + "->";
          loop_id_str += std::to_string(loop[i] - init_input_dff_size) + "->";
        }
        loop_name_str += read_circuit_string.gate_list_string[loop[0]
            - init_input_dff_size].output;
        loop_id_str += std::to_string(loop[0] - init_input_dff_size);
        LOG(ERROR) << "Loop name: " << loop_name_str << endl << "Loop ID: "
                   << loop_id_str << endl;
        return FAILURE;
      }
    } else {
      break;  // there is no unmarked gate.
    }
  }

  CHECK_EXPR_MSG(
      sorted_list->size() == init_input_dff_size + read_circuit.gate_size,
      "Some gates are not sorted.");

  return SUCCESS;
}

#ifdef HW_ACLRTR
int Reorder(const ReadCircuit &read_circuit, vector<int64_t>* sorted_list, vector<int64_t>* reordered_list, uint64_t pipe_stg){
	
	int64_t init_input_size = read_circuit.get_init_input_size();
	int64_t init_input_dff_size = init_input_size + read_circuit.dff_size;
	
	//first check current completion time
	vector<uint64_t> marks(init_input_dff_size + read_circuit.gate_size, 0);
	for (int64_t i = 0; i < init_input_dff_size; i++) {
		marks[i] = pipe_stg;
	}
	
	int64_t sorted_index, input0_index, input1_index;
	short type;
	bool XOR_placed, nonXOR_placed;
	uint64_t placed = 0;
	uint64_t cycles = 0;
	
	while (true) {
		cycles++;	
		XOR_placed = false;
		nonXOR_placed = false;
		
		for (uint64_t i = 0; i < init_input_dff_size + read_circuit.gate_size; i++) {
			if(marks[i]) marks[i]++;
		}
		
		if (placed == read_circuit.gate_size) {// Everything is placed.
			auto marks_ = min_element(marks.begin(), marks.end());
			if(marks_[0] >= pipe_stg) break;
		}
		else{			
			for (uint64_t j = 0; j < 2; j++){
				sorted_index = sorted_list->at(init_input_dff_size + placed);			
				input0_index = read_circuit.gate_list[sorted_index - init_input_dff_size].input[0];
				input1_index = read_circuit.gate_list[sorted_index - init_input_dff_size].input[1];
				type = read_circuit.gate_list[sorted_index - init_input_dff_size].type;
				
				if(((input0_index < 0)||(marks[input0_index] >= pipe_stg)) && ((input1_index < 0)||(marks[input1_index] >= pipe_stg))){
					if((XOR_placed == false)&&((type == XORGATE)||(type == XNORGATE)||(type == NOTGATE))){
						placed++;
						marks[sorted_index] = pipe_stg;					
						XOR_placed = true;
						if (placed == read_circuit.gate_size) break;
					}
					else if(nonXOR_placed == false){
						placed++;
						marks[sorted_index] = 1;				
						nonXOR_placed = true;
						if (placed == read_circuit.gate_size) break;
					}
				}
			}
		}
	}
	
	uint64_t cycles_before = cycles;
	
	//now reorder
	memset(&marks[0], 0, marks.size() * sizeof marks[0]);
	reordered_list->clear();
	
	for (int64_t i = 0; i < init_input_dff_size; i++) {
		marks[i] = pipe_stg;
		reordered_list->push_back(i);
	}
	
	cycles = 0;
	while (true) {
		cycles++;	
		XOR_placed = false;
		nonXOR_placed = false;
		
		for (uint64_t i = 0; i < init_input_dff_size + read_circuit.gate_size; i++) {
			if(marks[i]) marks[i]++;
		}
		
		if (reordered_list->size() == init_input_dff_size + read_circuit.gate_size) {// Everything is sorted.
			auto marks_ = min_element(marks.begin(), marks.end());
			if(marks_[0] >= pipe_stg) break;
		}		
		else{			
			for (uint64_t i = 0; i < read_circuit.gate_size; i++) {
				sorted_index = sorted_list->at(init_input_dff_size + i);			 
				if(marks[sorted_index] >= 1) continue;
				
				input0_index = read_circuit.gate_list[sorted_index - init_input_dff_size].input[0];
				input1_index = read_circuit.gate_list[sorted_index - init_input_dff_size].input[1];
				type = read_circuit.gate_list[sorted_index - init_input_dff_size].type;
				
				if(((input0_index < 0)||(marks[input0_index] >= pipe_stg)) && ((input1_index < 0)||(marks[input1_index] >= pipe_stg))){
					if((XOR_placed == false)&&((type == XORGATE)||(type == XNORGATE)||(type == NOTGATE))){
						reordered_list->push_back(sorted_index);
						marks[sorted_index] = pipe_stg;					
						XOR_placed = true;
					}
					else if(nonXOR_placed == false){
						reordered_list->push_back(sorted_index);
						marks[sorted_index] = 1;				
						nonXOR_placed = true;
					}
				}
				if(XOR_placed & nonXOR_placed) break;
			}
		}			
	}
	
	uint64_t cycles_after = cycles;
	
	LOG(INFO) << "Completion time: before reordering " << cycles_before 
				<< ", after reordering " << cycles_after << ". Improvement: " << (double)cycles_before/(double)cycles_after << endl;

	return SUCCESS; 	
}
#endif

int SortNetlist(ReadCircuit *read_circuit,
                const ReadCircuitString& read_circuit_string
#ifdef HW_ACLRTR
			  , uint64_t pipe_stg
#endif
) {

  int64_t init_input_size = read_circuit->get_init_input_size();
  int64_t init_input_dff_size = init_input_size + read_circuit->dff_size;

  vector<int64_t> sorted_list;
  if (TopologicalSort(*read_circuit, &sorted_list, read_circuit_string) == FAILURE)
    return FAILURE;

#ifdef HW_ACLRTR
	if(pipe_stg > 1){
		vector<int64_t> reordered_list;
		if (Reorder(*read_circuit, &sorted_list, &reordered_list, pipe_stg) == FAILURE)
			return FAILURE;	
		sorted_list = reordered_list;
	}
#endif

  vector<int64_t> sorted_list_1(sorted_list.size());  //reverse sorted list

  for (int64_t i = 0; i < (int64_t) sorted_list.size(); i++) {
    sorted_list_1[sorted_list[i]] = i;
  }

  read_circuit->task_schedule.clear();
  read_circuit->task_schedule.resize(read_circuit->gate_list.size(), 0);
  for (int64_t i = init_input_dff_size; i < (int64_t) sorted_list.size(); i++) {
    read_circuit->task_schedule[i - init_input_dff_size] = sorted_list[i]
        - init_input_dff_size;  // align index
  }

  for (int64_t i = 0; i < (int64_t) read_circuit->dff_size; i++) {
    if (read_circuit->dff_list[i].input[0] > 0) {  // Constant values are negative
      read_circuit->dff_list[i].input[0] = sorted_list_1[read_circuit->dff_list[i]
          .input[0]];
    }
    if (read_circuit->dff_list[i].input[1] > 0) {  // Constant values are negative
      read_circuit->dff_list[i].input[1] = sorted_list_1[read_circuit->dff_list[i]
          .input[1]];
    }
    read_circuit->dff_list[i].output = sorted_list_1[read_circuit->dff_list[i]
        .output];
  }
  for (int64_t i = 0; i < (int64_t) read_circuit->gate_size; i++) {
    read_circuit->gate_list[i].output = sorted_list_1[init_input_dff_size + i];
    if (read_circuit->gate_list[i].input[0] > 0) {  // Constant values are negative
      read_circuit->gate_list[i].input[0] =
          sorted_list_1[read_circuit->gate_list[i].input[0]];
      CHECK_EXPR_MSG(
          read_circuit->gate_list[i].input[0] < read_circuit->gate_list[i].output,
          "input 0 is larger than gate id");
    }
    if (read_circuit->gate_list[i].input[1] > 0) {  // Constant values are negative
      read_circuit->gate_list[i].input[1] =
          sorted_list_1[read_circuit->gate_list[i].input[1]];
      CHECK_EXPR_MSG(
          read_circuit->gate_list[i].input[1] < read_circuit->gate_list[i].output,
          "input 1 is larger than gate id");
    }
  }
  for (int64_t i = 0; i < (int64_t) read_circuit->output_size; i++) {
    read_circuit->output_list[i] = sorted_list_1[read_circuit->output_list[i]];
  }

  if (read_circuit->terminate_id != 0) {
    read_circuit->terminate_id = sorted_list_1[read_circuit->terminate_id];
  }

  return SUCCESS;
}

int WriteMapping(const ReadCircuitString& read_circuit_string,
                 const ReadCircuit &read_circuit,
                 const string& out_mapping_filename) {

  std::ofstream f(out_mapping_filename, std::ios::out);
  if (!f.is_open()) {
    LOG(ERROR) << "can't open " << out_mapping_filename << endl;
    return FAILURE;
  }

  for (int64_t i = 0; i < (int64_t) read_circuit.dff_size; i++) {
    f << read_circuit_string.dff_list_string[i].output << " "
        << read_circuit.dff_list[i].output << endl;
  }

  for (int64_t i = 0; i < (int64_t) read_circuit.gate_size; i++) {
    int64_t gid = read_circuit.task_schedule[i];
    f << read_circuit_string.gate_list_string[gid].output << " "
        << read_circuit.gate_list[gid].output << endl;
  }

  f << "terminate " << read_circuit.terminate_id << endl;

  f.close();

  return SUCCESS;
}
