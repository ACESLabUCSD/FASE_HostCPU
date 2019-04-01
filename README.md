# Host CPU of FASE powered by TinyGarble

This is the implementation of the host CPU of [`FASE`](https://github.com/siamumar/FASE) [1], 
an FPGA accelerator for Secure Function Evaluation (SFE) 
by employing the well-known cryptographic protocol named 
[Yao’s Garbled Circuit (GC)](https://en.wikipedia.org/wiki/Garbled_circuit).
The FPGA acclerator generates the garbled tables for Alice, the garbler.
The host CPU reads the generated garbled tables and executes Yao's protocol with Bob, the evaluator. 

This is based on the [TinyGarble](https://github.com/esonghori/TinyGarble) framework. Please read the README of TinyGarble for more details, especially on the dependencies, compilation, and general flow. A few of the features of TinyGarble (most notably <i>Skipgate</i>) is not available in this implementation. 

## Set up
- Install the dependencies. 
- Configure TinyGarble and then compile it in `bin` directory.
```
  $ ./configure
  $ cd bin
  $ make
```

- Set up the FPGA acclerator.

## Circuit Generation
The function to be executed through Yao's GC needs to be represented as a <i>netlist</i> of Boolean gates.
This repo includes the compiled Verilog netlists of a number of well-known benchmark functions.
They are extracted to the bin/scd/netlists directory during set up.
The steps to compile the behavioral Verilog code of any generic function to the netlist (also in Verilog)
along with the GC optimized synthesis library can be found at 
[`TinyGarbleCircuitSynthesis`](https://github.com/siamumar/TinyGarbleCircuitSynthesis).

The host CPU accepts the netlist in SCD format and the FPGA acclerator accepts the netlist in HSCD format.
Details of these formats are presented inside the [scd](/scd) directory.
To convert Verilog netlists of all the benchmark functions to SCD and HSCD formats run
```
$ cd bin/scd/
$ V2SCD_ALL.sh
```

To convert Verilog netlists of a specific benchmark function to SCD and HSCD formats run `bin/scd/V2SCD_Main`
```
  -h [ --help ]               produce help message.
  -i [ --netlist ] arg        Input netlist (verilog .v) file address.
  -o [ --scd ] arg            Output simple circuit description (scd) file
                              address.
  -w [ --hscd ] arg           Output hardware simple circuit description (hscd)
                              file address.
  -p [ --pipe_stg ] arg (=10) Number of pipelined stages for non-XOR gates.
```

#### Example:
```
$ cd bin/scd/
$ mkdir -p hw_aclrtr
$ mkdir -p hw_aclrtr/hamming_32bit_32cc
$ ./V2SCD_Main -i netlists/hamming_32bit_32cc.v -o netlists/hamming_32bit_32cc.scd -w  hw_aclrtr/hamming_32bit_32cc/Netlist.hscd --log2std
```

## Steps to securely evaluate a function through FASE (through Vivado Simulation)
1. Generate the garbled tables by following the intructions in `FASE`.

2. To execute Yao's protocol between Alice and Bob, run TinyGarble: `bin/garbled_circuit/TinyGarble`.
```
  -h [ --help ]                         produce help message
  -a [ --alice ]                        Run as Alice (server).
  -b [ --bob ]                          Run as Bob (client).
  -i [ --scd_file ] arg (=../scd/netlists/hamming_32bit_1cc.scd)
                                        Simple circuit description (.scd) file
                                        address.
  -p [ --port ] arg (=1234)             socket port
  -s [ --server_ip ] arg (=127.0.0.1)   Server's (Alice's) IP, required when
                                        running as Bob.
  --init arg (=0)                       Hexadecimal init for initializing DFFs.
  --input arg (=0)                      Hexadecimal input.
  -c [ --clock_cycles ] arg (=1)        Number of clock cycles to evaluate the
                                        circuit.
  --dump_directory arg                  Directory for dumping memory hex files.
  --disable_OT                          Disables Oblivious Transfer (OT) for
                                        transferring labels. WARNING: OT is
                                        crucial for GC security.
  --low_mem_foot                        Enables low memory footprint mode for
                                        circuits with multiple clock cycles. In
                                        this mode, OT is called at each clock
                                        cycle which degrades the performance.
  --output_mask arg (=0)                Hexadecimal mask for output. 0
                                        indicates that output belongs to Bob,
                                        and 1 belongs to Alice.
  -w [ --acc ]                          There is a HW accelerator generating
                                        the garbled tables.
  -d [ --acc_dir ] arg (=/hw_aclrtr)    Directory of HW accelerator generated
                                        garbled tables.
  --output_mode arg (=0)                0: normal, 1:separated by clock 2:last
                                        clock.
```
For generating the reference files to test the HW acclerator, run `TinyGarble` without the `-w` flag but with the `-d` flag. 

#### Example:

Alice's terminal
```
$ cd bin/garbled_circuit
$ ./TinyGarble -a -i ../scd/netlists/hamming_32bit_32cc.scd -w -d ../scd/hw_aclrtr/hamming_32bit_32cc/ --input AA -c 32 --log2std
```

Bob's terminal
```
$ cd bin/garbled_circuit
$ ./TinyGarble -b -i ../scd/netlists/hamming_32bit_32cc.scd --input F5 -c 32 --output_mode 2 --log2std
```

## Other binaries
By default, the generated garbled tables are stored in text format, which is easier to debug.
For faster operation, binary format should be used. 
To turn on the binary mode undefine `HW_ACLRTR_TEXT_IO` in `garbled_circuit/garbled_circuit.h`.
To convert text files to binary, run `util/txt2bin`
```
  -h [ --help ]                         produce help message
  -t [ --text ] arg (=/home/siam/git/hostCPU_TG/hw_aclrtr/Labels.txt)
                                        text file
  -b [ --bin ] arg (=/home/siam/git/hostCPU_TG/hw_aclrtr/Labels.bin)
                                        binary file
  -n [ --num_block ] arg (=512)         number of blocks
  -r [ --bin2text ]                     binary to text conversion, the default
                                        is text to binary
```

## References
[1] Siam U. Hussain and Farinaz Koushanfar, 
"[FASE: FPGA Acceleration of Secure Function Evaluation](http://aceslab.org/sites/default/files/FASE.pdf)",
<i>Field-Programmable Custom Computing Machines (FCCM)</i>, April, 2019.
