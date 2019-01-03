Host CPU-TinyGarble
=======
This is the host CPU of the [FPGA accelerator for Garbled Circuit](https://github.com/siamumar/MAXelerator).

This is based on the [TinyGarble](https://github.com/esonghori/TinyGarble) framework. Please read the README of TinyGarble for more details, especially on the dependencies, compilation, and general flow. A few of the features of TinyGarble (most notably <i>Skipgate</i>) is not available in this implementation. In addition, it does not include the GC optimized circuit synthesis library. Please go to the original TinyGarble repo for the libaray. 

TinyGarble main binary: `garbled_circuit/TinyGarble`
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
                                        the garbled tables.
  --output_mode arg (=0)                0: normal, 1:separated by clock 2:last
                                        clock.
```
For generating the reference files to test the HW acclerator, run `TinyGarble` without the `-w` flag. 

To generate the netlist from the synthesized Verilog files run `scd/V2SCD_Main`

```
  -h [ --help ]               produce help message.
  -i [ --netlist ] arg        Input netlist (verilog .v) file address.
  -b [ --brist_netlist ] arg  Input netlist (.txt) file address (in the format
                              given by www.cs.bris.ac.uk/Research/CryptographyS
                              ecurity/MPC/).
  -o [ --scd ] arg            Output simple circuit description (scd) file
                              address.
  -w [ --hscd ] arg           Output hardware simple circuit description (hscd)
                              file address.
  -p [ --pipe_stg ] arg (=10) Number of pipelined stages for non-XOR gates.
```

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
To turn on the binary mode undefine `HW_ACLRTR_TEXT_IO` in `garbled_circuit/garbled_circuit.h`
