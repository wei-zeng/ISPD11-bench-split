# ISPD11-bench-split
Split manufacturing challenge generator for ISPD 2011 benchmark (a variation of bookshelf format).

## Installation
    make

## Usage
    ./split-gen -design <design_name> -auxFile <aux_file> -rtFile <routing_solution_file> -layer -<split_layer> -outputRt <output_rt_file> -outputNets <output_nets_file> -outputKey <output_key_file> [-brokenNetsOnly] [-twoCutNetsOnly] [-floatingVpins] [-excludeNI]

The .aux file should contain a list of the .nodes, .nets, .wts, .pl, .scl, .shapes, and .route files for the benchmark. These files should all be located in the same directory as the aux file. 

Note: If you want to use a placement solution, you will have to replace the .pl file with the solution file and update the file name in the .aux file.

Note: The .route file and -rtFile are not the same; the .route file is provided in the benchmark and gives information for routing while the -rtFile is the actual global routing solution.

In our research paper [ObfusX (ASPDAC'21)](https://ieeexplore.ieee.org/document/9371556), we selected publicly available placement solutions and global router as described in https://github.com/wei-zeng/obfusX#readme, so that the design layouts can be overflow free. Routed designs are available at https://github.com/wei-zeng/obfusX/tree/master/circuits. Note that this tool only works for the superblue designs in ISPD'11 format.


This tool generates three files:

1) a modified .rt file with only the routes in public part, i.e., vpins and routes on or below the split layer. If a net is broken by the split, it will show as multiple nets in this .rt file. Note that (a) these are globe routes so the coordinates are the center of g-cells. (b) Not all wires end at a pin due to the split.

2) a modified .nets file where net pins are re-associated into new nets. Note that only pins in the original .nets are included. Each broken net will be shown as multiple new nets in this new .nets file. Each new net is connected to a v-pin. However, v-pins are not shown in this new .nets file for format compatibility reasons.

3) a .key file which states the correspondence of the new nets and the original nets. Each line is in the format of "new_net_name : original_net_name". If two or more new nets correspond to the same original net, they are connected in the original design.


There are four optional arguments that controls the split and output behaviors. All of them have a value of "false" if not specified. If specified, it will be set to "true".

`-brokenNetsOnly`: only nets that are broken by the split will be shown in .rt and .nets files.

`-twoCutNetsOnly`: only nets that are broken into exactly two parts by the split will be shown in .rt and .nets files.

`-floatingVpins`: include floating v-pins/part of net in .rt/.nets/.key files. A floating v-pin/broken part of net does not connect to any pin in the original design. It is "floating" as the wires above the split layer are removed.

`-excludeNI`: exclude NI terminals in the broken part of net. A NI terminal is a special fixed node as defined in ISPD 2011 contest. It is conjectured that it is not a cell pin.

We used options `-twoCutNetOnly` and `-excludeNI` in our research papers ([SplitMan TVLSI'19](https://ieeexplore.ieee.org/document/8789523) and [ObfusX ASP-DAC'21](https://ieeexplore.ieee.org/document/9371556)).

The ISPD 2011 benchmarks in this format can be found at http://www.ispd.cc/contests/11/ispd2011_contest.html#head-designs.

A description of the format of these files can be found at http://www.ispd.cc/contests/11/other_files/Benchmark_Format.pdf

The format of rt file is the same as the output format in ISPD 2008 contest and the output format of global routing solvers like NCTU-GR 2.0. A description can be found at http://www.ispd.cc/contests/08/ispd08rc.html#head-form.

## Troubleshooting
Problem 1: The program prints usage instead of running normally.

Solution 1: Please check all arguments (either required or optional) are correctly spelt. All of -auxFile, -rtFile, -layer, -outputRt, -outputNets, -outputKey are required.

Problem 2: File not found?

Solution 2: Please double-check the path to the file is correct; the path to the output file exists, and you have permission to write.

Problem 3: Who should I contact in case of further question?

Solution 3: Please direct questions to Wei Zeng and/or Prof. Azadeh Davoodi. E-mails are shown in our papers.

## Disclaimer
This tool is still in development. Bugs/issues may be found and fixed at any time.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
