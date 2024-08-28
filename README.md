# KDE: **K**ernel **D**ensity **E**stimation for Structural Distributions of Molecular Systems

**KDE** is a program designed to calculate the structural distributions of a system using Kernel Density Estimation (KDE).
The program focuses on computing Radial Distribution Function (RDF), Bond Distribution Function (BDF), Angle Distribution Function (ADF), Torsion Distribution Function (TDF), and Improper Distribution Function (IDF) from molecular dynamics simulations.

## Overview

This application processes data from LAMMPS simulations, allowing users to analyze structural distributions in a material or molecular system.
The program supports both serial and parallel processing to handle large datasets efficiently.

### Key Features

- **RDF Calculation:** Analyze the spatial distribution of particles in the system.
- **BDF Calculation:** Measure the distribution of bond lengths between particles.
- **ADF Calculation:** Evaluate the distribution of angles formed by triplets of particles.
- **TDF Calculation:** Assess the distribution of torsion angles formed by quadruplets of particles.
- **IDF Calculation:** Compute the distribution of improper dihedral angles in the system.
- **Parallel Processing:** Leverage multi-threading for efficient computation on modern processors.

## Prerequisites

Ensure the following dependencies are installed before deploying the code:

- `libeigen3-dev`: A C++ template library for linear algebra.
- `cmake`: A cross-platform build system.
- `gcc-9` or newer: The GNU Compiler Collection.

## Installation

Follow these steps to set up and compile the KDE program:

1. **Clone the Repository:**

   ```bash
   git clone https://github.com/your-repo/kde.git
   cd kde
   ```

2. **Install Required Libraries:**

   Navigate to the `lib` directory and run the script to download necessary libraries:

   ```bash
   cd lib
   ./get_fmt.sh
   ```

3. **Create a Build Directory:**

   Set up a separate build directory to keep things organized:

   ```bash
   mkdir release
   cd release
   ```

4. **Compile the Program:**

   Use `cmake` to generate the build files, and then compile the program:

   ```bash
   cmake -DCMAKE_BUILD_TYPE=Release ..
   make -j4
   ```

5. **Run the Program:**

   Execute the program with an input file:

   ```bash
   ./kde ../data/kde.ini
   ```

## Usage

The program requires an input `.ini` file that specifies the paths to the LAMMPS data and dump files, as well as the settings for the different distribution calculations.

### Sample Input File (`kde.ini`)

Below is an example of what your input file should look like:

```ini
###################################################
#             Input for the kde code              #
###################################################

################### Inputs ########################
# Path to LAMMPS data file (*.lammps)
# lammps_data <path to LAMMPS format data file>
lammps_data     ../tests/alpha_a1b1c1.lammps

# Path to LAMMPS trajectory (dump) file (*.lammpstrj)
# lammps_dump <path to LAMMPS format trajectory file>
lammps_dump     ../tests/alpha_a1b1c1.lammpstrj

# Outputs tag ("CG" is default)
# tag <output tag string>
tag     cg

################ Type Definition ###################
# System phase (crystal, amorphous(default))
# phase <crystal/amorphous>
phase crystal

# Input file atom types definition
# type  <type number in data file> <type name>
type 	1	c1
type	2	c2
type	3	c3
type	4	h

# Bead types definition
# bead <bead name> <atom 1 type name> ...
# bead <bead name> <atom 1 weight> ...
bead	A	  c3 	c1 	  c2
weight	A	0.333 0.334 0.333

############## Distribution Settings ###############
# Processing method (parallel(default), serial)
# process serial
# process parallel <number of threads>
# NOTE: If number of threads left bland used max available
process parallel

# Radial Distribution Function (RDF)
# rdf <min> <cutoff> <steps> <KDE_bandwidth>
rdf 0.0 15.0 1000 0.015

# Special bonds (0 = false, 1 = true)
# special_bonds <1st_neighbor> <2nd_neighbor> <3rd_neighbor>
special_bonds 0 0 1

# Bond Distribution Function (BDF)
# bdf <min> <max> <steps> <KDE_bandwidth>
bdf 0.0 10.0 1000 0.01

# Angle Distribution Function (ADF)
# adf <min> <max> <steps> <KDE_bandwidth>
adf 0.0 180 1000 0.18

# Torsion Distribution Function (TDF)
# tdf <min> <max> <steps> <KDE_bandwidth>
tdf -180.0 180.0 1000 0.36

# Improper Distribution Function (IDF)
# adf <min> <max> <steps> <KDE_bandwidth>
idf 0.0 180 1000 0.18
```

### Running Tests

Sample data and input files are provided in the `data` folder for testing purposes.
To run a test adjust the input file `kde.ini` and then run:

```bash
./kde ../data/kde.ini
```

This will process the sample data and output the computed distributions to files.

## Documentation

Comprehensive documentation is available and can be generated using Doxygen. The documentation covers all classes, methods, and functionalities of the KDE program.

### Generating Documentation

To generate the documentation:

1. **Install Doxygen:**

   To generate the documentation using Doxygen, you'll need to have `Graphviz` installed.
   This is required for generating diagrams in the documentation.
   You can install it using the following commands:

- **Ubuntu/Debian**:
   ```bash
   sudo apt-get install doxygen
   sudo apt-get install graphviz
   ```

2. **Run Doxygen:**

   In the root directory of the project, execute:

   ```bash
   ./generate_docs.sh
   ```

3. **View Documentation:**

   The documentation will be generated in the `root` directory.
   Open the `index.html` file to view it in your browser.

## Contributing

If you wish to contribute to the project, please fork the repository and submit a pull request with your changes. Contributions are welcome!

## License

This project is licensed under the MIT License. See the `LICENSE` file for details.

## Acknowledgments

Special thanks to the open-source community for providing the tools and libraries that made this project possible.

