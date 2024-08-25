# KDE

**K**ernel **D**ensity **E**stimation of structural distributions of the system.

This program calculate the structural distributions (RDF, BDF, ADF, and TDF) 
of the given system using Kernel Density Estimation (KDE).

## Prerequisites

The following are required to deploy the code:

1) libeigen3-dev \n

2) cmake

3) gcc-9 or newer

## Deployment Guideline

1) Get the libraries

    - cd lib
    - ./get_fmt.sh

2) Create a build folder

    - mkdir release
    - cd release

3) Compile

    - cmake -DCMAKE_BUILD_TYPE=Release ..
    - make -j4

4) Deploy (inside the release folder)

    - ./kde [input file]

### Inputs

    - **input file** : Path to the input file (.ini) comprising settings.
    

