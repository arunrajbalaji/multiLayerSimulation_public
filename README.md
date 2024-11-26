# multiLayerSimulation

## Purpose

A general-use solver for multi-layered electrochemical cells. The physical, transport, and chemical-kinetic properties of each “layer” are fully user-specifiable via a simple input file, allowing researchers to construct simulation domains with arbitrary arrangements of electrochemical components and to run efficient 1D simulations for a wide variety of applications without writing new code.

## What is included

This is a brief snapshot of parts of the code, pulled from the active research branch (a private repository belonging to the Mani Group). The intent is to provide a representative sample of my personal approach to scientific software development, as part of a professional portfolio. Here is a brief description of what you will find in each of the main folders/files:

- *src/*: main program files, with a wide variety of different sub-functions (including different versions that were used for experimentation with different schemes)
  - The most important file here is *main.m*, which drives the overall simulation
  - Several separately-defined functions are required to run *main.m*, many of which may be found in this folder
- *bsh/*: A variety of different bash files and SLURM scripts to facilitate batch-mode parallel simulations, sweeping large parameter spaces

## What is not included

Since publications related to this code have not been submitted yet, I have elected to withhold all of the post-processing and plotting components, in addition to certain components of the main simulation. In truth, post-processing and plotting code accounts for as much of the code in the active research repository as the actual simulation components. Additionally, several versions of the *main.m* file were used to prototype and develop different numerical methods and iteration schemes. I have only included intermediate versions of the file in this publicly-accessible repository.

## Running the simulation

I typically run the simulation in batch mode on HPC clusters, permitting quick sweeps of 1D or multi-D paramter spaces. The files to look at are located in *bsh/*. Check out [SLURM documentation](https://slurm.schedmd.com/documentation.html) to learn more about how these files work.

Also be sure to check out *src/runParallelJob.m*, where you can see exactly how parallel jobs are launched. This is also where you can learn details of SLURM scripts in order to sweep multiple paramters simultaneously, instead of just a single parameter.

## Input file formatting

The input for the simulation is provided in .JSON format. Here is a quick overview of the different sections of the input file. This will also help you get a quick sense of how flexible and adaptable the software is.

### Directory, time step, voltage, and structure

![Input file sample image 1!](/inputFileImages/input_1.png "Sample input file, figure 1.")

- Here, you can select the directory where the input file is located and where ouput files will be generated. You can also specify whether the simulation will be restarted from previous output files.
- It is possible to run using high-precision variables, if you anticipate large concentration spikes (i.e., with blocking electrode boundary conditions)
- Select the target time step, in addition to the time-intervals for output generation and plotting (mostly useful for real-time debugging purposes when running on an interactive node).
- You may choose to sweep the voltage to produce a voltammogram, but I usually run in potentiostatic mode.
- Specify the interface locations (one more interface than the total number of layers)

### Boundary conditions and physical constants

![Input file sample image 2!](/inputFileImages/input_2.png "Sample input file, figure 2.")

- You can modify this section to implement reactive (or blocking) electrode boundary conditions for each of the species.
  - Otherwise, each boundary condition will be a Dirichlet condition, set to the initial value for each species in the layer adjacent to the boundary
- Physical constants are set to standard values
- You can specify the voltage drop across the entire multi-layered structure (potentiostatic operation mode)

### Layer physical properties and aqueous chemistry

![Input file sample image 3!](/inputFileImages/input_3.png "Sample input file, figure 3.")

Here is an example of a single layer. Remember, the code allows you to stitch together as many layers as you would like, each with different physical properties and chemical compositions.

- Porosity and tortuosity are geometric characterizations of the porous structure
- "Background charge" is used to represent ion-selective membranes
- Each layer may have its own max and min mesh sizes, with smooth exponential stretching in between
  - This is necessary to resolve large gradients that occur near interfaces between layers with different properties
  - The gridSymmetry parameter permits users to specify whether mesh refinement is require near the left boundary alone, the right boundary alone, or both boundaries
- The list of species (with associated properties and initial values) may vary from layer to layer, as necessary for specific cases.

### Reaction mechanism

![Input file sample image 4!](/inputFileImages/input_4.png "Sample input file, figure 4.")

- The reaction mechanism in each layer may be uniquely specified, as shown here
- The species names must match those used in definitions above
- Forward and backward reactions are implemented as separate elementary reactions
- The Wien Coefficient was implemented to allow study of bipolar membranes, in which the large electric field at the junction may increase the effective rate coefficient for water splitting

## Mesh generation

## Boundary conditions

## Numerical methods

## Output format and processing

## Method of Manufactured Solutions

## Notes files

##
