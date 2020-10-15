# SMOKE [![DOI](https://zenodo.org/badge/39790736.svg)](https://zenodo.org/badge/latestdoi/39790736)

The Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling System was originally developed at MCNC to integrate emissions data processing with high-performance computing (HPC) sparse-matrix algorithms. SMOKE is now under active development at the Institute for the Environment and is partially supported by the Community Modeling and Analysis System (CMAS).

SMOKE is primarily an emissions processing system designed to create gridded, speciated, hourly emissions for input into a variety of air quality models such as CMAQ, REMSAD, CAMX and UAM. SMOKE supports area, biogenic, mobile (both onroad and nonroad), and point source emissions processing for criteria, particulate, and toxic pollutants. For biogenic emissions modeling, SMOKE uses the Biogenic Emission Inventory System, version 3.7 (BEIS3) and version 3.09, 3.14 and 3.61. SMOKE is also integrated with the on-road emissions model MOBILE6 and MOVES.

The sparse matrix approach used throughout SMOKE permits rapid and flexible processing of emissions data. Rapid processing is possible because SMOKE uses a series of matrix calculations rather than a less-efficient sequential approach used by previous systems. Flexible processing comes from splitting the processing steps of inventory growth, controls, chemical speciation, temporal allocation, and spatial allocation into independent steps whenever possible. The results from these steps are merged together in the final stage of processing using vector-matrix multiplication. This means that individual steps (such as adding a new control strategy, or processing for a different grid) can be performed and merged without having to redo all of the other processing steps.

The original SMOKE concept was envisioned in the early 1990's at MCNC by Dr. Carlie Coats. Marc Houyoux managed the development of SMOKE until his departure to the U.S. EPA Office of Air Quality Planning and Standards in 2002. With active-development continuing under the support from  U.S. EPA, lead SMOKE development was passed from Catherine Seppanen to Dr. B.H. Baek in 2005. The primary line of development has been managed by Dr. Baek under funding from the U.S. EPA since 2005.

# SMOKE-ready Input Data Files

SMOKE input data consist of emissions inventories, temporal and chemical speciation profiles, spatial surrogates, gridded meteorology and land use data, and other ancillary files for specifying the timing, location, and chemical nature of emissions. SMOKE is distributed with example data for getting started with the model. The example files distributed with SMOKE are for demonstration purposes only, they are not meant for real-world modeling applications.

The primary source for non-meteorology SMOKE input data is the U.S. EPA Clearinghouse for Inventories and Emissions Factors (CHIEF). The U.S. EPA Office of Air Quality Planning and Standards (OAQPS) Emissions Inventory and Analysis Group (EIAG) provides SMOKE inputs for different rule-making modeling platforms. These platforms include not only the NEI for both criteria air pollutants (CAPs) and hazardous air pollutants (HAPs), but also all of the SMOKE ancillary data files created by EPA for use in SMOKE. EPA uses CHIEF to provide these data.

Meteorology data must be generated for specific SMOKE applications using either MM5, WRF, or a similar model. The output data from meteorology models must be formatted for SMOKE using a program like MCIP.

