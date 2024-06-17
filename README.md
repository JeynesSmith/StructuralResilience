This repository contains the computer code (written in Matlab R2019b on Windows) for all mathematical analysis, and for generating all figures, in the publication "Identifying and Explaining Resilience in Ecological Networks", submitted to Ecology Letters.

Authors: Cailan Jeynes-Smith, Michael Bode, and Robyn P. Araujo

An exhaustive search of three-node networks under the generalised Lotka-Volterra equations is performed by running CommandCentre.m. This function will generate the networks and data used for most of the analysis. This function calls upon dependents in a Data folder and scripts: AutomaticsRPAIdentification.m, ExaminingEcosystemGBLoop.m, MonomialAnalysis.m, PerturbedTimeSeries.m, StabilityFeasibilityTest.m, TestingValidityOfSystems.m.

Structural resilience can be tested for a single network by using SingleNetworkTesting.m. This is currently coded to examine a network with structural resilience from the manuscript and the Supplementary Methods. This function can be adjusted to test networks of any size and with any modelling framework, as long as the modelling framework matches the requirements of the Groebner basis calculation (must be written in a polynomial form at steady state). 

The networks with structural resilience are visualised and studied using clustering.m. This function is used for analysis and generates Supplementary Figure 2. This function requires data generated from PerturbedTimeSeriesForTestingAllSpecies.m and AllPossibleStructures.mat. 

Figures 1, 4,6, and Supplementary Figures 2 and 4 are generated using SimulatedSystems.m. This function uses parameters from Data.
Figure 5 is generated using Simulated_Actions.m. This function uses parameters from Data.

For any further questions, please contact cailan.jeynessmith@hdr.qut.edu.au.
