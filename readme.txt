The attached codes generate part of the results presented in "Spatial Heterogeneity in Myelin Sheathing Impacts Signaling Reliability and Susceptibility to Injury, Talidou, A. & Lefebvre, J., eNeuro (2025), 12 (2) ENEURO.0402-24.2025; https://doi.org/10.1523/ENEURO.0402-24.2025".

-- The code "cortical.m" computes the membrane potential along axons of cortical myelin motifs with percentage of myelin coverage equal to 70%. Once the code runs, it returns the conduction delays and the corresponding conduction velocities of each axon. The axons are constructed using the Poisson distribution. For the numerical scheme, second order finite differences are used in space and first order backward Euler in time.

-- The code "callosal.m" computes the membrane potential along axons of callosal myelin motifs with percentage of myelin coverage equal to 70%. The axons are constructed using the Gamma distribution. The numerical scheme used is the same as in "cortical.m".

-- The code "corticalDemyelination.m" computes the membrane potential along healthy and demyelinated axons of cortical myelin motifs. The axons are constructed in the same way as in "cortical.m".

-- The code "callosalDemyelination.m" computes the membrane potential along healthy and demyelinated axons of callosal myelin motifs. The axons are constructed in the same way as in "callosal.m".

