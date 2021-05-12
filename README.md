# Drinking water exposure to PFAS in Nurses Health Study

## Description
Source code for paper:
Tap Water Contributions to Plasma Concentrations of Poly- and Perfluoroalkyl Substances (PFAS) in a Nationwide Prospective Cohort of U.S. Women Environmental Health Perspectives (2019) 127:6 CID: 067006 [https://doi.org/10.1289/EHP4093](https://ehp.niehs.nih.gov/doi/full/10.1289/EHP4093)

## Contents
All the data management and modeling code can be found in the /src folder
#### src
> `main.R`
    Creates the analytical file.
> Step_1_Descriptive.R
    Generates the summary statistics of the data and produces Table 1 - 3 in the paper
> Step_2_GAM_model.R
    Fits a generalized additive model between serum PFAS and drinking water PFAS, while controling for confounder variables. Produces Table 4, Figure 2 in the model
> Step_3_TKmodel_1990.R
    Estimates plasma PFAS using drinking water PFAS measured in samples collected in 1989/1990 and a one-compartment toxicokinetic model. Produces Figure 3 in the paper.
> Step_4_MC_simualtion.R
    Monte Carlo simulation simulations to investigate how variability in input parameters affected the modeled relative source contribution of tap water. Produces Table S8, Figure S2 of the paper
> Step_5_UCMR_TK_2016.R
    Estimates plasma PFAS using drinking water PFAS measured in sample collected in 2013-2015 (UCMR3) and a one-compartment toxicokinetic model. Non-detect was replaced by detection limit divided by the square root of two. Produce Figure 4 in the paper.
> Step_6_UCMR_TK_2016.sens.R
    Sensitivity analysis of Step_5_UCMR_TK_2016.R where non-detect was replaced by zero.
> TOC_makemap.R
    Makes the map included in Figure 1.

## Authors
* [Xindi C. Hu](https://www.mathematica.org/our-people/staff/cindy-hu), Mathematica, Inc.; Department of Environmental Health, Harvard T.H. Chan School of Public Health, Harvard University; Harvard John A. Paulson School of Engineering and Applied Sciences, Harvard University.
* Andrea K. Tokranov, Harvard John A. Paulson School of Engineering and Applied Sciences, Harvard University.
* Jahred Liddie, Harvard John A. Paulson School of Engineering and Applied Sciences, Harvard University.
* Xianming Zhang, Harvard John A. Paulson School of Engineering and Applied Sciences, Harvard University.
* Philippe Grandjean, Department of Environmental Health, Harvard T.H. Chan School of Public Health, Harvard University; Institute of Public Health, University of Southern Denmark.
* Jaime E. Hart, Department of Environmental Health, Harvard T.H. Chan School of Public Health, Harvard University
* Francine Laden, Department of Environmental Health, Harvard T.H. Chan School of Public Health, Harvard University
* Qi Sun, Department of Epidemiology, Harvard T.H. Chan School of Public Health, Harvard University.
* Leo W. Y. Yeung, MTM Research Centre, School of Science and Technology, Örebro University.
* [Elsie M. Sunderland](https://bgc.seas.harvard.edu/), Department of Environmental Health, Harvard T.H. Chan School of Public Health; Harvard John A. Paulson School of Engineering and Applied Sciences, Harvard University.



