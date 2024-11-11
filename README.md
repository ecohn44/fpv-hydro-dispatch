Welcome to the codebase for 'Lagrangian Informed Long-Term Dispatch Policy for Coupled Hydropower and Photovoltaic Systems'  

To access the input data used for simulations, refer to the folder data/fullsim.  

Notes on pre-processing:
1) Electricity Price Data: The historical price of electricity
was pulled from the CAISO regulated substation located at
Lake Mead (MEADS 2 N101). Hourly data is available for
download from June 2021 to present [12]. Therefore, 2022
and 2023 are the two full years of data available and is
the determining factor for the time frame of the simulations
presented in this case study.
2) Solar Radiation Data: The coordinates for Lake Mead
were used with NREL’s solar radiation database to pull the
historical hourly Direct Normal Irradation (DNI) from 1998
to 2022 [13]. The average DNI daily profile was calculated
for each month and normalized to create the ratio of available
solar radiation.
3) Hydrology Data: The daily average historical inflow rate
to Lake Powell is used as the system inflow. The daily average
historical outflow rate from Lake Mead is used as a proxy for
the water contract. The daily release was summed up over
each month to determine the monthly water contract. The
initial storage conditions used in the mass balance equation
combine the volume of both Lake Mead and Lake Powell
at the start of each month. All three datasets can be found
from the US Bureau of Reclamation Lower Colorado River
Operations Database [14].
4) Transmission Rating: The rating of the feeder capacity
is derived from 1,320 MW nameplate capacity of the Glen
Powell Dam [15]. Exploring the system behavior under a range
of feeder capacity ratings is explored in the sensitivity analysis.
5) FPV Field Rating: The proposed capacity of the FPV
field is 1 GW [16]. The average capacity density of floating
solar in the US is 10,000 m2/MW [17] and the total surface
area of both reservoirs is 500 square miles [18]- [15], trans-
lating to 0.75% coverage.
