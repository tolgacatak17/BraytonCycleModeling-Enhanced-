# Thermodynamic Model of Brayton Cycle (Enhanced)

This repository contains a project for the **Thermodynamics II (ME466)** course at **Bogazici University**, 
Mechanical Engineering Department. This project contains an engineering model to simulate a Brayton cycle 
with regeneration and intercooling using various methods in MATLAB.

---

## Problem Description

In this project, a thermodynamic model of a Brayton cycle with enhancements, specifically with intercooling, is modeled. 
The contribution of regeneration and the effects of ambient conditions on specific power output and energetic efficiency
are assessed. The schematic of the analyzed 
system is shown below:

<p align="center">
  <img src="https://github.com/user-attachments/assets/52474981-99f1-466a-b218-2957e48af1e4" alt="System Schematic" width="1000">
</p>

In this thermodynamic model:
- Intercooling is modeled as liquid water injection with varying mass flowrates.
- Combustion chamber is modeled as fuel combustion with corresponding psychrometry and reacting species.

---

## Assumptions

- **Working Fluid**:
  
  - The working fluid is an ideal gas mixture. The composition of the mixture at the first
compressor inlet is assumed as a mixture of N2, O2 and H2O (no CO2 and Ar).
  - Since atmospheric air includes these gaseous (as shown in "*Atomic or Molecular Weights and
  Critical Properties of Selected Elements and Compounds*" table) and yN2/yO2=3.76, dry air properties
  cannot be used directly throughout the engine model, especially for the psychometric analysis; instead,
  cp,j according to "*Variation of (c<sub>p</sub>) with Temperature for Selected Ideal Gases*" table should be used.

The following assumptions are made to simplify the modeling:

- Each component of the cycle is analyzed as a control volume at steady state.
- Each compressor stage is adiabatic.
- In the intercooler, liquid water is injected at P2 and all injected water evaporates. 
- Intercooler is modeled as adiabatic.
- Neglect effects of potential and kinetic energies.
- Pure methane (CH4) is injected to the combustion chamber at P4, and the reaction is complete.
- Combustion occurs with 100% theoretical air.
- Heat transfer from the combustion chamber occurs at Tmax.
- Turbine stage is adiabatic.
- No dissociation or accumulation occurs in the combustion chamber. 


## Model Inputs

The base case model inputs are provided below. These values will be varied parametrically for analysis:

<p align="center">
  <img src="https://github.com/user-attachments/assets/064c4199-f2b3-401f-ae9f-9529a08f19a0" alt="Model Inputs" width="700">
</p>


---

## Property Modeling

- **Reference Text**:

  All property values are based on "*Moran's Principles of Engineering Thermodynamics*".
  
- **Molecular Weights**:

  Taken from the table "*Atomic or Molecular Weights and Critical Properties of Selected Elements and Compounds*".
  
- **Specific Heat (c<sub>p,i</sub>):**

  Modeled as a function of temperature for N2, O2, CO2, and H2O using the table "*Variation of (c<sub>p</sub>) with Temperature for Selected Ideal Gases*".

- The isentropic turbine efficiency is defined as:

<p align="center">
  <img src="https://github.com/user-attachments/assets/1bdc7144-107c-4cd5-846e-bd6adf5744ab" alt="Isentropic Turbine Efficiency" width="450">
</p>

- The isentropic compressor efficiency is defined as:

<p align="center">
  <img src="https://github.com/user-attachments/assets/45fc8ed3-6ccc-419a-8803-de58cb68e134" alt="Isentropic Compressor Efficiency" width="450">
</p>


- **Enthalpy (h) and Entropy (s) Calculations**:

  Enthalpy and entropy are not taken from the tables but are calculated as follows:

  <p align="center">
    <img src="https://github.com/user-attachments/assets/18184d66-b8df-4b70-ab26-5e7d0562e32d" alt="Enthalpy Equation"   width="450">
  </p>

  where h = 0  at T<sub>ref</sub> = 0 K.


  <p align="center">
    <img src="https://github.com/user-attachments/assets/15081e5b-1524-4b42-b6cb-63400917025f" alt="Entropy Equation"   width="450">
  </p>

  where s(T<sub>ref</sub>, P<sub>ref</sub> = 1.70203 kJ/kg.K at T<sub>ref</sub> = 300 K and P<sub>ref</sub> = 1 atm.


- **Constants and Conversions**:

  The constants and conversions are given below:

  <p align="center">
    <img src="https://github.com/user-attachments/assets/22ba788e-e00c-48b5-a209-21bc374f26df" alt="Constants and       Conversions" width="300">
  </p>


---
