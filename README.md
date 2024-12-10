# Thermodynamic Model of Brayton Cycle (Enhanced)

This repository contains a project for the **Thermodynamics II (ME466)** course at **Bogazici University**, 
Mechanical Engineering Department. This project contains an engineering model to simulate a Brayton cycle 
with regeneration and intercooling using various methods in MATLAB. 

---

## Problem Description

In this project, a thermodynamic model of a Brayton cycle with enhancements, specifically with intercooling, is modeled.
Different from the basic version, this project considers the humidity in the air and the reaction at the 
combustion chamber. The effects of liquid water injection mass flow rate at the intercooler and the
fuel mass flow rate at the combustion chamber on the specific power output, energy, and exergy
efficiencies are assessed. The schematic of the analyzed system is shown below:

<p align="center">
  <img src="https://github.com/user-attachments/assets/00eaeb75-d1a2-4692-b60e-99f146a72f9c" alt="System Schematic" width="1000">
</p>

In this thermodynamic model:
- Intercooling is modeled as liquid water injection with varying mass flowrates.
- The combustion chamber is modeled as fuel combustion with corresponding psychrometry and reacting species.

---

## Assumptions

- **Working Fluid**:
  - The working fluid is an ideal gas mixture. The composition of the mixture at the first compressor inlet is assumed as a mixture of N<sub>2</sub>, O<sub>2</sub>, and H<sub>2</sub>O (no CO<sub>2</sub> or Ar).
  - Since atmospheric air includes these gases (as shown in the "*Atomic or Molecular Weights and Critical Properties of Selected Elements and Compounds*" table) and y<sub>N2</sub>/y<sub>O2</sub> = 3.76, dry air properties
  cannot be used directly throughout the engine model, especially for the psychometric analysis. Instead,
  c<sub>p,j</sub> values are taken from the "*Variation of (c<sub>p</sub>) with Temperature for Selected Ideal Gases*" table.

The following assumptions are made to simplify the modeling:

- Each component of the cycle is analyzed as a control volume at steady state.
- Each compressor stage is adiabatic.
- In the intercooler, liquid water is injected at P<sub>2</sub> and all injected water evaporates. 
- The intercooler is modeled as adiabatic.
- Neglect effects of potential and kinetic energies.
- Pure methane (CH<sub>4</sub>) is injected into the combustion chamber at P<sub>4</sub>, and the reaction is complete.
- Combustion occurs with 100% theoretical air.
- Heat transfer from the combustion chamber occurs at T<sub>max</sub>.
- The turbine stage is adiabatic.
- No dissociation or accumulation occurs in the combustion chamber. 

---

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

  Modeled as a function of temperature for N<sub>2</sub>, O<sub>2</sub>, CO<sub>2</sub>, and H<sub>2</sub>O using the table "*Variation of (c<sub>p</sub>) with Temperature for Selected Ideal Gases*".

- The isentropic turbine efficiency is defined as:

<p align="center">
  <img src="https://github.com/user-attachments/assets/1bdc7144-107c-4cd5-846e-bd6adf5744ab" alt="Isentropic Turbine Efficiency" width="450">
</p>

- The isentropic compressor efficiency is defined as:

<p align="center">
  <img src="https://github.com/user-attachments/assets/6285d0af-013f-4f2a-bb61-15104e3cc096" alt="Isentropic Compressor Efficiency" width="430">
</p>

- **Liquid and Vapor Water**:

  Liquid water is assumed to be an incompressible liquid with a constant specific heat of 4.179 kJ/kg.K.
  
  The saturation water vapor pressure is calculated as follows:
  <p align="center">
      <img src="https://github.com/user-attachments/assets/ba1f7541-ea8f-4aab-a2c3-c89f6c0b1c87" alt="Saturation Vapor Pressure" width="450">
  </p>

- **Enthalpy (h) and Entropy (s) Calculations**:

  Enthalpy and entropy are not taken from the tables but are calculated as follows:

  <p align="center">
    <img src="https://github.com/user-attachments/assets/19f027de-b2f3-44f2-ae0f-dd8868b2a4f9" alt="Enthalpy Equation" width="450">
  </p>

  where T<sub>ref</sub> = 298 K.

  <p align="center">
    <img src="https://github.com/user-attachments/assets/ff8e98be-c3a6-4b32-b9fb-9388acd3c0e0" alt="Entropy Equation" width="450">
  </p>

  where T<sub>ref</sub> = 298 K and P<sub>ref</sub> = 1 atm. The standard entropy values are taken from the table "*Thermochemical Properties of Selected Substances at 298K and 1 atm*".

- **Constants and Conversions**:

  The constants and conversions are given below:

  <p align="center">
    <img src="https://github.com/user-attachments/assets/22ba788e-e00c-48b5-a209-21bc374f26df" alt="Constants and Conversions" width="300">
  </p>

---


## Task

The following tasks are to be completed as part of this project:

### Analysis

1. **Base Case**:
   - Analyze the base case conditions and generate state and process tables for verification.

2. **Parametric Analysis 1**:
   - Vary the intercooler outlet relative humidity (ϕ) over a suitable range, considering at least 10 values.
   - Analyze and plot the system's performance metrics, including second law efficiency and backwork ratio, against liquid water injection mass flow rate.

3. **Parametric Analysis 2**:
   - Vary the maximum cycle temperature (T<sub>max</sub>) over a suitable range, considering at least 10 values.
   - Analyze and plot the system's performance metrics, including second law efficiency and backwork ratio, against fuel mass flow rate.

### Deliverables

1. **Base Case State and Process Tables**:
   - Produce state and process tables for the base case condition.
   - Verify the accuracy of the model using these tables.

2. **Parametric Analysis 1**:
   - Plot the following metrics against liquid water injection mass flow rate:
     - Second law efficiency (η<sub>2nd</sub>).
     - Backwork ratio (bwr).
   - Provide **two figures** for this analysis.

3. **Parametric Analysis 2**:
   - Plot the following metrics against fuel mass flow rate:
     - Second law efficiency (η<sub>2nd</sub>).
     - Backwork ratio (bwr).
   - Provide **two figures** for this analysis.

### Discussion

- Explain how the boundaries for relative humidity (ϕ) and maximum cycle temperature (T<sub>max</sub>) were determined. What parameters were considered?
- Describe the methodology used to calculate the exergetic efficiencies of the compressors, the combustor, the turbine, and the cycle.
- If the thermal efficiency of the system were to be calculated, how could it be defined?
