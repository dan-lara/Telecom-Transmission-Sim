# Telecom Transmission System & WebLab PA Linearization

## Project Overview
This project implements a complete telecommunication transmission chain using Matlab and integrates a hardware-in-the-loop component using the WebLab remote Power Amplifier (PA) testbench. The goal is to analyze signal transmission, noise impact, and PA nonlinearities, followed by Digital Predistortion (DPD) linearization.

**Institution:** Sorbonne Université  
**Course:** Project Télécom S9  
**Deadline:** January 5, 2026

## Project Structure

### Task I: Matlab Transmission Chain
Implementation of a TX/RX system including:
* **TX:** Random binary source, Modulation (BPSK/QPSK/16QAM), Upsampling/Filtering.
* **Channel:** Ideal channel initially.
* **RX:** Demodulation and BER calculation.

### Task II: Noise & Robustness
* Addition of AWGN (Additive White Gaussian Noise).
* Analysis of BER vs. SNR curves.
* Comparison of different modulation schemes.


### Task III: Power Amplifier (PA) & WebLab
* Hardware-in-the-loop testing via [WebLab](https://dpdcompetition.com/rfweblab/measurement-setup/).
* Measurement of EVM, ACPR, and Power Efficiency.
* **Optional:** Implementation of DPD (Digital Predistortion) to correct nonlinearities.

## Prerequisites
* **MATLAB:** (Ensure Communications and DSP System Toolboxes are installed).
* **WebLab Credentials:** Provided by the course instructors.

## Usage
1. **Simulation:** Run `main_simulation.m` to see the ideal Matlab chain.
2. **WebLab:** Open the `WebLab_Client` folder and run `main.m` to interact with the remote PA.

## Authors
* Daniel FERREIRA LARA
* Xiaochen YANG