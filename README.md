# Slow Growth Method for Modeling Activation Barriers Using the Trapezoidal Rule

This repository provides a computational implementation of the **slow-growth method** for modeling the activation barriers of chemical reactions. The method employs the **trapezoidal rule** to numerically integrate along the reaction coordinate, offering an efficient approach to compute free energy changes.

---

## Table of Contents
1. [Introduction](#introduction)
2. [Background](#background)
3. [Theoretical Framework](#theoretical-framework)
4. [Methodology](#methodology)
    - [Slow Growth Method](#slow-growth-method)
    - [Trapezoidal Rule](#trapezoidal-rule)
5. [Installation](#installation)

---

## Introduction
Chemical reactions often proceed through a transition state, requiring additional energy to overcome an **activation barrier**. The **slow-growth method** offers a systematic way to compute free energy changes by gradually perturbing the system along the reaction coordinate. By combining this approach with the **trapezoidal rule**, we can achieve accurate numerical integration over discrete points.

---

## Background
- **Slow Growth Method**: A free energy calculation technique where the system is perturbed incrementally, assuming it remains quasi-equilibrated at each step.
- **Trapezoidal Rule**: A numerical integration method that approximates the area under a curve by dividing it into trapezoidal segments.

---

## Theoretical Framework

To calculate activation barriers, we applied the **slow-growth method**, which involves gradually altering a **collective variable** — a geometric or structural property characterizing the system's configuration.

A constraint is applied to ensure the collective variable changes linearly with each simulation time step, allowing the system to maintain a **quasi-equilibrium state**. For example, the collective variable might represent the distance between two atoms or molecules, incrementally adjusted at a fixed rate.

The work required to transition the system from one state to another is determined by integrating the force acting on the system as the collective variable changes:

\[
\Delta \Omega = \int_{\xi_1}^{\xi_2} \left( \frac{\partial \Omega}{\partial \xi} \right) \cdot \dot{\xi} \, d\xi
\]

Where:
- \( \Omega(x) \): The grand-canonical free energy as a function of the reaction coordinate \( \xi \).
- \( \frac{\partial \Omega}{\partial \xi} \): The force acting along the reaction coordinate.
- \( \dot{\xi} \): The time derivative of the reaction coordinate (rate of change).  

This method provides valuable insights into energy barriers and reaction pathways by simulating system behavior under controlled, gradual changes, making it a powerful tool for studying molecular interactions and chemical reactions.

---

## Methodology

### Slow Growth Method
The reaction coordinate is divided into small steps, and the free energy difference \( \Delta G \) is computed iteratively. The system's potential energy surface is sampled at each step, allowing for the evaluation of the incremental change in energy.

### Trapezoidal Rule
To integrate the force \( F(\lambda) \) along the reaction coordinate \( \lambda \), we apply the trapezoidal rule:

\[
\Delta G \approx \sum \left(\frac{\Delta \lambda}{2} \right) \cdot \left[F(\lambda_i) + F(\lambda_{i+1})\right]
\]

Where:
- \( \lambda \): The reaction coordinate.
- \( F(\lambda) \): The force acting along the reaction coordinate.
- \( \Delta \lambda \): The step size between points.

---

## Installation
Clone the repository:
```bash
git clone https://github.com/theodorosP/HER_AU.git
cd HER_AU
