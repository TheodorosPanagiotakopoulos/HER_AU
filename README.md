# Slow Growth Method for Modeling Activation Barriers Using the Trapezoidal Rule

This repository provides a computational implementation of the **slow growth method** for modeling the activation barriers of chemical reactions. The method employs the **trapezoidal rule** to numerically integrate along the reaction coordinate, providing an efficient approach to compute free energy changes.

---

## Table of Contents
1. [Introduction](#introduction)
2. [Background](#background)
3. [Methodology](#methodology)
    - [Slow Growth Method](#slow-growth-method)
    - [Trapezoidal Rule](#trapezoidal-rule)
4. [Installation](#installation)
5. [Usage](#usage)
6. [Example](#example)
7. [Contributing](#contributing)
8. [License](#license)

---

## Introduction
Chemical reactions often proceed through a transition state, where the system requires additional energy to overcome an **activation barrier**. The **slow growth method** offers a systematic way to compute free energy changes by gradually perturbing the system along the reaction coordinate. By combining this approach with the **trapezoidal rule**, we can achieve accurate numerical integration over discrete points.

---

## Background
- **Slow Growth Method**: A free energy calculation technique where the system is perturbed incrementally, with the assumption that the system remains quasi-equilibrated at each step.
- **Trapezoidal Rule**: A numerical integration method that approximates the area under a curve by dividing it into trapezoidal segments.

---

## Methodology

### Slow Growth Method
The reaction coordinate is divided into small steps, and the free energy difference (\( \Delta G \)) is computed iteratively. The system's potential energy surface is sampled at each step, allowing us to evaluate the incremental change in energy.

### Trapezoidal Rule
To integrate the force (\( F \)) along the reaction coordinate (\( \lambda \)), we use the trapezoidal rule:

    `ΔG = ∫ F(λ) dλ ≈ Σ (Δλ / 2) * [F(λ_i) + F(λ_i+1)]`
  Where:
- `λ` represents the reaction coordinate.
- `F(λ)` is the force acting along the reaction coordinate.
- `Δλ` is the step size.

---


---

---

## Installation
Clone the repository:
```bash
git clone https://github.com/<your-username>/<repository-name>.git
cd <repository-name>
