<div align="center">

# Chainbox Robot — Sampling-Based Motion Planner

[![C++](https://img.shields.io/badge/C++-00599C?style=for-the-badge&logo=cplusplus&logoColor=white)](https://isocpp.org)
[![OMPL](https://img.shields.io/badge/OMPL-Motion_Planning-orange?style=for-the-badge)](https://ompl.kavrakilab.org)

*Planning collision-free motion for a 7-DOF chain-box robot in ℝ²×S¹×S⁴ configuration space — with self-collision detection and PRM\* clearance optimisation.*

</div>

---

## The robot

A **chain-box robot** — a chain of rigid links with a box at the end — operating in a 7-dimensional configuration space:

```
State: (x, y, θ, q₀, q₁, q₂, q₃)  ∈  ℝ² × S¹ × S⁴
```

The mixed topology (Euclidean + circular) makes standard Euclidean planners insufficient — requiring careful distance metrics and boundary handling.

---

## Implementation

### Custom validity checker
- Full **self-collision detection** between all link pairs
- Constraint: no two links of the chain may intersect each other
- Handles the compound geometry of boxes + connecting segments

### Scenarios

| Scenario | Planner | Config | Solve time |
|---|---|---|---|
| Narrow passage | RRTConnect | range = 0.9728 | ~20 s |
| Maximum clearance | PRM* | Clearance objective | — |

---

## Sampling strategy benchmark (PRM)

| Strategy | Behaviour |
|---|---|
| Uniform | Baseline, struggles in narrow passages |
| Gaussian | Biases samples near obstacle boundaries |
| Obstacle-Based | Concentrates in difficult regions |
| Bridge Test | Best for narrow corridors |

---

## Tech stack

`C++` · `OMPL` · `CMake`

---

<div align="center">
WPI Motion Planning (RBE 550) · <a href="https://github.com/Yami1106">Ashish Sukumar</a>
</div>
<!-- -->
