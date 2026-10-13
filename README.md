# Yarrowia_lipolytica_hypoxia_GSMM: Enzyme-Constrained Metabolic Modeling of *Yarrowia lipolytica* Across Shake-Flask and Bioreactor Cultivations

A reproducible workflow for enzyme-constrained genome-scale metabolic modeling (ecGEM) of *Yarrowia lipolytica* using `eciYali5` in Python (COBRApy). This repository investigates the bioprocess scale-up gap between uninstrumented shake flasks and controlled stirred-tank bioreactors across exponential growth and nitrogen-limited production phases. 

By integrating liquid-phase carbon exchange fluxes (glycerol, citrate, isocitrate, and polyols via HPLC) and total cellular protein measurements (BCA), the pipeline models infers oxygen consumption in shake flasks, models proteome reallocation under hypoxia, and evaluates metabolic control (FCC) and genetic engineering targets (FSEOF) for organic acid overproduction.

---
## Repository structure
Work in progress...

---

## References, Model Source, and Background

This project builds upon the *eciYali5* GECKO 3.0 reconstruction of *Yarrowia lipolytica*:

* **Foundational Paper:** De Biaggi et al. (2026) – *Enzyme-constrained genome-scale model of Yarrowia lipolytica predicts growth-phase specific metabolic engineering targets*. Applied Microbiology and Biotechnology, 110:123. [doi:10.1007/s00253-026-13791-4](https://doi.org/10.1007/s00253-026-13791-4)
* **Model Repository:** Base YAML models (`ecPARe_pooled.yml` and `ecPARn_pooled.yml`) obtained from the [bioengtaltech/eciYali5-GEM](https://github.com/bioengtaltech/eciYali5-GEM) repository (Zenodo: [10.5281/zenodo.17434109](https://doi.org/10.5281/zenodo.17434109)).
* **Biological Context:** Prototrophic chassis *Y. lipolytica* PAR (JMY2900) / wild-type W29 cultivated on defined glycerol minimal medium ($2\%\text{ w/v}$, $\text{C/N} = 86$).

---
