# Products and Cohomology of Homogeneous Vector Bundles on Grassmannians  

**Companion repository to the paper**  
*Window categories for a simple $9$-fold flop of Grassmannian type*  
by Will Donovan, Wahei Hara, Michał Kapustka, and Marco Rampazzo  

---

## Overview  

This repository collects the Python code used to perform cohomological computations of homogeneous vector bundles on Grassmannians.  

The implementation provides:  
- Classes representing partitions and homogeneous vector bundles.  
- Tensor products and duals of homogeneous bundles.  
- Borel–Weil–Bott computations of cohomology.  
- Ext-groups on Grassmannians and their total spaces.  
- Examples and vanishing checks for tilting bundles.  

The code is written in Python and relies only on the external package **[LRCalc](https://sites.math.rutgers.edu/~asbuch/lrcalc/)** by Anders S. Buch, which implements the Littlewood–Richardson rule.  

---

## Installation  

Requirements:  
- **Python 3.12** (other versions may work but are untested)  

Install the required dependency via pip:  

```bash
pip install lrcalc
```


## Main Functions

### `utility.py`

This file contains utility functions.

* **`has_repeats(lst)`** → checks if a list has repeated entries.
* **`is_decreasing(lst)`** → checks if a list is strictly decreasing.
* **`swap_first_increase(lst)`** → applies one step of the Borel–Weil–Bott algorithm.
* **`add_rho(a)` / `subtract_rho(a)`** → adds or subtracts the `ρ`-shift for weights.
* **`normalize_partition(a)`** → normalizes a partition so the smallest entry is zero.
* **`weyl_dim(partition)`** → computes the dimension of a GL$_n$ representation via the Weyl formula.
---

### `homogeneous_bundles.py`

Here we gather the main classes.

* **`Partition`**
  Represents a partition with multiplicity. Overloads the `*` operator to compute the Littlewood–Richardson product using `lrcalc`.

  * `__mul__(other)` → list of `Partition`
  * `__str__`, `__repr__` → clean display of partitions with multiplicity.

* **`HomogeneousIrreducible(k, n, first_partition, second_partition, multiplicity=1)`**
  Represents an irreducible homogeneous vector bundle on the Grassmannian `G(k, n)`.

  * `__mul__(other)` → tensor product of two bundles, using Littlewood–Richardson.
  * `rank()` → rank of the bundle via Weyl’s dimension formula.
  * `dual()` → computes the dual bundle.
  * `cohomology()` → applies Borel–Weil–Bott to compute cohomology groups.

* **`HomogeneousDirectSum`**
  A wrapper class representing direct sums of `HomogeneousIrreducible` bundles.

  * Distributes multiplication (`*`) over direct sums.
  * `dual()` → dual direct sum.
  * `cohomology()` → computes cohomology of all summands.
---

### `exts_on_tilting_bundles.py`

Running this file provides the proof of Theorem 5.8.

* **`ext_on_grassmannian(E, F)`**
  Computes the Ext groups between two bundles on the Grassmannian using Borel–Weil–Bott.

* **`ext_on_total_space(k, n, E, F, bundle, cutoff=8)`**
  Computes Ext groups on the total space of either `U*(-2)` or `Q(-2)` by expanding the pushforward of the structure sheaf up to a cutoff.

* **`higher_ext_on_total_space(k, n, E, F, bundle, cutoff=8)`**
  Returns `True` if there are no higher Ext groups (i.e. only degree 0 remains).

* **`higher_cohomlogy_on_total_space(k, n, E, bundle, cutoff=8)`**
  Special case of `higher_ext_on_total_space` for checking higher cohomology of a bundle.

* **`is_it_tilting(k, n, E, bundle, cutoff)`**
  Tests whether a bundle is (partially) tilting by checking vanishing of self-Exts on the total space.

* **`FirstWindowG35`, `SecondWindowG35`**
  Implementation of the window categories ... as `HomogeneousDirectSum` objects.
---

### `vanishings.py`

A collection of computations of higher cohomology and Ext groups, reproducing results from the paper.
Running this file provides the proof for Propositions 5.9, 5.10, 5.22, 5.25.
---

