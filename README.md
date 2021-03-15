
# Papers

## Structural Basis for Calcium Sensing by GCaMP2

- Fluorophore: Ser-Tyr-Gly
- GFP (class 1)
  - A = dim = protonated = neutral: major
    - excitation 395 nm, emission 503 nm
  - B = bright = deprotonated = ionized (anionic)
    - excitation 475 nm, emission 508 nm
- EGFP (class 2) = GFP-S65T, Ser-65 to Thr
  - deprotonated
    - excitation 489 nm
- cpEGFP = C-EGFP + N-EGFP
  - protonated (body pH)
    - excitation 399 nm
- GCaMP2 (GECI) = (RSET +) M13 + cpEGFP + CaM
  - No Ca = protonated
    - excitation 399 nm
  - With Ca = deprotonated
    - excitation 488 nm
- Ca-free state is compact, predocked => rapid kinetics

## Structural mechanism of voltage-dependent gating in an isolated voltage-sensing domain


# Simulation results

- All simulations are at 5 ns per frame.
- Load in simulations by typing `vmd -e load.tcl` on the command line from within the directory.
- You can also do `source load.tcl` in the VMD GUI TkConsole (note that this has to be done from within the simulation directory).
- Starting structures for `ASAP3b` and `ASAP4c` can be found in the `models` folder. These models were prepared by Siri. 
- Note that the residue numbering may be shifted as compared to that used in literature. 
- For both ASAP3b and ASAP4c, I simulated both a neutral (A-GFP) and anionic (B-GFP) GFP chromophore.

# Notes

- amino acid mass: S 105, A 89, G 75, M 149, F 165
- amino acid hydropathy: S −0.8, A 1.8, G −0.4, M 1.9, F 2.8
- `413S`, `F413T`
- `Q414K`, `aa414`
- `R`
- F = PHE, Q = GLN

# Messages

- Goal: identify kinetic barriers to the S4 helix movement
- Questions to answer 

1. Why `413S` is faster than `A`, `G`, `M`, or the parental `F` in `ASAP4`. It is faster for both activation and inactivation. It's not just because it's smaller, as the even smaller `A` and `G` aren't fast. One hypothesis is its polarity prevents tight packing as `S4` moves up through the VSD pore, especially when `413S` encounters the hydrophobic plug. 
MD simulations to either support or reject this model would be useful. BTW `413S` also shows a wider distribution at rest to brighter fluorescence, which would suggest some of the VSD population is already depolarized, i.e. there's a change in the relative energies of depolarized and polarized conformations as well. 

2. Recently we added the observation that `F413T` `Q414K` is faster in the downward-going `ASAP3` which we continue to improve. `aa414` is the last gating charge position, i.e. some VSDs have an `R` in that position, but it is completely intracellular at resting (polarized) membrane voltage. So `Q414K` might provide force in the context of a positive transmembrane potential to overcome a kinetic barrier, but only after the `S4` has already moved halfway up, when `414K` is at the hydrophobic plug and able to sense the transmembrane field. Or maybe it's `413T` that's important. We can test the single mutants to figure this out experimentally, but MD simulations to figure out the atomic mechanism of the observed differences would be insightful I think, and demonstrate the utility of dynamic simulations in fast indicator development.

3. Even more useful would be if MD simulations can propose mutants we haven't tested yet, in the context of `F413` `Q414` (parent) for `ASAP3` and `ASAP4`.

4. In addition we have the baffling observation in `ASAP3` that mutation of `T399`, which is outside the membrane, to `R` has no change to kinetics, but when combined with the fast mutant `Q414K`, which is cytosolically 胞质地 located, the combination is slower than parent (`T399` `Q414`). I would be extremely surprised if MD can explain this, but if one is up for a challenge, that is certainly a challenging question.

- corrections: 
  - 143S is a mistake and should be 413S, consistent with F413 and F413T. 
  - Both 413 and 414 are numbered according to the ASAP aa sequence. Since you don't have a sequence map of ASAP, you can also identify them in PyMol. They should be 365 and 366 in PyMol numbering. 

