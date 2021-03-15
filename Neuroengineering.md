
# Neuroengineering 

## Notes

- [PDB](./notes/PDB.md)

## Thoughts

- To understand GEVI, [Molecular dynamics (MD)](https://en.wikipedia.org/wiki/Molecular_dynamics) and [Time-dependent density-functional theory (TDDFT)](https://en.wikipedia.org/wiki/Time-dependent_density_functional_theory) needs to be combined. Because:
  - `MD` only simulates atoms whereas electrons are ignored, this is known as [Born-Oppenheimer (BO) approximation](https://en.wikipedia.org/wiki/Born%E2%80%93Oppenheimer_approximation)
  - Whereas `TDDFT` is a theory of calculating electronic structures and energy levels (**excited states**), note that `DFT` can only calculate **ground state** electron structure 
  - The fluorescence level of `GEVI` changes when the external electric field changes, this is a photon absorption-emission process, which must be due to the jump of electrons between different energy levels. My guess for the mechanism is:
    - voltage change => conformational change of the `VSD` => conformation change of the `GFP` => energy levels change => spectrum shift => fluorescence level change
  - Therefore, understanding `GEVI` requires both `MD` and `TDDFT` calculations
  - Alternatively, we may explore a recent work from [Klaus-Robert Müller (h=125)](https://scholar.google.com/citations?user=jplQac8AAAAJ&hl=en): [Towards exact molecular dynamics simulations with machine-learned force fields](https://www.nature.com/articles/s41467-018-06169-2)

  
## 1. Groups

**! order by increasing invasiveness**

### A. Stanford

- Input (stimulation):
  - **sonomagnetic** [Remote and localized neural activation using sonomagnetic stimulation](https://neuroscience.stanford.edu/research/funded-research/remote-and-localized-neural-activation-using-sonomagnetic-stimulation)
    - Amin Arbabian `arbabian@stanford.edu` (2011, 24):
    - Stephen Baccus `baccus@stanford.edu`:
  - **ultrasonic drug uncaging** [Ultrasonic neural control and neuroimaging in the awake, mobile, and behaving small rodent](https://neuroscience.stanford.edu/research/funded-research/ultrasonic-neural-control-and-neuroimaging-awake-mobile-and-behaving-small)
    - Raag Airan `rairan@stanford.edu` (2010, 23):
      - 2019: Polymeric perfluorocarbon nanoemulsions are ultrasound-activated wireless drug infusion catheters
      - 2018: Noninvasive Ultrasonic Drug Uncaging Maps Whole-Brain Networks
      - 2018: Hearing out Ultrasound Neuromodulation
      - 2017: Neuromodulation with Nanoparticles
      - 2017: Noninvasive Targeted Transcranial Neuromodulation via Focused Ultrasound Gated Drug Release from Nanoemulsions
    - Jeremy Dahl `jeremy.dahl@stanford.edu`:
    - Butrus Khuri-Yakub `khuri-yakub@stanford.edu`:
      - 2020: Spike frequency–dependent inhibition and excitation of neural activity by high-frequency ultrasound
      - 2019: Radiation Force as a Physical Mechanism for Ultrasonic Neurostimulation of the Ex Vivo Retina
  - **photovoltaic polymer** [Injectable photovoltaics for a wireless, gliosis-free neural stimulation interface](https://neuroscience.stanford.edu/research/funded-research/injectable-photovoltaics-wireless-gliosis-free-neural-stimulation-interface)
    - Guosong Hong `guosongh@stanford.edu` (2014, 46):
      - 2020: Conjugated Polymers Enable a Liquid Retinal Prosthesis
    - Marion Buckwalter `marion.buckwalter@stanford.edu`:
    - Alberto Salleo `asalleo@stanford.edu`: 
  - **NeuroRoots** [NeuroRoots, brain/computer interface solution for paralysis](https://neuroscience.stanford.edu/research/funded-research/neuroroots-braincomputer-interface-solution-paralysis)
    - Nicholas Melosh `nmelosh@stanford.edu` (2001, 44):
    - Jaimie Henderson `unknown`:
  - **microwire** [Massively parallel microwire arrays for deep brain stimulation](https://neuroscience.stanford.edu/research/funded-research/massively-parallel-microwire-arrays-deep-brain-stimulation)
    - Jun Ding `dingjun@stanford.edu` (2007, 23):
      - 2020: Massively Parallel Microwire Arrays Integrated with CMOS chips for Neural Recording
    - Nicholas Melosh `nmelosh@stanford.edu`:
    
- Output (recording):
  - **GEVI** [Enabling faster and more responsive voltage imaging through computational biophysics](https://neuroscience.stanford.edu/research/funded-research/enabling-faster-and-more-responsive-voltage-imaging-through-computational)
    - Ron Dror `ron.dror@stanford.edu` (2002, 64):
  - **biosensor** [Real-time biosensors for measuring multiple neuromodulators](https://neuroscience.stanford.edu/research/funded-research/real-time-biosensors-measuring-multiple-neuromodulators)
    - Hyongsok Soh `tsoh@stanford.edu` (1999, 66):
      - 2020: Combinatorial Polyacrylamide Hydrogels for Preventing Biofouling on Implantable Biosensors
    - Karen Parker `karen.parker@stanford.edu`:

- Others:
  - **Neurogrid** [lab](https://web.stanford.edu/group/brainsinsilicon/index.html)
    - Kwabena Boahen `boahen@stanford.edu`:

## 3. Methods

### A. [GEVI: Genetically encoded voltage indicator](https://en.wikipedia.org/wiki/Genetically_encoded_voltage_indicator)

- It is thought to be superior to conventional voltage detecting methods like **electrode-based electrophysiological recordings**, **calcium imaging**, or **voltage sensitive dyes**. It can show neuron signals with **subcellular spatial resolution**. It has fast temporal resolution (**sub-millisecond**), matching or surpassing that of the electrode recordings, and about one magnitude **faster than calcium imaging**. Researchers have used it to probe neural communications of an intact brain (of Drosophila[41] or mouse[42]), electrical spiking of bacteria (E. coli[PROPS]), and human stem-cell derived cardiomyocyte.[43][44]
- [41, Michael Nitabach, Yale](https://www.cell.com/cell/fulltext/S0092-8674(13)00898-2?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867413008982%3Fshowall%3Dtrue), 
- [42, Thomas Knopfel, ICL](https://www.sciencedirect.com/science/article/abs/pii/S1367593115000691), 
- [44, Joseph Wu, Stanford](https://www.sciencedirect.com/science/article/pii/S1934590919300670#!)
- Michael Lin `mzlin@stanford.edu`: 
  - `read` 2020: Two-photon voltage imaging of spontaneous activity from multiple neurons reveals network activity in brain tissue
  - 2019: **[ASAP3]** Ultrafast two-photon imaging of a high-gain voltage indicator in awake behaving mice
    - Differs from **[ASAP2f]** by the **mutations L146G S147T N149R S150G H151D R414Q**
  - 2017: **[ASAP2s]** Fast two-photon imaging of subcellular voltage dynamics in neuronal tissue with genetically encoded indicators
    - Differs from **[ASAP1]** by the **mutation of Arg-415 to Gln (R415Q)**, more responsive to depolarization (= more positive), but reduced in speed
  - 2016: **[ASAP2f]** Subcellular Imaging of Voltage and Calcium Signals Reveals Neural Processing In Vivo
    - Differs from **[ASAP1]** by a **shorter S3-cpGFP linker**, more responsive to hyperpolarization (= more negative)
  - 2014: **[ASAP1]** High-fidelity optical reporting of neuronal electrical activity with an ultrafast fluorescent voltage sensor
    - cpsfGFP-OPT was inserted between residues 147 and 148 of GgVSD
    - E field -> VSD -> GFP -> decrease of intensity or shift of spectrum
  - 2013: **[VSFP-CR]** Improving FRET Dynamic Range with Bright Green and Red Fluorescent Proteins
- Mark Schnitzer `mschnitz@stanford.edu`:
  - 2015: **[Ace GEVIs, Ace2N-mNeon]** High-speed recording of neural spikes in awake mice and flies with a fluorescent voltage sensor
  - 2014: **[Mac GEVIs]** Imaging neural spiking in brain tissue using FRET-opsin protein voltage sensors
- Viviana Gradinaru `viviana@caltech.edu`:
  - 2014: **[Archer]** Archaerhodopsin variants with enhanced voltage-sensitive fluorescence in mammalian and Caenorhabditis elegans neurons
- Ehud Isacoff `ehud@berkeley.edu`:
  - 1997: **[FlaSh]** A Genetically Encoded Optical Probe of Membrane Voltage
- Adam Cohen `cohen@chemistry.harvard.edu`: 
  - 2019: **[(pa)QuasAr3(-s)]** Voltage imaging and optogenetics reveal behaviour-dependent changes in hippocampal dynamics
  - 2016: **[FlicR1]** A Bright and Fast Red Fluorescent Protein Voltage Indicator That Reports Neuronal Activity in Organotypic Brain Slices
  - 2014: **[QuasAr1, QuasAr2]** All-optical electrophysiology in mammalian neurons using engineered microbial rhodopsins
  - 2012: **[Arch]** Optical recording of action potentials in mammalian neurons using a microbial rhodopsin
  - 2011: **[PROPS]** Electrical Spiking in Escherichia coli Probed with a Fluorescent Voltage-Indicating Protein
- Francisco Bezanilla `fbezanilla@UChicago.edu`:
  - 2017: **[ASAP-Y]** Biophysical Characterization of Genetically Encoded Voltage Sensor ASAP1: Dynamic Range Improvement
- Vincent Pieribone `vpieribo@jbpierce.org`:
  - 2012: **[ArcLight]** Single Action Potentials and Subthreshold Electrical Events Imaged in Neurons with a Fluorescent Protein Voltage Probe
  - 2012: **[Zahra, Zahra 2]** Genetically encoded fluorescent voltage sensors using the voltage-sensing domain of Nematostella and Danio phosphatases exhibit fast kinetics
  - 2002: **[SPARC]** A Genetically Targetable Fluorescent Probe of Channel Gating with Rapid Kinetics
- Gero Miesenbock `gero.miesenboeck@cncb.ox.ac.uk`:
  - 2008: **[hVOS]** Rational Optimization and Imaging In Vivo of a Genetically Encoded Optical Voltage Reporter
- Thomas Knopfel `t.knopfel@imperial.ac.uk`:
  - 2012: **[VSFP-Butterfly]** Imaging neural circuit dynamics with a voltage-sensitive fluorescent protein
  - 2009: **[Red-shifted VSFP's]** Red-Shifted Voltage-Sensitive Fluorescent Proteins
  - 2008: **[VSFP3.1]** Engineering of a Genetically Encodable Fluorescent Voltage Sensor Exploiting Fast Ci-VSP Voltage-Sensing Movements
  - 2007: **[Flare]** Three fluorescent protein voltage sensors exhibit low plasma membrane expression in mammalian cells
  - 2007: **[VSFP2's]** Engineering and Characterization of an Enhanced Fluorescent Protein Voltage Sensor
  - 2001: **[VSFP1]** Design and characterization of a DNA‐encoded, voltage‐sensitive fluorescent protein

## 4. MD

- [https://github.com/google/jax-md] and [https://github.com/prody/ProDy]

## 5. Abbreviations

- ASAP : Accelerated Sensor of Action Potentials
- CaM  : Calmodulin
- GCaMP: GFP + CaM + M13
- GECI : Genetically encoded calcium indicator
- GEVI : Genetically encoded voltage indicator
- GFP  : Green fluorescent protein, cpsf~: circularly permuted superfolder, E~: enhanced
- IUE  : In utero electroporation
- M13  : A peptide sequence from MLCK
- MLCK : Myosin light-chain kinase
- RAMP : Random-access multiphoton microscopy
- REST : RE1-Silencing Transcription factor
- TPEF : Two-photon excitation microscopy
- ULoVE: Ultrafast local volume excitation
- VSD  : Voltage-sensing domain, Gg~: Gallus gallus, Ci~: Ciona intestinalis
- VSP  : Voltage-sensitive phosphatase

## js 

```javascript
const projects = [
  { title:'Sonomagnetic neural stimulation' },
  { title:'Ultrasonic neural control and neuroimaging' },
  { title:'Photovoltaics neural stimulation' },
  { title:'NeuroRoots' },
  { title:'Microwires neural stimulation' },
  { title:'Voltage imaging' },
  { title:'Biosensors for measuring neuromodulators' },
]
const profs = [
  { name: 'Amin Arbabian', email: 'arbabian@stanford.edu', projects: [0] },
  { name: 'Stephen Baccus', email: 'baccus@stanford.edu', projects: [0] },

  { name: 'Raag Airan', email: 'rairan@stanford.edu', projects: [1] },
  { name: 'Jeremy Dahl', email: 'jeremy.dahl@stanford.edu', projects: [1] },
  { name: 'Butrus Khuri-Yakub', email: 'khuri-yakub@stanford.edu', projects: [1] },

  { name: 'Guosong Hong', email: 'guosongh@stanford.edu', projects: [2] },
  { name: 'Marion Buckwalter', email: 'marion.buckwalter@stanford.edu', projects: [2] },
  { name: 'Alberto Salleo', email: 'asalleo@stanford.edu', projects: [2] },

  { name: 'Nicholas Melosh', email: 'nmelosh@stanford.edu', projects: [3, 4] },

  { name: 'Jun Ding', email: 'dingjun@stanford.edu', projects: [4] },

  { name: 'Ron Dror', email: 'ron.dror@stanford.edu', projects: [5] },
  { name: 'Michael Lin', email: 'mzlin@stanford.edu', projects: [5] },

  { name: 'Hyongsok Tom Soh', email: 'tsoh@stanford.edu', projects: [6] },
  { name: 'Karen Parker', email: 'karen.parker@stanford.edu', projects: [6] },
]
```

