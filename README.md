# High-Throughput Multi-Engine Protein Mutation & Structural Dynamics Pipeline

A distributed bioinformatics workflow built to evaluate the structural and chemical impacts of single-point mutations using four deep-learning protein folding engines: AlphaFold2, ESMFold, RoseTTAFold2, and OmegaFold.

This pipeline replicates high-performance computing (HPC) workflows locally using Nextflow (DSL2) orchestration, hybrid API/local inference failovers, and automated PyMOL rendering tasks.

## Architecture Overview

### High-Level Data Flow

```
Plaintext    [ data/mutations.fasta ]
               │
               ▼ (splitFasta)
     Channel: (id, sequence)
               │
      ┌────────┼────────────────┬────────────────┐
      ▼        ▼                ▼                ▼
  [get_msa] [run_esmfold] [run_omegafold] [run_rosettafold]
      │        │                │                │
  (id, msa)    │                │                │
      │        │                │                │
  (join on id) │                │                │
      ▼        │                │                │
[run_alphafold]│                │                │
      │        │                │                │
      ▼        ▼                ▼                ▼
   (id, pdb) (id, pdb)       (id, pdb)        (id, pdb)
      │        │                │                │
      └────────┴────────┬───────┴────────────────┘
                        ▼ (.mix | .groupTuple)
               All PDBs grouped by ID
                        │
                        ▼
                [publish_results]
                        │
                        ▼
             [ results/ PDB Output Tree ]
                        │
                        ▼ (stage_inputs.py Manifest Step)
             [ comparison_metrics.nf ] ────► [ fetch_uniprot_annotations ]
                        │                                  │
                        ▼                                  ▼ (Value Channel Map)
             [ analyze_mutations.py ] ◄────────────────────┘
               (Biophysical Evaluation Engine)
                        │
                        ▼
             [ render_mutations.py ] (Automated PyMOL Wrapper)
                        │
                        ▼
             [ comparison_reports/ Data & Images ]
```

### Directory Structure

```
protein-folding/
├── folding.nf                 # Heavy structural prediction pipeline
├── comparison_metrics.nf      # Comparison analysis workflow 
├── nextflow.config.template   # Template Nextflow config for Docker and resources
├── Dockerfile                 # Container environment definition

Note: `nextflow.config` is not committed by default. Copy `nextflow.config.template` to `nextflow.config` and update the API secret before running the pipeline.
├── preprocess.py              # Data preprocessing utilities
├── staged_pairs.json          # Dynamic pairing data (generated)
│
├── bin/                       # Executable scripts
│   ├── stage_inputs.py        # Parse and organize PDB pairs for comparison
│   ├── get_msa.py             # Fetch MSA from NeuroSnap API (MMseqs2)
│   ├── get_active_sites.py    # Identify active sites via UniProt API
│   ├── analyze_mutations.py   # Compute mutation impact metrics
│   ├── render_mutations.py    # Visualize structural changes via PyMOL API
│   ├── run_alphafold.py       # AlphaFold2 inference wrapper
│   ├── run_esmfold.py         # ESMFold inference wrapper
│   ├── run_omegafold.py       # OmegaFold inference wrapper
│   └── run_rosettafold.py     # RosettaFold inference wrapper
│
├── data/                      # Input data
│   ├── test.fasta             # Test sequences
│   └── mutations.fasta        # Mutation sequences
│
├── results/                   # Output structures (generated)
│   ├── {PROTEIN_ID}_WT/       # Wild-type predictions
│   └── {PROTEIN_ID}_{MUTATION}/ # Mutant predictions
│
├── comparison_reports/        # Quantified physical chemistry metrics (.json)
├── comparison_images/         # Ray-traced structural alignment overlays (.png)
└── work/                      # Nextflow intermediate files (generated)
```

## Design Choices

### 1. Nextflow as Orchestrator
- Why: Enables reproducible, scalable execution across heterogeneous compute environments.
- Benefits:
  - DSL2 for modular process composition.
  - Built-in parallelization across multiple models.
  - Automatic error handling and retry logic.
  - Works with local, HPC, and cloud infrastructure.

### 2. Decoupled Model Drivers
- Each structure prediction engine operates via an independent Python wrapper instead of sharing a single script with complex parameter options.
- This encapsulates runtime dependencies, simplifies target tracking, and allows individual engines to be modified or swapped without changing the main workflow loops.

### 3. Local-to-API Topology Failover
- The pipeline checks the system execution space for local database footprints and model weights.
- If absent, it redirects tasks to cloud inference API endpoints.
- This enables rapid software development on lightweight environments while retaining bare-metal optimization for GPU clusters.

### 4. Separation of Concerns (Orchestration)
- Workflows are explicitly separated into two independent runs: `folding.nf` (heavy inference generation) and `comparison_metrics.nf` (downstream geometric calculations).
- This optimization isolates heavy computing steps so that altering downstream analysis thresholds, active-site targets, or visualization views does not re-trigger structural modeling.

### 5. Non-Blocking Value Channel Synchronization
- To resolve a 1-to-many stream pairing issue where linear Queue Channels tracking UniProt REST API queries would cause early pipeline termination, the pipeline uses a native Groovy `.reduce([:])` matrix accumulator step.
- This transforms single-use data streams into a static dictionary broadcasted to all 12 multiplexed computation tasks simultaneously.

## Docker Containerization
- The repository includes a `Dockerfile` that can build a runtime image with Python, Nextflow, and the required scientific dependencies.
- Nextflow is enabled for Docker with `docker.enabled = true`, and GPU access is requested using `containerOptions = '--gpus all'` for GPU-labeled processes.
- This image is a convenience helper for creating an isolated environment, but it is not required for host-based execution.
- If you are not fully containerizing the workflow, the Dockerfile can be treated as a dependency manifest rather than the only execution path.

### Docker Build & Usage
```bash
docker build -t protein-folding:latest .
# Use the image if you want an isolated runtime environment
docker run --rm -v "$PWD":"$PWD" -w "$PWD" protein-folding:latest nextflow run folding.nf --fasta data/mutations.fasta --outdir results
```

If you prefer to run Nextflow on the host, you can also skip the Docker image and run the pipeline directly once dependencies are installed.

## Model Personalities & Consensus Philosophy

Because deep-learning structural models are optimized to predict stable, baseline wild-type architectures, they exhibit distinct biases when evaluating single-point mutations. This pipeline leverages a multi-engine voting topology to counterbalance individual system limitations:

- **AlphaFold2 (The Evolutionary Anchor)**: Deeply dependent on co-evolutionary signals derived from Multiple Sequence Alignments (MSAs). It resolves conserved structural spaces with immense fidelity, but can over-respond to localized variant exceptions.
- **RoseTTAFold2 (The Physical Coordinate Anchor)**: Incorporates an explicit 3D structural track tracking continuous coordinate space throughout its network blocks. It is optimized to respect physical atom packing constraints and stereochemistry.
- **ESMFold (The Semantic Baseline)**: Uses a large language model trained on billions of unaligned protein sequences. It reads sequence context globally, making it highly resistant to local noise and excellent at preserving overall native fold topology.
- **OmegaFold (The High-Resolution Transformer)**: Combines a single-sequence language model with a GeoFormer geometry module. It skips MSA dependencies while displaying high sensitivity to localized loop patterns and dynamic shifts.

## Variant Structural Analysis Framework

To bridge the gap between pure structure prediction and functional variant effects, `analyze_mutations.py` evaluates the structural models using an automated biophysical decision matrix. The script performs a 3D Kabsch coordinate alignment between wild-type and mutant pairs, executing specific atomic screening steps:

- **Steric Clash Indexing**: Conducts all-atom spatial queries via a NeighborSearch tree to track non-bonded atoms that violate standard Van der Waals clearance boundaries (`(R_{vdwM} + R_{vdwW}) - 0.4Å`).
- **Chemical Property Classification**: Identifies electrostatic and polarity disruptions across sequence boundaries (e.g., Charged → Hydrophobic shifts).
- **Catalytic Proximity Mapping**: Maps coordinates retrieved from the UniProt database to the 3D model frame, calculating whether a point mutation sits within a strict `4.5Å` window of any critical active-site residue.

### Decision Flow

```
[Calculate Global & Local CA RMSD] ──► [Identify VdW Steric Clashes] ──► [Map Distance to UniProt Active Sites] ──► [Assign Final Classification]
```

## Benchmarked Case Studies

| Target Variant | Experimental Baseline | AlphaFold2 Behavior | ESMFold / OmegaFold Behavior | Pipeline Resolution Verdict |
|---|---|---|---|---|
| `1BTL_TEM1_E166A` (Beta-Lactamase) | Loss of Function. Truncates a critical catalytic residue; the 3D backbone remains perfectly rigid. | Neutral (0Å RMSD). Backbone ribbon tracks wild type perfectly; misses functional loss. | Neutral (0.016Å RMSD). Language model attention layers enforce extreme backbone rigidity. | Loss of Function (Catalytic Pocket Disruption). Custom decision logic flagged the <4.5Å spatial proximity to active site residues S70/K73 combined with the chemical property swap. |
| `2DN2_HBB_E6V` (Hemoglobin Beta) | Sickle Cell. Surface mutation driving intermolecular clumping; single-chain structure remains unaltered. | Backbone Distortion. Artificially unwinds a stable outer α-helix into an irregular beta-turn/loop. | Neutral. Language models correctly recognize that overall sequence context preserves the helix framework. | Neutral Variant. The mutation-site local RMSD calculation isolates AlphaFold's over-correction against the language model's stable baseline framework. |
| `2LZM_T4L_M99A` (T4 Lysozyme) | Core Cavity. Leaves an empty atomic void in the hydrophobic core; the outer framework remains folded. | Core Collapse. Catastrophically compresses surrounding helices inward to eliminate empty physical space. | Blind. Misses the internal structural cavity due to sequence-level generalizations. | Loss of Function (Severe Steric Hindrance / Core Collapse). The internal all-atom clash check and local loop tracking highlight the structural distortion introduced by the model's geometric optimization. |

## Workflow Summary

### 1. Input Stage
- Read FASTA file specified in `params.fasta` (default: `data/test.fasta`).
- Split into individual records via `splitFasta` to emit `(id, sequence)` tuples.

### 2. MSA Generation
- **Process**: `get_msa`
- Input: Sequence ID and protein sequence.
- Output: A3M alignment file.
- **Note**: Requires NeuroSnap API key configured in `nextflow.config`, typically by copying `nextflow.config.template` to `nextflow.config` and updating the secret.
- **Custom provider**: This pipeline is currently built for Neurosnap inference; if you use another provider, update the API key and endpoint URLs in the wrapper scripts accordingly.

### 3. Parallel Model Execution
- **AlphaFold2**: Consumes MSA and sequence; produces PDB coordinates.
  - Sanitizes MSA inputs to strip carriage returns and leading spaces.
- **ESMFold/OmegaFold/RosettaFold**: Accept sequence strings directly, bypassing MSA generation queues for fast initialization.

### 4. Result Publishing
- **Process**: `publish_results`
- Collects all generated PDB files per sequence structure using `.mix()` and `.groupTuple()`.
- Formats and copies files to the target `results/{ID}/` output directories alongside an execution summary log.

### 5. Downstream Analysis
- Run via `comparison_metrics.nf` configuration profiles.
- `comparison_metrics.nf` consumes `staged_pairs.json` and does not generate it automatically.
- Use `bin/stage_inputs.py` first if you need to create or refresh the manifest file, or supply an existing `staged_pairs.json`.
- The workflow then calls `analyze_mutations.py` to score spatial variables and `render_mutations.py` to execute targeted PyMOL ray-tracing routines.

## Key Configuration Parameters

### `nextflow.config`
- This file is created from `nextflow.config.template` and is not included in the repository by default.
- Copy the template first:
  ```bash
  cp nextflow.config.template nextflow.config
  ```
- Then update the API key and any local settings.
- **Docker**: Enabled at the config level, but container images are not explicitly declared in the repo.
- **GPU**: Configured with `--gpus all` for GPU-accelerated process tracks.
- **Resource Bounds**: Defaults to 6 CPUs, 6GB RAM per process task block (user-tunable).
- **API**: NeuroSnap key is required for MSA generation and any API-based inference fallback.

## Running the Pipeline

### Basic Execution
```bash
nextflow run folding.nf
```

### With Custom Input
```bash
nextflow run folding.nf --fasta data/mutations.fasta --outdir custom_results
```

### Resume Failed Runs
```bash
nextflow run folding.nf -resume
```

### Generate Report
```bash
nextflow run folding.nf -with-report execution_report.html
```

### Downstream Comparison Workflow
```bash
nextflow run comparison_metrics.nf --staged_json staged_pairs.json --outdir comparison_reports
```

If you need a new manifest, run `bin/stage_inputs.py` first or provide an updated `staged_pairs.json`.

## Output Structure

After execution:
```
results/
├── 1BTL_TEM1_WT/
│   ├── 1BTL_TEM1_WT_alphafold_result.pdb
│   ├── 1BTL_TEM1_WT_esmfold_result.pdb
│   ├── 1BTL_TEM1_WT_omegafold_result.pdb
│   ├── 1BTL_TEM1_WT_rosettafold_result.pdb
│   └── 1BTL_TEM1_WT_summary.txt
├── 1BTL_TEM1_E166A/
│   ├── 1BTL_TEM1_E166A_alphafold_result.pdb
│   ├── 1BTL_TEM1_E166A_esmfold_result.pdb
│   └── ...
```

## Troubleshooting

### MSA Generation Fails
- **Check**: NeuroSnap API key in `nextflow.config`.
- **Alternative**: Provide local MSA files or use local AlphaFold2 database.

### GPU Not Available
- **Check**: Docker daemon and NVIDIA runtime setup.
- **Fallback**: Remove GPU labels from processes in `nextflow.config`.

### Out of Memory
- **Adjust**: `memory` and `cpus` in `nextflow.config` based on hardware.
- **Reduce**: Batch size or use less demanding models.

## Dependencies

- **Nextflow**: >=21.0
- **Docker**: Enabled in config, but requires external images or runtime profile definitions.
- **NVIDIA GPU** (optional): For accelerated inference.
- **Python 3**: For preprocessing and analysis scripts.
- **Models**: AlphaFold2, ESMFold, OmegaFold, RosettaFold (Docker images)

## Future Enhancements

- [ ] Add confidence score extraction and reporting
- [ ] Implement structure quality metrics (Ramachandran, etc.)
- [ ] Support for multi-chain complexes
- [ ] Integration with AlphaFold-Multimer
- [ ] Automated model versioning and benchmarking
