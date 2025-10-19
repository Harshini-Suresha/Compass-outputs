# COMPASS: Constraint-Based Metabolic Activity Scoring Framework

## 1. Project Overview

The **COMPASS** framework (COmputational Metabolic Pathway Activity Scoring System) is a computational platform designed to infer **metabolic reaction activities** directly from gene expression data at the single-cell or bulk resolution. It bridges the gap between **transcriptomics** and **metabolic modeling** by translating RNA-seq measurements into quantitative estimates of pathway utilization and reaction flux feasibility.  

At its core, COMPASS formulates an optimization problem based on **constraint-based modeling (CBM)** principles. It integrates **gene expression profiles**, **stoichiometric matrices**, and **genome-scale metabolic reconstructions (GEMs)** such as **RECON2**, **Human1**, and **Mouse1**, and solves a series of **linear programming (LP)** problems using the **Gurobi optimizer**. This approach identifies the minimal "penalty" or cost required for each reaction to carry flux under the observed transcriptional constraints.

By applying this process to thousands of reactions and cells, COMPASS produces a **metabolic activity score matrix** that quantitatively reflects the predicted metabolic state of each individual cell or sample.

---

## 2. Scientific Context and Motivation

Single-cell RNA sequencing (scRNA-seq) has revolutionized our understanding of cellular heterogeneity, but it remains limited in its ability to directly capture **metabolic fluxes** or **pathway activities**. Metabolism operates through enzyme-catalyzed reactions and pathway interconnections that cannot be inferred purely from gene-level expression data.

Traditional metabolic modeling approaches like **Flux Balance Analysis (FBA)** require defined objective functions and flux constraints, making them unsuitable for high-dimensional, noisy single-cell data. COMPASS overcomes this limitation by:

1. Reformulating metabolic inference as a **penalty minimization problem** rather than explicit flux optimization.  
2. Employing **reaction-specific gene-reaction rules (GPRs)** to map transcript abundance to reaction feasibility.  
3. Using **Gurobi-based LP solvers** to estimate the minimal resistance (penalty) required to sustain each reaction.  
4. Applying **λ-diffusion smoothing** to account for transcriptional stochasticity and biological similarity among neighboring cells.  

The result is a biologically interpretable **reaction-by-cell matrix** that reveals which pathways are active, repressed, or reprogrammed under specific conditions. This makes COMPASS particularly useful for studying **cellular metabolism in immunology, oncology, and developmental biology**, where metabolic rewiring plays a critical regulatory role.

---

## 3. Purpose of the Repository

This repository serves as a **complete implementation and documentation hub** for the COMPASS workflow, including installation instructions, environment configuration, model dependencies, execution commands, and downstream analysis scripts.

It is designed for computational biologists, bioinformaticians, and systems biologists interested in:

- Performing **metabolic activity profiling** of single-cell or bulk RNA-seq datasets.  
- Comparing metabolic pathways between different cellular conditions (e.g., pathogenic vs. non-pathogenic Th17 cells).  
- Integrating **Gurobi-based optimization** with **transcriptomic datasets**.  
- Running **Turbo-Compass** (accelerated implementation) for large-scale datasets.  
- Generating statistical summaries and visualizations such as volcano plots and reaction consistency heatmaps.  

The repository encapsulates not only the COMPASS executable framework but also extensive documentation, test data, example notebooks, and reference models, enabling full end-to-end reproducibility.

---

## 4. Core Workflow Summary

The COMPASS pipeline operates through five structured computational stages:

1. **Input Definition**  
   - Accepts gene expression data in `.tsv`, `.mtx`, or `.h5ad` formats.  
   - Maps gene identifiers to metabolic reactions via model-specific GPR rules.  
   - Supports human (`homo_sapiens`) and mouse (`mus_musculus`) species with built-in ortholog mapping.

2. **Model Selection and Initialization**  
   - Loads genome-scale metabolic reconstructions (e.g., RECON2.2, Human1, Mouse1).  
   - Builds a stoichiometric matrix and defines metabolite balance constraints.

3. **Penalty Scoring Optimization**  
   - Translates expression levels into reaction-specific penalty functions.  
   - Minimizes the penalty required to enable reaction flux using the **Gurobi solver**.  
   - Outputs per-cell penalty matrices (`penalties.txt.gz`).

4. **Diffusion-Based Smoothing (λ)**  
   - Applies penalty diffusion across transcriptionally similar cells.  
   - Produces the final **Compass score matrix (`reactions.tsv`)**, representing smoothed metabolic activity.

5. **Postprocessing and Statistical Analysis**  
   - Computes **activity scores** using the transformation:  
     \( Activity_{r,c} = -\log(Penalty_{r,c} + 1) \)  
   - Performs Wilcoxon rank-sum tests and calculates effect sizes (Cohen’s d).  
   - Integrates reaction metadata for pathway-level interpretation.

---

## 5. Technical Highlights

| Component | Description |
|------------|-------------|
| **Optimization Engine** | Gurobi Linear Programming solver used to minimize flux penalties. |
| **Model Basis** | Genome-scale metabolic models (Human1, Mouse1, RECON2, RECON2.2). |
| **Programming Language** | Python (compatible with versions ≤ 3.9 for TensorFlow and Gurobi compatibility). |
| **Parallelization** | Multi-process execution with user-defined cores (`--num-processes`). |
| **Turbo Mode** | Accelerated matrix-based LP approximation for large datasets. |
| **Penalty Diffusion** | λ-controlled smoothing using k-nearest neighbor (KNN) graph structures. |
| **Statistical Evaluation** | Nonparametric Wilcoxon tests and Cohen’s d computation. |
| **Output Format** | Reaction-by-cell matrix (TSV or compressed), metadata-linked outputs. |

---

## 6. Relevance and Application

COMPASS has been successfully applied in multiple single-cell studies, most notably in the **analysis of Th17 immune cell metabolism** (Wagner et al., *Cell*, 2021). It enables the identification of metabolic programs underlying cell differentiation, immune activation, and disease-specific phenotypes.  

Through its integration with genome-scale metabolic models and high-performance solvers, COMPASS facilitates:
- Identification of **pathway-level metabolic shifts**.
- Comparative metabolic profiling across experimental conditions.
- Generation of interpretable, model-based metrics from transcriptomic data.
---
## 7. Installation Requirements and Environment Configuration

The COMPASS framework requires a **controlled Python environment** with a specific dependency structure, due to the compatibility requirements of TensorFlow, Gurobi, and certain scientific computing libraries. It is **strongly recommended** to use **Anaconda (or Miniconda)** for managing environments and installing dependencies.

All installation steps should be executed from an **Anaconda Prompt (Windows)** or **terminal (macOS/Linux)** with **administrator privileges**.

---

### 7.1. System Prerequisites

Before proceeding with installation, ensure the following software and system configurations are present:

| Requirement | Minimum Version | Notes |
|--------------|----------------|-------|
| **Operating System** | Windows 10 / macOS / Linux | Tested on Windows and Linux systems |
| **Python** | 3.8.x | TensorFlow and Gurobi are compatible only with Python ≤ 3.9 |
| **Anaconda / Miniconda** | Latest stable | Used for isolated environment management |
| **Gurobi Optimizer** | ≥ 11.0.0 | Required for linear programming optimization |
| **Git** | Latest stable | Required for cloning COMPASS from GitHub |

---

### 7.2. Setting Up a Dedicated Environment

To ensure reproducibility and dependency isolation, create a new Conda environment specifically for COMPASS.

```bash
# Create new environment with compatible Python version
conda create -n compass_env python=3.8

# Activate the environment
conda activate compass_env

```
---
### 7.3. Verifying the Active Environment

You can verify the active Conda environment using:
conda env list
If you need to deactivate at any point:
```bash
conda deactivate
```
---
### 7.4. Installing Core Computational Libraries
The following packages are essential for data manipulation, optimization, and numerical analysis.

```bash
conda install pandas numpy scipy matplotlib seaborn scikit-learn
```

Library Functions:

### Core Computational Libraries

| Library | Description |
| :--- | :--- |
| **pandas** | Provides data structures for loading and handling gene expression matrices. |
| **numpy** | Enables numerical operations on high-dimensional arrays. |
| **scipy** | Adds scientific computing functions for matrix algebra and optimization. |
| **matplotlib / seaborn** | Used for data visualization (volcano plots, reaction maps). |
| **scikit-learn** | Supplies basic machine learning and statistical models. |

---
### 7.5. Installing Deep Learning and Optimization Libraries
COMPASS uses machine learning and optimization modules for specific preprocessing tasks and solver integration.

```bash
conda install tensorflow=2.18 keras
conda install -c gurobi gurobi
```

### Deep Learning and Optimization Libraries

| Library | Description |
| :--- | :--- |
| **TensorFlow 2.18** | Required version for Python 3.8 compatibility; maintains GPU and solver compatibility. |
| **Keras** | High-level neural network API for TensorFlow, required for penalty diffusion and dimensionality reduction modules. |
| **Gurobi** | Core **LP/MILP solver** for metabolic flux optimization. **Requires license configuration (see Section 8).** |

---

### 7.6. Installing Supplementary Libraries
Several additional packages are used for visualization, data I/O, and algorithmic acceleration.

```bash
conda install imageio ipykernel jupyter theano cmake ninja
conda install -c conda-forge libuv=1.39
```


### Supplementary Libraries

| Library | Description |
| :--- | :--- |
| **imageio** | Reads and writes image/video data (useful for visualization output). |
| **ipykernel** | Required for launching Jupyter Notebooks inside the Conda environment. |
| **theano** | Backend engine for symbolic computation; occasionally required by Keras. |
| **cmake / ninja** | Lightweight build systems for compiling C++ extensions used by scientific libraries. |
| **libuv** | Dependency for asynchronous I/O operations (used indirectly through Gurobi and TensorFlow). |
---


### 7.7. Installing Advanced Machine Learning Packages
For large-scale data modeling, statistical acceleration, and integrated benchmarking, install the following libraries:

```Bash
conda install pycaret xgboost lightgbm catboost eli5
```

| Package | Description |
| :--- | :--- |
| **PyCaret** | Automated ML library for model comparison and tuning. |
| **XGBoost** | Gradient boosting library optimized for structured data. |
| **LightGBM** | High-performance tree-based boosting for fast classification. |
| **CatBoost** | Categorical feature boosting optimized for GPU acceleration. |
| **Eli5** | Provides model interpretation and visualization tools. |

These libraries assist in exploring data-driven interpretations of COMPASS scores and validating gene-reaction correlations through auxiliary ML analyses.

---

### 7.8. Installing Visualization and Data Exploration Tools
For high-quality, interactive data visualization and graphical postprocessing, install:

```Bash
conda install plotly bokeh opencv pattern fuel polars
```

### Visualization and Data Exploration Tools

| Package | Purpose |
| :--- | :--- |
| **Plotly / Bokeh** | Interactive charting and dashboard generation. |
| **OpenCV** | Image manipulation and visualization (useful for metabolic map overlays). |
| **Pattern** | Text-mining and web-scraping utility for biological annotations. |
| **Fuel** | Data pipeline management tool for ML model training. |
| **Polars** | High-performance DataFrame library supporting parallel computation. |

---
### 7.9. Installing the COMPASS Package
Finally, install the COMPASS framework directly from the official Wagner Lab GitHub repository.

```Bash
pip install git+https://github.com/wagnerlab-berkeley/Compass.git --upgrade
```

This command:

Clones the repository directly from GitHub.

Installs the COMPASS package and its submodules.

Registers the command-line interface (compass) globally within your Conda environment.

To verify installation:

```Bash
compass -h
```

If an error occurs due to Python version incompatibility, recheck your active environment and ensure Python ≤ 3.9 is installed.

---

### 7.10. Exporting and Reusing the Environment
To replicate the same setup on another system, export your environment to a .yml file:

```Bash
conda env export > environment.yml
```

Then, on another system, import it using:

```Bash
conda env create -f environment.yml
```

This ensures consistent software versions across collaborators and analysis platforms.

---
### 7.11. Summary of Key Version Requirements
Summary of Key Version Requirements

| Dependency | Minimum Version | Compatibility Note |
| :--- | :--- | :--- |
| **Python** | 3.8 | Required for TensorFlow $\le$ 2.18 |
| **TensorFlow** | 2.18 | Older versions ($< 2.19$) maintain backward support |
| **Gurobi** | $\ge 11.0.0$ | Requires academic or commercial license |
| **Pandas** | $\ge 1.0$ | For large-scale matrix handling |
| **NumPy** | $\ge 1.12$ | Core numerical array operations |
| **SciPy** | $\ge 1.0$ | Scientific computation backend |
| **scikit-learn** | $\ge 0.19$ | Required for statistical modeling |
| **Keras** | 3.0 | Compatible with TensorFlow 2.x API |
| **Matplotlib / Seaborn** | Latest | Visualization and plotting support |
---

### 7.12. Verifying Successful Installation
After all installations are complete, verify the core components:

```Bash
python -c "import numpy, pandas, tensorflow, gurobipy; print('All core libraries imported successfully.')"
```

If you see no errors and the message appears, your COMPASS computational environment has been configured correctly.

---
## 8. Gurobi License Configuration and Model Setup

The **Gurobi Optimizer** is the computational backbone of COMPASS, responsible for solving the **linear programming (LP)** problems that estimate metabolic reaction activity.  
Because Gurobi is a **commercial optimization solver**, it requires a valid **license file** (`gurobi.lic`) to run successfully.  

This section explains how to install, configure, and verify your Gurobi license and ensure it is properly linked with the COMPASS framework.

---

### 8.1. Purpose of Gurobi in COMPASS

COMPASS reformulates the metabolic activity inference as a **penalty minimization problem**, where each reaction’s activation cost is optimized subject to stoichiometric and thermodynamic constraints.  
This formulation requires a high-performance linear programming solver capable of handling thousands of constraints and variables efficiently.

**Gurobi** is used for:
- Solving LP problems to estimate minimal reaction penalties.
- Handling large-scale genome-scale models such as Human1 and RECON2.
- Supporting multi-core computation with deterministic numerical stability.
- Managing optimization accuracy through floating-point precision controls.

Without a valid Gurobi installation and license, COMPASS will fail at the optimization step, returning solver initialization errors.

---

### 8.2. Installing Gurobi (via Conda)

If Gurobi was not installed during environment setup, install it now using:

```bash
conda install -c gurobi gurobi
```
Once installed, you can verify the version using:

```Bash
python -c "import gurobipy; print(gurobipy.gurobi.version())"
```
This should print the version (e.g., (11, 0, 0)) if the installation was successful.

---
### 8.3. Requesting a Free Academic License
Gurobi provides free academic licenses for students, researchers, and faculty members at accredited institutions.
To obtain one:

Go to the official Gurobi website:
https://www.gurobi.com/downloads/end-user-license-agreement-academic/

Create a Gurobi account using your institutional email ID (e.g., .edu, .ac.in).

Once logged in, navigate to your account dashboard.

Click “Request Academic License”.

Follow the prompts to generate a license key specific to your device.

The system will generate a file named gurobi.lic.

---
### 8.4. Locating and Placing the License File
After you have downloaded the license file (gurobi.lic), it must be placed in a directory accessible to both the Gurobi engine and Python interpreter.

Recommended paths (depending on your OS):

Operating System	License Path Example
Windows	C:\Users\<username>\gurobi.lic
macOS/Linux	/home/<username>/gurobi.lic

You can confirm the file’s existence using a terminal or command prompt:

```Bash
dir C:\Users\<username>\gurobi.lic
```
or

```Bash
ls /home/<username>/gurobi.lic
```
---
### 8.5. Setting the License Path Manually
If Gurobi does not automatically detect your license file, you can manually define its location using:

```Bash
setx GUROBI_LICENSE "C:\Users\<username>\gurobi.lic"
```

For Linux or macOS:

```Bash
export GRB_LICENSE_FILE=/home/<username>/gurobi.lic
```

To verify the environment variable has been set:

```Bash
echo %GUROBI_LICENSE%
```

or on Linux/macOS:

```Bash
echo $GRB_LICENSE_FILE
```
---

### 8.6. Testing Gurobi Installation
After setting up the license, test your Gurobi installation using Python:

```Bash
python -c "import gurobipy; m = gurobipy.Model(); print('Gurobi license successfully validated')"
```

If no error appears and the message is printed, the license is active and the solver is ready for use within COMPASS.

In case of an error such as:

vbnet
GurobiError: License not found (error code 10009)
Double-check:

The license file path (gurobi.lic)

The environment variable configuration (GUROBI_LICENSE or GRB_LICENSE_FILE)

That your license has not expired or exceeded usage limits.

---
### 8.7. Gurobi License Renewal (Academic Licenses)
Academic licenses are valid for one year and must be renewed manually.
To renew:

Log in to your Gurobi account.

Navigate to “Academic License Management”.

Click Renew License.

Replace the existing gurobi.lic file with the new one and restart your terminal session.

---
### 8.8. Integrating Gurobi with COMPASS
Once Gurobi is installed and the license is validated, COMPASS can automatically detect and utilize the solver.

To verify this within the COMPASS environment, run:

```Bash
compass --check-license
```

If installed correctly, the output will display the solver’s version, license type, and license expiration date.

Example output:

```yaml

Gurobi Optimizer version 11.0.0
License type: Academic
License expires: 2026-02-14
License successfully linked to COMPASS.
```
---

### 8.9. Model Setup in COMPASS
COMPASS requires a valid metabolic reconstruction model to define reactions, metabolites, and stoichiometric relationships.
These models are the structural backbone of constraint-based modeling and determine how expression data is mapped to reaction penalties.

Supported models include:

### Supported Metabolic Models

| Model Name | Species | Source / Format |
| :--- | :--- | :--- |
| **Human1** | *Homo sapiens* | BiGG Models (JSON / MATLAB) |
| **RECON2 / RECON2.2** | *Homo sapiens* | RECON consortium |
| **Mouse1** | *Mus musculus* | BiGG Models / MMRR database |
Default models can be downloaded directly from the COMPASS GitHub repository or loaded locally.

Example (download manually):

Human1 Model (JSON)

RECON2.2 Model (Matlab format)

To specify the model during COMPASS execution:

```Bash
compass --data expression.tsv --model RECON2_mat --species homo_sapiens --num-processes 10
```

Explanation of Parameters:

Parameter	Description
--data	Path to the input gene expression matrix (.tsv, .mtx, or .h5ad).
--model	Selected metabolic model (e.g., Human1_mat, RECON2_mat, Mouse1_mat).
--species	Defines the biological organism (homo_sapiens or mus_musculus).
--num-processes	Number of CPU cores used for parallelization.

---

### 8.10. Troubleshooting License and Solver Issues
Error Message	Cause	Solution
GurobiError: License not found	License file missing or unreadable	Confirm gurobi.lic file path and re-export environment variable.
GurobiError: Feature not licensed	License type restricted	Verify you are using an Academic or Full license.
Segmentation fault	Incorrect Python-Gurobi binding	Reinstall Gurobi within the same Python environment.
Model failed to optimize	Invalid or incomplete metabolic model	Check that the selected model file is readable and in supported format.

Once Gurobi is correctly configured and the model is validated, COMPASS is ready for full execution using your transcriptomic dataset.

---
## 9. Input Data Preparation, Execution Commands, and Output Interpretation

The COMPASS workflow operates on gene expression data obtained from **bulk RNA-seq** or **single-cell RNA-seq (scRNA-seq)** experiments.  
To ensure reproducibility and accuracy, input files must be correctly formatted and preprocessed before being passed to the COMPASS command-line interface.

---

### 9.1. Accepted Input Data Formats

COMPASS accepts **three standard data formats** depending on the source and preprocessing pipeline used. Each format represents the same underlying data: gene identifiers (rows) and samples or cells (columns).

| Format | Files Required | Description | Typical Source |
|---------|----------------|-------------|----------------|
| `.tsv` | `expression.tsv` | Tab-separated matrix containing raw or normalized gene expression values. | Bulk RNA-seq / scRNA-seq post-quantification. |
| `.mtx` | `expression.mtx`, `genes.tsv`, `sample_names.tsv` | Sparse matrix format following the 10x Genomics standard. | 10x Genomics CellRanger output. |
| `.h5ad` | `anndata_object.h5ad` | Hierarchical data format for AnnData objects (Scanpy-compatible). | Python-based single-cell pipelines. |

**Example structure for each format:**

1. **TSV Format**

| Gene | Sample1 | Sample2 | Sample3 |
| :--- | :--- | :--- | :--- |
| **GAPDH** | 10.2 | 11.5 | 9.8 |
| **LDHA** | 5.3 | 6.1 | 7.0 |
| **HK2** | 4.5 | 4.8 | 3.9 |


2. **MTX Format**

matrix.mtx → Expression matrix in sparse Matrix Market format
genes.tsv → Gene identifiers (one per line)
sample_names.tsv → Sample or cell barcodes (one per line)


Ensure filenames are renamed to match COMPASS expectations:
matrix.mtx → expression.mtx
barcodes.tsv → sample_names.tsv
H5AD Format

If your data is in AnnData format (.h5ad), it can be used directly without conversion:

```bash
compass --data your_dataset.h5ad --num-processes 8 --species homo_sapiens
```
---
### 9.2. Preprocessing Recommendations
Before running COMPASS, it is recommended to perform the following preprocessing steps using Scanpy, Seurat, or any equivalent tool:

Normalization – Apply library size normalization and log-transformation.

Filtering – Remove low-quality cells and lowly expressed genes.

Gene Naming Consistency – Use consistent identifiers (e.g., Ensembl IDs or official gene symbols).

Batch Correction – If applicable, apply integration or batch correction across datasets.

Export – Save the processed matrix in .tsv, .mtx, or .h5ad format for COMPASS input.

---
### 9.3. Basic Command-Line Execution
The standard COMPASS execution command is as follows:

```Bash
compass --data expression.tsv --model RECON2_mat --species homo_sapiens --num-processes 10
```

Parameter Definitions:
Argument	Description
--data	Path to input gene expression file (.tsv, .mtx, or .h5ad).
--model	Selected genome-scale metabolic model (Human1_mat, RECON2_mat, Mouse1_mat).
--species	Defines the biological species (homo_sapiens or mus_musculus).
--num-processes	Number of processor cores to be used for parallel computation.

Example (Mouse dataset):

```Bash
compass --data expression.tsv --model Mouse1_mat --species mus_musculus --num-processes 12
```
---

### 9.4. Turbo Mode for Accelerated Computation
For large-scale datasets (e.g., >10,000 cells), COMPASS provides a Turbo Mode that implements vectorized penalty optimization and memory-efficient solver calls.
Enable Turbo Mode with:

```Bash
compass --data expression.tsv --turbo 1.0 --num-processes 16
```

Turbo Mode Notes:

The --turbo argument accepts a float value (e.g., 1.0, 0.5) controlling the trade-off between accuracy and speed.

When --turbo is enabled, COMPASS uses approximation-based penalty aggregation, which significantly reduces runtime but may slightly smooth fine-grained reaction variability.

Turbo Mode is ideal for exploratory or large-scale screening analyses.

---
### 9.5. Species Specification and Ortholog Mapping
When working with mouse data, COMPASS automatically performs ortholog mapping between Mus musculus and Homo sapiens models.
This ensures reaction penalties are computed on homologous metabolic reactions even if gene IDs differ.

```Bash
compass --data expression.tsv --species mus_musculus --model Mouse1_mat --num-processes 8
```

During runtime, a log file will indicate successful mapping:

```pgsql
[INFO] Ortholog mapping applied: 96.8% of genes mapped to Human1 model.
```
---
### 9.6. Output File Structure
Upon successful execution, COMPASS generates several structured output files. All outputs are written to the working directory unless specified otherwise using --output-dir.

File Name	Description
reactions.tsv	Smoothed penalty-based metabolic activity scores for all reactions across cells.
penalties.txt.gz	Raw unsmoothed reaction penalties before λ-diffusion.
model.json / model.json.gz	The metabolic model representation used in the analysis.
reaction_consistencies.csv	Derived Compass activity scores (-log transformed).
wilcox_results.csv	Statistical results from post-analysis (e.g., Wilcoxon rank-sum test, effect sizes).
_tmp/	Temporary computation files generated during solver execution.

---
### 9.7. Mathematical Computation of Compass Scores

Each reaction–cell pair receives a Compass penalty value representing the minimal cost required for that reaction to carry flux. To obtain interpretable activity measures, these penalties are transformed using the following formula:

$$\text{Activity}_{r,c} = - \log(\text{Penalty}_{r,c} + 1)$$

Where:
* $\text{Activity}_{r,c}$: Compass activity score for reaction $r$ in cell $c$.
* $\text{Penalty}_{r,c}$: Gurobi-optimized penalty value.

The log-transformation ensures diminishing sensitivity to large penalty values and normalizes the scale for comparison.

---

### 9.8. $\lambda$-Diffusion Smoothing

To account for **biological noise** and stochastic expression variability in single-cell data, COMPASS applies a penalty diffusion mechanism controlled by the parameter $\lambda$.

This process involves propagating penalty values across a cell similarity graph defined by $k$-nearest neighbors (KNN) based on transcriptomic distance. The **smoothed penalties** are recalculated as:

$$\text{Penalty}_{r,c}' = (1 - \lambda) \cdot \text{Penalty}_{r,c} + \lambda \cdot \sum_{i \in N(c)} w_{ci} \cdot \text{Penalty}_{r,i}$$

Where:
* $N(c)$ = neighborhood set of cell $c$.
* $w_{ci}$ = similarity weight between cells $c$ and $i$.
* $\lambda \in [0, 1]$ controls the diffusion strength.

Typical $\lambda$ values range from **0.1 to 0.3**, balancing signal preservation with noise reduction.

---
### 9.9. Example Execution Workflow
Below is a complete example pipeline demonstrating the end-to-end execution of COMPASS on a human dataset using the RECON2 model.

```Bash
# Step 1: Activate environment
conda activate compass_env

# Step 2: Run Compass on normalized gene expression data
compass --data expression.tsv \
        --model RECON2_mat \
        --species homo_sapiens \
        --num-processes 10 \
        --turbo 1.0

# Step 3: Verify successful run
ls -l | grep reactions.tsv
```

If successful, you should see files such as:

```pgsql
penalties.txt.gz
reactions.tsv
model.json.gz
wilcox_results.csv
```
---
### 9.10. Output Interpretation
The reactions.tsv file represents the key Compass output.
Each row corresponds to a metabolic reaction, and each column corresponds to a sample or cell.

Example (excerpt):

| Reaction | Cell_1 | Cell_2 | Cell_3 |
| :--- | :--- | :--- | :--- |
| **RXN_001** | 2.54 | 2.71 | 1.93 |
| **RXN_002** | 0.85 | 0.76 | 0.90 |
| **RXN_003** | 1.12 | 1.30 | 1.04 |

High Compass activity values correspond to reactions predicted to be biochemically active (low penalty), whereas low activity values indicate repressed or inactive metabolic pathways.

---
### 9.11. Post-Processing Overview
Downstream statistical and visualization steps are typically executed in Jupyter notebooks provided with the COMPASS package. These include:

Differential Compass score analysis between conditions (e.g., pathogenic vs. non-pathogenic Th17 cells).

Wilcoxon rank-sum tests for per-reaction significance.

Volcano and strip plots for pathway visualization.

Subsystem-level summarization (glycolysis, TCA cycle, lipid metabolism, etc.).

Example transformation:

```python
import numpy as np
activity_score = -np.log(penalty_score + 1)
```

At this stage, Compass has generated a complete set of reaction-level metabolic activity profiles ready for statistical evaluation and visualization.

---
## 10. Statistical Postprocessing, Visualization, and Biological Interpretation

Once COMPASS produces the reaction-level activity scores, downstream **statistical and visualization analyses** are used to identify biologically meaningful differences in metabolic pathway activity between experimental conditions.  
This postprocessing step transforms raw Compass outputs into interpretable insights, such as **upregulated pathways, altered flux patterns**, and **metabolic reprogramming events**.

---

### 10.1. Overview of Statistical Analysis Workflow

The postprocessing workflow involves several key analytical steps:

1. **Loading Compass Output Files**
   - Import `reactions.tsv` and metadata (e.g., sample annotations).
   - Align Compass scores with experimental group labels.

2. **Group-Based Statistical Testing**
   - Apply nonparametric tests (e.g., Wilcoxon rank-sum test) to identify significantly different reactions between conditions.

3. **Effect Size Estimation**
   - Compute **Cohen’s d** to assess the magnitude of metabolic activity changes.

4. **Multiple Testing Correction**
   - Adjust p-values for multiple comparisons using methods such as **Benjamini–Hochberg (FDR)**.

5. **Subsystem-Level Aggregation**
   - Group reactions into pathways or subsystems for higher-level biological interpretation.

6. **Visualization and Interpretation**
   - Generate volcano plots, violin plots, and subsystem heatmaps to visualize differential metabolism.

---

### 10.2. Required Python Libraries for Analysis

All analyses can be performed in Jupyter notebooks using the following libraries (already installed in previous steps):

```Bash
conda install pandas numpy scipy seaborn matplotlib plotly statsmodels
```

Imports Example:

```python
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import ranksums
from statsmodels.stats.multitest import multipletests
```
---
### 10.3. Loading Compass Output and Metadata

```python
# Load Compass reaction activity scores
compass_data = pd.read_csv("reactions.tsv", sep="\t", index_col=0)

# Load sample annotations
metadata = pd.read_csv("metadata.csv")  # Contains sample names and conditions

# Align Compass data with sample metadata
compass_data = compass_data.loc[:, metadata["SampleID"]]
metadata.head()
```

Expected metadata format:

| SampleID | Condition |
| :--- | :--- |
| **Cell_1** | Pathogenic |
| **Cell_2** | NonPathogenic |
| **Cell_3** | Pathogenic |

---
### 10.4. Wilcoxon Rank-Sum Test
The Wilcoxon rank-sum test is used to compare Compass activity scores between two conditions for each reaction.

```python
p_values = []
for rxn in compass_data.index:
    group1 = compass_data.loc[rxn, metadata["Condition"] == "Pathogenic"]
    group2 = compass_data.loc[rxn, metadata["Condition"] == "NonPathogenic"]
    stat, p = ranksums(group1, group2)
    p_values.append(p)

# Adjust for multiple testing
_, padj, _, _ = multipletests(p_values, method='fdr_bh')
```

The results can then be stored as:

```python
results = pd.DataFrame({
    "Reaction": compass_data.index,
    "p_value": p_values,
    "padj": padj
})
results.to_csv("wilcox_results.csv", index=False)
```
---
### 10.5. Effect Size Calculation (Cohen’s $d$)

Effect size quantifies the magnitude of metabolic difference between groups. It complements statistical significance by revealing biologically large but possibly subtle effects.

$$d = \frac{\bar{x}_1 - \bar{x}_2}{s_p}$$

Where:
* $\bar{x}_1, \bar{x}_2$ = mean Compass scores for groups 1 and 2.
* $s_p$ = pooled standard deviation.

* **Positive $d$**: Reaction more active in the first group (e.g., pathogenic cells).
**Negative $d$**: Reaction more active in the second group.
Implementation:

```python
def cohens_d(group1, group2):
    n1, n2 = len(group1), len(group2)
    s1, s2 = np.var(group1, ddof=1), np.var(group2, ddof=1)
    s_pooled = np.sqrt(((n1 - 1)*s1 + (n2 - 1)*s2) / (n1 + n2 - 2))
    return (np.mean(group1) - np.mean(group2)) / s_pooled

results["cohen_d"] = [
    cohens_d(compass_data.loc[rxn, metadata["Condition"] == "Pathogenic"],
             compass_data.loc[rxn, metadata["Condition"] == "NonPathogenic"])
    for rxn in compass_data.index
]
```
---
### 10.6. Volcano Plot Visualization
A volcano plot displays significance (−log10 p-value) versus effect size (Cohen’s d), helping identify key metabolic reactions that are significantly different between conditions.

```python
results["log_p"] = -np.log10(results["padj"])
sns.scatterplot(data=results, x="cohen_d", y="log_p", alpha=0.7)
plt.axhline(-np.log10(0.05), color="red", linestyle="--", label="FDR=0.05")
plt.xlabel("Effect Size (Cohen's d)")
plt.ylabel("-log10 Adjusted p-value")
plt.title("Volcano Plot: Differential Metabolic Activity")
plt.legend()
plt.show()
```

Interpretation:

Positive d → Reaction more active in the first group (e.g., pathogenic cells).

Negative d → Reaction more active in the second group.

Higher log p → Stronger statistical significance.

---
### 10.7. Subsystem-Level Aggregation
Reactions can be grouped into metabolic subsystems (e.g., glycolysis, TCA cycle, fatty acid synthesis) for pathway-level summarization.

If the model metadata includes subsystem information:

```python
model_meta = pd.read_csv("model_metadata.csv")  # Contains columns: Reaction, Subsystem
merged = results.merge(model_meta, on="Reaction")
# Compute subsystem mean effect sizes
subsystem_summary = merged.groupby("Subsystem")[["cohen_d", "padj"]].mean()
subsystem_summary.sort_values("cohen_d", ascending=False).head(10)
```

This highlights which metabolic subsystems are globally upregulated or repressed in one condition relative to the other.

---
### 10.8. Visualizing Subsystem Differences
Generate a bar plot or heatmap of differential subsystem activity:

```python
plt.figure(figsize=(10, 6))
sns.barplot(x=subsystem_summary.index, y=subsystem_summary["cohen_d"])
plt.xticks(rotation=90)
plt.ylabel("Mean Effect Size (Cohen's d)")
plt.title("Subsystem-Level Metabolic Reprogramming")
plt.tight_layout()
plt.show()
```

Interpretation:

Bars above zero → Subsystems upregulated in condition 1.

Bars below zero → Subsystems upregulated in condition 2.

---
### 10.9. Heatmap Visualization of Reaction Consistency
Reaction consistency can be visualized to assess reproducibility of Compass scores across cells or conditions.

```python
import seaborn as sns
sns.clustermap(compass_data.corr(), cmap="vlag", figsize=(12, 10))
plt.title("Reaction Activity Correlation Across Samples")
plt.show()
```

Such visualizations help identify clusters of co-active metabolic reactions and reveal pathway modules with coordinated expression.

---
### 10.10. Example Integrated Analysis Workflow
Below is a typical postprocessing pipeline integrating all the above steps:

```python
# Step 1: Load outputs
reactions = pd.read_csv("reactions.tsv", sep="\t", index_col=0)
metadata = pd.read_csv("metadata.csv")

# Step 2: Differential analysis
from scipy.stats import ranksums
p_values, cohen_d_vals = [], []

for rxn in reactions.index:
    g1 = reactions.loc[rxn, metadata["Condition"] == "Pathogenic"]
    g2 = reactions.loc[rxn, metadata["Condition"] == "NonPathogenic"]
    stat, p = ranksums(g1, g2)
    p_values.append(p)
    cohen_d_vals.append(cohens_d(g1, g2))

# Step 3: Multiple testing correction
_, padj, _, _ = multipletests(p_values, method="fdr_bh")

# Step 4: Compile results
df = pd.DataFrame({
    "Reaction": reactions.index,
    "p_value": p_values,
    "padj": padj,
    "Cohen_d": cohen_d_vals
})
df.to_csv("compass_statistical_results.csv", index=False)

# Step 5: Volcano plot visualization
df["log_p"] = -np.log10(df["padj"])
sns.scatterplot(data=df, x="Cohen_d", y="log_p", hue=df["padj"] < 0.05)
plt.title("Compass Differential Reaction Analysis")
plt.show()
```
---
### 10.11. Biological Interpretation
Once significant reactions and pathways are identified, interpret them in the biological context of your study.
For example:

Elevated glycolytic flux may suggest Warburg-like metabolic activation in immune or cancer cells.

Reduced oxidative phosphorylation activity could indicate hypoxia-induced adaptation.

Alterations in lipid metabolism may reflect membrane remodeling or signaling pathway changes.

To extend the analysis:

Perform Gene Ontology (GO) or KEGG pathway enrichment using significantly active reactions.

Integrate Compass outputs with Flux Balance Analysis (FBA) to simulate flux distributions.

Overlay Compass results onto metabolic maps for intuitive visualization.

---
### 10.12. Summary of Postprocessing Outputs
Summary of Postprocessing Outputs

| File Name | Description |
| :--- | :--- |
| **`wilcox_results.csv`** | Contains $p$-values and adjusted $p$-values for each reaction. |
| **`reaction_consistencies.csv`** | COMPASS-derived activity consistency metrics. |
| **`compass_statistical_results.csv`** | Consolidated statistical analysis (Wilcoxon Rank-Sum Test + Cohen’s $d$). |
| **`volcano_plot.png`** | Visualization of significant metabolic reactions. |
| **`subsystem_summary.csv`** | Aggregated mean effect sizes for each metabolic subsystem. |
---

### 10.13. Example Output Summary
Example Output Summary: Differential Subsystem Activity

| Subsystem | Mean Cohen’s $d$ | Interpretation |
| :--- | :--- | :--- |
| **Glycolysis / Gluconeogenesis** | +1.12 | Upregulated in pathogenic cells |
| **Fatty Acid Metabolism** | -0.84 | Reduced activity in pathogenic cells |
| **TCA Cycle** | +0.95 | Enhanced oxidative metabolism |
| **Purine Metabolism** | +0.60 | Increased nucleotide synthesis |
| **Pentose Phosphate Pathway** | +0.71 | Elevated NADPH production |
---
### 11. References
Wagner DE et al. (2021). Single-cell metabolic modeling identifies pathways underlying Th17 cell pathogenicity. Cell, 184(16), 4168–4186.

Brunk E et al. (2018). Recon3D enables a three-dimensional view of gene variation in human metabolism. Nature Biotechnology, 36(3), 272–281.

King ZA et al. (2016). BiGG Models: A platform for integrating, standardizing, and sharing genome-scale models. Nucleic Acids Research, 44(D1), D515–D522.

