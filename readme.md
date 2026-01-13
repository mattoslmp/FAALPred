# FAALPred: Fatty Acyl-AMP Ligases (FAAL) Prediction Tool

FAALPred is a reproducible, open-source workflow to predict the fatty acyl chain-length specificity of Fatty Acyl-AMP Ligases (FAALs), in the range from C4 to C18.  
It combines **protein-domain extraction**, **MAFFT-based alignment**, **Word2Vec protein embeddings**, **oversampling strategies**, and **Random Forest models with probability calibration**, all wrapped in an interactive **Streamlit** interface.

FAALPred is available as:

- 🌐 **Public web server**: https://faalpred.ciimar.up.pt/  
- 💻 **Local Streamlit app** (this repository)

![FAALPred overview](image.png)

---

## Table of Contents

- [Overview](#overview)
- [Method Summary](#method-summary)
- [Code Structure](#code-structure)
- [Requirements](#requirements)
- [Installation](#installation)
  - [1. Install Git and Conda](#1-install-git-and-conda)
  - [2. Clone this repository](#2-clone-this-repository)
  - [3. Create and activate the `faalpred` environment](#3-create-and-activate-the-faalpred-environment)
- [Running FAALPred](#running-faalpred)
  - [Default training data vs custom training](#default-training-data-vs-custom-training)
  - [Input formats](#input-formats)
  - [Main outputs](#main-outputs)
- [Running FAALPred with Docker](#running-faalpred-with-docker)
- [Auxiliary Tool: Automatic FAAL Domain Extraction](#auxiliary-tool-automatic-faal-domain-extraction)
- [Reproducibility, Logging and Directories](#reproducibility-logging-and-directories)
- [Supplementary Methodology and Internal API](#supplementary-methodology-and-internal-api)
- [Citation](#citation)
- [Contact and Acknowledgements](#contact-and-acknowledgements)

---

## Overview

Fatty Acyl-AMP Ligases (FAALs) activate fatty acids of different chain lengths for incorporation into natural product biosynthetic pathways.  
FAALPred predicts **chain-length specificity profiles** (C4–C18) of FAAL domains based only on their amino acid sequences.

The core model was developed and tested on FAAL domains identified with conserved domain annotations (e.g. **CDD cd05931: FAAL**), and is intended to be used **on FAAL domains**, not full-length multidomain proteins.  

The graphical interface provides:

- Guided training and prediction workflow (with default training data or user-supplied data),
- FAAL domain extraction helper using InterProScan + bedtools,
- Multiple quality-control plots (UMAP, oversampling diagnostics, F1 per class, calibration, etc.),
- Tabular output with a **single “best specificity block”** plus a **continuous confidence score (0–1)**.

---

## Method Summary

At a high level, FAALPred performs:

1. **FAAL domain extraction (optional but recommended)**  
   - Via an auxiliary pipeline that uses the InterProScan REST API and `bedtools getfasta` to extract FAAL domains from full-length proteins.

2. **Sequence alignment and preprocessing**
   - If sequences are not aligned, FAALPred calls **MAFFT** (localpair / maxiterate 1000) to generate a multiple sequence alignment.
   - The model operates on the aligned sequences.

3. **Protein sequence embedding**
   - Sequences are tokenized into overlapping k-mers (default `k=3`, user-configurable `step_size`).
   - A **Word2Vec** model (`gensim`) is trained (or reused if already present) on k-mer “sentences”.
   - Per-protein embeddings are built by:
     - generating k-mer embeddings,
     - standardizing the number of k-mers per sequence (`min_kmers`),
     - aggregating k-mer embeddings (currently using the **mean** aggregation by default).

4. **Feature scaling**
   - Embeddings are standardized using `StandardScaler`.
   - The scaler is saved (`scaler_associated.pkl`) for later reuse in prediction.

5. **Class balancing and oversampling**
   - Class distributions are often highly imbalanced.
   - FAALPred applies:
     - **RandomOverSampler** (with class-specific minimum counts ≥ CV folds + 1),
     - followed by **SMOTE** oversampling.

6. **Model training and calibration**
   - A **Random Forest** classifier is trained on the oversampled embeddings.
   - Hyperparameters are tuned via **GridSearchCV** with stratified cross-validation.
   - The best model is then calibrated with **isotonic regression** (`CalibratedClassifierCV`).
   - Evaluation includes:
     - F1 scores (global and per class),
     - ROC AUC (binary or multi-class OVO),
     - Precision–Recall AUC,

7. **Prediction on new sequences**
   - New sequences are aligned (if required), embedded, scaled with the trained scaler, and fed into the calibrated RF model.
   - For each query:
     - full ranked probabilities for all classes are computed,
     - FAALPred maps single-chain-length labels into **chain-length “blocks”** (e.g. `C4–C6–C8`, `C10–C12–C14`, …),
     - a **single main block** is returned together with a **continuous confidence score in [0,1]**, based on normalized probabilities across top blocks.

---

## Code Structure

The implementation is contained in a single main Streamlit script (`faalpred.py`) plus supporting data and image folders.

### Main Classes

- **`Support`**
  - Random Forest training with oversampling and cross-validation
  - Grid search for hyperparameter tuning
  - Probability calibration with `CalibratedClassifierCV`
  - Learning-curve plotting

- **`ProteinEmbeddingGenerator`**
  - Handles sequence alignment (MAFFT if needed),
  - Generates k-mer based Word2Vec embeddings,
  - Manages the `min_kmers` logic and aggregation strategy,
  - Produces standardized embedding matrices and associated labels (`associated_variable`).

### Key Standalone Functions

- `are_sequences_aligned`  
- `create_unique_model_directory`  
- `realign_sequences_with_mafft`  
- `plot_roc_curve_global`  
- `plot_precision_recall_curve_global`  
- `get_class_rankings_global`  
- `adjust_predictions_global`  
- `format_and_sum_probabilities`  
- `plot_predictions_scatterplot_custom`  
- `plot_prediction_confidence_bar`

### Optional functions
- `visualize_latent_space_with_similarity`  
- `plot_umap_3d_combined`  
- `plot_oversampling_quality`  
- `plot_confusion_and_calibration`  


### Auxiliary Tools (FAAL domain extraction)

- InterProScan REST API workflow:
  - `submit_job`, `poll_status`, `retrieve_result`
- FAAL-specific domain extraction:
  - `extract_faal`, `create_bed`, `run_bedtools`, `process_faal_domain`

A more detailed, method-oriented description is provided in the file  
**`Supplementary_methodology.dox`**, which serves as an extended “Supplementary Methods” reference for the accompanying manuscript.

---

## Requirements

- **Python**: 3.9 (Python ≥ 3.8 should work, but the environment is built for 3.9).
- **Conda** (Anaconda or Miniconda)
- **Operating system**: Linux, macOS, or Windows (64-bit recommended)
- **MAFFT** (for sequence alignment; must be available in `PATH`)
- **bedtools** (for the optional FAAL domain extraction tool; see below)

All Python dependencies are specified in the `faalpred_env.yml` file.

---

## Installation

### 1. Install Git and Conda

If not already installed:

- **Git**: https://git-scm.com/downloads  
- **Anaconda**: https://www.anaconda.com/download  
  or  
- **Miniconda** (lightweight): https://docs.conda.io/en/latest/miniconda.html  

After installation, open a **terminal** (or **Anaconda Prompt** on Windows).

### 2. Clone this repository

```bash
git clone https://github.com/CNP-CIIMAR/FAALPred.git
cd FAALPred
```

After this step, the file **`faalpred_env.yml`** is available in the current directory.

### 3. Create and activate the `faalpred` environment

Create the Conda environment from the YAML file:

```bash
conda env create -f faalpred_environment.yml
conda activate faalpred
```

Every time you want to use FAALPred locally:

```bash
conda activate faalpred
```

To leave the environment:

```bash
conda deactivate
```

---

## Running FAALPred

With the environment activated (`conda activate faalpred`) and standing in the repository root:

```bash
streamlit run faalpred.py
```

Streamlit will print a local URL such as:

```text
Local URL: http://localhost:8501
```

Open this URL in your browser to access the FAALPred interface.

If running on a remote server, you can configure Streamlit via `~/.streamlit/config.toml` (for example):

```toml
[server]
headless = true
enableCORS = false
enableXsrfProtection = false
address = "0.0.0.0"
port = 8501
```

---

### Default training data vs custom training

On the **left sidebar**, the app provides:

- **Use Default Training Data** (checkbox)  
  - If enabled, FAALPred uses the internal training set:
    - `data/train.fasta`
    - `data/train_table.tsv`
- If you uncheck this option, you can upload your own:
  - **Training FASTA** (aligned or unaligned protein sequences),
  - **Training table (TSV)** with at least the columns:
    - `Protein.accession`
    - `Target variable`
    - `Associated variable`  
    (`Associated variable` is the chain-length specificity label used for training.)

For **prediction**, you must upload a **FASTA file** containing the **query FAAL domains**.

Other tunable parameters in the sidebar:

- `K-mer Size` (default: 3)  
- `Step Size` (default: 1)  
- `Aggregation Method` (currently `mean`)  
- Optional Word2Vec settings: window size, number of workers, number of epochs.

---

### Input formats

#### 1. Training FASTA

- Multi-FASTA protein file with one sequence per FAAL domain.
- Sequences may be aligned or unaligned:
  - If unaligned, FAALPred automatically runs MAFFT and creates an `_aligned.fasta` file.

#### 2. Training table (TSV)

A tab-delimited file where each row corresponds to a FAAL domain in the training FASTA.  
Required columns:

- `Protein.accession` – identifier matching sequence headers (up to first whitespace),
- `Target variable` – optional or unused in some runs,
- `Associated variable` – chain-length specificity label used as the **class**.

#### 3. Prediction FASTA

- Multi-FASTA file with the **FAAL domains to be predicted**.
- Recommended: use the **Auxiliary FAAL domain tool** to extract domains first (see below).

---

### Main outputs

By default, outputs are created under:

```text
results/models_<aggregation_method>/
```

For example, with the default aggregation method `mean`:

```text
results/models_mean/
```

Key files include:

- **`predictions.tsv`**  
  Tab-separated file with, for each protein:
  - Protein ID,
  - Predicted specificity label (`Associated_Prediction`),
  - Full probability ranking across classes.

- **`results.xlsx`**  
  Excel version of the formatted results table:

  - **Query Name**
  - **SS Prediction Specificity** (e.g. `C10-C12-C14`)
  - **Prediction Confidence ( range: 0 - 1 )**  
    (continuous score derived from normalized probabilities of the top specificity blocks)

- **`formatted_results.txt`**  
  The same table as ASCII/Markdown-style grid (via `tabulate`).

- **`scatterplot_predictions.png`**  
  Publication-ready scatter plot:  
  - Y-axis: protein IDs  
  - X-axis: chain-length positions (C4 to C18)  
  - Lines mark the predicted block for each protein.


- **Model & scaler files**
  - `word2vec_model.bin`
  - `scaler_associated.pkl`
  - `model_best_associated.pkl`
  - `calibrated_model_associated.pkl`

- **Performance plots**
  - `learning_curve_associated.png`
  - `roc_curve_associated.png`

From the Streamlit interface, you can also:

- Download **CSV** and **Excel** tables directly,
- Download a **`results.zip`** archive containing all outputs in the chosen run directory.

---

## 2. Running FAALPred with Docker (Docker Hub image: `mattoslmp/faalpred`)

FAALPred is available as a ready-to-use Docker image on Docker Hub: **`mattoslmp/faalpred`**.  
This allows you to run the full Streamlit app (with all dependencies) without manually installing the Conda environment.

The recommended workflow is:

1. Clone the FAALPred project from GitHub  
2. Ensure local permissions for output directories (optional but recommended)  
3. Install Docker  
4. Pull and run the FAALPred Docker image **`mattoslmp/faalpred`** from Docker Hub

---

### 2.1. Clone the FAALPred repository and set permissions

First, clone this repository and move into the project folder:

```bash
git clone https://github.com/mattoslmp/FAALPred.git
cd FAALPred
```

If you plan to **mount local folders** for results/logs (recommended for persistence), create them and set permissions:

```bash
mkdir -p results logs
chmod 775 results logs
```

> 💡 If you are using WSL2, it is recommended to clone the repository inside the Linux filesystem  
> (e.g. `/home/<user>/FAALPred`) instead of under `/mnt/c` for better performance and fewer I/O issues.

---

### 2.2. Install Docker on Ubuntu or WSL2

On a fresh Ubuntu (or WSL2 Ubuntu) system, you can install Docker using the official convenience script:

```bash
sudo apt-get update
sudo apt-get install -y ca-certificates curl gnupg lsb-release

curl -fsSL https://get.docker.com -o get-docker.sh
sudo sh get-docker.sh
```

After installation, add your user to the `docker` group so you can run Docker without `sudo`:

```bash
sudo usermod -aG docker $USER
```

**Log out and log back in** (or fully restart your session) so group changes take effect.  
Then test:

```bash
docker ps
```

You should see an empty list of containers (no “permission denied” error).

> If you prefer the long-form installation following Ubuntu’s documentation, see:  
> https://docs.docker.com/engine/install/ubuntu/

---

### 2.3. Pull the FAALPred Docker image from Docker Hub

Once Docker is installed and working, you can **download (pull)** the FAALPred image from Docker Hub:

```bash
docker pull mattoslmp/faalpred:latest
```

This will download the FAALPred image `mattoslmp/faalpred` to your local Docker installation.

You can verify that the image is available with:

```bash
docker images | grep faalpred
```

You should see a line similar to:

```text
mattoslmp/faalpred   latest   <IMAGE_ID>   <CREATED>   <SIZE>
```

---

### 2.4. Run the FAALPred container (simple mode)

To start the Streamlit app in a container and expose it on port **8501** of your machine:

```bash
docker run --rm -p 8501:8501 mattoslmp/faalpred:latest
```

- `--rm` removes the container when it stops.
- `-p 8501:8501` maps container port 8501 (Streamlit) to host port 8501.
- `mattoslmp/faalpred:latest` is the image name on Docker Hub.

After the container starts, open in your browser:

```text
http://localhost:8501
```

If you are running Docker on a remote server, replace `localhost` with the server IP or hostname, and ensure that port 8501 is open in the firewall.

To stop FAALPred in this mode, go to the terminal where `docker run` is running and press `Ctrl + C`.  
Because we used `--rm`, the container will be automatically removed.

---

### 2.5. Run with persistent outputs and local data (recommended)

By default, any result files generated inside the container (e.g. under `results/`) will be **lost when the container stops**.  
To keep outputs and optionally use your local `data/` and `validation/` directories, you can mount them as volumes.

From inside the cloned repository (where `results/` and `logs/` exist):

```bash
docker run --rm -it \
  -p 8501:8501 \
  -v "$(pwd)/results:/app/results" \
  -v "$(pwd)/logs:/app/logs" \
  mattoslmp/faalpred:latest
```

- `results/` and `logs/` on the host will store all outputs and logs created by FAALPred.
- The container still runs the embedded app and environment.

If you also want to mount your own training/validation data:

```bash
docker run --rm -it \
  -p 8501:8501 \
  -v "$(pwd)/data:/app/data" \
  -v "$(pwd)/validation:/app/validation" \
  -v "$(pwd)/results:/app/results" \
  -v "$(pwd)/logs:/app/logs" \
  mattoslmp/faalpred:latest
```

In this setup:

- The code inside the container continues to refer to `data/`, `validation/`, `results/`, and `logs/`,
- But the actual files live in your cloned GitHub project on the host, making it easy to inspect, back up, or version-control them.

---

## 3. (Optional) Building the FAALPred Docker image from source

If you prefer to build the Docker image yourself (instead of pulling from Docker Hub), make sure Docker is installed and then run:

```bash
git clone https://github.com/mattoslmp/FAALPred.git
cd FAALPred

docker build -t mattoslmp/faalpred:latest .
```

After the build completes, run as before:

```bash
docker run --rm -p 8501:8501 mattoslmp/faalpred:latest
```

---

## 4. Development notes

- The main entry point of the app is `faal_pred.py`.
- Environment dependencies are defined in `faalpred_env.yml` (Python 3.9 + scientific stack).
- Most helper functions and reusable code are under the `utilities/` directory.

If you plan to modify the app or extend the models:

1. Create a feature branch in your fork.  
2. Adjust or add utilities under `utilities/`.  
3. Update `faalpred_env.yml` if you add new dependencies.  
4. Rebuild the Docker image (if needed) with:

   ```bash
   docker build -t mattoslmp/faalpred:latest .
   ```

------

## Auxiliary Tool: Automatic FAAL Domain Extraction

FAALPred includes an **Auxiliary Tool** in the sidebar: **“Get your FAAL domain”**.

This workflow:

1. Accepts a **full-length FASTA** file uploaded by the user.
2. Submits the file to the **EBI InterProScan REST API** (`iprscan5`).
3. Monitors the job until completion.
4. Parses the TSV output to find hits containing `"FAAL"`.
5. Builds a **BED** file with FAAL coordinates.
6. Uses `bedtools getfasta` to extract the corresponding FAAL regions.
7. Returns a processed FASTA file containing only the FAAL domains.

The final FASTA can then be used as:

- **Training FASTA** (if you are building your own model), or
- **Prediction FASTA** (for direct specificity prediction).

> **Important:**  
> The model in this repository was trained **only on FAAL domains**, identified by conserved-domain signatures (e.g. **CDD cd05931: FAAL**).  
> For best performance and interpretability, **always extract the FAAL domain** before running predictions.

### Additional requirement for this tool

- `bedtools` must be installed and available in your `PATH`. For example:

```bash
conda install -c bioconda bedtools
```

The InterProScan calls require Internet access and comply with EBI’s usage policies.

---

## Reproducibility, Logging and Directories

- Outputs are organized under `results/` by **aggregation method**:
  - e.g. `results/models_mean/`
- Oversampling statistics are recorded in:
  - `oversampling_counts.txt`
  - `training_sample_counts_after_oversampling.txt`
- A simple visit counter for the web interface uses:
  - `logs/visit_count.txt` (created automatically if missing)
- Progress bars and status messages within Streamlit track the main pipeline steps.

Random seeds (`SEED = 42`) are fixed in the code for NumPy and Python’s random module to improve reproducibility.

---

## Supplementary Methodology and Internal API

A detailed description of the FAALPred workflow, including:

- k-mer generation strategy,
- embedding dimensionality and Word2Vec training hyperparameters,
- oversampling and cross-validation setup,
- ROC / PR AUC evaluation and calibration,
- Others.

is provided in the file:

- **`Supplementary_methodology.dox`**

which can be cited as Supplementary Methods in the Protein Science article.

---

## Citation

If you use FAALPred in your research, please cite:

**Diversity of FAAL enzymes and prediction of their substrate specificity using FAALPred**  
Anne Liong¹˒²,†, Leandro de Mattos Pereira¹,† and Pedro N. Leão¹,*  
¹ Interdisciplinary Centre of Marine and Environmental Research (CIIMAR/CIMAR), University of Porto, Matosinhos, 4450-208, Portugal  
² ICBAS – School of Medicine and Biomedical Sciences, University of Porto, Porto, 4050-313, Portugal  

**To whom correspondence should be addressed.* Email: [pleao@ciimar.up.pt](mailto:pleao@ciimar.up.pt)  
† Anne Liong and Leandro de Mattos Pereira contributed equally to this work. *Artigo em revisão.**

---

## Contact and Acknowledgements

The entire codebase was implemented by Leandro de Mattos Pereira during his postdoctoral appointment in the group led by Dr. Pedro N. Leão.

- Leandro de Mattos Pereira  
- Anne Liong  
- Pedro N. Leão  

For questions, bug reports or feature requests, please open an issue in this repository or contact Leandro de Mattos Pereira (contact: mattoslmp@gmail.com).
The web server and code are supported by CIIMAR/CNP and associated funding agencies, as shown in the app footer (lab and institutional logos).
