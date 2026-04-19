🧠 MOXY-bids: Raw Acquisition Data (BIDS Format)

This directory contains the raw data acquired during the experimental phase of the OxyMove project. The data is organized according to the BIDS (Brain Imaging Data Structure) standard, specifically adapted for NIRS (Near-Infrared Spectroscopy).
📌 BIDS Structure

The BIDS format ensures standardization, making the dataset interoperable and easy to process with automated tools such as mne-nirs.
Plaintext

/MOXY-bids
├── dataset_description.json  # Global project metadata
├── participants.tsv          # Demographic information (age, sex, etc.)
├── participants.json         # Data dictionary for participants
├── sub-01/                   # Data for Subject 01
│   └── nirs/
│       ├── sub-01_task-contraction_nirs.snirf  # Raw NIRS data
│       └── sub-01_task-contraction_nirs.json   # Acquisition parameters
├── sub-02/                   # Data for Subject 02
│   └── ...
└── ... (N=23)

📊 Acquisition Details

The data stored here reflects the experimental protocol described in the main analysis pipeline:

    Modality: NIRS (Near-Infrared Spectroscopy).

    Task: Muscle contractions (Concentric vs. Eccentric) at 30% and 50% of MVC (Maximum Voluntary Contraction).

    Equipment: Semaxone and Portalite devices (stored in .snirf format).

🛠️ Usage

The .snirf files in this directory are designed to be processed by the Python scripts located in the parent directory:

    Processing Script: main_scripts/NIRS/NIRS_Treatment/mNIRS_V0.0.1.ipynb.

    Recommended Libraries: mne-bids-nirs or the Snirf library for Python.

📝 Naming Conventions

    Subjects: sub-XX (e.g., sub-01 to sub-23).

    Task: task-contraction.

    File Format: Standard .snirf to ensure cross-platform compatibility for analysis.

⚠️ Privacy & Anonymization

All data has been anonymized in accordance with health data protection guidelines. No real names or personally identifiable information (PII) are stored in this repository.
