# CRC Gut Microbiome–Metabolome Analysis

Extracted from `CRC_Gut_Analysis_Preprocessing.ipynb`.

## Directory Structure

```
01_preprocessing/          Raw TSV files → Preprocessed Data
02_association_maps/       Preprocessed Data → Association Maps
03_validated_associations/ Association Maps → Validated Associations
04_source_attribution/     Validated Associations → Outputs
utils.py                   Shared utility functions
```

## Dependencies

```
pandas numpy scipy scikit-learn matplotlib seaborn
statsmodels pingouin networkx xgboost shap
```

## Data

Place raw TSV files in the path expected by `utils.DATA_DIR`.
