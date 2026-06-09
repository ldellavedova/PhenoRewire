# data/

Place your input files here before running PhenoRewire.

PhenoRewire requires two comma-separated files:

| File | Description |
|------|-------------|
| `feature_table.csv` | One row per metabolite feature. Include identifier columns such as feature ID, m/z, retention time, optional name, and one intensity column per sample. |
| `metadata.csv` | One row per sample. Include the sample identifier plus the variables used to define biological groups and/or timepoints. |

The root config templates (`config_phenotype.yaml` and `config_temporal.yaml`) already point to this folder. Replace the example files with your own data and then update the column names and group labels in the config.

The files currently in this folder are synthetic and are included so you can test the package immediately after installation.

For the optional 2D molecular network annotation workflow, this folder also includes:

| File | Description |
|------|-------------|
| `example_molecular_network.graphml` | Small demo molecular/chemical similarity network whose node IDs match the example `feature_table.csv` feature IDs. It is intended only for testing `phenorewire --config ... --run-2dnetwork`. |
