# guide-assignment

This module is part of the [CRISPRa Analysis Pipeline](https://github.com/mangochiral/CRISPRa_Analysis_pipeline) and provides tools for assigning guide RNAs (gRNAs) to genetic targets for CRISPR activation (CRISPRa) screening experiments.

## Features

- Assigns candidate guide RNAs to genes or other genomic features.
- Supports input from custom design files (e.g., CSV, Excel, or BED).
- Provides filtering based on sequence quality, location, and user-defined constraints.
- Outputs assignment tables compatible with downstream analysis and visualization.
- Easily integrates with the overall CRISPRa pipeline.

## Usage

1. **Input Preparation**
   - Prepare a file containing guide RNA sequences and their intended targets.
   - Supported formats: `.csv`, `.tsv`, `.xlsx`, or customized formats as described in the code.

2. **Running guide-assignment**
   - Launch the corresponding Jupyter Notebook or Python script.
   - Example command (if script-based):

     ```bash
     python assign_guides.py --input guides.csv --output assigned_guides.csv
     ```

   - Or open the notebook in JupyterLab and proceed through the documented cells.

3. **Parameters and Options**
   - Filtering options (e.g., minimum on-target score).
   - Custom assignments by region or transcript.
   - Aggregated guide statistics.

Refer to the provided notebooks or scripts for in-depth usage and parameter settings.

## Example

```python
# Example function call
from assign_guides import assign_guides_to_targets

assign_guides_to_targets('guides.csv', 'targets.csv', output='assigned_guides.csv')
```

Check the `example/` subdirectory for demonstration files and outputs.

## Dependencies

- Python 3.7+
- pandas
- numpy
- Jupyter Notebook (optional, for interactive runs)

Install with:

```bash
pip install -r requirements.txt
```

## License

This module is released under the MIT License. See the [main repository](https://github.com/mangochiral/CRISPRa_Analysis_pipeline) for details.

## Acknowledgements

Inspired by published CRISPRa libraries and community-developed analysis pipelines.

---

For questions or contributions, please open an Issue or Pull Request on the [main repository](https://github.com/mangochiral/CRISPRa_Analysis_pipeline).
