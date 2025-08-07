
# Analysis and Visualization of Synaptic Morphometric Data

This Python script is designed to automate statistical analysis and generate publication-ready plots from morphometric data. It performs comparisons between different experimental groups, automatically selects the appropriate statistical test, and generates violin plots with the analysis results.

-----

## 🧬 Key Features

  * **Data Loading**: Loads data from a `.csv` file.
  * **Automatic Data Processing**: Cleans missing values (`"none"`, `"destroyed"`, empty cells) and converts fractional formats (e.g., "1/2") into numbers.
  * **Statistical Analysis**:
      * Checks for data normality using the **Shapiro-Wilk test**.
      * Automatically chooses between **Welch's t-test** (for normal distributions) and the **Mann-Whitney U test** (for non-normal distributions).
      * Calculates **p-values** and **Cohen's d effect size** for parametric tests.
  * **Publication-Style Visualization**:
      * Generates violin plots overlaid with individual data points (jitter plot).
      * Adheres to formatting requirements: Arial font, specified figure dimensions.
      * Automatically adds brackets with statistical significance stars (`*`, `**`, `***`).
      * Saves plots as high-resolution (900 DPI) **TIFF** files.

-----

## 🛠️ Requirements

  * Python 3.6+
  * Python Libraries:
      * `pandas`
      * `matplotlib`
      * `numpy`
      * `scipy`

### Installing Dependencies

You can install all the necessary libraries with a single command using `pip`:

```bash
pip install pandas matplotlib numpy scipy
```

-----

## 🚀 Usage

1.  **Prepare your data**. Ensure you have an `amygdala_data.csv` file in the same directory as the script.
2.  **Configure parameters** within the script (see the "Customization" section).
3.  **Run the script** from your terminal or command line:
    ```bash
    python main.py
    ```
4.  **Get the results**. Statistical summaries will be printed to the console, and the plots will be saved as `.tif` files in the script's directory.

### Input Data Format

The `amygdala_data.csv` file must contain at least the following columns:

  * `group`: The name of the experimental group (e.g., "control", "PPA 5 days").
  * Columns with the measured parameters, e.g., `Length of AZ(nm)`, `Width of AZ(nm)`, etc.

-----

## ⚙️ Customization

You can easily adapt the script to your needs by modifying the variables in its upper section.

  * **Changing Groups and Comparisons**:

    ```python
    # Define the comparison pairs
    comparisons = [("control", "PPA 5 days"), ("control", "PPA 10 days"), ("PPA 5 days", "PPA 10 days")]
    ```

  * **Changing the List of Analyzed Parameters**:

    ```python
    # List of parameters to compare
    parameters = ['Length of AZ(nm)', 'Width of AZ(nm)', 'Distance bw presyn mitochon & AZ', ...]
    ```

  * **Customizing Plot Labels**:
    To replace abbreviated parameter names with full names (e.g., "bw" with "between") on the plots only, use the `display_names` dictionary:

    ```python
    display_names = {
        'Distance bw presyn mitochon & AZ': 'Distance between presyn mitochon & AZ',
        'Distance bw postsyn mitoch & AZ(nm)': 'Distance between postsyn mitoch & AZ(nm)'
    }
    ```

  * **Customizing Plot Appearance**:
    All visualization settings are located within the `plot_param()` function. You can change:

      * **Figure dimensions**: the `width_mm` and `height_mm` variables.
      * **Title display**: the `ax.set_title(...)` line is commented out.
      * **Letter labels (A, B, C...)**: the `ax.text(...)` line is commented out.

-----

## 📤 Output

  * **Console Output**: Detailed statistical information for each comparison, including the test used, p-value, mean values, and standard errors.
  * **Image Files**: For each parameter in the `parameters_for_plotting` list, a separate high-resolution TIFF file is created (e.g., `Fig1_82mm.tif`, `Fig2_82mm.tif`).
