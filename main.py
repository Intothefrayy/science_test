import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import re
import string
from scipy import stats
import textwrap

# --- Journal Guideline Settings ---
# Set font to Arial, per journal requirements
plt.rcParams['font.family'] = 'Arial'
# Ensure fonts are embedded correctly in vector formats (like EPS)
plt.rcParams['pdf.fonttype'] = 42

# Figure dimensions in mm, converted to inches for matplotlib
width_mm = 82
height_mm = 120 # A reasonable height, less than the 235mm max
width_in = width_mm / 25.4
height_in = height_mm / 25.4
# --- End of Settings ---


def safe_convert(x):
    if isinstance(x, str):
        x = x.strip()  # removing empty spaces at the edges
        # If the cleaned string is empty or equals "none", return NaN
        if x == "" or x.lower() == "none":
            return np.nan
        # If the string contains a fraction in the form 'numerator/denominator'
        if '/' in x:
            try:
                numerator, denominator = x.split('/')
                return float(numerator) / float(denominator)
            except Exception:
                return np.nan
        else:
            try:
                return float(x)
            except Exception:
                return np.nan
    else:
        return x

# Loading data from CSV-file.
# File should contain following columns: Group, ActiveZoneLength, ActiveZoneWidth, PreDistance, PreArea, PreDiameter, PostDistance, PostDiameter, PostArea
data = pd.read_csv("amygdala_data.csv", na_values=["", " ", "none", "destroyed"], skipinitialspace=True)

# List of parameters to compare
parameters = ['Length of AZ(nm)', 'Width of AZ(nm)', 'Distance bw presyn mitochon & AZ', 'Area of pre Mitoch',
              'D of pre Mitoch', 'Area of post Mitoch', 'D of post Mitoch', 'Distance bw postsyn mitoch & AZ(nm)']

#print (data)
#print (parameters)


for param in parameters:
    # Apply the safe_convert function to all elements of the column.
    data[param] = data[param].apply(safe_convert)

# Groups found checking
groups_found = data['group'].unique()
print("Found groups:", groups_found)

# Define comparison pairs: control vs exp5, control vs exp10, exp5 vs exp10
#comparisons = [("control", "PPA 5 days"), ("control", "PPA 10 days"), ("PPA 10 days", "PPA 5 days")]
comparisons = [("control", "PPA 5 days"), ("control", "PPA 10 days"), ("PPA 5 days", "PPA 10 days")]


def compare_groups(data1, data2):
    """
    This function takes two samples (as pandas Series) for a given parameter,
    checks for normality of distribution, and selects the appropriate statistical test.
    Returns:
      - the name of the test used,
      - the p-value,
      - the effect size (Cohen's d for t-test) or None for a non-parametric test,
      - normality test results (Shapiro–Wilk p-values for each sample).
"""
    # Excluding empty spaces
    data1 = data1.dropna()
    data2 = data2.dropna()

    # If there is insufficient data, return a message
    if len(data1) < 3 or len(data2) < 3:
        return None, "Not enough data", None, (None, None)

    # Normality check using the Shapiro–Wilk test
    stat1, p_norm1 = stats.shapiro(data1)
    stat2, p_norm2 = stats.shapiro(data2)
    normal1 = (p_norm1 > 0.05)
    normal2 = (p_norm2 > 0.05)

    # If both samples are normally distributed, the Welch’s t-test is used.
    if normal1 and normal2:
        t_stat, p_value = stats.ttest_ind(data1, data2, equal_var=False)
        test_used = "Welch's t-test"
        # Calculating Cohen’s d
        diff = data1.mean() - data2.mean()
        pooled_std = np.sqrt((data1.std() ** 2 + data2.std() ** 2) / 2)
        effect_size = diff / pooled_std
    else:
        # When the distribution is non-normal, the Mann–Whitney U test is used.
        u_stat, p_value = stats.mannwhitneyu(data1, data2, alternative='two-sided')
        test_used = "Mann–Whitney U test"
        effect_size = None  # Effect size calculation for non-parametric tests can be added separately

    return test_used, p_value, effect_size, (p_norm1, p_norm2)


# A dictionary to store results for each parameter and comparison pair
results = {}

for param in parameters:
    results[param] = {}
    for group1, group2 in comparisons:
        group_data1 = data[data['group'] == group1][param]
        group_data2 = data[data['group'] == group2][param]
        test_used, p_value, effect_size, norm_results = compare_groups(group_data1, group_data2)
        results[param][f"{group1} vs {group2}"] = (test_used, p_value, effect_size, norm_results)

# Reporting results with interpretation
alpha = 0.05  # significance threshold
for param, comp in results.items():
    print(f"\nParameter: {param}")
    for comparison_name, res in comp.items():
        test_used, p_value, effect_size, (p_norm1, p_norm2) = res

        group1, group2 = comparison_name.split(" vs ")

        # Calculating medium values and standard errors
        group1_data = data[data['group'] == group1][param].dropna()
        group2_data = data[data['group'] == group2][param].dropna()
        mean1 = group1_data.mean()
        se1 = group1_data.std(ddof=1) / np.sqrt(len(group1_data)) if len(group1_data) > 0 else np.nan
        mean2 = group2_data.mean()
        se2 = group2_data.std(ddof=1) / np.sqrt(len(group2_data)) if len(group2_data) > 0 else np.nan

        # If there is no data
        if test_used is None:
            print(f"  {comparison_name}: {p_value}")
            continue

        # Interpreting normality
        group1, group2 = comparison_name.split(" vs ")
        norm_msg = (f"{group1}: {'normal' if p_norm1 > 0.05 else 'ubnormal'}, "
                    f"{group2}: {'normal' if p_norm2 > 0.05 else 'ubnormal'}")

        # Interpreting test results
        significance = "statistically relevant" if p_value < alpha else "not statistically relevant"
        result_msg = (f"Differences between groups {group1} and {group2} {significance} (p = {p_value:.4f}).")

        print(f"  Comparison {comparison_name}:")
        print(f"    Normality test: {norm_msg}")
        print(f"    Test used: {test_used}")
        print(f"    p-value: {p_value:.4f}")
        if effect_size is not None:
            # If the effect size has been calculated, the direction of the difference can be specified.
            direction = "more" if effect_size > 0 else "less"
            print(f"    Effect size (Cohen's d): {effect_size:.4f}")
            print(f"    Interpretation: on average values in group {group1} {direction} values group {group2}.")
        else:
            print("    Effect size not calculated for non-parametric test.")
        print(f"    Mean Values + Std. Er., Group: {group1}: {mean1:.2f}+{se1:.2f}")
        print(f"    Mean Values + Std. Er., Group: {group2}: {mean2:.2f}+{se2:.2f}")
        print(f"    Conclusion: {result_msg}")

#visuals

group_names = ["control", "PPA 5 days", "PPA 10 days"]

# The dictionary of measurement units is used ONLY for plotting

units = {
    'Length of AZ(nm)': 'nm', 'Width of AZ(nm)': 'nm', 'Distance bw presyn mitochon & AZ': 'nm',
    'Area of pre Mitoch': 'nm²', 'Area of post Mitoch': 'nm²', 'Distance bw postsyn mitoch & AZ(nm)': 'nm',
}

display_names = {
    'Distance bw presyn mitochon & AZ': 'Distance between presyn mitochon & AZ',
    'Distance bw postsyn mitoch & AZ(nm)': 'Distance between postsyn mitoch & AZ(nm)'
}

# Select parameters for plotting (those that have defined units)
parameters_for_plotting = list(units.keys())

def get_asterisks(p):
    """Returns significance stars based on the p-value."""
    if p is None or p > 0.05: return ""
    if p < 0.001: return "***"
    if p < 0.01: return "**"
    return "*"


def plot_param(ax, param, letter_label, units_dict, display_names):
    """Plots a single chart with smart wrapping of long titles and extra padding on the right."""
    arrs = [data[data['group'] == g][param].dropna().values for g in group_names]
    x_pos = np.arange(len(group_names))

    display_param = display_names.get(param, param)

    parts = ax.violinplot(arrs, positions=x_pos, showmedians=False, showextrema=False)
    for pc in parts['bodies']:
        pc.set_facecolor('white')
        pc.set_edgecolor('black')
        pc.set_linewidth(1)

    for xi, arr in zip(x_pos, arrs):
        jitter = np.random.normal(0, 0.04, size=len(arr))
        ax.scatter(xi + jitter, arr, s=15, facecolors='none', edgecolors='black', linewidth=0.7, zorder=2)

    medians = [np.median(a) for a in arrs if len(a) > 0]
    ax.scatter(x_pos, medians, marker='_', color='black', s=200, zorder=3, linewidth=1.5)

    label = display_param.split('(')[0].strip()
    unit = units_dict.get(param, '')
    title_text = f"{label} ({unit})" if unit else label
    wrapped_title = "\n".join(textwrap.wrap(title_text, width=35))

    #title management

    #ax.set_title(wrapped_title, fontsize=12, pad=15)

    ax.set_ylabel(f"{label} ({unit})", fontsize=10)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(group_names, rotation=30, ha='right')

    #letter labeling making

    #ax.text(-0.2, 1.17, letter_label, transform=ax.transAxes, fontsize=14, va='top', ha='left')
    ax.tick_params(axis='both', which='major', labelsize=8)
    ax.spines[['right', 'top']].set_visible(True)
    ax.grid(axis='y', linestyle=':', color='gray', linewidth=0.5)

    y_max = data[param].max()
    bracket_configs = {
        ("control", "PPA 5 days"): {'pos': [0, 1], 'height_factor': 1.05},
        ("PPA 5 days", "PPA 10 days"): {'pos': [1, 2], 'height_factor': 1.15},
        ("control", "PPA 10 days"): {'pos': [0, 2], 'height_factor': 1.25}
    }
    for (g1, g2), config in bracket_configs.items():
        p_value = results[param].get(f"{g1} vs {g2}", (None, 1.0, None, None))[1]
        stars = get_asterisks(p_value)
        if stars:
            x1, x2 = config['pos']
            height = config['height_factor'] * y_max
            bar_height = height * 0.02
            ax.plot([x1, x1, x2, x2], [height, height + bar_height, height + bar_height, height], c='black', lw=1)
            ax.text((x1 + x2) / 2, height + bar_height, stars, ha='center', va='bottom', fontsize=10)

    # Set axis limits
    ax.set_ylim(bottom=0, top=y_max * 1.4)
    # --- Adding white space on the right ---
    # The data points are located at 0, 1, 2 We're extending the right axis limit to 2.5 to create extra padding.
    ax.set_xlim(left=-0.5, right=2.5)


# Constructing and saving figures
for i, param in enumerate(parameters_for_plotting):
    fig, ax = plt.subplots(figsize=(width_in, height_in))
    letter_label = chr(ord('A') + i)

    plot_param(ax, param, letter_label, units, display_names)

    # --- Controlling spacing by adjusting figure margins ---
    # Replacing fig.tight_layout() with plt.subplots_adjust() for precise control.
    # The values are fractions of the figure's width and height.
    # right=0.9 means the right edge of the plot will be at 90% of the canvas width,
    # leaving 10% white space on the right.
    plt.subplots_adjust(left=0.2, right=0.9, top=0.85, bottom=0.2)

    filename = f"Fig{i + 1}_{width_mm}mm.tiff"
    plt.savefig(filename, format='tiff', dpi=900)
    plt.close(fig)
    print(f"Graph has been saves as: {filename}")

