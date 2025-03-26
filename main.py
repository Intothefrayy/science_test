import pandas as pd
import numpy as np
from scipy import stats

def safe_convert(x):
    if isinstance(x, str):
        x = x.strip()  # удаляем пробелы по краям
        # Если после очистки строка пуста или равна "none", возвращаем NaN
        if x == "" or x.lower() == "none":
            return np.nan
        # Если строка содержит дробь вида "числитель/знаменатель"
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

# Загрузка данных из CSV-файла.
# Файл должен содержать колонки: Group, ActiveZoneLength, ActiveZoneWidth, PreDistance, PreArea, PreDiameter, PostDistance, PostDiameter, PostArea
data = pd.read_csv("put_your_filepath", na_values=["", " ", "none", "destroyed"], skipinitialspace=True)

# Список параметров для сравнения
parameters = ['Length of AZ(nm)', 'Width of AZ(nm)', 'Distance bw presyn mitochon & AZ', 'Area of pre Mitoch',
              'D of pre Mitoch', 'Area of post Mitoch', 'D of post Mitoch', 'Distance bw postsyn mitoch & AZ(nm)']

#print (data)
#print (parameters)


for param in parameters:
    # Применяем функцию safe_convert ко всем элементам колонки
    data[param] = data[param].apply(safe_convert)

# Проверка найденных групп
groups_found = data['group'].unique()
print("Found groups:", groups_found)

# Определяем пары сравнений: контроль с exp5, контроль с exp10, exp5 с exp10
#comparisons = [("control", "PPA 5 days"), ("control", "PPA 10 days"), ("PPA 10 days", "PPA 5 days")]
comparisons = [("control", "PPA 5 days"), ("control", "PPA 10 days"), ("PPA 5 days", "PPA 10 days")]


def compare_groups(data1, data2):
    """
    Функция принимает две выборки (pandas Series) для заданного параметра,
    проверяет нормальность распределения и выбирает подходящий статистический тест.
    Возвращает:
      - название теста,
      - p-значение,
      - размер эффекта (Cohen's d для t-теста) или None для непараметрического теста,
      - результаты теста нормальности (p-значения Shapiro–Wilk для каждой выборки).
    """
    # Исключаем пропуски
    data1 = data1.dropna()
    data2 = data2.dropna()

    # Если данных недостаточно, возвращаем сообщение
    if len(data1) < 3 or len(data2) < 3:
        return None, "Not enough data", None, (None, None)

    # Проверка нормальности распределения (Shapiro–Wilk test)
    stat1, p_norm1 = stats.shapiro(data1)
    stat2, p_norm2 = stats.shapiro(data2)
    normal1 = (p_norm1 > 0.05)
    normal2 = (p_norm2 > 0.05)

    # Если обе выборки нормально распределены – используем Welch-тест
    if normal1 and normal2:
        t_stat, p_value = stats.ttest_ind(data1, data2, equal_var=False)
        test_used = "Welch's t-test"
        # Расчёт Cohen’s d
        diff = data1.mean() - data2.mean()
        pooled_std = np.sqrt((data1.std() ** 2 + data2.std() ** 2) / 2)
        effect_size = diff / pooled_std
    else:
        # При ненормальном распределении – используем Mann–Whitney U test
        u_stat, p_value = stats.mannwhitneyu(data1, data2, alternative='two-sided')
        test_used = "Mann–Whitney U test"
        effect_size = None  # Расчёт размера эффекта для непараметрических тестов можно добавить отдельно

    return test_used, p_value, effect_size, (p_norm1, p_norm2)


# Словарь для хранения результатов по каждому параметру и паре сравнения
results = {}

for param in parameters:
    results[param] = {}
    for group1, group2 in comparisons:
        group_data1 = data[data['group'] == group1][param]
        group_data2 = data[data['group'] == group2][param]
        test_used, p_value, effect_size, norm_results = compare_groups(group_data1, group_data2)
        results[param][f"{group1} vs {group2}"] = (test_used, p_value, effect_size, norm_results)

# Вывод результатов с интерпретацией
alpha = 0.05  # порог значимости
for param, comp in results.items():
    print(f"\nParameter: {param}")
    for comparison_name, res in comp.items():
        test_used, p_value, effect_size, (p_norm1, p_norm2) = res

        group1, group2 = comparison_name.split(" vs ")

        # Расчёт средних значений и стандартной ошибки для каждой группы
        group1_data = data[data['group'] == group1][param].dropna()
        group2_data = data[data['group'] == group2][param].dropna()
        mean1 = group1_data.mean()
        se1 = group1_data.std(ddof=1) / np.sqrt(len(group1_data)) if len(group1_data) > 0 else np.nan
        mean2 = group2_data.mean()
        se2 = group2_data.std(ddof=1) / np.sqrt(len(group2_data)) if len(group2_data) > 0 else np.nan

        # Если данных недостаточно
        if test_used is None:
            print(f"  {comparison_name}: {p_value}")
            continue

        # Интерпретация нормальности
        group1, group2 = comparison_name.split(" vs ")
        norm_msg = (f"{group1}: {'normal' if p_norm1 > 0.05 else 'ubnormal'}, "
                    f"{group2}: {'normal' if p_norm2 > 0.05 else 'ubnormal'}")

        # Интерпретация результата теста
        significance = "statistically relevant" if p_value < alpha else "not statistically relevant"
        result_msg = (f"Differences between groups {group1} and {group2} {significance} (p = {p_value:.4f}).")

        print(f"  Comparison {comparison_name}:")
        print(f"    Normality test: {norm_msg}")
        print(f"    Test used: {test_used}")
        print(f"    p-value: {p_value:.4f}")
        if effect_size is not None:
            # Если размер эффекта рассчитан, можно указать направление разницы
            direction = "more" if effect_size > 0 else "less"
            print(f"    Effect size (Cohen's d): {effect_size:.4f}")
            print(f"    Interpretation: on average values in group {group1} {direction} values group {group2}.")
        else:
            print("    Effect size not calculated for non-parametric test.")
        print(f"    Mean Values + Std. Er., Group: {group1}: {mean1:.2f}+{se1:.2f}")
        print(f"    Mean Values + Std. Er., Group: {group2}: {mean2:.2f}+{se2:.2f}")
        print(f"    Conclusion: {result_msg}")

import matplotlib.pyplot as plt
import seaborn as sns

# Предполагается, что DataFrame data и список параметров parameters уже созданы и обработаны
for param in parameters:
    plt.figure(figsize=(8, 6))
    sns.boxplot(x='group', y=param, data=data)
    plt.title(f"Boxplot for {param}")
    plt.xlabel("Group")
    plt.ylabel(param)
    plt.show()

import matplotlib.pyplot as plt
import numpy as np

# Список групп, если он фиксирован
group_names = ["control", "PPA 5 days", "PPA 10 days"]

for param in parameters:
    means = []
    errors = []
    for grp in group_names:
        grp_data = data[data['group'] == grp][param].dropna()
        means.append(grp_data.mean())
        errors.append(grp_data.std(ddof=1) / np.sqrt(len(grp_data)) if len(grp_data) > 0 else 0)
    plt.figure(figsize=(8, 6))
    plt.bar(group_names, means, yerr=errors, capsize=5, color='skyblue', edgecolor='black')
    plt.title(f"Bar plot for {param}")
    plt.xlabel("Group")
    plt.ylabel(param)
    plt.show()
