import pandas as pd
from typing import Union, Tuple, Literal
from scipy.stats import ttest_ind, mannwhitneyu
import statsmodels.api as sm
import statsmodels.formula.api as smf


def compare_signal_groups(df: pd.DataFrame,
                          group_col: str,
                          value_col: str,
                          test: Literal["t-test", "mannwhitney"] = "t-test") -> float:
    """
    Compare signal values between two groups using a statistical test.

    Parameters:
    - df: DataFrame containing the data.
    - group_col: Column name with group labels (must have exactly two unique values).
    - value_col: Column name with signal values.
    - test: Which statistical test to use ('t-test' or 'mannwhitney').

    Returns:
    - float: p-value of the test.

    Raises:
    - ValueError: If the group column does not contain exactly two groups.
    """
    groups = df[group_col].unique()
    if len(groups) != 2:
        raise ValueError(f"Expected exactly 2 groups in '{group_col}', found {len(groups)}")

    g1 = df[df[group_col] == groups[0]][value_col]
    g2 = df[df[group_col] == groups[1]][value_col]

    if test == "t-test":
        _, p = ttest_ind(g1, g2, equal_var=False)
    elif test == "mannwhitney":
        _, p = mannwhitneyu(g1, g2, alternative="two-sided")
    else:
        raise ValueError("test must be 't-test' or 'mannwhitney'")

    return p

def run_one_way_anova(df: pd.DataFrame,
                      factor: str,
                      response: str,
                      return_model: bool = False) -> Union[pd.DataFrame, Tuple[pd.DataFrame, sm.regression.linear_model.RegressionResultsWrapper]]:
    """
    Run one-way ANOVA on a numeric response grouped by a single categorical factor.

    Parameters:
    - df (pd.DataFrame): Input DataFrame.
    - factor (str): Name of the categorical grouping column.
    - response (str): Name of the numeric response column.
    - return_model (bool): If True, also return the fitted statsmodels OLS model.

    Returns:
    - pd.DataFrame: ANOVA summary table.
    - RegressionResultsWrapper (optional): The fitted model if return_model=True.

    Raises:
    - ValueError: If the specified columns are missing from the DataFrame.
    """
    if factor not in df.columns or response not in df.columns:
        raise ValueError(f"Missing columns: {factor}, {response}")

    df_temp = df[[factor, response]].copy().rename(columns={response: "y", factor: "x"})
    model = smf.ols("y ~ C(x)", data=df_temp).fit()
    anova_table = sm.stats.anova_lm(model, typ=2)

    return (anova_table, model) if return_model else anova_table

def run_two_way_anova(df: pd.DataFrame,
                      factor1: str,
                      factor2: str,
                      response: str,
                      include_interaction: bool = True,
                      return_model: bool = False) -> Union[pd.DataFrame, Tuple[pd.DataFrame, sm.regression.linear_model.RegressionResultsWrapper]]:
    """
    Run two-way ANOVA on a numeric response using two categorical factors.

    Parameters:
    - df (pd.DataFrame): Input DataFrame.
    - factor1 (str): First categorical factor.
    - factor2 (str): Second categorical factor.
    - response (str): Name of the numeric response column.
    - include_interaction (bool): Whether to include the interaction term C(factor1):C(factor2).
    - return_model (bool): If True, also return the fitted statsmodels OLS model.

    Returns:
    - pd.DataFrame: ANOVA summary table.
    - RegressionResultsWrapper (optional): The fitted model if return_model=True.

    Raises:
    - ValueError: If any of the specified columns are missing from the DataFrame.
    """
    for col in [factor1, factor2, response]:
        if col not in df.columns:
            raise ValueError(f"Missing column: {col}")

    df_temp = df[[factor1, factor2, response]].copy().rename(
        columns={factor1: "a", factor2: "b", response: "y"}
    )

    formula = "y ~ C(a) + C(b)"
    if include_interaction:
        formula += " + C(a):C(b)"

    model = smf.ols(formula, data=df_temp).fit()
    anova_table = sm.stats.anova_lm(model, typ=2)

    return (anova_table, model) if return_model else anova_table
