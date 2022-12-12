from lifelines.datasets import load_rossi
from lifelines import CoxPHFitter
import pandas as pd
import matplotlib.pyplot as plt
import logging
import numpy as np
import streamlit as st
from datetime import date
# 1. Basic Setting
today = date.today()
logger = logging.getLogger()
logging.basicConfig(filename = 'output_cox.log', filemode = 'a',
            format="%(asctime)s--%(filename)s--%(message)s", level=logging.INFO)
summary_col = ['coef', 'exp(coef)', 'se(coef)', 'coef lower 95%', 'coef upper 95%',
       'exp(coef) lower 95%', 'exp(coef) upper 95%', 'cmp to', 'z', 'p',
       '-log2(p)']
dur_eve = ['week','arrest']


    # return final_tb[final_tb['p']<0.05] #filtered by p value
# 2. Get the summary with p, HR for multicox regression
def get_summary(df,covar, var):
    # Get the hear of summary matrix
    global summary_col, dur_eve
    final_tb  = pd.DataFrame(columns=summary_col)
    for i in range(len(var)):
        covar.append(var[i])
        final_col = covar + dur_eve
        sub_df = df[final_col]
        cph = CoxPHFitter(penalizer=0.0001)
        cph.fit(sub_df,dur_eve[0], dur_eve[1])
        final_tb.loc[var[i]] = cph.summary.loc[var[i]].tolist()

    return final_tb.sort_values(by='coef', ascending=False)

# 3. Start ploting

def plot_cox(df,ax=None, hazard_ratios = False,**errorbar_kwargs):

    if ax is None:
        fig,ax = plt.subplots()

    errorbar_kwargs.setdefault("c", "k")
    errorbar_kwargs.setdefault("fmt", "s")
    errorbar_kwargs.setdefault("markerfacecolor", "white")
    errorbar_kwargs.setdefault("markeredgewidth", 1.25)
    errorbar_kwargs.setdefault("elinewidth", 1.25)
    errorbar_kwargs.setdefault("capsize", 3)

    columns = df.index
    yaxis_locations = list(range(len(columns)))
    
    if hazard_ratios: #exp(coef)
        exp_log_hazards = df['exp(coef)'].tolist() 
        errors = np.subtract(df['exp(coef) upper 95%'].tolist(),df['exp(coef)'].tolist())
        ax.errorbar(
            exp_log_hazards,
            yaxis_locations,
            xerr=errors,
            **errorbar_kwargs)
        ax.set_xlabel("HR (%g%% CI)" % ((1) * 100))

    else: #coef
        log_hazards = df['coef'].tolist() 
        errors = np.subtract(df['coef upper 95%'].tolist(),df['coef'].tolist()) 
        ax.errorbar(log_hazards, yaxis_locations, 
            xerr=errors, 
            **errorbar_kwargs)
        ax.set_xlabel("log(HR) (%g%% CI)" % ((1) * 100))



    best_ylim = ax.get_ylim()
    ax.vlines(1 if hazard_ratios else 0, -2, len(columns) + 1, linestyles="dashed", linewidths=1, alpha=0.65, color="k")
    ax.set_ylim(best_ylim)

    tick_labels = columns

    ax.set_yticks(yaxis_locations)
    ax.set_yticklabels(tick_labels)

    return fig

df = load_rossi()



count_matrices = ('tpm',
'tumor',
'Tcell',
'CD8Tex')


# 4. Page display
st.sidebar.text(f"Edited by Xiyuan\nModifyied on {today.strftime('%b %d, %Y')}")
st.write('This is the test for cox regresssion using rossi dataset')
options = st.multiselect(
    'Select covariants:',
    ['fin','age'])
final_tb = get_summary(df,covar = options, var = ['race','mar'])
cox_fig = plot_cox(final_tb)
st.selectbox("Select Dataset: ", count_matrices)
cox = st.button('Cox regression')
if cox:
    st.pyplot(cox_fig)
