import numpy as np
import matplotlib.pyplot as plt
import argparse
from lifelines.datasets import load_rossi
import os
from lifelines import CoxPHFitter
import pandas as pd

summary_col = ['coef', 'exp(coef)', 'se(coef)', 'coef lower 95%', 'coef upper 95%',
       'exp(coef) lower 95%', 'exp(coef) upper 95%', 'cmp to', 'z', 'p',
       '-log2(p)']
dur_eve = ['week','arrest']

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
        tmp = cph.summary
        final_tb.loc[var[i]] = tmp.loc[var[i]].tolist()

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



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot harzard ratoio(HR) for uni variable analysis with assigned covariates')
    parser.add_argument('-c', '--covariates', action='store', dest='covar',
                    type=str, nargs='*', default=['race'] ,
                    help="Examples: -covariates item1 item2 item3")
    parser.add_argument('-v','--variates',action='store',dest='var',
                    type=str, nargs='*',default=['mar'],help = 'Vaiants for uni-variable cox regression, seperated by space')
    
    parser.add_argument('--output',type = str,required=True,help = 'Dictionary to store the cox fig')
    parser.add_argument('--hr',type = bool,help = 'True: exp(coef), False: coef')
    
    args = parser.parse_args()
    df = load_rossi()
    df = get_summary(df,covar = args.covar, var = args.var)
    fig = plot_cox(df, hazard_ratios= args.hr)
    if not os.path.exists():
        os.makedirs(f'{args.output}')
    plt.savefig(f"{args.output}/cox.png")
    

