import pandas as pd
import statsmodels.api as sm

def load_phenotype():
    phenotype_path = "~/projects/uosa/Silvia_Paracchini/20240813_Grip_PRS_bayesian_Filippo/quant.hand.10102021.full.sample.with.pcs"
    phenotype = pd.read_table(phenotype_path)
    phenotype.dropna(subset=['grip.best.hand', 'grip.worst.hand', 'sex', 'age_grip',
                         'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10'], inplace=True)
    phenotype_subset = phenotype[['grip.best.hand', 'grip.worst.hand', 'FID', 'IID', 'sex', 'age_grip',
                              'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10']]
    phenotype_subset.rename(columns={'grip.best.hand': 'GSD', 'grip.worst.hand': 'GSND'}, inplace=True)
    return phenotype_subset

def load_plink_scores(base_gwas):
    plink_scores = f"~/scratch/private/PRScs/results/plink_score_{base_gwas}.profile"
    plink_scores = pd.read_table(plink_scores, sep='\s+')
    return plink_scores

base_GWAS = ['AD', 'ADHD', 'ASD', 'BIP', 'GEOS', 'HBD.L', 'HBD.R','SCZ', 'UKB_CAD']
phenotypes = ['GSD', 'GSND']

phenotype_file = load_phenotype()

all_dfs = []
for pheno in phenotypes:
    for gwas in base_GWAS:
        plink_scores = load_plink_scores(gwas)
        merged = phenotype_file.merge(plink_scores, on=["FID", "IID"])
        X = merged[[pheno, 'SCORESUM', 'sex', 'age_grip',
                              'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10']]
        X.dropna(inplace=True)
        y = X.pop(pheno)    
        X = sm.add_constant(X) # Add a constant to the independent variable (required for statsmodels)
        model = sm.OLS(y, X).fit()
        
        r2 = model.rsquared
        r2_adjusted = model.rsquared_adj
        coef = model.params['SCORESUM']
        std_err = model.bse['SCORESUM']
        p_value = model.pvalues['SCORESUM']
        
        # Creating a DataFrame
        results_df = pd.DataFrame({
            'Target trait': [pheno],
            'Base GWAS': [gwas],
            'R2': [r2],
            'Adjusted R2': [r2_adjusted],
            'Coefficient': [coef],
            'Std Error': [std_err],
            'P>|t|': [p_value]
        })
        all_dfs.append(results_df)

df = pd.concat(all_dfs, ignore_index=True)

df.to_csv("PRScs_grip_results.csv", index=False)