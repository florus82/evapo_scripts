from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import train_test_split
from sklearn.ensemble import GradientBoostingRegressor
import seaborn as sns
import joblib
import time
from sklearn import metrics
import sys
sys.path.append('/home/potzschf/repos/')
from helperToolz.helpsters import *
from scipy.stats import gaussian_kde
workhorse = True
if workhorse:
    origin = 'Aldhani/eoagritwin/'
else:
    origin = ''



#################################################################### build dictionary for results, functions and load pixel subsets for agri thresh
# build result container
def makeRES():
    keys = ['Iteration', 'Metric', 'Year', 'Month', 'Thresh', 'Total_sample_size', 'Train_size', 'Test_size', 
            'r2', 'mse', 'rmse', 'y_true', 'y_pred', 'max_depth', 'learning_rate', 'min_samples_leaf', 'max_features']
    vals = [list() for i in range(len(keys))]
    res  = dict(zip(keys, vals))
    return res

def Model(yValues, predictors, CVout, nCores):
    param_grid = {'learning_rate': [0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001],
                  'max_depth': [4, 6],
                  'min_samples_leaf': [3, 5, 9, 17],
                  'max_features': [1.0, 0.5, 0.3, 0.1]}
    est = GradientBoostingRegressor(n_estimators=7500)
    gs_cv = GridSearchCV(est, param_grid, cv=5, refit=True, n_jobs=nCores).fit(predictors, yValues)
    # Write outputs to disk and return elements from function
    joblib.dump(gs_cv, CVout)
    return (gs_cv)

# model performance
def ModPerfor(cv_results, yData, xData, res):
    # Calculate Predictions of the true values
    y_true, y_pred = yData, cv_results.predict(xData)

    res['r2'].append(metrics.r2_score(y_true, y_pred))
    res['mse'].append(metrics.mean_squared_error(y_true, y_pred))
    res['rmse'].append((metrics.mean_squared_error(y_true, y_pred))**0.5)
    res['y_true'].append(list(y_true))
    res['y_pred'].append(list(y_pred))

    # Get parameters from the best estimator
    res['max_depth'].append(cv_results.best_params_['max_depth'])
    res['learning_rate'].append(cv_results.best_params_['learning_rate'])
    res['min_samples_leaf'].append(cv_results.best_params_['min_samples_leaf'])
    res['max_features'].append(cv_results.best_params_['max_features'])
    #return res

# plot function
def hist_with_deciles(data, metric='median', **kwargs):
    sns.histplot(data[metric], kde=True, bins=10, **kwargs)
    deciles = data[metric].quantile([0.1 * i for i in range(1, 10)])
    for p in deciles:
        plt.axvline(p, color='red', linestyle='--', alpha=0.5)

# get weights per decile vor value range
def getKDEweights(values):
    #  Fit KDE once
    kde = gaussian_kde(values)

    # Calculate decile edges
    decile_edges = np.percentile(values, np.arange(0, 110, 10))

    # Calculate midpoints of each decile bin
    decile_midpoints = (decile_edges[:-1] + decile_edges[1:]) / 2

    # Evaluate KDE only at decile midpoints
    kde_weights_per_decile = kde.evaluate(decile_midpoints)

    # Combine decile labels and weights into a DataFrame
    decile_labels_low = [decile_edges[i] for i in range(len(decile_edges)-1)]
    decile_labels_up = [decile_edges[i+1] for i in range(len(decile_edges)-1)]
    df_decile_weights = pd.DataFrame({
        'decile_low': decile_labels_low,
        'decile_up': decile_labels_up,
        'kde_weight': kde_weights_per_decile,
        'kde_weight_norm': kde_weights_per_decile / kde_weights_per_decile.sum()
    })

    return df_decile_weights

def subsetONThresh(rows, cols, df):
    '''
    very customized function, that subsets a df with the columns row & col based on 2 lists containing valid (pairwise) values on which will be subsetted
    '''
    valid_pairs = set(zip(rows, cols)) # use a set, although list() would gain the same number of combis --> set much faster to lookUp
    pairs_in_df = list(zip(df['row'], df['col']))
    mask = [pair in valid_pairs for pair in pairs_in_df]
    return df[mask]

Month_Int_to_Str = {'1': 'January', 
                    '2': 'February',
                    '3': 'March',
                    '4': 'April',
                    '5': 'May',
                    '6': 'June',
                    '7': 'July',
                    '8': 'August',
                    '9': 'September',
                    '10': 'October',
                    '11': 'November',
                    '12': 'December'}

# read csv for valid row_cols for samples to draw. They are based on the share of agriculture (HR Landcover maps) within a S3 pixel
thresh_csv = pd.read_csv(f'/data/{origin}et/Auxiliary/landcover/csv/row_cols.csv')
rowsL, colsL, indL = [], [], []
for col in thresh_csv.columns:
    #print(f'finding indices for {col}')
    nested = [entry.split('_') for entry in thresh_csv[col] if type(entry) == str]
    rows, cols = zip(*nested)
    rowsL.append([int(row) for row in rows])
    colsL.append([int(col) for col in cols])
    indL.append(col)


################################################ build function (per year per month ML-model)
def ModelRun():
    # load data for respective run
    raw_files = getFilelist(f'/data/{origin}et/Training_ML/training_data/combined_extracts_monthly/', 'parquet')
    # subset by year
    for i in range(2017,2025,1):
        # subset by month
        for j in range(3,10,1):
            annual_files = pd.concat([pd.read_parquet(raw_file) for raw_file in raw_files if f'/{i}_' in raw_file and f'_{j}' in raw_file], ignore_index=True)
            # subset by agri-threshold
            for row, col, ind in zip(rowsL, colsL, indL):
                if ind != 'Thresh100':
                    continue
                sub_annual_files = subsetONThresh(row, col, annual_files)
                # subset on metric
                for metric in ['median']:#, 'mean']:
                    print(f'Working on year: {i}, month: {Month_Int_to_Str[str(j)]}, {ind} and metric: {metric}')
                    if metric == 'median':
                        dropper = 'mean'
                    else:
                        dropper = 'median'
                    
                    # set sample size & test_ratio & number of cores to use
                    samp_size = 3000
                    test_ratio = 0.66
                    number_cores = 15

                    # get sample size / 10 samples per decile of data
                    weights = getKDEweights(sub_annual_files[metric])
                    subs = []
                    for binni in range(weights.shape[0]):
                        subs.append(sub_annual_files[(sub_annual_files[metric] >= weights['decile_low'][binni]) & (sub_annual_files[metric] < weights['decile_up'][binni])].sample(int(samp_size/10)))
                    sub_df = pd.concat(subs, ignore_index=True)
                    
                    # draw sample randomly from kde weighted distribution of LST values per iteration
                    for iteri in range(0,5):
                        res = makeRES()
                        res['Iteration'].append(str(iteri))
                        res['Metric'].append(metric)
                        res['Year'].append(str(i))
                        res['Month'].append(Month_Int_to_Str[str(j)])
                        res['Thresh'].append(ind)
                        res['Total_sample_size'].append(samp_size)
                        sub = sub_df.drop(['year', 'doy', 'row', 'col', dropper], axis = 1)
                        sub.reset_index()
            
                        x_Train, x_Test, y_Train, y_Test = train_test_split(
                                    sub.iloc[:, np.where((sub.columns.values == metric) == False)[0]],
                                    sub[metric], random_state = iteri, test_size = test_ratio)
                        
                        res['Train_size'].append(x_Train.shape[0])
                        res['Test_size'].append(x_Test.shape[0])

                        stor = f'/data/Aldhani/eoagritwin/et/Training_ML/models/Iteration_{str(iteri)}_{metric}_{str(i)}_{Month_Int_to_Str[str(j)]}_{ind}.sav'
                        ModPerfor(Model(y_Train, x_Train, stor, number_cores),
                                y_Test, x_Test, res)
                        df = pd.DataFrame(data = res)
                        df.to_csv(f'/data/Aldhani/eoagritwin/et/Training_ML/models/Iteration_{str(iteri)}_{metric}_{str(i)}_{Month_Int_to_Str[str(j)]}_{ind}.csv', sep=',', index=False)
                        print(iteri)


# ##### run gbr 
if __name__ == '__main__':
    starttime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    print("--------------------------------------------------------")
    print("Starting process, time: " + starttime)
    print("")

    # run model and store performances in results-container
    ModelRun()

    print("")
    endtime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    print("--------------------------------------------------------")
    print("--------------------------------------------------------")
    print("start: " + starttime)
    print("end: " + endtime)
    print("")

# df = pd.DataFrame(data = res)
# df.to_csv(f'/data/Aldhani/eoagritwin/et/Training_ML/output/Per_year_per_month.csv', sep=',', index=False)