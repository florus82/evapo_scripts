from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import train_test_split
from sklearn.ensemble import GradientBoostingRegressor
import joblib
import time
from sklearn import metrics


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
def ModPerfor(cv_results, yData, xData):
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

# make results container
# build result container
keys = ['Iteration','Input', 'DbH_Set', 'no_obs','r2', 'mse', 'rmse', 'y_true', 'y_pred', 'max_depth', 'learning_rate', 'min_samples_leaf', 'max_features']
vals = [list() for i in range(len(keys))]
res  = dict(zip(keys, vals))

# access test files created with Prepare_S1_points.R
files = getFilelist('K:/Seafile/Uni_Life/CarbonPaper/2nd_round/Modelling/S1_predictor_tests/input_files_for_python/with_Dante_old_correct', 'csv')

# before dante; with dbh differentiation
# def ModelRun():
#     # iterate over input files
#     for iteri in range(0,100):
#         for file in files:
#             dat = pd.read_csv(file)
#             y_mod = dat.DbH_year.unique()
#             for setti in y_mod:
#                 sub = dat[dat['DbH_year'] == setti].drop(['Point_ID', 'DbH_year'], axis = 1)
#                 sub.reset_index()
#                 res['no_obs'].append(sub.shape[0])
#                 x_Train, x_Test, y_Train, y_Test = train_test_split(
#                             sub.iloc[:, np.where((sub.columns.values == 'AGB') == False)[0]],
#                             sub['AGB'], random_state = iteri, test_size = 0.3)
#                 res['Iteration'].append(str(iteri))
#                 res['Input'].append(file.split('/')[-1].split('.')[0])
#                 res['DbH_Set'].append(setti)
#
#
#                 stor = 'K:/Seafile/Uni_Life/CarbonPaper/2nd_round/Modelling/S1_predictor_tests/savs/' + \
#                        file.split('/')[-1].split('.')[0] + '_' + setti + '_' + str(iteri) + '.sav'
#                 ModPerfor(Model(y_Train, x_Train, stor, 54),
#                           y_Test, x_Test)
#                 print(setti)
#             print(file)
#         print(iteri)

def ModelRun():
    # iterate over input files
    for iteri in range(0,100):
        for file in files:
            dat = pd.read_csv(file)
            sub = dat.drop(['POINT_ID', 'YEAR', 'DbH'], axis = 1)
            sub.reset_index()
            res['no_obs'].append(sub.shape[0])
            x_Train, x_Test, y_Train, y_Test = train_test_split(
                        sub.iloc[:, np.where((sub.columns.values == 'AGB') == False)[0]],
                        sub['AGB'], random_state = iteri, test_size = 0.3)
            res['Iteration'].append(str(iteri))
            res['Input'].append(file.split('/')[-1].split('.')[0])
            res['DbH_Set'].append(file.split('_')[-1].split('.')[0])

            stor = 'K:/MSc_outside_Seafile\savs_2nd_round\withDante_old_BM_correct/' + file.split('/')[-1].split('.')[0] + str(iteri) + '.sav'
            ModPerfor(Model(y_Train, x_Train, stor, 15),
                      y_Test, x_Test)
            print(file)
        print(iteri)


       # ##### run gbr once
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

df = pd.DataFrame(data = res)
df.to_csv('K:/Seafile/Uni_Life/CarbonPaper/2nd_round/Modelling/S1_predictor_tests/S1_combis_100_iter_Dante_and_Old_BM_correct.csv', sep=',', index=False)