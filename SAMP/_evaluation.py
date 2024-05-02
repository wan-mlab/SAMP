## This file is not a part of the original codebase, but it is a useful addition to the codebase.
from sklearn import metrics
def eval_metrics(y_true,y_pred,prob):
    print('Accuracy:', round(metrics.accuracy_score(y_true,y_pred),3))
    print('Precision:', round(metrics.precision_score(y_true, y_pred),3))
    print('Sensitivity:', round(metrics.recall_score(y_true,y_pred),3))
    temp_tn, temp_fp, temp_fn, temp_tp = metrics.confusion_matrix(y_true,y_pred).ravel()
    temp_specificity = temp_tn / (temp_tn + temp_fp)
    print('Specificity:', round(temp_specificity,3))
    print('F1-score:',round(metrics.f1_score(y_true,y_pred),3))
    print('AUC: ', round(metrics.roc_auc_score(y_true, prob),3))
    print('MCC: ', round(metrics.matthews_corrcoef(y_true,y_pred),3))