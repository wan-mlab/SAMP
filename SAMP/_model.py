import numpy as np
import pandas as pd
from sklearn.random_projection import GaussianRandomProjection, SparseRandomProjection
from sklearn.model_selection import KFold, cross_val_predict, GridSearchCV
from sklearn.svm import SVC
from sklearn import metrics

class SAMP:
    def __init__(self, dim=50, method='Gaussian', runs=10):
        self.models = []
        self.dim = dim
        self.method = method
        self.runs = runs

    def calculate_metrics(self, y_true, y_pred, scores):
        tn, fp, fn, tp = metrics.confusion_matrix(y_true, y_pred).ravel()
        return {
            "accuracy": metrics.accuracy_score(y_true, y_pred),
            "sensitivity": metrics.recall_score(y_true, y_pred),  # Sensitivity is recall
            "specificity": tn / (tn + fp),
            "precision": metrics.precision_score(y_true, y_pred),
            "auc": metrics.roc_auc_score(y_true, scores),
            "f1": metrics.f1_score(y_true, y_pred),
            "mcc": metrics.matthews_corrcoef(y_true, y_pred)
        }

    def score_ensemble(self, scores, y_tr):
        # Average the score across all runs of RP
        mean_scores = np.mean(scores, axis=1)
        y_pred = (mean_scores > 0).astype(int)
        
        # Evaluate metrics for average score
        ensemble_metrics = self.calculate_metrics(y_tr, y_pred, mean_scores)

        # Calculate metrics for each run
        run_metrics = {metric: [] for metric in ensemble_metrics}
        for i in range(self.runs):
            y_pred_per_run = (scores[f'Run_{i}'] > 0).astype(int)
            metrics_per_run = self.calculate_metrics(y_tr, y_pred_per_run, scores[f'Run_{i}'])
            for metric in run_metrics:
                run_metrics[metric].append(metrics_per_run[metric])

        # Calculate standard deviations
        std_metrics = {metric: np.std(values) for metric, values in run_metrics.items()}

        # Print results
        print(f"Dimension {self.dim} =>>")
        for metric in ensemble_metrics:
            print(f"{metric.capitalize()} : {ensemble_metrics[metric] * 100:.2f} ± {std_metrics[metric]:.3f}")


    def random_projection(self, X, random_state):
        if self.method == 'Gaussian':
            transformer = GaussianRandomProjection(n_components=self.dim, random_state=random_state)
        elif self.method == 'Sparse':
            transformer = SparseRandomProjection(n_components=self.dim, random_state=random_state)
        else:
            raise ValueError("Method should be 'Gaussian' or 'Sparse'")
        return transformer.fit_transform(X)

    def train_svm(self, X, y, param_grid, cv_splits=10, random_state=1):
        print("Starting SVM training using GridSearchCV...")
        svm_model = SVC(random_state=random_state)
        grid_search = GridSearchCV(svm_model, param_grid, cv=cv_splits, n_jobs=-1, scoring='accuracy')
        grid_search.fit(X, y)
        print("Grid search completed. Best parameters found.")
        best_model = SVC(**grid_search.best_params_, random_state=random_state)
        best_model.fit(X, y)
        return best_model

    def fit(self, X, y):
        scores = pd.DataFrame()
        param_grid = {'C': [0.1, 1, 10, 100], 'gamma': [1, 0.1, 0.01, 0.001, 0.0001], 'kernel': ['rbf']}
        for i in range(self.runs):
            print(f"Starting run {i+1}/{self.runs}...")  # Progress of runs
            X_transformed = self.random_projection(X, i)
            print(f"Random projection {i+1} completed.")
            # Train model and store it
            model = self.train_svm(X_transformed, y, param_grid)
            self.models.append(model)
            
            # Get scores for each run
            score = cross_val_predict(model, X_transformed, y, cv=KFold(n_splits=10, shuffle=True, random_state=i), method='decision_function')
            scores[f'Run_{i}'] = score

        self.score_ensemble(scores, y)
        print("All runs completed. Ensemble scoring finished.")

    def predict(self, X, y):
        scores = pd.DataFrame()
        # Use the trained models to predict
        for i, model in enumerate(self.models):
            X_test_transformed = self.random_projection(X, i)
            score = model.decision_function(X_test_transformed)
            scores[f'Run_{i}'] = score

        self.score_ensemble(scores, y)
