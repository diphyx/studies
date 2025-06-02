import cuml
import cupy as cp
from cuml.datasets.classification import make_classification
from cuml.ensemble import RandomForestClassifier
from cuml.metrics import accuracy_score, roc_auc_score
from sklearn.model_selection import train_test_split

def gpu_predict_molecular_activity():
    # Generate synthetic dataset for molecular activity (on GPU)
    X, y = make_classification(n_samples=10000, n_features=10, 
                               n_informative=6, n_classes=2, random_state=42)
    
    # Ensure GPU arrays
    X = cp.array(X)
    y = cp.array(y)

    # Split using sklearn but convert back to cupy for GPU
    X_train, X_test, y_train, y_test = train_test_split(X.get(), y.get(), test_size=0.2, random_state=42)
    X_train, X_test = cp.array(X_train), cp.array(X_test)
    y_train, y_test = cp.array(y_train), cp.array(y_test)

    # Train a cuML Random Forest
    clf = RandomForestClassifier(n_estimators=100, max_depth=10)
    clf.fit(X_train, y_train)

    # Predict
    y_pred = clf.predict(X_test)
    y_proba = clf.predict_proba(X_test)[:, 1]

    # Evaluate
    acc = accuracy_score(y_test, y_pred)
    roc = roc_auc_score(y_test, y_proba)

    print(f"Accuracy: {acc:.2f}")
    print(f"ROC AUC Score: {roc:.2f}")

if __name__ == "__main__":
    gpu_predict_molecular_activity()