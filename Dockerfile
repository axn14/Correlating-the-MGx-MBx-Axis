FROM python:3.12-slim-bookworm

# System dependencies: build tools (scikit-bio C extensions) + R (mofapy2)
RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc g++ build-essential \
    r-base r-base-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Install Python dependencies in layers (cache-friendly)
RUN pip install --no-cache-dir \
    numpy>=1.26 pandas>=2.1 scipy>=1.11 statsmodels>=0.14

RUN pip install --no-cache-dir \
    scikit-learn>=1.4 xgboost>=2.0 lightgbm>=4.1

RUN pip install --no-cache-dir \
    shap>=0.44 networkx>=3.2 biopython>=1.83 scikit-bio>=0.6.0

RUN pip install --no-cache-dir \
    mofapy2==0.7.0 tqdm>=4.66 joblib>=1.3

RUN pip install --no-cache-dir \
    streamlit>=1.35 plotly>=5.18 Pillow>=10.0 matplotlib>=3.8 seaborn>=0.13

# Copy repository
COPY . .

# Create directories for user data and results
RUN mkdir -p /app/user_data /app/results/intermediate /app/results/figures /app/results/tables

EXPOSE 8501

HEALTHCHECK CMD curl --fail http://localhost:8501/_stcore/health || exit 1

ENTRYPOINT ["streamlit", "run", "streamlit_app/app.py", \
            "--server.port=8501", "--server.address=0.0.0.0", \
            "--server.headless=true"]
