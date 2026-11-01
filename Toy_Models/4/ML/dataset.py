"""
1_dataset.py
===============================================================================
Pipeline Machine Learning - Step 1: Data Ingestion, Windowing & PyTorch Dataset
===============================================================================

Modulo responsabile di:
1. Caricamento file .mat HDF5 generati da MATLAB (sani e danneggiati).
2. Data splitting a livello di CORSE (Run-level) per evitare data leakage:
   Train (70%), Validation (15%), Test Sano (15%).
3. Normalizzazione Sensor-wise (StandardScaler: media e std stimate SOLO sul train sano).
4. Windowing temporale (finestre scorrevoli con overlap su serie a 200 Hz).
5. Costruzione della classe PyTorch Dataset e DataLoader per Autoencoder (AE).
6. Definizione del grafo e coordinate dei 9 sensori (preparazione per GNN).
"""

from pathlib import Path
from typing import Dict, Optional, Tuple, Union
import h5py
import numpy as np
import torch
from torch.utils.data import DataLoader, Dataset

# Percorsi di default del progetto
SCRIPT_DIR = Path(__file__).resolve().parent
DATA_DIR = SCRIPT_DIR.parent / "data"

FILE_HEALTHY = DATA_DIR / "raw_healthy_runs.mat"
FILE_DAMAGED = DATA_DIR / "raw_damaged_runs_1.mat"


# =============================================================================
# 1. CARICAMENTO DATI RAW MAT
# =============================================================================
def load_raw_mat(path_file: Union[str, Path]) -> Tuple[np.ndarray, np.ndarray]:
    """
    Carica i dati accelerometrici e i metadati da file MATLAB v7.3 (.mat / HDF5).
    
    Output:
        acc: np.ndarray con shape (N_runs, N_sensors=9, N_time_steps=3000)
        meta: np.ndarray con shape (N_runs, 7)
              Colonne: [RunID, Load_N, Speed_kmh, Speed_ms, IsDamaged, Scenario/Element, Severity]
    """
    path_file = Path(path_file)
    if not path_file.exists():
        raise FileNotFoundError(f"File non trovato: {path_file}")

    with h5py.File(path_file, "r") as f:
        if "acc_healthy" in f and "meta_healthy" in f:
            # In MATLAB v7.3 HDF5 l'ordine delle dimensioni è invertito: (3000, 9, N)
            acc = np.array(f["acc_healthy"]).transpose(2, 1, 0)  # -> (N, 9, 3000)
            meta = np.array(f["meta_healthy"]).T                  # -> (N, 7)
        elif "acc_damaged" in f and "meta_damaged" in f:
            acc = np.array(f["acc_damaged"]).transpose(2, 1, 0)
            meta = np.array(f["meta_damaged"]).T
        else:
            raise KeyError(f"Struttura non riconosciuta all'interno di {path_file.name}")

    return acc.astype(np.float32), meta.astype(np.float32)


# =============================================================================
# 2. DATA SPLITTING A LIVELLO DI CORSA (RUN-LEVEL)
# =============================================================================
def train_val_test_split_runs(
    acc: np.ndarray,
    meta: np.ndarray,
    val_size: float = 0.15,
    test_size: float = 0.15,
    random_state: int = 42,
) -> Tuple[Tuple[np.ndarray, np.ndarray], Tuple[np.ndarray, np.ndarray], Tuple[np.ndarray, np.ndarray]]:
    """
    Esegue la suddivisione casuale a livello di CORSA intera per evitare data leakage
    temporale tra le finestre di addestramento e di test.
    """
    np.random.seed(random_state)
    n_runs = acc.shape[0]
    indices = np.arange(n_runs)
    np.random.shuffle(indices)

    n_test = int(n_runs * test_size)
    n_val = int(n_runs * val_size)
    n_train = n_runs - n_val - n_test

    train_idx = indices[:n_train]
    val_idx = indices[n_train : n_train + n_val]
    test_idx = indices[n_train + n_val :]

    train_set = (acc[train_idx], meta[train_idx])
    val_set = (acc[val_idx], meta[val_idx])
    test_set = (acc[test_idx], meta[test_idx])

    return train_set, val_set, test_set


# =============================================================================
# 3. NORMALIZZAZIONE SENSOR-WISE (STANDARD SCALER)
# =============================================================================
class SensorStandardScaler:
    """
    StandardScaler per canali accelerometrici (9 sensori).
    Media e deviazione standard vengono stimate ESCLUSIVAMENTE sui dati di train sani
    e applicate a tutte le partizioni (train, val, test, damaged).
    """

    def __init__(self, eps: float = 1e-8):
        self.eps = eps
        self.mean: Optional[np.ndarray] = None  # shape (1, N_sensors, 1)
        self.std: Optional[np.ndarray] = None   # shape (1, N_sensors, 1)

    def fit(self, x: np.ndarray) -> "SensorStandardScaler":
        """
        Calcola media e std su asse 0 (runs) e asse 2 (campioni temporali) per ciascun sensore.
        x ha shape (N_runs, N_sensors, N_pts).
        """
        self.mean = np.mean(x, axis=(0, 2), keepdims=True)
        self.std = np.std(x, axis=(0, 2), keepdims=True)
        # Evita divisione per zero nel caso di varianza nulla
        self.std = np.where(self.std < self.eps, 1.0, self.std)
        return self

    def transform(self, x: np.ndarray) -> np.ndarray:
        if self.mean is None or self.std is None:
            raise RuntimeError("Lo scaler deve essere addestrato con fit() prima di trasformare i dati.")
        return (x - self.mean) / self.std

    def fit_transform(self, x: np.ndarray) -> np.ndarray:
        return self.fit(x).transform(x)

    def inverse_transform(self, x: np.ndarray) -> np.ndarray:
        if self.mean is None or self.std is None:
            raise RuntimeError("Lo scaler deve essere addestrato con fit() prima di inverse_transform.")
        return (x * self.std) + self.mean


# =============================================================================
# 4. WINDOWING TEMPORALE (FINESTRE SCORREVOLI)
# =============================================================================
def create_windows(
    acc: np.ndarray,
    meta: np.ndarray,
    window_size: int = 200,
    stride: int = 100,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Estrae finestre temporali scorrevoli con sovrapposizione dalle corse.
    
    Parametri:
        acc: array shape (N_runs, N_sensors, N_pts)
        meta: array shape (N_runs, 7)
        window_size: lunghezza della finestra (es. 200 punti = 1.0 s @ 200 Hz)
        stride: passo di avanzamento (es. 100 punti = 50% overlap)
        
    Output:
        windows_acc: array shape (N_windows, N_sensors, window_size)
        windows_meta: array shape (N_windows, 8)
                      Colonne: [RunID, Load_N, Speed_kmh, Speed_ms, IsDamaged, Scenario/Element, Severity, WindowID]
    """
    n_runs, n_sensors, n_pts = acc.shape
    if window_size > n_pts:
        raise ValueError(f"window_size ({window_size}) non può essere maggiore della durata della corsa ({n_pts}).")

    windows_acc_list = []
    windows_meta_list = []

    start_indices = list(range(0, n_pts - window_size + 1, stride))

    for run_i in range(n_runs):
        run_acc = acc[run_i]      # shape (N_sensors, N_pts)
        run_meta = meta[run_i]    # shape (7,)

        for win_id, start in enumerate(start_indices):
            end = start + window_size
            window = run_acc[:, start:end]  # shape (N_sensors, window_size)
            windows_acc_list.append(window)

            # Metadata esteso con ID finestra
            extended_meta = np.append(run_meta, win_id)
            windows_meta_list.append(extended_meta)

    windows_acc = np.stack(windows_acc_list, axis=0).astype(np.float32)
    windows_meta = np.stack(windows_meta_list, axis=0).astype(np.float32)

    return windows_acc, windows_meta


# =============================================================================
# 5. CLASSE PYTORCH DATASET PER AUTOENCODER
# =============================================================================
class BridgeAEDataset(Dataset):
    """
    Dataset PyTorch per Autoencoder su dati di accelerazione multi-canale del ponte.
    
    Ogni campione restituito da __getitem__ è una tupla:
        (x, meta)
    dove x è un torch.FloatTensor con shape (N_sensors=9, window_size).
    Per un Autoencoder, il target di ricostruzione è x stesso.
    """

    def __init__(self, windows: np.ndarray, meta: Optional[np.ndarray] = None):
        super().__init__()
        self.windows = torch.from_numpy(windows).float()
        if meta is not None:
            self.meta = torch.from_numpy(meta).float()
        else:
            self.meta = torch.zeros((windows.shape[0], 8), dtype=torch.float32)

    def __len__(self) -> int:
        return self.windows.shape[0]

    def __getitem__(self, idx: int) -> Tuple[torch.Tensor, torch.Tensor]:
        return self.windows[idx], self.meta[idx]


# =============================================================================
# 6. TOPOLOGIA DEI SENSORI (GRAFO DEL PONTE PER FUTURE GNN)
# =============================================================================
def get_sensor_graph(sigma_dist: float = 6.0) -> Dict[str, Union[np.ndarray, torch.Tensor]]:
    """
    Restituisce coordinate fisiche, matrice di adiacenza e archi del grafo dei 9 sensori.
    I sensori sono posizionati a x = [3, 6, 9, 12, 15, 18, 21, 24, 27] m.
    
    I pesi degli archi sono calcolati tramite Kernel Gaussiano:
        W_ij = exp(- (d_ij)^2 / (2 * sigma_dist^2))
    """
    # Coordinate x dei 9 sensori virtuali (in metri)
    sensor_x = np.array([3.0, 6.0, 9.0, 12.0, 15.0, 18.0, 21.0, 24.0, 27.0], dtype=np.float32)
    n_sensors = len(sensor_x)

    # Matrice delle distanze euclidee
    dist_matrix = np.abs(sensor_x[:, None] - sensor_x[None, :])

    # Pesi gaussiani
    adj_matrix = np.exp(-(dist_matrix ** 2) / (2.0 * sigma_dist ** 2))
    np.fill_diagonal(adj_matrix, 1.0)

    # Line graph adiacenza binaria (connessione solo tra vicini immediati)
    line_adj = np.zeros((n_sensors, n_sensors), dtype=np.float32)
    for i in range(n_sensors - 1):
        line_adj[i, i + 1] = 1.0
        line_adj[i + 1, i] = 1.0
    np.fill_diagonal(line_adj, 1.0)

    return {
        "sensor_nodes_fem": np.array([4, 7, 10, 13, 16, 19, 22, 25, 28]),
        "sensor_x_m": sensor_x,
        "dist_matrix": dist_matrix,
        "adj_gaussian": torch.from_numpy(adj_matrix).float(),
        "adj_line": torch.from_numpy(line_adj).float(),
    }


# =============================================================================
# 7. FACTORY ORCHESTRATORE PER PYTORCH DATALOADERS
# =============================================================================
def get_ae_dataloaders(
    healthy_file: Union[str, Path] = FILE_HEALTHY,
    damaged_file: Union[str, Path] = FILE_DAMAGED,
    batch_size: int = 64,
    window_size: int = 200,
    stride: int = 100,
    val_size: float = 0.15,
    test_size: float = 0.15,
    random_state: int = 42,
    num_workers: int = 0,
) -> Dict[str, Union[DataLoader, SensorStandardScaler]]:
    """
    Pipeline completa per l'Autoencoder:
    1. Carica i dati grezzi.
    2. Divide i dati sani in Train, Val e Test a livello di corsa.
    3. Normalizza con StandardScaler fit SOLO su Train.
    4. Crea le finestre temporali per ciascun set.
    5. Restituisce i DataLoader PyTorch e lo scaler.
    """
    # 1. Caricamento
    acc_h, meta_h = load_raw_mat(healthy_file)
    acc_d, meta_d = load_raw_mat(damaged_file)

    # 2. Split corse sane
    (acc_train, meta_train), (acc_val, meta_val), (acc_test_h, meta_test_h) = train_val_test_split_runs(
        acc_h, meta_h, val_size=val_size, test_size=test_size, random_state=random_state
    )

    # 3. Fit scaler su train sano e trasformazione
    scaler = SensorStandardScaler()
    scaler.fit(acc_train)

    acc_train_norm = scaler.transform(acc_train)
    acc_val_norm = scaler.transform(acc_val)
    acc_test_h_norm = scaler.transform(acc_test_h)
    acc_d_norm = scaler.transform(acc_d)

    # 4. Windowing temporale
    win_train, meta_win_train = create_windows(acc_train_norm, meta_train, window_size=window_size, stride=stride)
    win_val, meta_win_val = create_windows(acc_val_norm, meta_val, window_size=window_size, stride=stride)
    win_test_h, meta_win_test_h = create_windows(acc_test_h_norm, meta_test_h, window_size=window_size, stride=stride)
    win_d, meta_win_d = create_windows(acc_d_norm, meta_d, window_size=window_size, stride=stride)

    # 5. Dataset PyTorch
    ds_train = BridgeAEDataset(win_train, meta_win_train)
    ds_val = BridgeAEDataset(win_val, meta_win_val)
    ds_test_h = BridgeAEDataset(win_test_h, meta_win_test_h)
    ds_test_d = BridgeAEDataset(win_d, meta_win_d)

    # 6. DataLoader
    loader_train = DataLoader(ds_train, batch_size=batch_size, shuffle=True, num_workers=num_workers)
    loader_val = DataLoader(ds_val, batch_size=batch_size, shuffle=False, num_workers=num_workers)
    loader_test_h = DataLoader(ds_test_h, batch_size=batch_size, shuffle=False, num_workers=num_workers)
    loader_test_d = DataLoader(ds_test_d, batch_size=batch_size, shuffle=False, num_workers=num_workers)

    return {
        "train_loader": loader_train,
        "val_loader": loader_val,
        "test_healthy_loader": loader_test_h,
        "test_damaged_loader": loader_test_d,
        "scaler": scaler,
    }


# =============================================================================
# VERIFICA E DEBUG
# =============================================================================
if __name__ == "__main__":
    print("=" * 70)
    print("Verifica Pipeline Dataset Autoencoder (AE)...")
    print("=" * 70)

    loaders = get_ae_dataloaders(
        healthy_file=FILE_HEALTHY,
        damaged_file=FILE_DAMAGED,
        batch_size=64,
        window_size=200,  # 1.0 s a 200 Hz
        stride=100,       # 50% overlap -> 0.5 s step
    )

    train_loader = loaders["train_loader"]
    val_loader = loaders["val_loader"]
    test_h_loader = loaders["test_healthy_loader"]
    test_d_loader = loaders["test_damaged_loader"]
    scaler = loaders["scaler"]

    print(f"[OK] Train set finestre:        {len(train_loader.dataset)} campioni ({len(train_loader)} batch)")
    print(f"[OK] Validation set finestre:   {len(val_loader.dataset)} campioni ({len(val_loader)} batch)")
    print(f"[OK] Test Sano finestre:        {len(test_h_loader.dataset)} campioni ({len(test_h_loader)} batch)")
    print(f"[OK] Test Danneggiato finestre: {len(test_d_loader.dataset)} campioni ({len(test_d_loader)} batch)")

    # Controllo primo batch
    x_batch, meta_batch = next(iter(train_loader))
    print(f"\nShape del Batch x:    {x_batch.shape} (Batch_size, N_sensors, Window_size)")
    print(f"Shape del Batch meta: {meta_batch.shape} (Batch_size, Metadata_cols)")
    print(f"Tipo dato x:          {x_batch.dtype}")
    print(f"Media batch:          {x_batch.mean().item():.4f} (attesa ~0)")
    print(f"Std batch:            {x_batch.std().item():.4f} (attesa ~1)")

    # Verifica Grafo
    graph = get_sensor_graph()
    print(f"\nMatrice adiacenza Gaussiana dei sensori: shape {graph['adj_gaussian'].shape}")
    print("=" * 70)
    print("Pipeline Dataset verificata con successo!")
    print("=" * 70)
