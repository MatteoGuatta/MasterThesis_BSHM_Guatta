"""
train.py
===============================================================================
Master Thesis in Bridge Structural Health Monitoring (BSHM)
Pipeline Machine Learning - Step 3: Addestramento ConvAutoencoder2D
===============================================================================

Modulo responsabile di:
1. Caricamento dei DataLoader (Train e Validation sano).
2. Istanziazione del modello ConvAutoencoder2D e spostamento su GPU/CPU.
3. Definizione della funzione di costo (MSE) e dell'ottimizzatore (Adam).
4. Esecuzione del loop di training con Early Stopping basato sulla val_loss.
5. Salvataggio del checkpoint dei pesi migliori (best_model.pth).
6. Salvataggio del grafico delle curve di apprendimento (loss_history.png).
"""

from pathlib import Path
import sys
import time
import matplotlib.pyplot as plt
import numpy as np
import torch
import torch.nn as nn
from torch.optim import Adam
from torch.optim.lr_scheduler import ReduceLROnPlateau

# Aggiunta della directory corrente al path di sistema
SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

# Import standard dai moduli locali del progetto
from dataset import get_ae_dataloaders
from models import ConvAutoencoder2D


# =============================================================================
# 1. FUNZIONE DI ADDESTRAMENTO PER SINGOLA EPOCA
# =============================================================================
def train_epoch(
    model: nn.Module,
    loader: torch.utils.data.DataLoader,
    criterion: nn.Module,
    optimizer: torch.optim.Optimizer,
    device: torch.device,
) -> float:
    """
    Esegue una singola epoca di addestramento su tutti i batch del train set.
    """
    model.train()  # Attiva la modalita training (attiva batchnorm e dropout se presenti)
    running_loss = 0.0

    for x, _ in loader:
        # Input dal dataset: shape [B, 9, 200]
        # Adattamento per Conv2d: shape [B, 1, 9, 200]
        x = x.unsqueeze(1).to(device)

        # 1. Azzera i gradienti residui
        optimizer.zero_grad()

        # 2. Forward pass (ricostruzione)
        x_recon = model(x)

        # 3. Calcolo della loss di ricostruzione (MSE)
        loss = criterion(x_recon, x)

        # 4. Backward pass (calcolo gradienti)
        loss.backward()

        # 5. Aggiornamento dei pesi
        optimizer.step()

        running_loss += loss.item() * x.size(0)

    # Restituisce l'MSE medio per campione
    return running_loss / len(loader.dataset)


# =============================================================================
# 2. FUNZIONE DI VALIDAZIONE
# =============================================================================
def validate_epoch(
    model: nn.Module,
    loader: torch.utils.data.DataLoader,
    criterion: nn.Module,
    device: torch.device,
) -> float:
    """
    Calcola la loss sui dati di validazione senza calcolare i gradienti.
    """
    model.eval()  # Attiva la modalita valutazione
    running_loss = 0.0

    with torch.no_grad():  # Disabilita autograd per risparmiare memoria e tempo
        for x, _ in loader:
            x = x.unsqueeze(1).to(device)
            x_recon = model(x)
            loss = criterion(x_recon, x)
            running_loss += loss.item() * x.size(0)

    return running_loss / len(loader.dataset)


def safe_save_checkpoint(payload: dict, target_path: Path, max_retries: int = 5, delay: float = 0.2):
    """
    Salvataggio atomico e robusto per Windows:
    Evita il crash "RuntimeError: File cannot be opened" causato dal blocco temporaneo
    del file .pth da parte dell'antivirus (Windows Defender) o dell'indicizzatore.
    """
    target_path = Path(target_path)
    temp_path = target_path.with_suffix(".tmp")
    
    # Salva prima su file temporaneo
    torch.save(payload, str(temp_path))
    
    # Rinomina atomica con tentativi ripetuti in caso di lock momentaneo
    for attempt in range(max_retries):
        try:
            temp_path.replace(target_path)
            return
        except PermissionError:
            time.sleep(delay)
    
    # Se i retry non bastano, prova una sostituzione diretta
    temp_path.replace(target_path)


# =============================================================================
# 3. PIPELINE COMPLETA DI ADDESTRAMENTO
# =============================================================================
def train_pipeline(
    latent_dim: int = 32,
    batch_size: int = 64,
    window_size: int = 200,
    stride: int = 100,
    lr: float = 1e-3,
    weight_decay: float = 1e-5,
    epochs: int = 80,
    patience: int = 10,
    checkpoint_name: str = "best_model.pth",
    plot_name: str = "loss_history.png",
):
    print("=" * 70)
    print("AVVIO ADDESTRAMENTO: ConvAutoencoder2D (Bridge SHM)")
    print("=" * 70)

    # 1. Selezione del Device
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"[Device] Dispositivo in uso: {device}")

    # 2. Caricamento dati tramite la factory
    print("\n[Data] Caricamento e preparazione dei dati (dataset.py)...")
    loaders = get_ae_dataloaders(
        batch_size=batch_size,
        window_size=window_size,
        stride=stride,
    )
    train_loader = loaders["train_loader"]
    val_loader = loaders["val_loader"]

    print(f"  - Train finestre:      {len(train_loader.dataset)} ({len(train_loader)} batch)")
    print(f"  - Validation finestre: {len(val_loader.dataset)} ({len(val_loader)} batch)")

    # 3. Inizializzazione del Modello
    model = ConvAutoencoder2D(latent_dim=latent_dim).to(device)
    total_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
    print(f"\n[Model] ConvAutoencoder2D inizializzato (Spazio Latente Z: {latent_dim})")
    print(f"  - Parametri allenabili: {total_params:,}")

    # 4. Definizione Loss, Ottimizzatore e Scheduler
    criterion = nn.MSELoss()
    optimizer = Adam(model.parameters(), lr=lr, weight_decay=weight_decay)
    scheduler = ReduceLROnPlateau(optimizer, mode="min", factor=0.5, patience=4)

    # 5. Variabili di controllo per Early Stopping
    best_val_loss = float("inf")
    patience_counter = 0
    checkpoint_path = SCRIPT_DIR / checkpoint_name
    plot_path = SCRIPT_DIR / plot_name

    train_history = []
    val_history = []

    start_time = time.time()
    print(f"\n[Training] Avvio ciclo su {epochs} epoche (Patience: {patience})...\n")

    for epoch in range(1, epochs + 1):
        t_epoch_start = time.time()

        train_loss = train_epoch(model, train_loader, criterion, optimizer, device)
        val_loss = validate_epoch(model, val_loader, criterion, device)

        # Aggiornamento learning rate basato sulla val_loss
        scheduler.step(val_loss)
        current_lr = optimizer.param_groups[0]["lr"]

        train_history.append(train_loss)
        val_history.append(val_loss)

        t_epoch = time.time() - t_epoch_start

        # Verifica miglioramento per Early Stopping e Checkpoint
        is_best = val_loss < best_val_loss
        if is_best:
            best_val_loss = val_loss
            patience_counter = 0
            safe_save_checkpoint(
                {
                    "epoch": epoch,
                    "model_state_dict": model.state_dict(),
                    "optimizer_state_dict": optimizer.state_dict(),
                    "val_loss": val_loss,
                    "latent_dim": latent_dim,
                },
                checkpoint_path,
            )
            marker = " -> [MIGLIORE SALVATO]"
        else:
            patience_counter += 1
            marker = f" (patience: {patience_counter}/{patience})"

        print(
            f"Epoca [{epoch:02d}/{epochs:02d}] "
            f"Train MSE: {train_loss:.6f} | "
            f"Val MSE: {val_loss:.6f} | "
            f"LR: {current_lr:.1e} | "
            f"Tempo: {t_epoch:.2f}s{marker}"
        )

        if patience_counter >= patience:
            print(f"\n[Early Stopping] Nessun miglioramento per {patience} epoche consecutive.")
            print(f"Addestramento interrotto all'epoca {epoch}.")
            break

    total_time = time.time() - start_time
    print("=" * 70)
    print(f"Addestramento completato in {total_time:.2f} s")
    print(f"Miglior Validation MSE: {best_val_loss:.6f}")
    print(f"Pesi salvati in:        {checkpoint_path}")
    print("=" * 70)

    # 6. Generazione e salvataggio del grafico delle curve di apprendimento
    try:
        plt.figure(figsize=(9, 5))
        epochs_range = range(1, len(train_history) + 1)
        plt.plot(epochs_range, train_history, "b-o", markersize=3, label="Train Loss (MSE)")
        plt.plot(epochs_range, val_history, "r-s", markersize=3, label="Validation Loss (MSE)")

        # Evidenzia il punto di minimo della validation loss
        best_epoch = np.argmin(val_history) + 1
        plt.axvline(best_epoch, color="g", linestyle="--", alpha=0.7, label=f"Best Epoch ({best_epoch})")

        plt.xlabel("Epoca", fontsize=11)
        plt.ylabel("Reconstruction Error (MSE)", fontsize=11)
        plt.title("Curva di Apprendimento - ConvAutoencoder2D (Bridge Baseline Sana)", fontsize=12)
        plt.yscale("log")
        plt.grid(True, which="both", linestyle="--", alpha=0.5)
        plt.legend(fontsize=10)
        plt.tight_layout()
        plt.savefig(plot_path, dpi=300)
        plt.close()
        print(f"Grafico salvato con successo: {plot_path}")
    except Exception as e:
        print(f"[Warning] Impossibile generare il grafico: {e}")

    return {
        "model": model,
        "train_history": train_history,
        "val_history": val_history,
        "best_val_loss": best_val_loss,
    }


# =============================================================================
# ENTRY POINT
# =============================================================================
if __name__ == "__main__":
    train_pipeline(
        latent_dim=32,
        batch_size=64,
        window_size=200,
        stride=100,
        lr=1e-3,
        weight_decay=1e-5,
        epochs=80,
        patience=10,
    )
