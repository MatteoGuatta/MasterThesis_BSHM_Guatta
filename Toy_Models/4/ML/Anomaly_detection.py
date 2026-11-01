"""
Anomaly_detection.py
===============================================================================
Master Thesis in Bridge Structural Health Monitoring (BSHM)
Pipeline Machine Learning - Step 4: Anomaly Detection & Damage Localization
===============================================================================

Orchestratore principale che:
1. Carica il modello addestrato (best_model.pth).
2. Carica i dati di test (Sani e Danneggiati).
3. Calcola i Damage Indices (globali e locali per sensore) tramite metrics.py.
4. Fissa la soglia di anomalia statistica (Dual-Threshold Strategy di Giorgi et al.).
5. Genera e salva tutti i grafici di Rilevamento e Localizzazione in figures/ (plots.py).
"""

from pathlib import Path
import sys
import numpy as np
import torch

# Gestione percorsi per importare in modo pulito e affidabile
SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from dataset import get_ae_dataloaders
from metrics import compute_reconstruction_errors, compute_thresholds
from models import load_trained_model
import plots

# Percorsi di riferimento
CHECKPOINT_PATH = SCRIPT_DIR / "best_model.pth"
FIGURES_DIR = SCRIPT_DIR / "figures"
FIGURES_DIR.mkdir(exist_ok=True)

# Selezione Device
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")


def run_anomaly_detection():
    print("=" * 75)
    print("AVVIO ANOMALY DETECTION & DAMAGE LOCALIZATION (Bridge SHM)")
    print("=" * 75)
    print(f"[Device] In uso: {device}")

    # 1. Caricamento del Modello Addestrato
    print("\n1. Caricamento modello da checkpoint...")
    model = load_trained_model(checkpoint_path=CHECKPOINT_PATH, device=device)

    # 2. Caricamento dei Dataloader di Test
    print("\n2. Caricamento DataLoader di test (dataset.py)...")
    loaders = get_ae_dataloaders()
    test_healthy_loader = loaders["test_healthy_loader"]
    test_damaged_loader = loaders["test_damaged_loader"]

    print(f"  - Test Sano finestre:        {len(test_healthy_loader.dataset)}")
    print(f"  - Test Danneggiato finestre: {len(test_damaged_loader.dataset)}")

    # 3. Calcolo errori sul set Sano e definizione delle Soglie
    print("\n3. Calcolo errori sul Test Set Sano...")
    di_sensors_h, di_global_h, meta_h = compute_reconstruction_errors(model, test_healthy_loader, device)

    # Calcolo Dual-Threshold (3-sigma)
    tau_sensors, tau_global = compute_thresholds(di_sensors_h, di_global_h, method="3sigma")

    print(f"  [Soglia Globale]  tau_global = {tau_global.item():.6f}")
  # print(f"  [Soglie Locali]   tau_sensors (9 canali) = {np.round(tau_sensors.numpy(), 5)}")

    # 4. Calcolo errori sul set Danneggiato
    print("\n4. Calcolo errori sul Test Set Danneggiato...")
    di_sensors_d, di_global_d, meta_d = compute_reconstruction_errors(model, test_damaged_loader, device)

    # 5. Valutazione Metriche di Rilevamento (Detection)
    alarms_h = di_global_h > tau_global
    fpr = alarms_h.float().mean().item() * 100.0  # False Positive Rate

    alarms_d = di_global_d > tau_global
    tpr = alarms_d.float().mean().item() * 100.0  # True Positive Rate (Sensitivity)

    print("\n" + "=" * 75)
    print("RISULTATI DI RILEVAMENTO (DETECTION PERFORMANCE)")
    print("=" * 75)
    print(f"  - Tasso Falsi Allarmi su Test Sano (FPR):    {fpr:.2f}% ({alarms_h.sum().item()}/{len(di_global_h)})")
    print(f"  - Tasso Rilevamento Danno Globale (TPR):    {tpr:.2f}% ({alarms_d.sum().item()}/{len(di_global_d)})")

    # Dettaglio per severita
    severities = np.unique(meta_d[:, 6].numpy())
    severities = severities[severities > 0]
    print("\n  Dettaglio Detection Rate per Severita di Danno:")
    for s in severities:
        mask = meta_d[:, 6].numpy() == s
        tpr_sev = np.mean(di_global_d.numpy()[mask] > tau_global.item()) * 100.0
        print(f"    * Severita -{int(s*100)}%:  {tpr_sev:.2f}% di finestre anomale rilevate")

    # 6. Generazione di tutti i Grafici (plots.py)
    print("\n" + "=" * 75)
    print("GENERAZIONE GRAFICI PER LA TESI (figures/)")
    print("=" * 75)

    # Grafico 1: Istogramma Sano vs Danneggiato con soglia tau
    plots.plot_detection_histogram(
        di_healthy=di_global_h.numpy(),
        di_damaged=di_global_d.numpy(),
        tau_global=tau_global.item(),
    )

    # Grafico 2: Serie temporale DI(t) per singola corsa (Run-Level)
    # Estraiamo la prima corsa sana (29 finestre) e una corsa danneggiata a severita 25% o 50%
    run_id_h = meta_h[0, 0].item()
    mask_run_h = (meta_h[:, 0] == run_id_h).numpy()
    healthy_run_di = di_global_h.numpy()[mask_run_h]

    # Corsa danneggiata rappresentativa (es. severita 25%)
    mask_run_d = (meta_d[:, 6] == 0.25).numpy()
    if np.any(mask_run_d):
        run_id_d = meta_d[mask_run_d][0, 0].item()
        mask_single_d = (meta_d[:, 0] == run_id_d).numpy()
        damaged_run_di = di_global_d.numpy()[mask_single_d]
    else:
        damaged_run_di = di_global_d.numpy()[:len(healthy_run_di)]

    plots.plot_transit_time_history(
        healthy_run_di=healthy_run_di,
        damaged_run_di=damaged_run_di,
        tau_global=tau_global.item(),
    )

    # Grafico 3: Detection Rate per Severita (Bar Chart)
    plots.plot_detection_rate_by_severity(
        di_damaged=di_global_d.numpy(),
        meta_damaged=meta_d.numpy(),
        tau_global=tau_global.item(),
    )

    # Grafico 4: Curva ROC e AUC
    plots.plot_roc_curve_metric(
        di_healthy=di_global_h.numpy(),
        di_damaged=di_global_d.numpy(),
    )

    # Grafico 5: Profilo di Localizzazione lungo i 9 Sensori
    plots.plot_damage_localization_profile(
        di_sensors_healthy=di_sensors_h.numpy(),
        di_sensors_damaged=di_sensors_d.numpy(),
        tau_sensors=tau_sensors.numpy(),
        damage_x_target=7.5,  # Elemento 8 del FEM -> x = 7.5 m
    )

    # Grafico 6: Mappa Spazio-Tempo (Heatmap)
    # Prendiamo una finestra danneggiata dal dataloader
    with torch.no_grad():
        x_sample, _ = next(iter(test_damaged_loader))
        x_in = x_sample[0:1].unsqueeze(1).to(device)
        x_out = model(x_in)
        plots.plot_space_time_heatmap(
            x_real=x_sample[0].numpy(),
            x_recon=x_out.squeeze(0).squeeze(0).cpu().numpy(),
            window_id=int(meta_d[0, 7].item()),
        )

    print("=" * 75)
    print("Tutti i grafici sono stati generati con successo in:")
    print(f"  {FIGURES_DIR.resolve()}")
    print("=" * 75)


if __name__ == "__main__":
    run_anomaly_detection()
