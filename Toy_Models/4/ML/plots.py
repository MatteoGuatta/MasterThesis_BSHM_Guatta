"""
plots.py
===============================================================================
Master Thesis in Bridge Structural Health Monitoring (BSHM)
Visualizzazione e Grafici di Rilevamento e Localizzazione del Danno
===============================================================================

Modulo dedicato alla generazione di grafici scientifici (ad alta risoluzione)
per la tesi magistrale:
1. RILEVAMENTO:
   - Istogramma comparativo degli errori (Sano vs Danneggiato) con soglia tau.
   - Serie temporale del transito DI(t) per singola corsa (sana vs danneggiata).
   - Detection Rate in funzione della severita del danno (10%, 25%, 50%, 75%).
   - Curva ROC (Receiver Operating Characteristic) con calcolo dell'AUC.
2. LOCALIZZAZIONE:
   - Profilo spaziale dell'errore di ricostruzione lungo i 9 sensori (x = 3..27 m)
     confrontato con la posizione reale del difetto (Elemento 8 a L/4, x ~ 7.5m).
   - Mappa di calore spazio-temporale (Heatmap 2D: Sensori vs Tempo).
"""

from pathlib import Path
import matplotlib
matplotlib.use("Agg")  # Backend headless per evitare conflitti con la GUI Windows
import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import auc, roc_curve
import torch

SCRIPT_DIR = Path(__file__).resolve().parent
FIGURES_DIR = SCRIPT_DIR / "figures"
FIGURES_DIR.mkdir(exist_ok=True)

# Coordinate reali dei 9 sensori virtuali sul ponte da 30 metri (da config_benchmark.m)
SENSOR_X = np.array([3.0, 6.0, 9.0, 12.0, 15.0, 18.0, 21.0, 24.0, 27.0])
SENSOR_LABELS = [f"S{i+1}\n({x:.0f}m)" for i, x in enumerate(SENSOR_X)]


# =============================================================================
# 1. RILEVAMENTO: ISTOGRAMMA DISTRIBUZIONE ERRORI (SANO VS DANNEGGIATO)
# =============================================================================
def plot_detection_histogram(
    di_healthy: np.ndarray,
    di_damaged: np.ndarray,
    tau_global: float,
    save_path: Path = FIGURES_DIR / "1_detection_histogram.png",
):
    """
    Istogramma della distribuzione del Damage Index globale per campioni sani e danneggiati.
    Mostra la separabilita statistica e la soglia di anomalia tau.
    """
    plt.figure(figsize=(9, 5))
    
    # Istogramma sano
    plt.hist(
        di_healthy, bins=35, density=True, alpha=0.6,
        color="royalblue", edgecolor="black", label=f"Test Sano (N={len(di_healthy)})"
    )
    
    # Istogramma danneggiato
    plt.hist(
        di_damaged, bins=45, density=True, alpha=0.6,
        color="crimson", edgecolor="black", label=f"Test Danneggiato (N={len(di_damaged)})"
    )
    
    # Linea verticale della soglia
    plt.axvline(
        tau_global, color="black", linestyle="--", linewidth=2.0,
        label=f"Soglia Globale $\\tau$ = {tau_global:.5f}"
    )

    plt.xlabel("Damage Index Globale (MSE)", fontsize=11)
    plt.ylabel("Densita di Probabilita", fontsize=11)
    plt.title("Rilevamento Anomalia: Distribuzione DI Sano vs Danneggiato", fontsize=12, fontweight="bold")
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.legend(fontsize=10, loc="upper right")
    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    plt.close()
    print(f"[Grafico salvato] {save_path.name}")


# =============================================================================
# 2. RILEVAMENTO: STORIA TEMPORALE DEL TRANSITO (RUN-LEVEL DI(t))
# =============================================================================
def plot_transit_time_history(
    healthy_run_di: np.ndarray,
    damaged_run_di: np.ndarray,
    tau_global: float,
    dt_stride: float = 0.5,
    save_path: Path = FIGURES_DIR / "2_transit_time_history.png",
):
    """
    Confronta l'evoluzione temporale del Damage Index durante i 15 secondi di una corsa:
    una corsa sana vs una corsa con danno.
    """
    plt.figure(figsize=(10, 5))
    time_axis = np.arange(len(healthy_run_di)) * dt_stride  # 29 finestre con step 0.5s -> 0..14s

    plt.plot(time_axis, healthy_run_di, "b-o", linewidth=1.8, markersize=4, label="Corsa Sana (Baseline)")
    plt.plot(time_axis, damaged_run_di, "r-s", linewidth=1.8, markersize=4, label="Corsa Danneggiata (Danno L/4)")

    # Linea di soglia
    plt.axhline(tau_global, color="black", linestyle="--", linewidth=2.0, label=f"Soglia $\\tau$ = {tau_global:.5f}")

    plt.xlabel("Tempo di transito [s]", fontsize=11)
    plt.ylabel("Damage Index Globale (MSE)", fontsize=11)
    plt.title("Evoluzione Temporale dell'Anomalia durante il Transito del Convoglio", fontsize=12, fontweight="bold")
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.legend(fontsize=10)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    plt.close()
    print(f"[Grafico salvato] {save_path.name}")


# =============================================================================
# 3. RILEVAMENTO: DETECTION RATE PER SEVERITA DEL DANNO (BAR CHART)
# =============================================================================
def plot_detection_rate_by_severity(
    di_damaged: np.ndarray,
    meta_damaged: np.ndarray,
    tau_global: float,
    save_path: Path = FIGURES_DIR / "3_detection_rate_by_severity.png",
):
    """
    Calcola e plotta la percentuale di rilevamento (True Positive Rate) per ciascuna
    severita di danno presente nel dataset (10%, 25%, 50%, 75%).
    Colonna 6 di meta_damaged = Severity.
    """
    severities = np.unique(meta_damaged[:, 6])
    severities = severities[severities > 0]  # Esclude lo 0 se presente
    detection_rates = []

    for sev in severities:
        mask = meta_damaged[:, 6] == sev
        di_sev = di_damaged[mask]
        tpr = np.mean(di_sev > tau_global) * 100.0
        detection_rates.append(tpr)

    plt.figure(figsize=(8, 5))
    bars = plt.bar(
        [f"-{int(s*100)}%" for s in severities],
        detection_rates,
        color="coral", edgecolor="black", width=0.5, alpha=0.85
    )

    # Scrittura della percentuale sopra ciascuna barra
    for bar in bars:
        h = bar.get_height()
        plt.text(
            bar.get_x() + bar.get_width() / 2.0, h + 1.5,
            f"{h:.1f}%", ha="center", va="bottom", fontsize=10, fontweight="bold"
        )

    plt.ylim(0, 115)
    plt.xlabel("Severita del Danno Flessionale (Riduzione EJ)", fontsize=11)
    plt.ylabel("Detection Rate (Sensitivity) [%]", fontsize=11)
    plt.title("Capacita di Rilevamento in Funzione della Gravita del Danno", fontsize=12, fontweight="bold")
    plt.grid(True, axis="y", linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    plt.close()
    print(f"[Grafico salvato] {save_path.name}")


# =============================================================================
# 4. RILEVAMENTO: CURVA ROC E CALCOLO AUC
# =============================================================================
def plot_roc_curve_metric(
    di_healthy: np.ndarray,
    di_damaged: np.ndarray,
    save_path: Path = FIGURES_DIR / "4_roc_curve.png",
):
    """
    Genera la Curva ROC (Receiver Operating Characteristic) e calcola l'Area Under Curve (AUC).
    """
    y_true = np.concatenate([np.zeros(len(di_healthy)), np.ones(len(di_damaged))])
    y_scores = np.concatenate([di_healthy, di_damaged])

    fpr, tpr, _ = roc_curve(y_true, y_scores)
    roc_auc = auc(fpr, tpr)

    plt.figure(figsize=(7, 6))
    plt.plot(fpr, tpr, color="darkorange", linewidth=2.5, label=f"ConvAutoencoder2D (AUC = {roc_auc:.4f})")
    plt.plot([0, 1], [0, 1], color="navy", linestyle="--", linewidth=1.5, label="Random Guess (AUC = 0.5000)")

    plt.xlim([-0.02, 1.0])
    plt.ylim([0.0, 1.02])
    plt.xlabel("False Positive Rate (Falsi Allarmi)", fontsize=11)
    plt.ylabel("True Positive Rate (Danni Rilevati)", fontsize=11)
    plt.title("Curva ROC - Prestazioni Complessive di Anomaly Detection", fontsize=12, fontweight="bold")
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.legend(fontsize=10, loc="lower right")
    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    plt.close()
    print(f"[Grafico salvato] {save_path.name}")


# =============================================================================
# 5. LOCALIZZAZIONE: PROFILO SPAZIALE LUNGO I 9 SENSORI DEL PONTE
# =============================================================================
def plot_damage_localization_profile(
    di_sensors_healthy: np.ndarray,
    di_sensors_damaged: np.ndarray,
    tau_sensors: np.ndarray,
    damage_x_target: float = 7.5,  # Elemento 8 -> ~7.5m (quarto campata L/4)
    save_path: Path = FIGURES_DIR / "5_damage_localization_profile.png",
):
    """
    Plotta l'errore medio di ricostruzione per ciascuno dei 9 sensori lungo l'asse del ponte (30 m).
    Evidenzia visivamente dove si concentra l'anomalia rispetto alla posizione vera del difetto.
    """
    mean_err_h = np.mean(di_sensors_healthy, axis=0)
    mean_err_d = np.mean(di_sensors_damaged, axis=0)

    plt.figure(figsize=(10, 5.5))

    # Zona reale del danno (Elemento 8: tra x=7m e x=8m)
    plt.axvspan(
        damage_x_target - 0.7, damage_x_target + 0.7,
        color="red", alpha=0.15, label=f"Posizione Reale Danno (L/4, x={damage_x_target:.1f}m)"
    )

    # Profilo sano
    plt.plot(SENSOR_X, mean_err_h, "b-o", linewidth=2.0, markersize=6, label="Test Sano (Profilo Medio)")

    # Profilo danneggiato
    plt.plot(SENSOR_X, mean_err_d, "r-s", linewidth=2.2, markersize=6, label="Test Danneggiato (Profilo Medio)")

    # Soglie per sensore
    plt.plot(SENSOR_X, tau_sensors, "k--^", linewidth=1.5, markersize=5, label="Soglie Locali Sensore $\\tau_{c}$")

    plt.xticks(SENSOR_X, SENSOR_LABELS, fontsize=9)
    plt.xlim(1.0, 29.0)
    plt.xlabel("Posizione Spaziale dei Sensori lungo la Trave [m]", fontsize=11)
    plt.ylabel("Errore di Ricostruzione Medio (DI Locale)", fontsize=11)
    plt.title("Localizzazione del Danno: Profilo dell'Errore sui 9 Sensori", fontsize=12, fontweight="bold")
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.legend(fontsize=10)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    plt.close()
    print(f"[Grafico salvato] {save_path.name}")


# =============================================================================
# 6. LOCALIZZAZIONE: HEATMAP SPAZIO-TEMPO (SENSORI VS TEMPO)
# =============================================================================
def plot_space_time_heatmap(
    x_real: np.ndarray,
    x_recon: np.ndarray,
    window_id: int = 14,
    save_path: Path = FIGURES_DIR / "6_space_time_heatmap.png",
):
    """
    Mostra la mappa 2D dell'errore assoluto di ricostruzione |x - x_recon|
    (Sensori sull'asse Y vs Tempo sull'asse X) per visualizzare l'anomalia nel piano spazio-tempo.
    x_real e x_recon hanno shape [9, 200].
    """
    error_2d = np.abs(x_real - x_recon)  # [9, 200]

    plt.figure(figsize=(10, 5))
    im = plt.imshow(
        error_2d, aspect="auto", cmap="hot", origin="lower",
        extent=[0, 1.0, 0.5, 9.5]  # tempo da 0 a 1s, sensori da 1 a 9
    )
    plt.colorbar(im, label="Errore Assoluto di Ricostruzione $|x - \\hat{x}|$")

    plt.yticks(range(1, 10), [f"S{i} ({x:.0f}m)" for i, x in enumerate(SENSOR_X)])
    plt.xlabel("Tempo all'interno della Finestra [s]", fontsize=11)
    plt.ylabel("Sensori lungo la Travata", fontsize=11)
    plt.title(f"Mappa Spazio-Temporale dell'Errore (Finestra {window_id})", fontsize=12, fontweight="bold")
    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    plt.close()
    print(f"[Grafico salvato] {save_path.name}")
