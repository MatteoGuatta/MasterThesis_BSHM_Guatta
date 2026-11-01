"""
metrics.py
===============================================================================
Master Thesis in Bridge Structural Health Monitoring (BSHM)
Funzioni di Calcolo per Damage Indices (DI) e Soglie di Anomalia (Thresholds)
===============================================================================
"""

import torch
from torch.utils.data import DataLoader

def compute_reconstruction_errors(model, dataloader, device):
    """
    Calcola il Damage Index globale e per singolo sensore per tutte le finestre del dataloader.
    Restituisce:
      - all_di_sensors: shape [N_tot, 9] (DI locale per sensore)
      - all_di_global:  shape [N_tot]    (DI globale per finestra)
      - all_meta:       shape [N_tot, 8] (Metadati associati a ciascuna finestra)
    """
    model.eval()
    all_di_global = []
    all_di_sensors = []
    all_meta = []
    
    with torch.no_grad():
        for x, meta in dataloader:
            # 1. Adattamento per Conv2d: da [B, 9, 200] a [B, 1, 9, 200]
            inputs = x.unsqueeze(1).to(device)
            outputs = model(inputs)
            
            # 2. Errore quadratico [B, 9, 200] (rimuovendo il canale 1)
            diff = (inputs.squeeze(1) - outputs.squeeze(1)) ** 2
            
            # 3. DI Locale (media sul tempo dim=2) -> [B, 9]
            reconstruction_error_local = torch.mean(diff, dim=2)
            
            # 4. DI Globale (media su sensori e tempo dim=(1,2)) -> [B]
            reconstruction_error_global = torch.mean(diff, dim=(1, 2))
            
            all_di_global.append(reconstruction_error_global.cpu())
            all_di_sensors.append(reconstruction_error_local.cpu())
            all_meta.append(meta.cpu())
    
    # 5. Concatenazione finale di tutti i batch
    all_di_global = torch.cat(all_di_global, dim=0)    # shape: [N_tot]
    all_di_sensors = torch.cat(all_di_sensors, dim=0)  # shape: [N_tot, 9]
    all_meta = torch.cat(all_meta, dim=0)              # shape: [N_tot, 8]
    
    return all_di_sensors, all_di_global, all_meta


def compute_thresholds(all_di_sensors, all_di_global, method="3sigma", alpha=0.01):
    """
    Calcola le soglie di anomalia (1 globale e 9 locali per sensore).
    Supporta:
      - method="3sigma": regola della media + 3 * deviazione standard
      - method="percentile": percentile (es. 99° con alpha=0.01)
    """
    if method == "3sigma":
        threshold_global = torch.mean(all_di_global) + 3 * torch.std(all_di_global)
        thresholds_sensor = torch.mean(all_di_sensors, dim=0) + 3 * torch.std(all_di_sensors, dim=0)
    elif method == "percentile":
        q = (1.0 - alpha) * 100.0
        # Calcolo percentile con numpy
        q_glob = float(torch.quantile(all_di_global, 1.0 - alpha))
        threshold_global = torch.tensor(q_glob)
        thresholds_sensor = torch.quantile(all_di_sensors, 1.0 - alpha, dim=0)
    else:
        raise ValueError(f"Metodo non supportato: {method}")
    
    return thresholds_sensor, threshold_global
