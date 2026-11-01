import torch
import torch.nn as nn
import torch.nn.functional as F

class ConvAutoencoder2D(nn.Module):
    def __init__(self, latent_dim: int = 32):
        super().__init__()
        self.latent_dim = latent_dim

        # -------------------------------------------------------------
        # 1. ENCODER CONVOLUZIONALE (Spazio-Temporale)
        # -------------------------------------------------------------
        self.conv_encoder = nn.Sequential(
            # Layer 1: [B, 1, 9, 200] -> [B, 16, 9, 100]
            # Stride=(1, 2) riduce a metà SOLO il tempo, lasciando i 9 sensori intatti
            nn.Conv2d(in_channels=1, out_channels=16, kernel_size=3, stride=(1, 2), padding=1),
            nn.BatchNorm2d(16),
            nn.ReLU(),

            # Layer 2: [B, 16, 9, 100] -> [B, 32, 5, 50]
            # Stride=(2, 2) riduce sia i sensori (9 -> 5) sia il tempo (100 -> 50)
            nn.Conv2d(in_channels=16, out_channels=32, kernel_size=3, stride=(2, 2), padding=1),
            nn.BatchNorm2d(32),
            nn.ReLU(),
        )

        # -------------------------------------------------------------
        # 2. BOTTLENECK LINEARE (Spazio Latente Z)
        # -------------------------------------------------------------
        # 32 canali * 5 sensori * 50 passi temporali = 8000 feature
        self.flat_features = 32 * 5 * 50

        # Da mappa 2D a vettore latente compatto Z
        self.fc_enc = nn.Linear(self.flat_features, self.latent_dim)
        
        # Da vettore Z di nuovo a 8000 valori
        self.fc_dec = nn.Linear(self.latent_dim, self.flat_features)

        # -------------------------------------------------------------
        # 3. DECODER DE-CONVOLUZIONALE (Ricostruzione)
        # -------------------------------------------------------------
        self.conv_decoder = nn.Sequential(
            # Deconv 1: [B, 32, 5, 50] -> [B, 16, 9, 100]
            # output_padding=(0, 1) serve per compensare l'arrotondamento del tempo (50*2 = 100)
            nn.ConvTranspose2d(
                in_channels=32, out_channels=16,
                kernel_size=3, stride=(2, 2), padding=1, output_padding=(0, 1)
            ),
            nn.BatchNorm2d(16),
            nn.ReLU(),

            # Deconv 2 (Output): [B, 16, 9, 100] -> [B, 1, 9, 200]
            nn.ConvTranspose2d(
                in_channels=16, out_channels=1,
                kernel_size=3, stride=(1, 2), padding=1, output_padding=(0, 1)
            )
            # NOTA: Niente ReLU qui! Il segnale normalizzato ha anche valori negativi.
        )

    def encode(self, x: torch.Tensor) -> torch.Tensor:
        """Estrae il vettore latente compatto Z a partire dall'immagine dei sensori."""
        h = self.conv_encoder(x)
        h = h.view(h.size(0), -1)  # Appiattimento a vettore [B, 8000]
        z = self.fc_enc(h)          # Compressione a [B, latent_dim]
        return z

    def decode(self, z: torch.Tensor) -> torch.Tensor:
        """Ricostruisce la griglia del ponte a partire dal vettore latente Z."""
        h = self.fc_dec(z)                                  # [B, 8000]
        h = h.view(-1, 32, 5, 50)                           # Rimodellamento a griglia 2D [B, 32, 5, 50]
        reconstructed = self.conv_decoder(h)                # [B, 1, 9, 200]
        return reconstructed

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        z = self.encode(x)
        x_recon = self.decode(z)
        return x_recon
    
# ===============================================================================
def load_trained_model(checkpoint_path, device="cpu"):
    """
    Ricarica il modello ConvAutoencoder2D con i pesi addestrati e lo imposta in eval().
    """
    checkpoint = torch.load(checkpoint_path, map_location=device, weights_only=False)
    latent_dim = checkpoint.get("latent_dim", 32)
    
    model = ConvAutoencoder2D(latent_dim=latent_dim).to(device)
    model.load_state_dict(checkpoint["model_state_dict"])
    model.eval()
    return model
    
    
if __name__ == "__main__":
    model = ConvAutoencoder2D(latent_dim=32)
    
    # Batch casuale di 8 finestre da 9 sensori e 200 campioni
    x = torch.randn(8, 1, 9, 200)
    
    z = model.encode(x)
    x_out = model(x)
    
    print(f"Input shape:        {x.shape}")        # torch.Size([8, 1, 9, 200])
    print(f"Spazio Latente Z:   {z.shape}")        # torch.Size([8, 32])
    print(f"Output ricostruito: {x_out.shape}")   # torch.Size([8, 1, 9, 200])
    
    assert x.shape == x_out.shape, "Le dimensioni di input e output non coincidono!"
    print("Test superato con successo!")
