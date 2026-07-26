import torch.nn.functional as F
import torch.nn as nn

class BFactorHead(nn.Module):

    def __init__(self, hidden_channels):
        super().__init__()

        self.mlp = nn.Sequential(
            nn.Linear(hidden_channels, 128),
            nn.ReLU(),
            nn.Linear(128, 1)
        )

    def forward(self, x):
        return self.mlp(x)