import torch
import torch.nn.functional as F
import torch.nn as nn
from torch_geometric.nn import GATv2Conv
from torch_geometric.nn import MetaLayer

from ml.prediction_head import BFactorHead

class ProteinGNN(torch.nn.Module):
    def __init__(self, edge_dim, node_dim):
        super().__init__()

        self.layers = nn.ModuleList([
            MetaLayer(
                EdgeModel (edge_dim, node_dim),
                NodeModel(36, 32, 64, 4, 5)
            ),
            MetaLayer(
                EdgeModel (128, 64),
                NodeModel(64, 32, 64, 4, 128)
            ),
            MetaLayer(
                EdgeModel (128, 64),
                NodeModel(64, 32, 64, 4, 128)
            ),
            MetaLayer(
                EdgeModel (128, 64),
                NodeModel(64, 32, 64, 4, 128)
            ),
            MetaLayer(
                EdgeModel (128, 64),
                NodeModel(64, 32, 64, 4, 128)
            ),
        ])

        self.head = BFactorHead(hidden_channels=64)

    def forward(self, data):
        x = data.x
        edge_attr = data.edge_attr

        for layer in self.layers:
            x, edge_attr, _ = layer(x, data.edge_index, edge_attr, None, data.batch)

        return self.head(x)

class EdgeModel(torch.nn.Module):

    def __init__(self, edge_dim, node_dim):
        super().__init__()
        self.edge_mlp = nn.Sequential(
            nn.Linear(2 * node_dim + edge_dim, 256),
            nn.ReLU(),
            nn.Linear(256, 128)
        )

    def forward(self, src, dst, edge_attr, u=None, batch=None):
        out = torch.cat([src, dst, edge_attr], dim=1)
        return self.edge_mlp(out)

class NodeModel(torch.nn.Module):
    def __init__(self, in_channels, hidden_channels, out_channels, heads, edge_dim):
        super().__init__()
        self.gat1 = GATv2Conv(in_channels, hidden_channels, heads, edge_dim=edge_dim, dropout=0.3)
        self.gat2 = GATv2Conv(hidden_channels * heads, out_channels, heads=1, edge_dim=edge_dim, concat=False, dropout=0.3)

    def forward (self, x, edge_index, edge_attr, u=None, batch=None):

        x = F.dropout(x, p=0.3, training=self.training)
        x = F.elu(self.gat1(x, edge_index, edge_attr))
        x = F.dropout(x, p=0.3, training=self.training)
        x = self.gat2(x, edge_index, edge_attr)
        return x

