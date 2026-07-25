import torch
import torch.nn.functional as F
import torch.nn as nn
from torch_geometric.nn import GATConv

class EdgeModel(torch.nn.Module):
    def __init__(self, edge_dim, node_dim):
        super().__init__()
        self.edge_mlp = nn.Sequential(
            nn.Linear(2 * node_dim + edge_dim, 256),
            nn.ReLU(),
            nn.Linear(256, 128)
        )

    def forward(self, src, dst, edge_attr):
        out = torch.cat([src, dst, edge_attr], dim=1)
        return self.edge_mlp(out)

class NodeModel(torch.nn.Module):
    def __init__(self, in_channels, hidden_channels, out_channels, heads):
        super().__init__()
        self.gat1 = GATConv(in_channels, hidden_channels, heads, dropout=0.3)
        self.gat2 = GATConv(hidden_channels * heads, out_channels, heads=1, concat=False, dropout=0.3)

    def forward (self, x, edge_index):

        x = F.dropout(x, p=0.3, training=self.training)
        x = F.elu(self.gat1(x, edge_index))
        x = F.dropout(x, p=0.3, training=self.training)
        x = self.gat2(x, edge_index)
        return x

    



    

    
    