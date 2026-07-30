from pathlib import Path
import torch 
from torch_geometric.data import Dataset

class ProteinDataset(Dataset):

    def __init__(self, root):
        self.files = sorted(Path(root).glob("*.pt"))
        super().__init(root)

    def len(self):
        return len(self.files)

    def get(self, idx):
        return torch.load(self.files[idx])
