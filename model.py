from metalayer import ProteinGNN
import train as train
from dataset import ProteinDataset

import torch
from torch.optim import AdamW
from torch_geometric.loader import DataLoader
from torch.utils.data import random_split

def main(): 

    if torch.cuda.is_available():
        device = torch.device('cuda')
    else:
        device = torch.device('cpu')

    dataset = ProteinDataset(root="processed")

    generator = torch.Generator().manual_seed(42)

    train_size = int(0.8 * len(dataset))
    val_size = int(0.1 * len(dataset))
    test_size = len(dataset) - train_size - val_size

    train_dataset, val_dataset, test_dataset = random_split(
    dataset, [train_size, val_size, test_size], generator = generator)

    train_loader = DataLoader(train_dataset, batch_size=1, shuffle=True)

    model = ProteinGNN(5, 36).to(device)

    optimizer = AdamW(model.parameters(), lr=0.0007, betas=(0.9, 0.999), weight_decay=0.0001)

    epochs = 50

    for epoch in range(epochs):

        for data in train_loader: 
            data = data.to(device)
            train.train_epoch(epoch, model, optimizer, data)

if __name__ == "__main__":
    main()