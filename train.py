import time

import torch
from torch import nn
import torch.nn.functional as F
from torch.optim import AdamW
from scipy.stats import pearsonr

def train_iter (epoch, model, optimizer, input, target, print_every=10):
    start_t = time.time()
    model.train()
    optimizer.zero_grad()

    output = model(*input)
    loss = F.l1_loss(output, target)

    loss.backward()
    optimizer.step()

    if epoch % print_every == 0:
         loss_val = test(model, input, target)

         print(f"Epoch: {epoch:04d} " 
               f"time: ({(time.time() - start_t):.4f}s "
               f"loss_val: {loss_val:.4f} "
            )
         
    return loss.item

def test(model, input, target):
    model.eval()
    with torch.no_grad():
        output=model(*input)
        
        mae = F.l1_loss(output, target)

    return mae.item()




