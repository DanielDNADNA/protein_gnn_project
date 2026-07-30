import time

import torch
import torch.nn.functional as F

def train_epoch (epoch, model, optimizer, data, print_every=10):
    start_t = time.time()
    model.train()
    optimizer.zero_grad()

    output = model(data)
    loss = F.l1_loss(output, data.y)

    loss.backward()
    optimizer.step()

    if epoch % print_every == 0:
         loss_val = test(model, data)

         print(f"Epoch: {epoch:04d} " 
               f"time: {(time.time() - start_t):.4f}s "
               f"loss_val: {loss_val:.4f} "
            )
         
    return loss.item()

def test(model, data):
    model.eval()
    with torch.no_grad():
        output=model(data)
        
        mae = F.l1_loss(output, data.y)

    return mae.item()




