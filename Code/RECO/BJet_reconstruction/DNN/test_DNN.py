import configparser
import importlib
import numpy as np
import os
import pandas as pd
import pickle
import random
import sys
import torch
import torch.nn as nn

from torch.utils.data import random_split
from torch.optim.lr_scheduler import LambdaLR
from torch.utils.data import Dataset, DataLoader


# ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----

class TabularDataset(Dataset):
    """Torch container for a spiral dataset."""

    def __init__(self, X, y):
        """Initializes instance of class StudentsPerformanceDataset.
        Args:
            input_collection: sample randomly generated with generate_spiral_data
        """
        # Save target and predictors
        self.values = X
        self.labels = y

    def __len__ (self) :
        return len (self.values)

    def __getitem__ (self, idx) :
        return [self.values[idx], self.labels[idx]]


#------------------------------------------------------------------

def init_weights (model_element):
    
    classname = model_element.__class__.__name__
    if classname.find ('Linear') != -1: 
        torch.nn.init.xavier_uniform_(model_element.weight)
    if classname.find ('Conv') != -1:
        torch.nn.init.normal_ (model_element.weight, 0.0, 0.02)
    elif classname.find ('BatchNorm') != -1:
        torch.nn.init.normal_(model_element.weight, 1.0, 0.02)
        torch.nn.init.zeros_(model_element.bias)


#------------------------------------------------------------------


class EarlyStopping:

    def __init__(self, patience=10, verbose=False):
        self.patience = patience
        self.verbose = verbose
        self.counter = 0
        self.best_score = None
        self.early_stop = False
        self.best_epoch = None
        self.best_model_state_dict = None

    def __call__(self, val_loss, model, epoch):
        score = -val_loss  

        if self.best_score is None:
            self.best_score = score
            self.best_model_state_dict = model.state_dict()
            self.best_epoch = epoch
        elif score < self.best_score:
            self.counter += 1
            if self.counter >= self.patience:
                self.early_stop = True
        else:
            self.best_score = score
            self.best_model_state_dict = model.state_dict()
            self.best_epoch = epoch
            self.counter = 0
            if self.verbose:
                print(f"EarlyStopping: improvement at epoch {epoch}")


#------------------------------------------------------------------
'''
def expo_lr_scheduler (epoch, initial_lr, final_lr, delay = 0) :
  
  #returns a multiplicative number to the learning rate
  #so that the learning rate decrease is exponential 
  #and has an horizontal asymptote at final_lr
  

  if (epoch < delay) : return 1.

  lr_frac = final_lr / initial_lr
  return (lr_frac + np.exp (lr_frac - epoch + delay - 1))/(lr_frac + np.exp (lr_frac - epoch + delay))
'''    

def expo_lr_scheduler(epoch, initial_lr, final_lr, total_epochs, delay=0):
    if epoch < delay:
        return 1.0
    decay_rate = (final_lr / initial_lr) ** (1 / (total_epochs - delay))
    return decay_rate ** (epoch - delay)

#------------------------------------------------------------------


def train_loop (dataloader, model, loss_fn, optimizer, device, verbosity = 1):
    ''' 
    the function to be called to train the model for one epoch,
    looping over serval batches
    '''
    size = len (dataloader.dataset)
    
    model.train ()
    loss = -1.
    current = -1.
    
    for batch, (X, y) in enumerate (dataloader):      
        X, y = X.to (device), y.to (device)

        #Compute prediction error
        pred = model (X)
        loss = loss_fn (pred, y)       

        #Backpropagation
        loss.backward ()		
        optimizer.step ()
        optimizer.zero_grad ()

        if verbosity != 0 and batch % 100 == 0:
            loss, current = loss.item (), (batch + 1) * len (X)
            print (f"loss: {loss:>7f}  [{current:>5d}/{size:>5d}]")


# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

def test_loop (dataloader, model, loss_fn, verbosity = 1):
    ''' 
    test the model 
    '''
    
    model.eval ()
    size = len (dataloader.dataset)
    num_batches = len (dataloader)
    test_loss, correct = 0., 0.
    TP, TN, FP, FN = 0., 0., 0., 0.

    
    with torch.no_grad () :
        
        for X, y in dataloader :
            pred = model (X)

            test_loss += loss_fn (pred, y).item () * X.size (0)

            TP += ((pred.argmax (1) == 1) & (y == 1)).sum ().item () #true positive etc
            TN += ((pred.argmax (1) == 0) & (y == 0)).sum ().item ()
            FP += ((pred.argmax (1) == 1) & (y == 0)).sum ().item ()
            FN += ((pred.argmax (1) == 0) & (y == 1)).sum ().item ()

           
    correct = TP + TN
    test_loss /= size
    accuracy = correct / size				                    
    precision = TP / (TP + FP) if (TP + FP) > 0 else 0          
    recall = TP / (TP + FN) if (TP + FN) > 0 else 0	            
    F1 = 2 * (precision * recall) / (precision + recall) if TP > 0 else 0 

    if verbosity > 0 : 
        
        print (f"{dataloader} metrics:")
        print (f' - Accuracy : {(100*accuracy):>0.1f}%')
        print (f' - Avg loss : {test_loss:>8f}')
        print (f' - Precision: {precision:>8f}')
        print (f' - Recall   : {recall:>8f}')
        print (f' - F1       : {F1:>8f}')

    model.train()

    return test_loss, accuracy, precision, recall, F1
    
# ------------------------------------------------------------------------------------------

def prepare_samples(config, output_path):

    train_frac = config.getfloat('training', 'train_frac', fallback = 0.7)
    valid_frac = config.getfloat('training', 'valid_frac', fallback = 0.15)
    batch_size = config.getint('training', 'batch_size', fallback = 64)

    df = pd.read_pickle(config.get('samples', 'file'))
    df["Tag"] = df["Tag"].map({"b jet": 1, "non b jet": 0})     

    X = torch.tensor(df.iloc[:, 1:-1].values, dtype = torch.float32)
    Y = torch.tensor(df.iloc[:, -1].values, dtype = torch.long)

    print(f"feature shape = {X.shape}")
    print(f"tags shape = {Y.shape}")

    idxs = torch.randperm(X.shape[0])
    X = X[idxs]
    Y = Y[idxs]

    input_dataset = TabularDataset(X, Y)

    N_train = int(len(input_dataset) * train_frac)
    N_valid = int(len(input_dataset) * valid_frac)
    N_test  = (len(input_dataset) - N_train - N_valid)    #avoid improper sum to 1
    train_set, valid_set, test_set = random_split(input_dataset, [N_train, N_valid, N_test])

    var_names = list(df.columns[1:-1])
    print(f"List of input features = {var_names}")

    # build train and validation dataloader for training
    train_dataloader = DataLoader(train_set, batch_size = batch_size, shuffle = True)
    valid_dataloader = DataLoader(valid_set, batch_size = batch_size, shuffle = True)
    #test_dataloader = DataLoader(test_set, batch_size = batch_size, shuffle = True)

    # dump test set in a separate file -> use for evaluation
    print(var_names)
    X_test = []
    Y_test = []

    for X, y in test_set:
        X_test.append(X.numpy())
        Y_test.append(y.numpy())

    X_test = np.array(X_test)
    Y_test = np.array(Y_test)

    df_test = pd.DataFrame(X_test, columns = var_names) #build df with variables name
    df_test["Tag"] = Y_test    
    
    print(df)
    print(df_test)                         

    filename = "test_set_dataframe.pkl"
    df_test.to_pickle(os.path.join(output_path, filename))

    return train_dataloader, valid_dataloader


# ------- ------- ------- ------- ------- ------- ------- ------- ------- ------- -------


def model_training(config, train_dataloader, valid_dataloader, device, verbosity, model_class):

    epochs               = config.getint ('training','epochs', fallback = 1000)
    layer_norm           = config.getboolean ('training','layer_norm', fallback = True)
    learning_rate        = config.getfloat ('training','learning_rate', fallback = 1e-3)
    lr_annealing         = config.getfloat ('training','lr_annealing', fallback = 0.)
    lr_minimum           = config.getfloat ('training','lr_minimum', fallback = 1.)
    lr_delay             = config.getfloat ('training','lr_delay', fallback = 1.)
    weight_decay         = config.getfloat ('training','weight_decay', fallback = 1e-5)   
    dropout_prob         = config.getfloat ('training','dropout_prob', fallback = 0.)
    dropout_annealing    = config.getfloat ('training','dropout_annealing', fallback = 0.)
    if dropout_annealing >= 1. : dropout_annealing = -1
    out_dim              = config.getint ('model', 'classes_num')
    hidden_layers        = config.getint ('training','hidden_layers', fallback = 3)
    hidden_units         = config.getint ('training','hidden_units', fallback = 2)
    save_model           = config.getboolean ('saving','save_model', fallback = True)

    apply_early_stopping = config.getboolean ('training','apply_early_stopping', fallback = True)
    es_patience          = config.getint ('training','es_patience', fallback = 10)  
    es_starting          = config.getint ('training','es_starting', fallback = 10)
    
    #print parameters to screen

    params = [
              "epochs", "layer_norm", "learning_rate", "lr_annealing", "lr_minimum", "lr_delay",
              "weight_decay", "dropout_prob", "dropout_annealing", "out_dim", "hidden_layers",
              "save_model", "apply_early_stopping", "es_patience", "es_starting"
              ]

    params_values = [
                    epochs, layer_norm, learning_rate, lr_annealing, lr_minimum, lr_delay,
                    weight_decay, dropout_prob, dropout_annealing, out_dim, hidden_layers,
                    save_model, apply_early_stopping, es_patience, es_starting
                    ]

    for n, v in zip(params, params_values):
        print(f"{n} = {v}")
    
    out_path = config.get('saving', 'out_path', fallback = None)
    out_file_base        = config.get ('saving', 'out_file_base', fallback = None)

    loss_fn = nn.CrossEntropyLoss (reduction ='mean')

    input_dimension = 0
    for batch in train_dataloader:
        inp, targets = batch
        input_dimension = inp.shape[1]
        break

    hidden_units = input_dimension

    # ------------------ importing the model ------------------

    model = model_class(input_dimension, hidden_units, hidden_layers, out_dim, dropout_prob, layer_norm)
    model = model.to(device)
    print(str(model))

    # ------------------ define optimizer ------------------
    model.apply(init_weights)
    optimiser = torch.optim.Adam(model.parameters(), lr = learning_rate, weight_decay=weight_decay)

    early_stopping = EarlyStopping(patience = es_patience, verbose = True)
    lr_annealer = lambda e : expo_lr_scheduler (e, learning_rate, lr_minimum, epochs, lr_delay)
    scheduler = LambdaLR (optimiser, lr_lambda = lr_annealer)
    
    # ------------------ model training ------------------

    train_log = []
    valid_log = []
        
    for t in range(epochs):

        #if t % 10 == 0:
        print(f" - epoch {t+1}")
        print("LR:", scheduler.get_last_lr()[0])
        print("dropout:", model.do_layers[0].p)

        
        # ------------------------ train vs validation overview ------------------------
        train_loop(train_dataloader, model, loss_fn, optimiser, device)        
        train_log.append(test_loop(train_dataloader, model, loss_fn, verbosity=1))
        valid_log.append(test_loop(valid_dataloader, model, loss_fn, verbosity=1))


        if dropout_annealing > 0. : model.scale_do (dropout_annealing)
        if lr_annealing > 0       : scheduler.step ()


        partial_file = out_file_base + '_trend_partial.pickle'
        with open(os.path.join(out_path, partial_file), "wb") as f:
            pickle.dump((train_log, valid_log), f)


        # ------------------------ early stopping ------------------------
        if not apply_early_stopping or t < es_starting: 
            best_epoch = epochs
            continue

        early_stopping(valid_log[-1][0], model, t)

        # --- save partial best model at every epoch ---
        torch.save(early_stopping.best_model_state_dict,
            os.path.join(out_path, out_file_base + "_best_model_partial.pth"))


        if early_stopping.early_stop:
            print(f'EARLY STOPPING triggered at epoch {t} with best epoch {early_stopping.best_epoch}\n')
            best_epoch = early_stopping.best_epoch

            while len(train_log) < epochs:
                train_log.append(train_log[-1])
            while len(valid_log) < epochs:
                valid_log.append(valid_log[-1])

            break   

    print(f"Finished training at epoch {best_epoch} with train loss of {train_log[-1][0]} and validation loss of {valid_log[-1][0]}")

    # ------------------------ save trends ------------------------

    train_valid_trends = (train_log, valid_log)
    output_file = out_file_base + '_trend.pickle'
    with open(os.path.join(out_path, output_file), "wb") as f:
        pickle.dump(train_valid_trends, f)
    print(f"Training and validation trends saved in {os.path.join(out_path, output_file)}")

    # ------------------------ save model ------------------------
    if save_model:
        output_file = out_file_base + '_best_model.pth'
        torch.save(early_stopping.best_model_state_dict, os.path.join(out_path, output_file))
    print(f"Best DNN model saved in {os.path.join(out_path, output_file),}")

    return model


# ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----


def main():

    if len(sys.argv) < 2:
        print(f"usage: {sys.argv[0]} input.cfg")
        sys.exit(1)

    # ---------------------------------- read config file ----------------------------------
    
    config = configparser.ConfigParser()
    config.read('config.cfg')

    out_path = config.get('saving', 'out_path', fallback=None)
    
    if out_path == None:
        sys.exit(1)
    
    os.makedirs(out_path, exist_ok=True)


    print(f"output path = {out_path}")
    print(f"Sections in config file = {config.sections()}")

    with open(out_path + 'config.cfg', 'w') as configfile:
        config.write(configfile)


    verbosity = config.getint ('running','verbosity', fallback = 20)
    seed = config.getint('training', 'seed', fallback=-1)
    if seed > 0:
        random.seed(seed)
        np.random.seed(seed)
        torch.manual_seed(seed)

    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

    # ---------------------------------- prepare samples ----------------------------------

    train_dataloader, valid_dataloader = prepare_samples(config, out_path)

    # ---------------------------------- define model ----------------------------------
    
    library_name = config.get('model', 'library_name')
    model_name = config.get('model', 'model_name')

    print(f"Using library = {library_name}")
    print(f"Using model = {model_name}")

    model_library = importlib.import_module(library_name)
    model_class = getattr(model_library, model_name)

    # ---------------------------------- train model ----------------------------------

    model = model_training(config, train_dataloader, valid_dataloader, device, verbosity, model_class)

if __name__ == "__main__" : 
    main ()