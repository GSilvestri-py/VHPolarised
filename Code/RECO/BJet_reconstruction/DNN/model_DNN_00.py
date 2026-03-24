# choose import in config file

import torch
import torch.nn as nn

# ------- ------- ------- ------- ------- ------- ------- ------- ------- ------- -------

class ResidualBlock(nn.Module):
    def __init__(self, dim, dropout_prob=0.0, layer_norm=True):
        super().__init__()

        layers = []
        layers.append(nn.Linear(dim, dim))

        if layer_norm:
            layers.append(nn.LayerNorm(dim))

        layers.append(nn.GELU())

        self.block = nn.Sequential(*layers)

        self.dropout = nn.Dropout(dropout_prob) if dropout_prob > 0 else nn.Identity()

    def forward(self, x):
        out = self.block(x)
        out = self.dropout(out)
        return x + out   

# ------- ------- ------- ------- ------- ------- ------- ------- ------- ------- -------

class DNN_ResNet(nn.Module):
    def __init__(self, N_dims, N_hidden, N_layers, N_classes,
                 dropout_prob=0.1, layer_norm=True):
        super().__init__()

        self.do_layers = []   

        # input layer
        layers = []
        layers.append(nn.Linear(N_dims, N_hidden))
        if layer_norm:
            layers.append(nn.LayerNorm(N_hidden))
        layers.append(nn.GELU())

        # hidden layers
        self.blocks = nn.ModuleList()
        for _ in range(N_layers):
            block = ResidualBlock(N_hidden, dropout_prob, layer_norm)
            self.blocks.append(block)
            if dropout_prob > 0:
                self.do_layers.append(block.dropout)

        # output layer
        self.output = nn.Linear(N_hidden, N_classes)

        self.model = nn.Sequential(*layers)

        # Weight initialization
        self.apply(self._init_weights)

    def _init_weights(self, m):
        if isinstance(m, nn.Linear):
            nn.init.xavier_uniform_(m.weight)

    def forward(self, x):
        x = self.model(x)
        for block in self.blocks:
            x = block(x)
        return self.output(x)

    def scale_do(self, scale):
        for layer in self.do_layers:
            if isinstance(layer, nn.Dropout):
                layer.p *= scale


# ------- ------- ------- ------- ------- ------- ------- ------- ------- ------- -------


class DNN_00(nn.Module):
    def __init__(self, N_dims, N_hidden, N_layers, N_classes,
                 dropout_prob=0.1, layer_norm=True):
        super().__init__()

        layers = []
        self.do_layers = []   

        in_dim = N_dims
        hidden_dim = N_hidden

        for i in range(N_layers):
            layers.append(nn.Linear(in_dim, hidden_dim))

            if layer_norm:
                layers.append(nn.LayerNorm(hidden_dim))

            layers.append(nn.ReLU())

            if dropout_prob > 0:
                do = nn.Dropout(dropout_prob)
                layers.append(do)
                self.do_layers.append(do)   

            in_dim = hidden_dim

        layers.append(nn.Linear(hidden_dim, N_classes))

        self.model = nn.Sequential(*layers)

    def forward(self, x):
        return self.model(x)

    def scale_do(self, scale):
        for layer in self.do_layers:
            layer.p *= scale

'''
class DNN_00 (nn.Module) : 	
    def __init__ (self, N_dims, N_hidden, N_layers, N_classes, dropout_prob = 0., layer_norm = True):
        Set the model structure:
        - N_dims      = number of input variables
        - N_hidden    = number of nodes per hidden layer
        - N_layers    = number of layers betwen input and output
        - N_classes   = number of output categories
        - droput_prob = dropout probability (default 0)
        - layer_norm  = whether to apply layer normalisation (default True)
        

        super ().__init__ ()
        if (N_layers < 1) : N_layers = 1 

        layers = []
        self.do_layers = []

        #first layer: from input to hidden
        layers.append (torch.nn.Linear (N_dims, N_hidden))          
        if (dropout_prob > 0) : 
            layers.append (nn.Dropout (dropout_prob))
            self.do_layers.append (layers[-1])
        if (layer_norm) : layers.append (nn.LayerNorm (N_hidden))
        layers.append (torch.nn.ReLU ())

        #intermediate layers
        print("number of layers: ", N_layers)
        for i_layer in range (N_layers):
            next_hid_dim = max(1, int(N_hidden * 2/3))
            layers.append (torch.nn.Linear (N_hidden, next_hid_dim))
            if (dropout_prob > 0) : 
                layers.append (nn.Dropout (dropout_prob))
                self.do_layers.append (layers[-1])
            if (layer_norm) : layers.append (nn.LayerNorm (next_hid_dim))
            layers.append (torch.nn.ReLU ())
            
            N_hidden = next_hid_dim
        
        #last layer: from hidden to output
        layers.append (torch.nn.Linear (N_hidden, N_classes))

        self.model = nn.Sequential (*layers) 

    # def get_weights (self):
    #     return (self.model[0].weight.clone (), self.model[2].weight.clone ())

    def forward (self, x):
        logits = self.model (x)
        return logits

    def scale_do (self, scale):
        for layer in self.do_layers: layer.p *= scale
'''