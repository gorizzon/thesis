# The MoGCN scripts can only be run in the terminal. I copied the functions here to import them
# in other files. Nearly all of the code is not mine, it is the same as the original files.
# The original MoGCN scripts are: AE_run.py, SNF.py, and MoGCN_run.py.

import pandas as pd
import numpy as np
import argparse
from tqdm import tqdm
import autoencoder_model
import torch
import torch.utils.data as Data
import snf
import seaborn as sns
import glob
import os
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import f1_score
import torch.nn.functional as F
from gcn_model import GCN
from utils import load_data
from utils import accuracy

def setup_seed(seed):
    torch.manual_seed(seed)
    np.random.seed(seed)

##########################################################################################

def work(data, in_feas, lr=0.001, bs=32, epochs=100, device=torch.device('cpu'), a=0.4, b=0.3, c=0.3, mode=0, topn=100):
    
    bs = int(bs)
        
    #name of sample
    sample_name = data['Sample'].tolist()

    #change data to a Tensor
    X,Y = data.iloc[:,1:].values, np.zeros(data.shape[0])
    TX, TY = torch.tensor(X, dtype=torch.float, device=device), torch.tensor(Y, dtype=torch.float, device=device)
    #train a AE model
    if mode == 0 or mode == 1:
        print('Training model...')
        Tensor_data = Data.TensorDataset(TX, TY)
        train_loader = Data.DataLoader(Tensor_data, batch_size=bs, shuffle=True)

        #initialize a model
        mmae = autoencoder_model.MMAE(in_feas, latent_dim=100, a=a, b=b, c=c)
        mmae.to(device)
        mmae.train()
        mmae.train_MMAE(train_loader, learning_rate=lr, device=device, epochs=epochs)
        mmae.eval()       #before save and test, fix the variables
        torch.save(mmae, 'model/AE/MMAE_model.pkl')

    #load saved model, used for reducing dimensions
    if mode == 0 or mode == 2:
        print('Get the latent layer output...')
        mmae = torch.load('model/AE/MMAE_model.pkl')
        omics_1 = TX[:, :in_feas[0]]
        omics_2 = TX[:, in_feas[0]:in_feas[0]+in_feas[1]]
        omics_3 = TX[:, in_feas[0]+in_feas[1]:in_feas[0]+in_feas[1]+in_feas[2]]
        latent_data, decoded_omics_1, decoded_omics_2, decoded_omics_3 = mmae.forward(omics_1, omics_2, omics_3)
        latent_df = pd.DataFrame(latent_data.detach().cpu().numpy())
        latent_df.insert(0, 'Sample', sample_name)
        #save the integrated data(dim=100)
        latent_df.to_csv('result/latent_data.csv', header=True, index=False)

    print('Extract features...')
    extract_features(data, in_feas, epochs, topn)
    return

def extract_features(data, in_feas, epochs, topn=100):
    # extract features
    #get each omics data
    data_omics_1 = data.iloc[:, 1: 1+in_feas[0]]
    data_omics_2 = data.iloc[:, 1+in_feas[0]: 1+in_feas[0]+in_feas[1]]
    data_omics_3 = data.iloc[:, 1+in_feas[0]+in_feas[1]: 1+in_feas[0]+in_feas[1]+in_feas[2]]

    #get all features of each omics data
    feas_omics_1 = data_omics_1.columns.tolist()
    feas_omics_2 = data_omics_2.columns.tolist()
    feas_omics_3 = data_omics_3.columns.tolist()

    #calculate the standard deviation of each feature
    std_omics_1 = data_omics_1.std(axis=0)
    std_omics_2 = data_omics_2.std(axis=0)
    std_omics_3 = data_omics_3.std(axis=0)

    #record top N features every 10 epochs
    topn_omics_1 = pd.DataFrame()
    topn_omics_2 = pd.DataFrame()
    topn_omics_3 = pd.DataFrame()

    #used for feature extraction, epoch_ls = [10,20,...], if epochs % 10 != 0, add the last epoch
    epoch_ls = list(range(10, epochs+10,10))
    if epochs %10 != 0:
        epoch_ls.append(epochs)
    for epoch in tqdm(epoch_ls):
        #load model
        mmae = torch.load('model/AE/model_{}.pkl'.format(epoch))
        #get model variables
        model_dict = mmae.state_dict()

        #get the absolute value of weights, the shape of matrix is (n_features, latent_layer_dim)
        weight_omics1 = np.abs(model_dict['encoder_omics_1.0.weight'].detach().cpu().numpy().T)
        weight_omics2 = np.abs(model_dict['encoder_omics_2.0.weight'].detach().cpu().numpy().T)
        weight_omics3 = np.abs(model_dict['encoder_omics_3.0.weight'].detach().cpu().numpy().T)

        weight_omics1_df = pd.DataFrame(weight_omics1, index=feas_omics_1)
        weight_omics2_df = pd.DataFrame(weight_omics2, index=feas_omics_2)
        weight_omics3_df = pd.DataFrame(weight_omics3, index=feas_omics_3)

        #calculate the weight sum of each feature --> sum of each row
        weight_omics1_df['Weight_sum'] = weight_omics1_df.apply(lambda x:x.sum(), axis=1)
        weight_omics2_df['Weight_sum'] = weight_omics2_df.apply(lambda x:x.sum(), axis=1)
        weight_omics3_df['Weight_sum'] = weight_omics3_df.apply(lambda x:x.sum(), axis=1)
        weight_omics1_df['Std'] = std_omics_1
        weight_omics2_df['Std'] = std_omics_2
        weight_omics3_df['Std'] = std_omics_3

        #importance = Weight * Std
        weight_omics1_df['Importance'] = weight_omics1_df['Weight_sum']*weight_omics1_df['Std']
        weight_omics2_df['Importance'] = weight_omics2_df['Weight_sum']*weight_omics2_df['Std']
        weight_omics3_df['Importance'] = weight_omics3_df['Weight_sum']*weight_omics3_df['Std']

        #select top N features
        fea_omics_1_top = weight_omics1_df.nlargest(topn, 'Importance').index.tolist()
        fea_omics_2_top = weight_omics2_df.nlargest(topn, 'Importance').index.tolist()
        fea_omics_3_top = weight_omics3_df.nlargest(topn, 'Importance').index.tolist()

        #save top N features in a dataframe
        col_name = 'epoch_'+str(epoch)
        topn_omics_1[col_name] = fea_omics_1_top
        topn_omics_2[col_name] = fea_omics_2_top
        topn_omics_3[col_name] = fea_omics_3_top

    #all of top N features
    topn_omics_1.to_csv('result/topn_omics_1.csv', header=True, index=False)
    topn_omics_2.to_csv('result/topn_omics_2.csv', header=True, index=False)
    topn_omics_3.to_csv('result/topn_omics_3.csv', header=True, index=False)


def AE_run(path1, path2, path3, mode=0, seed=0, batchsize=32, learningrate=0.001, epoch=100, latent=100, a=0.4, b=0.3, c=0.3, topn=100):
        
        #read data
    omics_data1 = pd.read_csv(path1, header=0, index_col=None)
    omics_data2 = pd.read_csv(path2, header=0, index_col=None)
    omics_data3 = pd.read_csv(path3, header=0, index_col=None)

    #Check whether GPUs are available
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    #set random seed
    setup_seed(seed)

    if a + b + c != 1.0:
        print('The sum of weights must be 1.')
        exit(1)

    #dims of each omics data
    in_feas = [omics_data1.shape[1] - 1, omics_data2.shape[1] - 1, omics_data3.shape[1] - 1]
    omics_data1.rename(columns={omics_data1.columns.tolist()[0]: 'Sample'}, inplace=True)
    omics_data2.rename(columns={omics_data2.columns.tolist()[0]: 'Sample'}, inplace=True)
    omics_data3.rename(columns={omics_data3.columns.tolist()[0]: 'Sample'}, inplace=True)

    omics_data1.sort_values(by='Sample', ascending=True, inplace=True)
    omics_data2.sort_values(by='Sample', ascending=True, inplace=True)
    omics_data3.sort_values(by='Sample', ascending=True, inplace=True)

    #merge the multi-omics data, calculate on common samples
    Merge_data = pd.merge(omics_data1, omics_data2, on='Sample', how='inner')
    Merge_data = pd.merge(Merge_data, omics_data3, on='Sample', how='inner')
    Merge_data.sort_values(by='Sample', ascending=True, inplace=True)

    #train model, reduce dimensions and extract features
    work(Merge_data, in_feas, lr=learningrate, bs=batchsize, epochs=epoch, device=device, a=a, b=b, c=c, mode=mode, topn=topn)
    print('Success! Results can be seen in result file')


##########################################################################################

def SNF_run(path1, path2, path3, metric='sqeuclidean', K=20, mu=0.5):
    print('Load data files...')
    omics_data_1 = pd.read_csv(path1, header=0, index_col=None)
    omics_data_2 = pd.read_csv(path2, header=0, index_col=None)
    omics_data_3 = pd.read_csv(path3, header=0, index_col=None)
    print(omics_data_1.shape, omics_data_2.shape, omics_data_3.shape)

    if omics_data_1.shape[0] != omics_data_2.shape[0] or omics_data_1.shape[0] != omics_data_3.shape[0]:
        print('Input files must have same samples.')
        exit(1)

    omics_data_1.rename(columns={omics_data_1.columns.tolist()[0]: 'Sample'}, inplace=True)
    omics_data_2.rename(columns={omics_data_2.columns.tolist()[0]: 'Sample'}, inplace=True)
    omics_data_3.rename(columns={omics_data_3.columns.tolist()[0]: 'Sample'}, inplace=True)

    # align samples of different data
    omics_data_1.sort_values(by='Sample', ascending=True, inplace=True)
    omics_data_2.sort_values(by='Sample', ascending=True, inplace=True)
    omics_data_3.sort_values(by='Sample', ascending=True, inplace=True)

    print('Start similarity network fusion...')
    affinity_nets = snf.make_affinity([omics_data_1.iloc[:, 1:].values.astype(np.float), omics_data_2.iloc[:, 1:].values.astype(np.float), omics_data_3.iloc[:, 1:].values.astype(np.float)],
                                      metric=metric, K=K, mu=mu)

    fused_net =snf.snf(affinity_nets, K=K)

    print('Save fused adjacency matrix...')
    fused_df = pd.DataFrame(fused_net)
    fused_df.columns = omics_data_1['Sample'].tolist()
    fused_df.index = omics_data_1['Sample'].tolist()
    fused_df.to_csv('result/SNF_fused_matrix.csv', header=True, index=True)

    np.fill_diagonal(fused_df.values, 0)
    fig = sns.clustermap(fused_df.iloc[:, :], cmap='vlag', figsize=(8,8),)
    fig.savefig('result/SNF_fused_clustermap.png', dpi=300)
    print('Success! Results can be seen in result file')

##########################################################################################

def train(epoch, optimizer, features, adj, labels, idx_train, device, GCN_model):
    '''
    :param epoch: training epochs
    :param optimizer: training optimizer, Adam optimizer
    :param features: the omics features
    :param adj: the laplace adjacency matrix
    :param labels: sample labels
    :param idx_train: the index of trained samples
    '''
    labels.to(device)

    GCN_model.train()
    optimizer.zero_grad()
    output = GCN_model(features, adj)
    loss_train = F.cross_entropy(output[idx_train], labels[idx_train])
    acc_train = accuracy(output[idx_train], labels[idx_train])
    loss_train.backward()
    optimizer.step()
    if (epoch+1) % 10 ==0:
        print('Epoch: %.2f | loss train: %.4f | acc train: %.4f' %(epoch+1, loss_train.item(), acc_train.item()))
    return loss_train.data.item()

def test(features, adj, labels, idx_test, GCN_model):
    '''
    :param features: the omics features
    :param adj: the laplace adjacency matrix
    :param labels: sample labels
    :param idx_test: the index of tested samples
    '''
    GCN_model.eval()
    output = GCN_model(features, adj)
    loss_test = F.cross_entropy(output[idx_test], labels[idx_test])

    #calculate the accuracy
    acc_test = accuracy(output[idx_test], labels[idx_test])

    #output is the one-hot label
    ot = output[idx_test].detach().cpu().numpy()
    #change one-hot label to digit label
    ot = np.argmax(ot, axis=1)
    #original label
    lb = labels[idx_test].detach().cpu().numpy()
    print('predict label: ', ot)
    print('original label: ', lb)

    #calculate the f1 score
    f = f1_score(ot, lb, average='weighted')

    print("Test set results:",
          "loss= {:.4f}".format(loss_test.item()),
          "accuracy= {:.4f}".format(acc_test.item()))

    #return accuracy and f1 score
    return acc_test.item(), f


def GCN_cross_validation(featuredata, adjdata, labeldata, cv_splits, epochs=150, learningrate=0.001, weight_decay=0.01, hidden=64, dropout=0.5, nclass=4, patience=20):
    
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    
    # set random seed
    setup_seed(seed)

    # load input files
    adj, data, label = load_data(adjdata, featuredata, labeldata)

    # change dataframe to Tensor
    adj = torch.tensor(adj, dtype=torch.float, device=device)
    features = torch.tensor(data.iloc[:, 1:].values, dtype=torch.float, device=device)
    #labels = torch.tensor([label.items()["Tumor"]], dtype=torch.long, device=device)
    labels = torch.tensor(label.iloc[:, 1].values, dtype=torch.long, device=device)

    print('Begin training model...')
    
    skf = StratifiedKFold(n_splits=cv_splits, shuffle=True)
    acc_res, f1_res = [], []

    for idx_train, idx_test in skf.split(features, labels):
        GCN_model = GCN(n_in=features.shape[1], n_hid=hidden, n_out=nclass, dropout=dropout)
        GCN_model.to(device)
        optimizer = torch.optim.Adam(GCN_model.parameters(), lr=learningrate, weight_decay=weight_decay)

        idx_train, idx_test = torch.tensor(idx_train, dtype=torch.long, device=device), torch.tensor(idx_test, dtype=torch.long, device=device)

        for epoch in range(epochs):
            train(epoch, optimizer, features, adj, labels, idx_train, device, GCN_model)

        ac, f1 = test(features, adj, labels, idx_test, GCN_model)
        acc_res.append(ac)
        f1_res.append(f1)
        
    print('Finished!')

    return np.mean(acc_res), np.std(acc_res), np.mean(f1_res), np.std(f1_res)

def GCN_train_test(features, adj, labels, testsample, epochs=150, learningrate=0.001, weight_decay=0.01, hidden=64, dropout=0.5, nclass=4, patience=20):
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    
    # set random seed
    setup_seed(seed)

    # load input files
    adj, data, label = load_data(adjdata, featuredata, labeldata, threshold)

    # change dataframe to Tensor
    adj = torch.tensor(adj, dtype=torch.float, device=device)
    features = torch.tensor(data.iloc[:, 1:].values, dtype=torch.float, device=device)
    #labels = torch.tensor([label.items()["Tumor"]], dtype=torch.long, device=device)
    labels = torch.tensor(label.iloc[:, 1].values, dtype=torch.long, device=device)

    print('Begin training model...')
    
    test_sample_df = pd.read_csv(testsample, header=0, index_col=None)
    test_sample = test_sample_df.iloc[:, 0].tolist()
    all_sample = data['Sample'].tolist()
    train_sample = list(set(all_sample)-set(test_sample))
    train_idx = data[data['Sample'].isin(train_sample)].index.tolist()
    test_idx = data[data['Sample'].isin(test_sample)].index.tolist()

    GCN_model = GCN(n_in=features.shape[1], n_hid=hidden, n_out=nclass, dropout=dropout)
    GCN_model.to(device)
    optimizer = torch.optim.Adam(GCN_model.parameters(), lr=learningrate, weight_decay=weight_decay)
    idx_train, idx_test = torch.tensor(train_idx, dtype=torch.long, device=device), torch.tensor(test_idx, dtype=torch.long, device=device)

    loss_values = []
    bad_counter, best_epoch = 0, 0
    best = 1000

    for epoch in range(epochs):
        loss_values.append(train(epoch, optimizer, features, adj, labels, idx_train))
        if loss_values[-1] < best:
            best = loss_values[-1]
            best_epoch = epoch
            bad_counter = 0
        else:
            bad_counter += 1

        if bad_counter == patience:
            break

        torch.save(GCN_model.state_dict(), 'model/GCN/{}.pkl'.format(epoch))

        files = glob.glob('model/GCN/*.pkl')
        for file in files:
            name = file.split('/')[-1]
            epoch_nb = int(name.split('.')[0])
            if epoch_nb != best_epoch:
                os.remove(file)

    print('Training finished.')
    print('The best epoch model is ', best_epoch)
    GCN_model.load_state_dict(torch.load('model/GCN/{}.pkl'.format(best_epoch)))
    predict(features, adj, all_sample, test_idx)

    print('Finished!')

