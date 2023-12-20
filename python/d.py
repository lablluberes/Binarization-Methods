from dash import Dash, dash_table, dcc, html, Input, Output, callback, State
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import math
import matplotlib as plt

def Y_a_b(genes, a, b):
    #print("Mean ", np.mean(genes[a:b]), " from a ",a , " to b ", b)
    return np.mean(genes[a:b])

def C_a_b(genes, a, b):
    mean = Y_a_b(genes, a, b+1)
    return sum( (np.array(genes[a:b+1]) - mean) ** 2 )


def determine_h(P, i, j, genes):
    N = len(genes)

    if (i == 0 and j > 0):
        return Y_a_b(genes, P[i][j], P[i+1][j]) - Y_a_b(genes, 0, P[i][j]);
    elif (i == j and j > 0):
        return Y_a_b(genes, P[i][j], N) - Y_a_b(genes, P[i-1][j], P[i][j]);
    elif (i == 0 and j == 0):
        return Y_a_b(genes, P[i][j], N) - Y_a_b(genes, 0, P[i][j]);
    else:
        return Y_a_b(genes, P[i][j], P[i+1][j]) - Y_a_b(genes, P[i-1][j], P[i][j]);

def BASC_A(gene):
    gene_og = gene
    gene = np.sort(gene)
    N = len(gene)

    cost_matrix = [[0 for _ in range(N - 1)] for _ in range(N)]
    ind_matrix = [[0 for _ in range(N - 2)] for _ in range(N - 1)]
    P = [[0 for _ in range(N - 2)] for _ in range(N - 2)]

    # Step 1: Compute a Series of Step Function

    # initialization C_i_(0) = c_i_N
    # calculate first cost matrix column with no intermidiate break points
    for i in range(N):
        cost_matrix[i][0] = C_a_b(gene, i, N)

    # Algorithm 1: Calculate optimal step functions
    for j in range(N-2):
        for i in range(N-j-1):
            min_value = math.inf
            min_index = math.inf

            for d in range(N-j-1):
                if(i <= d):
                    curr_value = C_a_b(gene, i, d) + cost_matrix[d+1][j]

                if(curr_value < min_value):
                    min_value = curr_value
                    min_index = d

            cost_matrix[i][j+1] = min_value
            ind_matrix[i][j] = min_index + 1

    #  Algorithm 2: Compute the break points of all optimal step functions
    for j in range(N-2):
        z = j
        P[0][j] = ind_matrix[0][z]
        if(j > 0):
            z = z - 1
            for i in range(1, j+1):
                P[i][j] = ind_matrix[P[i-1][j]][z]
                z = z - 1

    # Step 2: Find Strongest Discontinuity in Each Step Function
    v = [0] * (N-2)

    for j in range(N-2):
        max_value = -math.inf
        max_index = j
        for i in range(j+1):
            h = determine_h(P, i, j, gene)
            z = (gene[P[i][j]] + gene[P[i][j]-1]) / 2
            e = sum( (np.array(gene) - z) ** 2 )
            q_score = h / e
            if(q_score > max_value):
                max_value = q_score
                max_index = i

        v[j] = P[max_index][j]

    # Step 3: Estimate Location and Variation of the Strongest Discontinuities
    thr = (gene[round(np.median(v))-1] + gene[round(np.median(v))]) / 2

    return thr, P

from sklearn.cluster import KMeans

def K_Means(genes):
    data = np.array(genes).reshape(-1, 1)
    kmeans = KMeans(n_clusters=2, n_init='auto')
    kmeans.fit(data)
    c=kmeans.labels_
    genes = np.array(genes)
    groupOne = genes[c==1]
    groupZero = genes[c==0]
    
    thr1 = np.mean(groupOne)
    thr2 = np.mean(groupZero)
    
    thr = (thr1 + thr2) / 2

    return thr

def getSSTOT(x, n, xmean):
    m = 0
    for i in range(n):
        m = m + (x[i] - xmean)**2
    return m


def onestep(x):
    
    n = len(x)
    #step = 0
    xmean = np.mean(x)
    SSTOT = getSSTOT(x, n, xmean)
    
    SSEmin = SSTOT
    
    for i in range(n-1):
        leftMean = np.mean(x[0:i+1])
    
        rightMean = np.mean(x[i+1:n])
        
        SSE = 0
        
        for j in range(n):
            if j < i+1:
                SSE = SSE + (x[j] - leftMean)**2
            else:
                SSE = SSE + (x[j] - rightMean)**2
                    
        
        if SSEmin > SSE:
            SSEmin = SSE
            #print("1:",SSEmin1)
                
            t = (leftMean + rightMean)/2
        
    
    return t

def shmulevich(x):
    
    n = len(x)
    s = np.sort(x)
    d = np.empty(n)
    
    for i in range(n-2):
        d[i] = s[i+1] - s[i]
    
    t = (s[n-1] - s[0])/(n-1)
    
    mn = s[n-1]
    index = 0
    
    for i in range(n-1):
        if d[i] > t and d[i] < mn:
            mn = d[i]
            index = i
            
    z = s[index + 1]
   
    
    return z

#threshold/displacement indexes are vectors with values in them
#for different algorithms

def binarizationVoting(gene, threshold, displacement):

    #2d array that stores algorithm results
    #x axis is gene
    #y axis is algorithm

    alg = len(threshold)
    n = len(gene)

    algos = np.zeros([alg, n])

    #get results for every algorithm
    for i in range(alg): 
        for j in range(n):
            if (gene[j] + displacement[i] < threshold[i]) and (abs(threshold[i] - gene[j]) > displacement[i]):
                algos[i][j] = 1
            elif(abs(threshold[i] - gene[j]) <= displacement[i]):
                algos[i][j] = 2
            else:
                algos[i][j] = 3

    #id 
    # 1 -> U
    # 2 -> N
    # 3 -> E

    # start counting votes
    majority = np.zeros(n)
       
    for i in range(n):

        #array of tally votes
        results = algos[:,i]
        U = np.count_nonzero(results == 1)
        N = np.count_nonzero(results == 2)
        E = np.count_nonzero(results == 3)

        #binarize tally votes
        tally = [U,N,E]
        for j in range(3):
            #if algorithms is an even number this will be half, else itll be a
            #number rounded up from half ex. 3 -> 1.5 becomes 2
            if tally[j] >= math.ceil(alg/2):
                tally[j] = 1
            else:
                tally[j] = 0

        #first check if final tally contradicts itself
        if sum(tally) != 1:
            majority[i] = np.nan
        #then check other cases
        else:
            if tally[0] == 1:
                majority[i] = 0
            elif tally[1] == 1:
                majority[i] = np.nan
            else:
                majority[i] = 1



     #return 2d array of votes + final votes
    return algos, majority


def interpolationConverge(vect, method, tolerance):

    thr = []
    converge = False
    n = (len(vect) - 1)
    newSize = len(vect) + (len(vect) - 1)
    gene = vect
    sample = []
    limit = 10000
    conver = 0
    
    while(~converge and newSize < limit):
        
        sample.append(newSize)
        indices = np.arange(len(vect))

        interpolated_values = np.interp(
            np.linspace(0, len(vect) - 1, len(vect) + n),
            indices,
            gene
        )
        
        #print(interpolated_values, "\n")
        
        gene = interpolated_values
        
        if(method == 'K-Means'):
            thr.append(K_Means(vect))
        elif(method == 'BASC A'):
            t, _ = BASC_A(vect)
            thr.append(t)
        elif(method == 'Onestep'):
            thr.append(onestep(vect))
        else:
            thr.append(shmulevich(vect))
    
        #print(interpolated_values)
        
        n_thr = len(thr)
        
        for i in range(n_thr):
            for j in range(n_thr):
                if(i != j):
                    difference = abs(thr[i] - thr[j])
                    
                    if(difference <= tolerance):
                        converge = True
                        conver = thr[j]
                        break
              
        if(~converge):
            n = newSize - 1
            newSize = newSize + (newSize-1)
            vect = interpolated_values
        


    tMin = min(thr)
    tMax = max(thr)
    
    return tMin, tMax, conver, newSize

def three_interpolation(vect, method):
    thr = []
    n = (len(vect) - 1)
    newSize = len(vect) + (len(vect) - 1)
    gene = vect
    sample = []
    for i in range(3):
        sample.append(newSize)
        indices = np.arange(len(gene))

        interpolated_values = np.interp(
            np.linspace(0, len(gene) - 1, len(gene) + n),
            indices,
            gene
        )
        
        #print(interpolated_values, "\n")
        
        gene = interpolated_values
        
        if(method == 'K-Means'):
            thr.append(K_Means(gene))
        elif(method == 'BASC A'):
            t, _ = BASC_A(gene)
            thr.append(t)
        elif(method == 'Onestep'):
            thr.append(onestep(gene))
        else:
            thr.append(shmulevich(gene))
        
        n = newSize - 1
        newSize = newSize + (newSize-1)
    
    return thr, sample

import networkx as nx

def A_state(state):
    return state[1]

def B_state(state):
    return int(not state[0])

def C_state(state):
    return int(state[0] or state[1])

def create_boolean_network(selected, method, data):
    
    gene_data = pd.DataFrame(data)
    
    genes = gene_data.iloc[selected].values.astype(float)
    print(genes)
    
    df = pd.DataFrame(columns=['A', 'B', 'C', "A'", "B'", "C'"])
    
    binarize_array = []
    binarize = []
    
    for gene in genes:
        if(method == "K-Means"):
            thr = K_Means(gene)
        elif(method == "Shmulevich"):
            thr = shmulevich(gene)
        elif(method == "Onestep"):
            thr = onestep(gene)
        else:
            thr, _ = BASC_A(gene)
            
        for exp in gene:
            if(exp <= thr):
                binarize.append(0)
            else:
                binarize.append(1)
                
        binarize_array.append(binarize)
        binarize = []
    
    
    df["A"] = binarize_array[0]
    df["B"] = binarize_array[1]
    df["C"] = binarize_array[2]
    
    #print(df)
    
    result = df[['A', 'B', 'C']]
    
    for index, row in result.iterrows():
        state = [row['A'], row['B'], row['C']]
        
        df.at[index, "A'"] = A_state(state)
        df.at[index, "B'"] = B_state(state)
        df.at[index, "C'"] = C_state(state)
        
    
    
    return df

import matplotlib.pyplot as plt

def create_boolean_network_graph(selected, data):
    
    KMeans = pd.DataFrame(create_boolean_network(selected, "K-Means", data))
    BASCA = pd.DataFrame(create_boolean_network(selected, "BASC A", data))
    Onestep = pd.DataFrame(create_boolean_network(selected, "Onestep", data))
    Shmul = pd.DataFrame(create_boolean_network(selected, "Shmulevich", data))
    
    G1 = nx.DiGraph()
    G2 = nx.DiGraph()
    G3 = nx.DiGraph()
    G4 = nx.DiGraph()
    
    for index, row in KMeans.iterrows():
        KMeans_curr_state = str(row["A"]) + str(row["B"]) + str(row["C"])
        KMeans_next_state = str(row["A'"]) + str(row["B'"]) + str(row["C'"])
        
        BASC_curr_state = str(BASCA.at[index ,"A"]) + str(BASCA.at[index, "B"]) + str(BASCA.at[index, "C"])
        BASC_next_state = str(BASCA.at[index ,"A'"]) + str(BASCA.at[index, "B'"]) + str(BASCA.at[index, "C'"])
        
        One_curr_state = str(Onestep.at[index ,"A"]) + str(Onestep.at[index, "B"]) + str(Onestep.at[index, "C"])
        One_next_state = str(Onestep.at[index ,"A'"]) + str(Onestep.at[index, "B'"]) + str(Onestep.at[index, "C'"])
        
        Shm_curr_state = str(Shmul.at[index ,"A"]) + str(Shmul.at[index, "B"]) + str(Shmul.at[index, "C"])
        Shm_next_state = str(Shmul.at[index ,"A'"]) + str(Shmul.at[index, "B'"]) + str(Shmul.at[index, "C'"])
        
        '''if(KMeans_curr_state == BASC_curr_state == One_curr_state == Shm_curr_state):
            
            G1.add_node(KMeans_curr_state, color = 'gray')

            G2.add_node(BASC_curr_state, color = 'gray')

            G3.add_node(One_curr_state, color = 'gray')

            G4.add_node(Shm_curr_state, color = 'gray')
            
        else:
            
            G1.add_node(KMeans_curr_state, color = 'red')

            G2.add_node(BASC_curr_state, color = 'red')

            G3.add_node(One_curr_state, color = 'red')

            G4.add_node(Shm_curr_state, color = 'red')
            
        if(KMeans_next_state == BASC_next_state == One_next_state == Shm_next_state):
            
            G1.add_node(KMeans_next_state, color = 'gray')

            G2.add_node(BASC_next_state, color = 'gray')
            
            G3.add_node(One_next_state, color = 'gray')

            G4.add_node(Shm_next_state, color = 'gray')
            
        else:
            
            G1.add_node(KMeans_next_state, color = 'red')

            G2.add_node(BASC_next_state, color = 'red')

            G3.add_node(One_next_state, color = 'red')

            G4.add_node(Shm_next_state, color = 'red')'''
        
        G1.add_node(KMeans_curr_state, color = 'gray')
        G2.add_node(BASC_curr_state, color = 'gray')
        G3.add_node(One_curr_state, color = 'gray')
        G4.add_node(Shm_curr_state, color = 'gray')
        
        G1.add_node(KMeans_next_state, color = 'gray')
        G2.add_node(BASC_next_state, color = 'gray')   
        G3.add_node(One_next_state, color = 'gray')
        G4.add_node(Shm_next_state, color = 'gray')
            
        G1.add_edge(KMeans_curr_state, KMeans_next_state)
        G2.add_edge(BASC_curr_state, BASC_next_state)
        G3.add_edge(One_curr_state, One_next_state)
        G4.add_edge(Shm_curr_state, Shm_next_state)
        
    fig = plt.figure()
    
    pos = nx.spring_layout(G1) 
    plt.subplot(221)
    plt.title("K-Means")
    nx.draw(G1, pos, with_labels=True, node_color=[G1.nodes[n]['color'] for n in G1.nodes])
    #plt.show()
    
    #plt.figure()
    
    plt.subplot(222)
    plt.title("BASCA")
    pos = nx.spring_layout(G2) 
    nx.draw(G2, pos, with_labels=True, node_color=[G2.nodes[n]['color'] for n in G2.nodes])
    #plt.show()
    
    #plt.figure()
    
    plt.subplot(223)
    plt.title("Onestep")
    pos = nx.spring_layout(G3) 
    nx.draw(G3, pos, with_labels=True, node_color=[G3.nodes[n]['color'] for n in G3.nodes])
    #plt.show()
    
    #plt.figure()
    
    plt.subplot(224)
    plt.title("Shmulevich")
    pos = nx.spring_layout(G4) 
    nx.draw(G4, pos, with_labels=True, node_color=[G4.nodes[n]['color'] for n in G4.nodes])
    plt.show()
    
    return fig
    

def graph_three_interpolation_network(selected, method, df):
    
    inter1Net = pd.DataFrame(columns=['A', 'B', 'C', "A'", "B'", "C'"])
    inter2Net = pd.DataFrame(columns=['A', 'B', 'C', "A'", "B'", "C'"])
    inter3Net = pd.DataFrame(columns=['A', 'B', 'C', "A'", "B'", "C'"])
    
    gene_data = pd.DataFrame(df)
    
    genes = gene_data.iloc[selected].values.astype(float)
    
    thr_arr_geneA = three_interpolation(genes[0], method)
    thr_arr_geneB = three_interpolation(genes[1], method)
    thr_arr_geneC = three_interpolation(genes[2], method)
    
    #print(thr_arr_geneA, thr_arr_geneB, thr_arr_geneC  )
    
    interA_Binarized = []
    interB_Binarized = []
    interC_Binarized = []
    
    binarized1 = []
    binarized2 = []
    binarized3 = []
    
    #print(genes[0], genes[1], genes[2])
    
    for i in range(3):
        
        for j in range(len(genes[0])):
            if(genes[0][j] <= thr_arr_geneA[0][i]):
                binarized1.append(0)
            else:
                binarized1.append(1)
                
            if(genes[1][j] <= thr_arr_geneB[0][i]):
                binarized2.append(0)
            else:
                binarized2.append(1)
                
            if(genes[2][j] <= thr_arr_geneC[0][i]):
                binarized3.append(0)
            else:
                binarized3.append(1)
                
        interA_Binarized.append(binarized1)
        interB_Binarized.append(binarized2)
        interC_Binarized.append(binarized3)
        
        binarized1 = []
        binarized2 = []
        binarized3 = []
                
    inter1Net["A"] = interA_Binarized[0]
    inter1Net["B"] = interB_Binarized[0]
    inter1Net["C"] = interC_Binarized[0]
    
    inter2Net["A"] = interA_Binarized[1]
    inter2Net["B"] = interB_Binarized[1]
    inter2Net["C"] = interC_Binarized[1]
    
    inter3Net["A"] = interA_Binarized[2]
    inter3Net["B"] = interB_Binarized[2]
    inter3Net["C"] = interC_Binarized[2]
    
    #print(interA_Binarized)
    
    #print(inter1Net)
    #print(inter2Net)
    #print(inter3Net)
    
    result = inter1Net[['A', 'B', 'C']]
    
    for index, row in result.iterrows():
        state1 = [row['A'], row['B'], row['C']]
        state2 = [inter2Net.at[index, 'A'], inter2Net.at[index, 'B'], inter2Net.at[index, 'C']]
        state3 = [inter3Net.at[index, 'A'], inter3Net.at[index, 'B'], inter3Net.at[index, 'C']]
        
        inter1Net.at[index, "A'"] = A_state(state1)
        inter1Net.at[index, "B'"] = B_state(state1)
        inter1Net.at[index, "C'"] = C_state(state1)
        
        inter2Net.at[index, "A'"] = A_state(state2)
        inter2Net.at[index, "B'"] = B_state(state2)
        inter2Net.at[index, "C'"] = C_state(state2)
        
        inter3Net.at[index, "A'"] = A_state(state3)
        inter3Net.at[index, "B'"] = B_state(state3)
        inter3Net.at[index, "C'"] = C_state(state3)
        
    #print(inter1Net)
    #print(inter2Net)
    #print(inter3Net)
    
    
    #plot1 = graph_network(inter1Net, thr_arr_geneA[1][0])
    
    #plot2 = graph_network(inter2Net, thr_arr_geneA[1][1])
    
    #plot3 = graph_network(inter3Net, thr_arr_geneA[1][2])
    
    G = nx.DiGraph()
    
    for index, row in inter1Net.iterrows():
        curr_state = str(row["A"]) + str(row["B"]) + str(row["C"])
        next_state = str(row["A'"]) + str(row["B'"]) + str(row["C'"])
        
        G.add_node(curr_state)
        G.add_node(next_state)
        G.add_edge(curr_state, next_state)
    
    fig = plt.figure()
    plt.subplot(221)
    pos = nx.spring_layout(G)
    plt.title("Interpolation of Size " + str(thr_arr_geneA[1][0]))
    nx.draw(G, pos, with_labels=True)
    
    G1 = nx.DiGraph()
    
    for index, row in inter2Net.iterrows():
        curr_state = str(row["A"]) + str(row["B"]) + str(row["C"])
        next_state = str(row["A'"]) + str(row["B'"]) + str(row["C'"])
        
        G1.add_node(curr_state)
        G1.add_node(next_state)
        G1.add_edge(curr_state, next_state)
    
    pos = nx.spring_layout(G1)
    plt.subplot(222)
    plt.title("Interpolation of Size " + str(thr_arr_geneA[1][1]))
    nx.draw(G1, pos, with_labels=True)
    
    G2 = nx.DiGraph()
    
    for index, row in inter3Net.iterrows():
        curr_state = str(row["A"]) + str(row["B"]) + str(row["C"])
        next_state = str(row["A'"]) + str(row["B'"]) + str(row["C'"])
        
        G2.add_node(curr_state)
        G2.add_node(next_state)
        G2.add_edge(curr_state, next_state)
    
    pos = nx.spring_layout(G2)
    plt.subplot(223)
    plt.title("Interpolation of Size " + str(thr_arr_geneA[1][2]))
    nx.draw(G2, pos, with_labels=True)

    return fig

import dash_bootstrap_components as dbc 
import base64
import io
from io import BytesIO

col_names = {'basc_thr':[], 'kmeans_thr':[], 'onestep_thr':[], 'shmulevich_thr':[]}
final_df = pd.DataFrame(col_names)

app = Dash(__name__)
server = app.server

app.layout = html.Div([
    html.Div([
                html.H4('Binarize Genes Visualization'), 
        
                dcc.Store(id='stored-data', storage_type='session'),
           
                dcc.Upload(
                    id='upload-data',
                    children=html.Button("Upload Gene Expression File"),
                    multiple=True
                ),
              
                html.Div(id='output-data-upload'),
            
                html.Div(id='select-all-binarize', style={"width": "25%"}),
               
                html.Div(id='binarize-all'),
                html.Hr()
               
          
    ]),
    
    html.Div([
        html.Div([html.Div(id='dropdown-methods', style={"width": "25%"}),
                  html.Div(id='heatmap-binarize'),
                  html.Div(id='interpolation-heatmap'),
                  html.Div(id='rules-dropdown'),
                  html.Img(id='graph_rules'),
                 ], style={'padding': 10, 'flex': 1}),
        html.Div([
            html.Div(id='select-gene-binarize', style={"width": "25%"}),
            html.Div(id='graph-gene-binarize'),
            html.Div(id='select-basc-discontinuity', style={"width": "25%"}),
            html.Div(id='graph-gene-discontinuity'),
            html.Div(id='tolerance-dropdown', style={"width": "25%"}),
            html.Div(id='interpolation-graph')
        ], style={'padding': 10, 'flex': 1})
    ], style={'display': 'flex', 'flexDirection': 'row'})
])
#, style={"width":"80%", "margin":"auto"}


# function to parse the contents of thje selected file
def parse_contents(contents, filename, date):
    content_type, content_string = contents.split(',')

    # decode the content
    decoded = base64.b64decode(content_string)
    
    # if it is a csv then read it 
    if 'csv' in filename:
        df = pd.read_csv(io.StringIO(decoded.decode('utf-8')), header=None)

        # return a new html div that has the file name and the datatable with selectable rows
        return html.Div([
           # name of file
           html.H5(filename),
            
           # dash datatable of data with rows that can be selected
           dash_table.DataTable(
                        id='datatable-interactivity',
                        columns=[
                            {"name": str(i), "id": str(i)} for i in df.columns
                        ],
                        data=df.to_dict('records'),
                        column_selectable="single",
                        row_selectable="multi",
                        selected_columns=[],
                        selected_rows=[],
                        page_action="native",
                        page_current= 0,
                        page_size= 10,
                    ),
            # store the read dataframe
            dcc.Store(id='stored-data', data=df.to_dict('records')),
            html.Hr(), 

        ])
    
    else:
        return "The file needs to be a csv."

# this callout receives the contents of the file and outputs the component
# output-data-upload
@callback(Output('output-data-upload', 'children'),
              Input('upload-data', 'contents'),
              State('upload-data', 'filename'),
              State('upload-data', 'last_modified'))
# function parses and update the output of the selected dataset
def update_output(list_of_contents, list_of_names, list_of_dates):
    # if there is a selected file
    if list_of_contents is not None:
        # parse the content
        children = [
            parse_contents(c, n, d) for c, n, d in
            zip(list_of_contents, list_of_names, list_of_dates)]
        return children

# callout to show the select all genes to binarize component
# the input is the stored data
@app.callback(
    Output('select-all-binarize', 'children'),
    Input('stored-data','data'),
     prevent_initial_call=True)
# function to show dropdowm to binarize all
def select_all_binarize(data):
    # if there is no data then return none
    if not data:
        return None
    
    # return the dropdown component to binarize all
    return [html.H5('Select option in order to show all the thresholds by each algorithm of each gene:'),
                dcc.Dropdown(
                            [{'label': 'Binarize All Genes', 'value':'all'}],
                            placeholder="Select to binarize all genes and get thresholds",
                            id="dropdown-binarize-all",
                            searchable=False)]


@app.callback(
    Output('binarize-all', 'children'),
    Input('dropdown-binarize-all', 'value'),
    Input('stored-data','data'))
def binarize_all(all_rows, data):
    if all_rows is None:
        return None
    
    df = pd.DataFrame(data)
    
    genes = df.values
    rows = df.shape[0]
    
    col_names = {'basc_thr':[], 'kmeans_thr':[], 'onestep_thr':[], 'shmulevich_thr':[]}
    
    for i in range(rows):
            k_means = K_Means(genes[i])
            basc_a, _ = BASC_A(genes[i])
            one_step = onestep(genes[i])
            shmulevich_ = shmulevich(genes[i])
            
            new_row = {'basc_thr':basc_a, 'kmeans_thr':k_means, 'onestep_thr':one_step, 'shmulevich_thr': shmulevich_}
            final_df.loc[len(final_df)] = new_row
        
    return [dash_table.DataTable(final_df.to_dict('records'), [{"name": i, "id": i} for i in final_df.columns],
                               page_size= 10), html.Button("Download CSV", id="btn_csv"),
        dcc.Download(id="download-dataframe-csv")]

@callback(
    Output("download-dataframe-csv", "data"),
    Input("btn_csv", "n_click"),
    prevent_initial_call=True,
)
def download_csv(n_click):
    return dcc.send_data_frame(final_df.to_csv, "thr.csv")
    
@app.callback(
    Output('dropdown-methods', 'children'),
    Input('datatable-interactivity', 'selected_rows'))
def display_selected_data(selected_rows):
    if not selected_rows:
        return None
    
    return  [html.H5('Select binarization method to binarize selected rows and visualize:'), 
            dcc.Dropdown(
                ['All','BASC A', 'K-Means', 'Onestep', 'Shmulevich'], 'All',
                placeholder="Select binarization method",
                id="dropdown-method",
                searchable=False)]

@app.callback(
    Output('heatmap-binarize', 'children'),
    Input('dropdown-method', 'value'),
    Input('datatable-interactivity', 'selected_rows'), 
    Input('stored-data','data'), 
    prevent_initial_call=True)
def heatmap_binarize(selected_method, selected_rows, data):
    if selected_rows is None:
        return None
    if selected_method is None:
        return None
    if selected_method == 'All':
        return None
    
    df = pd.DataFrame(data)
    
    binarize_vect = []
    labels = []
        
    # for each row binarize
    for row in selected_rows:
            selected = df.iloc[row]
            gene = selected.values
            sizeGene = len(gene)
            binarize = []
            
            # add gene label
            labels.append("Gene " + str(row+1))
            
            # get threshold of gene
            if(selected_method == 'K-Means'):
                thr = K_Means(gene)
            elif(selected_method == 'BASC A'):
                thr, _ = BASC_A(gene)
            elif(selected_method == 'Onestep'):
                thr = onestep(gene)
            else:
                thr = shmulevich(gene)

            # binarize based on thr
            for j in range(sizeGene):
                if(gene[j] <= thr):
                    binarize.append(0) 
                else:
                    binarize.append(1) 
                    
            # add binarization to array
            binarize_vect.append(binarize) 
        
    # get the selected rows values
    selected = df.iloc[selected_rows]
    genes = selected.values
        
    # create the figure heatmap with the genes, labels, and binarization
    # also add title text to the graph and return the figure
    data = go.Figure(data=go.Heatmap(
                    z=genes,
                    y = labels,
                    text=binarize_vect,
                    texttemplate="%{text}",
                    textfont={"size":20}))
    data.update_layout(
            title={
        'text': "Heatmap of Binarized Genes using " + selected_method,
            }
        ) 
    return dcc.Graph(figure=data)

@app.callback(
    Output('select-gene-binarize', 'children'),
    Input('dropdown-method', 'value'),
    Input('datatable-interactivity', 'selected_rows'), prevent_initial_call=True)
def select_gene_binarize(selected_method, selected_rows):
    if selected_method is None:
        return None

    if selected_rows is None:
        return None

    return [html.H5('Select gene to visualize:'), dcc.Dropdown(
        options=[{'label': 'Gene ' + str(row+1), 'value': row} for row in selected_rows], 
        value = selected_rows[0],
        placeholder="Select rows",
        id="dropdown-selected-rows")]

@app.callback(
    Output('graph-gene-binarize', 'children'),
    Input('dropdown-method', 'value'),
    Input('dropdown-selected-rows', 'value'), 
    Input('stored-data','data'),
    prevent_initial_call=True)
def graph_gene_algorithm(selected_method, selected_gene, data):
    if selected_gene is None:
        return None
    
    df = pd.DataFrame(data)
    
    selected = df.iloc[selected_gene]
    
    gene = selected.values
    sizeGene = len(gene)
    
    if(selected_method == 'BASC A'):
        thr , _ = BASC_A(gene)
        data = go.Figure(go.Scatter(x=np.arange(1,sizeGene+1), y=gene, name="Gene "+ str(selected_gene+1)))
        #data.add_hline(y=thr, line_width=3, line_dash="dash", line_color="green")
        data.add_trace(go.Scatter(x=np.arange(1,sizeGene+1), y=np.full(sizeGene, thr), line_dash="dash",
                                  mode="lines", name="BASC A Threshold"))
        
        data.update_layout(
            title={
        'text': "BASC A Threshold of Gene " + str(selected_gene+1),
            }
        )
        return dcc.Graph(figure=data)
    
    elif(selected_method == 'Onestep'):
        thr = onestep(gene)
        data = go.Figure(go.Scatter(x=np.arange(1,sizeGene+1), y=gene, name="Gene "+ str(selected_gene+1)))
        #data.add_hline(y=thr, line_width=3, line_dash="dash", line_color="green")
        data.add_trace(go.Scatter(x=np.arange(1,sizeGene+1), y=np.full(sizeGene, thr), line_dash="dash",
                                  mode="lines", name="Onestep Threshold"))
        
        data.update_layout(
            title={
        'text': "Onestep Threshold of Gene " + str(selected_gene+1),
            }
        )
        return dcc.Graph(figure=data)
    
    elif(selected_method == 'Shmulevich'):
        thr = shmulevich(gene)
        data = go.Figure(go.Scatter(x=np.arange(1,sizeGene+1), y=gene, name="Gene "+ str(selected_gene+1)))
        #data.add_hline(y=thr, line_width=3, line_dash="dash", line_color="green")
        data.add_trace(go.Scatter(x=np.arange(1,sizeGene+1), y=np.full(sizeGene, thr), line_dash="dash",
                                  mode="lines", name="Shmulevich Threshold"))
        
        data.update_layout(
            title={
        'text': "Shmulevich Threshold of Gene " + str(selected_gene+1),
            }
        )
        return dcc.Graph(figure=data)
        
    elif(selected_method == 'K-Means'):
        thr = K_Means(gene)
        data = go.Figure(go.Scatter(x=np.arange(1,sizeGene+1), y=gene, name="Gene "+ str(selected_gene+1)))
        #data.add_hline(y=thr, line_width=3, line_dash="dash", line_color="green")
        data.add_trace(go.Scatter(x=np.arange(1,sizeGene+1), y=np.full(sizeGene, thr), line_dash="dash",
                                  mode="lines", name="K-Means Threshold"))
        
        data.update_layout(
            title={
        'text': "K-Means Threshold of Gene " + str(selected_gene+1),
            }
        )
        
        data2 = go.Figure(go.Scatter(x=gene, y=np.ones(sizeGene, dtype=int),mode = 'markers', name="Gene"))
        data2.add_trace(go.Scatter(x=np.full(sizeGene+2,thr), y=np.arange(0,sizeGene+2)
                                  ,mode="lines", line_dash="dash", line_color="green", name="Threshold"))
        #data2.add_vline(x=thr, line_width=3, line_dash="dash", line_color="green")
        
        data2.update_layout(
            title={
        'text': "K-Means Clusters of Gene " + str(selected_gene+1),
            }
        )
        
        return [dcc.Graph(figure=data), dcc.Graph(figure=data2)]
    
    elif(selected_method == "All"):
        thrBasc, _ = BASC_A(gene)
        thrKmeans = K_Means(gene)
        thrShmu = shmulevich(gene)
        thrOnestep = onestep(gene)
        
        data = go.Figure(go.Scatter(x=np.arange(1,sizeGene+1), y=gene, name="Gene "+ str(selected_gene+1)))
        data.add_trace(go.Scatter(x=np.arange(1,sizeGene+1), y=np.full(sizeGene, thrBasc), line_dash="dash",
                                  mode="lines", name="BASC A Threshold"))
        data.add_trace(go.Scatter(x=np.arange(1,sizeGene+1), y=np.full(sizeGene, thrKmeans), line_dash="dash",
                                  mode="lines", name="K-Means Threshold"))
        data.add_trace(go.Scatter(x=np.arange(1,sizeGene+1), y=np.full(sizeGene, thrShmu), line_dash="dash",
                                  mode="lines", name="Shmulevich Threshold"))
        data.add_trace(go.Scatter(x=np.arange(1,sizeGene+1), y=np.full(sizeGene, thrOnestep), line_dash="dash",
                                  mode="lines", name="Onestep Threshold"))
        data.update_layout(
            title={
        'text': "Threshold of Gene " + str(selected_gene+1) + " using all methods",
            }
        )
            
        return dcc.Graph(figure=data)
    
@app.callback(
    Output('select-basc-discontinuity', 'children'),
    Input('dropdown-method', 'value'),
    Input('dropdown-selected-rows', 'value'), 
     Input('stored-data','data'), 
    prevent_initial_call=True)
def select_basc_discontinuity(selected_method, selected_gene, data):
    if selected_method is None:
        return None

    if selected_gene is None:
        return None
    
    df = pd.DataFrame(data)
    
    if(selected_method == 'BASC A'):
    
        selected = df.iloc[selected_gene]
        gene = selected.values
        sizeGene = len(gene)
        
        options = [{'label':'All','value':0}]
        for i in range(sizeGene-2):
            options.append({'label':'Discontinuity ' + str(i+1), 'value': i+1})
            
        return [html.H5('Select discontinuity to visualize:'), dcc.Dropdown(
            options=options,
            placeholder="Select discontinuity",
            id="dropdown-selected-discontinuity", value=0)]
    
@app.callback(
    Output('graph-gene-discontinuity', 'children'),
    Input('dropdown-selected-rows', 'value'),
    Input('dropdown-selected-discontinuity', 'value'), 
    Input('dropdown-method', 'value'), 
     Input('stored-data','data'), 
    prevent_initial_call=True)
def graph_discontinuity(selected_gene, selected_discontinuity, method, data):
    if selected_discontinuity is None:
        return None

    if selected_gene is None:
        return None
    if method != "BASC A":
        return None
    
    df = pd.DataFrame(data)
    
    selected = df.iloc[selected_gene]
    gene = np.sort(selected.values)
    sizeGene = len(gene)
    
    if(selected_discontinuity > 0):
            _ , P = BASC_A(gene) 
            
            sizeP = len(P)
            y = np.sort(gene)

            x_dis = []
            y_dis = []

            for i in range(sizeP):
                if(P[i][selected_discontinuity-1]>0):

                    index = P[i][selected_discontinuity-1]

                    if not x_dis:
                        x_dis.append(0)
                    else:
                        x_dis.append(P[i-1][selected_discontinuity-1])
                    x_dis.append(index)
                    x_dis.append(None)

                    y_dis.append(y[index-1])
                    y_dis.append(y[index-1])
                    y_dis.append(None)
                    #print(index)

            x_dis.append(index)
            x_dis.append(len(gene))
            x_dis.append(None)

            y_dis.append(y[index])
            y_dis.append(y[index])
            y_dis.append(None)

            data2 = go.Figure(go.Scatter(x=x_dis, y=y_dis))
            
            for i in range(sizeP):
                if(P[i][selected_discontinuity-1]>0):
                    index = P[i][selected_discontinuity-1]
                    data2.add_trace(go.Scatter(x=[index, index], y=[y[index], y[index-1]], line_dash="dash",
                                  mode="lines", line_color="green", name="Jump"))
            
            data2.update_layout(
            title={
                'text': "BASC A Discontinuity " + str(selected_discontinuity) + " of Gene " + str(selected_gene+1),
                    }
                )
            
            return dcc.Graph(figure=data2)
                
    else:
        x = np.arange(1,sizeGene+1)
        y = gene
        x_dis = []
        y_dis = []
        for i in range(len(x)):
                x_dis.append(i)
                x_dis.append(i+1)
                x_dis.append(None)
                
                y_dis.append(y[i])
                y_dis.append(y[i])
                y_dis.append(None)
            
        data2 = go.Figure(go.Scatter(x=x_dis, y=y_dis))
        
        data2.update_layout(
            title={
                'text': "BASC A Discontinuities of Gene " + str(selected_gene+1),
                    })
            
        return dcc.Graph(figure=data2)
    

@app.callback(
    Output('tolerance-dropdown', 'children'),
    Input('dropdown-selected-rows', 'value'),
    Input('dropdown-method', 'value'), 
    prevent_initial_call=True)
def interpolation_dropdown(selected_row, selected_method):
    if selected_row is None:
        return None
    if selected_method == 'All':
        return None
    return  [html.H5('Select tolerance for interpolation:'), 
            dcc.Dropdown(
                [0.1, 0.01, 0.001, 0.0001, 0.00001],
                placeholder="Select tolerance",
                id="dropdown-tolerance",
                searchable=False)]


@app.callback(
    Output('interpolation-graph', 'children'),
    Input('dropdown-selected-rows', 'value'),
    Input('dropdown-method', 'value'),
    Input('dropdown-tolerance', 'value'),
    Input('stored-data','data'),
    prevent_initial_call=True)
def interpolation_graph(selected_row, selected_method, tolerance, data):
    if selected_row is None:
        return None
    if tolerance is None:
        return None
    if selected_method is None:
        return None
    
    df = pd.DataFrame(data)
    
    if(selected_method != 'All'):
        
        thr = 0
        selected = df.iloc[selected_row]
        gene = selected.values
        sizeGene = len(gene)
        
        thrMin, thrMax, conver, samples = interpolationConverge(gene, selected_method, tolerance)
        
        #return thrMin, thrMax
        
        if(selected_method == 'K-Means'):
            thr = K_Means(gene)
            
        elif(selected_method == 'BASC A'):
            thr, _ = BASC_A(gene)
            
        elif(selected_method == 'Onestep'):
            thr = onestep(gene)
            
        else:
            thr = shmulevich(gene)

        data = go.Figure(go.Scatter(x=np.arange(1,sizeGene+1), y=gene, name="Gene"))
        data.add_trace(go.Scatter(x=np.arange(0,sizeGene+2), y=np.full(sizeGene+2,thr) 
                                  ,mode="lines", line_dash="dash", line_color="green", name="Threshold"))
        data.add_hline(y=thrMin, line_width=3, line_dash="dash", line_color="red")
        data.add_hline(y=thrMax, line_width=3, line_dash="dash", line_color="red")
        
        data.update_layout(
            title={
        'text': "Gene "+ str(selected_row+1) + 
                                    " Sample Size " + str(samples) 
                                    + " " + selected_method
            }
        )
        return dcc.Graph(figure=data)
    
    else:
        return None

@app.callback(
    Output('interpolation-heatmap', 'children'),
    Input('dropdown-selected-rows', 'value'),
    Input('dropdown-method', 'value'), 
    Input('stored-data','data'),
    prevent_initial_call=True)
def interpolation_heatmap(selected_row, method, data):
    if selected_row is None:
        return None
    if method is None or method == 'All':
        return None
    
    df = pd.DataFrame(data)
    
    selected = df.iloc[selected_row]
    gene = selected.values
    thr, sample = three_interpolation(gene , method)
    sizeGene = len(gene)
    inter = 1
    binarize_vect = []
    labels = []
    
    for i in range(len(thr)):
        binarize = []
        labels.append("Size " + str(sample[i]) + " Thr: " + str("{:.2f}".format(thr[i])))
        inter += 1
        
        for j in range(sizeGene):
                if(gene[j] <= thr[i]):
                    binarize.append(0) 
                else:
                    binarize.append(1)     
        binarize_vect.append(binarize) 
        
    genes = []
    genes.append(gene)
    genes.append(gene)
    genes.append(gene)
    
    data2 = go.Figure(data=go.Heatmap(
                        z=genes,
                        y = labels,
                        text=binarize_vect,
                        texttemplate="%{text}",
                        textfont={"size":20}))
    data2.update_layout(
                title={
            'text': method + " Heatmap of Three Samples of Gene " + str(selected_row+1),
                }
            )
    
    return dcc.Graph(figure=data2)


@app.callback(
    Output('rules-dropdown', 'children'),
    Input('datatable-interactivity', 'selected_rows'),
    Input('dropdown-method', 'value'), 
    prevent_initial_call=True)
def rules_dropdown(selected_rows, selected_method):
    if selected_rows is None:
        return None
    if selected_method is None:
        return None
    
    if len(selected_rows) >= 3:
        return [ html.H5('Select A rule is A = B:'), 
                dcc.Dropdown(
                options=[{'label': 'Gene ' + str(row+1), 'value': row} for row in selected_rows], 
                placeholder="Select Gene A",
                id="selected-A"),

                html.H5('Select B rule is B = ~A:'),
                dcc.Dropdown(
                options=[{'label': 'Gene ' + str(row+1), 'value': row} for row in selected_rows], 
                placeholder="Select Gene B",
                id="selected-B"),

                html.H5('Select C rule is C = A or B:'),
                dcc.Dropdown(
                options=[{'label': 'Gene ' + str(row+1), 'value': row} for row in selected_rows], 
                placeholder="Select Gene C",
                id="selected-C")]
    else:
        return None

@app.callback(
    Output('graph_rules', 'src'),
    Input('selected-A', 'value'),
    Input('selected-B', 'value'),
    Input('selected-C', 'value'),
    Input('dropdown-method', 'value'), 
    Input('stored-data','data'),
    prevent_initial_call=True)
def process_rules(A, B, C, selected_method, data):
    if A is None or B is None or C is None:
        return None
    if selected_method is None:
        return None
    
    selected = [A, B, C]
    
    if(selected_method == 'All'):
        fig = create_boolean_network_graph(selected, data)
        
        # Save it to a temporary buffer.
        buf = BytesIO()
        fig.savefig(buf, format="png")
        # Embed the result in the html output.
        fig_data = base64.b64encode(buf.getbuffer()).decode("ascii")
        fig_bar_matplotlib = f'data:image/png;base64,{fig_data}'

        return fig_bar_matplotlib
    
    else:
        fig = graph_three_interpolation_network(selected, selected_method, data)
        
        # Save it to a temporary buffer.
        buf = BytesIO()
        fig.savefig(buf, format="png")
        # Embed the result in the html output.
        fig_data = base64.b64encode(buf.getbuffer()).decode("ascii")
        fig_bar_matplotlib = f'data:image/png;base64,{fig_data}'

        return fig_bar_matplotlib
        
    
    return selected

if __name__ == '__main__':
    app.run_server(debug=False)