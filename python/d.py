from d import Dash, dash_table, dcc, html, Input, Output, callback
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



from scipy.interpolate import splev, splrep

def interpolation(gene, method, tolerance, limit):
    
    converge = False
    samples = 10
    thr_array = []
    tMin = 0
    tMax = 0
    thr = 0
    
    while(~converge and samples != limit):
        n = len(gene)
        x = np.arange(1, n + 1, dtype=int)
        y = gene

        xx = np.arange(1, n + 1, (n) / (samples))

        p = splrep(x, y, s=0)
        yy = splev(xx, p)
        
        if(method == 'K-Means'):
            thr = K_Means(yy)
            
        elif(method == 'BASC A'):
            thr, _ = BASC_A(yy)
            
        elif(method == 'Onestep'):
            thr = onestep(yy)
            
        else:
            thr = shmulevich(yy)
        
        thr_array.append(thr)
        
        n_thr = len(thr_array)
        
        for i in range(n_thr):
            for j in range(n_thr):
                if(i != j):
                    difference = abs(thr_array[i] - thr_array[j])
                    
                    if(difference <= tolerance):
                        converge = True
                        break
                        
            if(converge or samples != limit):
                tMin = min(thr_array)
                tMax = max(thr_array)
                break
                
        if(~converge):
            samples += 10
            
            
    
    return tMin, tMax
        

from scipy.interpolate import splev, splrep

def interpolationConverge(gene, method, tolerance):
    
    converge = False
    samples = 10
    thr_array = []
    tMin = 0
    tMax = 0
    thr = 0
    conver = 0
    #change later
    #limit = 10000
    limit = 10000
    
    while(~converge and samples != limit):
        n = len(gene)
        x = np.arange(1, n + 1, dtype=int)
        y = gene

        xx = np.arange(1, n + 1, (n) / (samples))

        p = splrep(x, y, s=0)
        yy = splev(xx, p)
        
        if(method == 'K-Means'):
            thr = K_Means(yy)
            
        elif(method == 'BASC A'):
            thr, _ = BASC_A(yy)
            
        elif(method == 'Onestep'):
            thr = onestep(yy)
            
        else:
            thr = shmulevich(yy)
        
        thr_array.append(thr)
        
        n_thr = len(thr_array)
        
        for i in range(n_thr):
            for j in range(n_thr):
                if(i != j):
                    difference = abs(thr_array[i] - thr_array[j])
                    
                    if(difference <= tolerance):
                        converge = True
                        conver = thr_array[j]
                        break
                
        if(~converge):
            samples += 10

        
            


    tMin = min(thr_array)
    tMax = max(thr_array)
    
    return tMin, tMax, conver, samples



import dash_bootstrap_components as dbc 

df = pd.read_csv('HIVIn(Matlab).csv')

col_names = {'basc_thr':[], 'kmeans_thr':[], 'onestep_thr':[], 'shmulevich_thr':[]}
final_df = pd.DataFrame(col_names)

app = Dash(__name__)
server = app.server

app.layout = html.Div([
    html.Div([
                html.H4('Binarize Genes Visualization'), 
        
                html.Div(dash_table.DataTable(
                    id='datatable-interactivity',
                    columns=[
                        {"name": i, "id": i} for i in df.columns
                    ],
                    data=df.to_dict('records'),
                    column_selectable="single",
                    row_selectable="multi",
                    selected_columns=[],
                    selected_rows=[],
                    page_action="native",
                    page_current= 0,
                    page_size= 10,
                )),
                html.H5('Select option in order to show all the thresholds by each algorithm of each gene:'),
                html.Div([
                dcc.Dropdown(
                            [{'label': 'Binarize All Genes', 'value':'all'}],
                            placeholder="Select to binarize all genes and get thresholds",
                            id="dropdown-binarize-all",
                            searchable=False)
                ], style={"width": "25%"}),
                html.Div(id='binarize-all'),
                html.Hr()
               
          
    ]),
    
    html.Div([
        html.Div([html.Div(id='dropdown-methods', style={"width": "25%"}),
                  html.Div(id='heatmap-binarize')
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

@app.callback(
    Output('binarize-all', 'children'),
    Input('dropdown-binarize-all', 'value'))
def binarize_all(all_rows):
    if all_rows is None:
        return None
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
    Input("btn_csv", "click"),
    prevent_initial_call=True,
)
def download_csv(click):
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
    Input('datatable-interactivity', 'selected_rows'), prevent_initial_call=True)
def heatmap_binarize(selected_method, selected_rows):
    if selected_rows is None:
        return None
    if selected_method is None:
        return None
    
    #selected = df.iloc[selected_rows]
    #gene = selected.values
    #sizeGene = len(gene)
    
    if(selected_method == "K-Means"):
        binarize_vect = []
        labels = []

        for row in selected_rows:
            selected = df.iloc[row]
            gene = selected.values
            sizeGene = len(gene)
            binarize = []
            labels.append("Gene " + str(row+1))
            thr = K_Means(gene)

            for j in range(sizeGene):
                if(gene[j] <= thr):
                    binarize.append(0) 
                else:
                    binarize.append(1)     

            binarize_vect.append(binarize) 
            
        selected = df.iloc[selected_rows]
        genes = selected.values
        
        data = go.Figure(data=go.Heatmap(
                    z=genes,
                    y = labels,
                    text=binarize_vect,
                    texttemplate="%{text}",
                    textfont={"size":20}))
        data.update_layout(
            title={
        'text': "Heatmap of Binarized Genes using K-Means",
            }
        ) 
        return dcc.Graph(figure=data)
        
    elif(selected_method == "BASC A"):   
        binarize_vect = []
        labels = []

        for row in selected_rows:
            selected = df.iloc[row]
            gene = selected.values
            sizeGene = len(gene)
            binarize = []
            labels.append("Gene " + str(row+1))
            thr, _ = BASC_A(gene)

            for j in range(sizeGene):
                if(gene[j] <= thr):
                    binarize.append(0) 
                else:
                    binarize.append(1)     

            binarize_vect.append(binarize) 
            
        selected = df.iloc[selected_rows]
        genes = selected.values
        
        data = go.Figure(data=go.Heatmap(
                    z=genes,
                    y = labels,
                    text=binarize_vect,
                    texttemplate="%{text}",
                    textfont={"size":20}))
        data.update_layout(
            title={
        'text': "Heatmap of Binarized Genes using BASC A",
            }
        )
        return dcc.Graph(figure=data)
    
    elif(selected_method == "Onestep"):
        binarize_vect = []
        labels = []

        for row in selected_rows:
            selected = df.iloc[row]
            gene = selected.values
            sizeGene = len(gene)
            binarize = []
            labels.append("Gene " + str(row+1))
            thr = onestep(gene)

            for j in range(sizeGene):
                if(gene[j] <= thr):
                    binarize.append(0) 
                else:
                    binarize.append(1)     

            binarize_vect.append(binarize) 
            
        selected = df.iloc[selected_rows]
        genes = selected.values
        
        data = go.Figure(data=go.Heatmap(
                    z=genes,
                    y = labels,
                    text=binarize_vect,
                    texttemplate="%{text}",
                    textfont={"size":20}))
        data.update_layout(
            title={
        'text': "Heatmap of Binarized Genes using Onestep",
            }
        )
        return dcc.Graph(figure=data)
    
    elif(selected_method == "Shmulevich"):
        binarize_vect = []
        labels = []

        for row in selected_rows:
            selected = df.iloc[row]
            gene = selected.values
            sizeGene = len(gene)
            binarize = []
            labels.append("Gene " + str(row+1))
            thr = shmulevich(gene)

            for j in range(sizeGene):
                if(gene[j] <= thr):
                    binarize.append(0) 
                else:
                    binarize.append(1)     

            binarize_vect.append(binarize) 
            selected = df.iloc[selected_rows]
            genes = selected.values

            data = go.Figure(data=go.Heatmap(
                        z=genes,
                        y = labels,
                        text=binarize_vect,
                        texttemplate="%{text}",
                        textfont={"size":20}))
            data.update_layout(
                title={
            'text': "Heatmap of Binarized Genes using Shmulevich",
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

    if not selected_rows:
        return None

    return [html.H5('Select gene to visualize:'), dcc.Dropdown(
        options=[{'label': 'Gene ' + str(row+1), 'value': row} for row in selected_rows], 
        value = selected_rows[0],
        placeholder="Select rows",
        id="dropdown-selected-rows")]

@app.callback(
    Output('graph-gene-binarize', 'children'),
    Input('dropdown-method', 'value'),
    Input('dropdown-selected-rows', 'value'), prevent_initial_call=True)
def graph_gene_algorithm(selected_method, selected_gene):
    if not selected_gene:
        return None
    
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
        #data.add_hline(y=thrBasc, line_width=3, line_dash="dash", line_color="green")
        #data.add_annotation(xref="paper", y=thrBasc, text="BASC A", showarrow=False)
        
        #data.add_hline(y=thrKmeans, line_width=3, line_dash="dash", line_color="red")
        #data.add_annotation(xref="paper", y=thrKmeans, text="KMeans", showarrow=False)
        
        #data.add_hline(y=thrOnestep, line_width=3, line_dash="dash", line_color="orange")
        #data.add_annotation(xref="paper", y=thrOnestep, text="OneStep", showarrow=False)
        
        #data.add_hline(y=thrShmu, line_width=3, line_dash="dash", line_color="purple")
        #data.add_annotation(xref="paper", y=thrShmu, text="Shmulevich", showarrow=False)
        
        data.update_layout(
            title={
        'text': "Threshold of Gene " + str(selected_gene+1) + " using all methods",
            }
        )
            
        return dcc.Graph(figure=data)
    
@app.callback(
    Output('select-basc-discontinuity', 'children'),
    Input('dropdown-method', 'value'),
    Input('dropdown-selected-rows', 'value'), prevent_initial_call=True)
def select_basc_discontinuity(selected_method, selected_gene):
    if selected_method is None:
        return None

    if not selected_gene:
        return None
    
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
    Input('dropdown-method', 'value'), prevent_initial_call=True)
def graph_discontinuity(selected_gene, selected_discontinuity, method):
    if selected_discontinuity is None:
        return None

    if not selected_gene:
        return None
    if method != "BASC A":
        return None
    
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
    if not selected_row:
        return None
    if selected_method == 'All':
        return None
    return  [html.H5('Select tolerance for interpolation:'), 
            dcc.Dropdown(
                [0.1, 0.01, 0.001, 0.0001, 0.00001], 0.1,
                placeholder="Select tolerance",
                id="dropdown-tolerance",
                searchable=False)]


@app.callback(
    Output('interpolation-graph', 'children'),
    Input('dropdown-selected-rows', 'value'),
    Input('dropdown-method', 'value'),
    Input('dropdown-tolerance', 'value'), 
    prevent_initial_call=True)
def interpolation_graph(selected_row, selected_method, tolerance):
    if not selected_row:
        return None
    if not tolerance:
        return None
    if not selected_method:
        return None
    
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
                                    " interpolation of sample size " + str(samples) 
                                    + " and tolerance " + str(tolerance),
            }
        )
        return dcc.Graph(figure=data)
    
    else:
        return None


if __name__ == '__main__':
    app.run_server(debug=True, port=8004)