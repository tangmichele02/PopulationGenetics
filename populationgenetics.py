# import web app modules
import streamlit as st
import plotly.graph_objs as go


def makePlot(sims, AA, Aa, aa, pop, gen, stA):
    import plotly.graph_objs as go
    import numpy as np
    from numpy import random as rd
    from plotly.graph_objs import Scatter, Layout, Figure, Histogram
    from plotly import tools
    from plotly.subplots import make_subplots
 

    # Allele frequencies
    initA = 0.50
    inita = 0.50

    #Fitnesses
    fAA = AA
    fAa = Aa
    faa = aa

    #Pop Size
    pop = pop

    #Number of generations
    generation = int(gen)

    #Number of simulations
    simulation = int(sims)

    # to store all iterations to graph 
    allgraphs = []
    # to store all final allelic frequency values for A 
    finalA = []
    # to graph generations 
    x_vals = [i for i in range(generation)]
    fig = make_subplots(rows = 2, cols = 1)

    for i in range(simulation): 
        # frequencies for A and a
        p = initA
        q = inita
        
        valsA = []
        for j in range(generation): 

            # expected frequency of reproducing adults 
            # wbar is for calculating the denominator later 
            wbar = ((p**2) * fAA) + (2*p*q * fAa) + ((q**2) * faa) 
            prime_AA_x = ((p**2) * fAA) / wbar
            prime_Aa_x = (2*p*q * fAa) / wbar
            prime_aa_x = ((q**2) * faa) / wbar

            # observed count of reproducing adults 
            prime_AA = rd.binomial(pop, prime_AA_x)
            prime_Aa = rd.binomial(pop, prime_Aa_x)
            prime_aa = rd.binomial(pop, prime_aa_x)

            # summing for actual population 
            actual_pop = prime_AA + prime_Aa + prime_aa

            # actual count of reproducing adults 
            prime_AA = prime_AA / actual_pop
            prime_Aa = prime_Aa / actual_pop
            prime_aa = prime_aa / actual_pop 

            # allelic frequency in reproducing adults 
            prime_A = prime_AA + 0.5 * prime_Aa 
            prime_a = prime_aa + 0.5 * prime_Aa
            
            # add allelic frequency for A into list 
            valsA.append(prime_A)

            # reinitialize frequencies for A and a to calculate next generation 
            p = prime_A
            q = prime_a
        
        finalA.append(prime_A)
        graph = Scatter(x = x_vals, y = valsA, mode = 'lines', name = "sim " + str(j + 1))
        fig.append_trace(graph, 1, 1)

    histo = Histogram(x=finalA, opacity=0.5, autobinx=False, xbins=dict(
        start=-0.01,
        end=1.01,
        size=0.01
    ))

    fig.append_trace(histo, 2, 1)
    fig['layout']['xaxis2'].update(range=[-0.01, 1.01])
    fig['layout']['yaxis1'].update(range=[0, 1])
    fig['layout'].update(height=600, width=800)
    fig['layout'].update(showlegend=False)

    
    return fig

with st.sidebar:
    nSim = st.slider('Number of Simulations:',min_value=1,max_value=100,step=1)
    AA = st.slider('Fitness of AA:',min_value=0.,max_value=1.,step=0.05,value=1.)
    Aa = st.slider('Fitness of Aa:',min_value=0.,max_value=1.,step=0.05,value=1.)
    aa = st.slider('Fitness of aa:',min_value=0.,max_value=1.,step=0.05,value=1.)
    pop = st.select_slider('Population Size:',[10,50,100,500,1000])
    gen = st.slider("Number of Generations:",min_value=100,max_value=1000,step=100)
    stA = st.slider("Starting frequency of A:", min_value=0.01,max_value=0.99,step=0.01, value=0.5)

plot = makePlot(nSim, AA, Aa, aa, pop, gen,stA)
st.write(plot)



